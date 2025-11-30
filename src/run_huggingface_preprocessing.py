import hashlib
import mimetypes
import os
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Tuple

import click
import pandas as pd
from datasets import Dataset


TABULAR_EXTENSIONS = {".csv", ".tsv", ".txt"}
DEFAULT_ROOT = Path("data").resolve()
DEFAULT_OUT = DEFAULT_ROOT / "hf"


@dataclass
class FileRecord:
    task_id: str
    relative_path: str
    absolute_path: str
    extension: str
    size_bytes: int
    is_gz: bool
    mime_type: Optional[str]
    sha256: Optional[str]


def iter_task_dirs(root_dir: Path) -> Iterable[Tuple[str, Path]]:
    for child in sorted(root_dir.iterdir()):
        if child.is_dir():
            yield child.name, child


def detect_delimiter(sample: str) -> Optional[str]:
    if "\t" in sample and "," in sample:
        return "\t"
    if "\t" in sample:
        return "\t"
    if "," in sample:
        return ","
    if "  " in sample:
        return r"\s+"
    return None


def is_tabular_file(path: Path) -> bool:
    ext = path.suffix.lower()
    if ext not in {'.csv', '.tsv'}:
        return False
    return True

def is_fastq_file(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".fastq") or name.endswith(".fq") or name.endswith(".fastq.gz") or name.endswith(".fq.gz"):
        return True
    return False

def is_fasta_file(path: Path) -> bool:
    name = path.name.lower()
    if (
        name.endswith(".fa")
        or name.endswith(".fasta")
        or name.endswith(".fna")
        or name.endswith(".fa.gz")
        or name.endswith(".fasta.gz")
        or name.endswith(".fna.gz")
    ):
        return True
    return False

def is_vcf_file(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".vcf") or name.endswith(".vcf.gz") or name.endswith(".eff.vcf") or name.endswith(".eff.vcf.gz"):
        return True
    return False

def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def build_file_index(
    task_id: str, task_dir: Path, compute_checksums: bool, checksum_bytes: Optional[int]
) -> Dataset:
    records: List[FileRecord] = []
    for file_path in sorted(task_dir.rglob("*")):
        if not file_path.is_file():
            continue
        rel = str(file_path.relative_to(DEFAULT_ROOT))
        abs_path = str(file_path.resolve())
        ext = file_path.suffix.lower()
        try:
            size = file_path.stat().st_size
        except FileNotFoundError:
            size = 0
        is_gz = abs_path.endswith(".gz")
        mime, _ = mimetypes.guess_type(abs_path)
        sha = compute_sha256(file_path, max_bytes=checksum_bytes) if compute_checksums else None
        records.append(
            FileRecord(
                task_id=task_id,
                relative_path=rel,
                absolute_path=abs_path,
                extension=ext,
                size_bytes=size,
                is_gz=is_gz,
                mime_type=mime,
                sha256=sha,
            )
        )
    df = pd.DataFrame([r.__dict__ for r in records])
    return Dataset.from_pandas(df, preserve_index=False)


def sanitize_for_dir(name: str) -> str:
    name = name.replace(os.sep, "__").replace("/", "__")
    return "".join(ch if (ch.isalnum() or ch in ("_", "-", ".")) else "_" for ch in name)


def fastq_reader(path: str, max_reads: Optional[int], source_sha256: str) -> Iterator[dict]:
    is_gz = path.endswith(".gz")
    opener = gzip.open if is_gz else open
    emitted = 0
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not seq or not plus or not qual:
                break
            hdr = header.strip()
            read_id = hdr[1:].split()[0] if hdr.startswith("@") and len(hdr) > 1 else hdr
            yield {
                "read_id": read_id,
                "sequence": seq.strip(),
                "quality": qual.strip(),
                "source_file": os.path.basename(path),
                "source_sha256": source_sha256,
            }
            emitted += 1
            if max_reads is not None and max_reads > 0 and emitted >= max_reads:
                break


def fasta_reader(path: str, max_seqs: Optional[int], source_sha256: str) -> Iterator[dict]:
    is_gz = path.endswith(".gz")
    opener = gzip.open if is_gz else open
    emitted = 0
    seq_id: Optional[str] = None
    desc: str = ""
    seq_chunks: List[str] = []
    def emit_current():
        if seq_id is None:
            return None
        sequence = "".join(seq_chunks)
        return {
            "seq_id": seq_id,
            "description": desc,
            "sequence": sequence,
            "source_file": os.path.basename(path),
            "source_sha256": source_sha256,
        }
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                # Emit previous
                record = emit_current()
                if record:
                    yield record
                    emitted += 1
                    if max_seqs is not None and max_seqs > 0 and emitted >= max_seqs:
                        return
                # Start new
                header = line[1:].rstrip("\n")
                parts = header.split(None, 1)
                seq_id = parts[0] if parts else ""
                desc = parts[1] if len(parts) > 1 else ""
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        # Emit last
        record = emit_current()
        if record:
            yield record


def vcf_reader(path: str, source_sha256: str) -> Iterator[dict]:
    is_gz = path.endswith(".gz")
    opener = gzip.open if is_gz else open
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        header_cols: Optional[List[str]] = None
        for raw in fh:
            if not raw:
                continue
            if raw.startswith("##"):
                continue
            if raw.startswith("#CHROM"):
                header_cols = raw.lstrip("#").strip().split("\t")
                continue
            if raw.startswith("#"):
                continue
            if header_cols is None:
                # Fallback to standard 8+ columns if header missing
                header_cols = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
            fields = raw.rstrip("\n").split("\t")
            # Pad fields if fewer than header
            if len(fields) < len(header_cols):
                fields = fields + [""] * (len(header_cols) - len(fields))
            # Map standard columns
            chrom = fields[0] if len(fields) > 0 else ""
            pos = fields[1] if len(fields) > 1 else ""
            vid = fields[2] if len(fields) > 2 else ""
            ref = fields[3] if len(fields) > 3 else ""
            alt = fields[4] if len(fields) > 4 else ""
            qual = fields[5] if len(fields) > 5 else ""
            flt = fields[6] if len(fields) > 6 else ""
            info = fields[7] if len(fields) > 7 else ""
            fmt = fields[8] if len(fields) > 8 else ""
            samples = fields[9:] if len(fields) > 9 else []
            yield {
                "chrom": chrom,
                "pos": int(pos) if pos.isdigit() else None,
                "id": vid,
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": flt,
                "info": info,
                "format": fmt,
                "samples": samples,
                "source_file": os.path.basename(path),
                "source_sha256": source_sha256,
            }


def read_tabular(path: Path, max_rows: Optional[int]) -> pd.DataFrame:
    ext = path.suffix.lower()
    sep: Optional[str] = None
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            head = "".join([fh.readline() for _ in range(3)])
    except Exception:
        head = ""
    if ext == ".tsv":
        sep = "\t"
    elif ext == ".csv":
        sep = ","
    else:
        sep = detect_delimiter(head) or ","
    nrows = max_rows if max_rows and max_rows > 0 else None
    df = pd.read_csv(
        path,
        sep=sep,
        engine="python",
        on_bad_lines="skip",
        nrows=nrows,
    )
    return df


def save_dataset(ds: Dataset, out_dir: Path, overwrite: bool) -> None:
    if out_dir.exists():
        if not overwrite:
            return
        for p in out_dir.iterdir():
            if p.is_file():
                p.unlink()
            elif p.is_dir():
                for sub in p.rglob("*"):
                    if sub.is_file():
                        sub.unlink()
                try:
                    p.rmdir()
                except OSError:
                    pass
    out_dir.mkdir(parents=True, exist_ok=True)
    ds.save_to_disk(str(out_dir))


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--root-dir",
    type=click.Path(path_type=Path),
    default=str(DEFAULT_ROOT),
    show_default=True,
    help="Root directory containing dataset subfolders.",
)
@click.option(
    "--out-dir",
    type=click.Path(path_type=Path),
    default=str(DEFAULT_OUT),
    show_default=True,
    help="Directory to write Hugging Face datasets.",
)

@click.option(
    "--max-rows",
    type=int,
    default=0,
    show_default=True,
    help="Limit rows when converting tabular files (0 means all rows).",
)
@click.option(
    "--overwrite/--no-overwrite",
    default=True,
    show_default=True,
    help="Overwrite existing converted datasets.",
)
@click.option(
    "--flat-output/--nested-output",
    default=True,
    show_default=True,
    help="Save each dataset as its own sibling directory directly under out-dir.",
)
@click.option(
    "--max-fastq-reads",
    type=int,
    default=0,
    show_default=True,
    help="Limit number of reads per FASTQ file (0 means all reads).",
)
@click.option(
    "--max-fasta-seqs",
    type=int,
    default=0,
    show_default=True,
    help="Limit number of sequences per FASTA file (0 means all sequences).",
)
def main(
    root_dir: Path,
    out_dir: Path,
    max_rows: int,
    overwrite: bool,
    flat_output: bool,
    max_fastq_reads: int,
    max_fasta_seqs: int,
) -> None:
    root_dir = root_dir.resolve()
    out_dir = out_dir.resolve()
    if not root_dir.exists():
        raise SystemExit(f"Root directory not found: {root_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)

    for task_id, task_path in iter_task_dirs(root_dir):
        try:
            if task_path.resolve() == out_dir:
                continue
        except Exception:
            pass
        click.echo(f"[{task_id}] Processing directory: {task_path}")


        tables_out_base = out_dir / task_id / "tables"
        fastq_out_base = out_dir / task_id / "fastq"
        fasta_out_base = out_dir / task_id / "fasta"
        vcf_out_base = out_dir / task_id / "vcf"
        converted_tables = 0
        converted_fastq = 0
        converted_fasta = 0
        converted_vcf = 0
        for file_path in sorted(task_path.rglob("*")):
            if not file_path.is_file():
                continue
            # Tabular conversion
            if is_tabular_file(file_path):
                try:
                    sha = compute_sha256(file_path)
                    df = read_tabular(file_path, max_rows=max_rows if max_rows > 0 else None)
                    if df.empty:
                        continue
                    df.insert(0, "source_sha256", sha)
                    ds = Dataset.from_pandas(df, preserve_index=False)
                    rel_to_task = file_path.relative_to(task_path)
                    if flat_output:
                        rel_str = str(rel_to_task)
                        if rel_str.endswith(".gz"):
                            rel_str = rel_str[:-3]
                        rel_no_ext = str(Path(rel_str).with_suffix(""))
                        flat_name = f"{task_id}__tables__{sanitize_for_dir(rel_no_ext)}"
                        out_subdir = out_dir / flat_name
                    else:
                        stem = file_path.stem
                        subdir_parts = list(rel_to_task.parts[:-1])  # parent folders under task
                        safe_subdir = Path(*subdir_parts) if subdir_parts else Path()
                        out_subdir = tables_out_base / safe_subdir / stem
                    save_dataset(ds, out_subdir, overwrite=overwrite)
                    converted_tables += 1
                except Exception as e:
                    click.echo(f"[{task_id}] Skipping table {file_path.name}: {e}", err=True)
                continue

            if is_fastq_file(file_path):
                try:
                    sha = compute_sha256(file_path)
                    gen_kwargs = {
                        "path": str(file_path),
                        "max_reads": None if max_fastq_reads <= 0 else max_fastq_reads,
                        "source_sha256": sha,
                    }
                    ds = Dataset.from_generator(fastq_reader, gen_kwargs=gen_kwargs)
                    rel_to_task = file_path.relative_to(task_path)
                    if flat_output:
                        rel_str = str(rel_to_task)
                        if rel_str.endswith(".gz"):
                            rel_str = rel_str[:-3]
                        rel_no_ext = str(Path(rel_str).with_suffix(""))
                        flat_name = f"{task_id}__fastq__{sanitize_for_dir(rel_no_ext)}"
                        out_subdir = out_dir / flat_name
                    else:
                        name = file_path.name[:-3] if file_path.name.endswith(".gz") else file_path.name
                        stem = str(Path(name).with_suffix(""))
                        subdir_parts = list(rel_to_task.parts[:-1])
                        safe_subdir = Path(*subdir_parts) if subdir_parts else Path()
                        out_subdir = fastq_out_base / safe_subdir / stem
                    save_dataset(ds, out_subdir, overwrite=overwrite)
                    converted_fastq += 1
                except Exception as e:
                    click.echo(f"[{task_id}] Skipping FASTQ {file_path.name}: {e}", err=True)
                    continue
            
            if is_fasta_file(file_path):
                try:
                    sha = compute_sha256(file_path)
                    gen_kwargs = {
                        "path": str(file_path),
                        "max_seqs": None if max_fasta_seqs <= 0 else max_fasta_seqs,
                        "source_sha256": sha,
                    }
                    ds = Dataset.from_generator(fasta_reader, gen_kwargs=gen_kwargs)
                    rel_to_task = file_path.relative_to(task_path)
                    if flat_output:
                        rel_str = str(rel_to_task)
                        if rel_str.endswith(".gz"):
                            rel_str = rel_str[:-3]
                        rel_no_ext = str(Path(rel_str).with_suffix(""))
                        flat_name = f"{task_id}__fasta__{sanitize_for_dir(rel_no_ext)}"
                        out_subdir = out_dir / flat_name
                    else:
                        name = file_path.name[:-3] if file_path.name.endswith(".gz") else file_path.name
                        stem = str(Path(name).with_suffix(""))
                        subdir_parts = list(rel_to_task.parts[:-1])
                        safe_subdir = Path(*subdir_parts) if subdir_parts else Path()
                        out_subdir = fasta_out_base / safe_subdir / stem
                    save_dataset(ds, out_subdir, overwrite=overwrite)
                    converted_fasta += 1
                except Exception as e:
                    click.echo(f"[{task_id}] Skipping FASTA {file_path.name}: {e}", err=True)
                    continue

            if is_vcf_file(file_path):
                try:
                    sha = compute_sha256(file_path)
                    ds = Dataset.from_generator(
                        vcf_reader,
                        gen_kwargs={"path": str(file_path), "source_sha256": sha},
                    )
                    rel_to_task = file_path.relative_to(task_path)
                    if flat_output:
                        rel_str = str(rel_to_task)
                        if rel_str.endswith(".gz"):
                            rel_str = rel_str[:-3]
                        rel_no_ext = str(Path(rel_str).with_suffix(""))
                        flat_name = f"{task_id}__vcf__{sanitize_for_dir(rel_no_ext)}"
                        out_subdir = out_dir / flat_name
                    else:
                        name = file_path.name[:-3] if file_path.name.endswith(".gz") else file_path.name
                        stem = str(Path(name).with_suffix(""))
                        subdir_parts = list(rel_to_task.parts[:-1])
                        safe_subdir = Path(*subdir_parts) if subdir_parts else Path()
                        out_subdir = vcf_out_base / safe_subdir / stem
                    save_dataset(ds, out_subdir, overwrite=overwrite)
                    converted_vcf += 1
                except Exception as e:
                    click.echo(f"[{task_id}] Skipping VCF {file_path.name}: {e}", err=True)
                    continue
        click.echo(f"[{task_id}] Converted tabular files: {converted_tables}")
        click.echo(f"[{task_id}] Converted FASTQ files: {converted_fastq}")
        click.echo(f"[{task_id}] Converted FASTA files: {converted_fasta}")
        click.echo(f"[{task_id}] Converted VCF files: {converted_vcf}")
    click.echo(f"Done. Outputs in: {out_dir}")

if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parents[1])
    main()

