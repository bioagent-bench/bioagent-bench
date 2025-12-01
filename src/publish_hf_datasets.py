import os
from pathlib import Path
from typing import List

import click
from datasets import load_from_disk, Dataset
from huggingface_hub import HfApi


def find_dataset_dirs(root: Path) -> List[Path]:
    if not root.exists():
        return []
    candidates: List[Path] = []
    for info_file in root.rglob("dataset_info.json"):
        ds_dir = info_file.parent
        if (ds_dir / "state.json").exists():
            candidates.append(ds_dir)
    return sorted({p.resolve() for p in candidates})


def get_published_dataset_name(rel: Path) -> str:
    rel_str = rel.as_posix()
    rel_parts = rel_str.split("__")
    task_name = rel_parts[0]
    dataset_name = rel_parts[-1]

    return f"{task_name}.{dataset_name}"


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--out-dir",
    type=click.Path(path_type=Path),
    default=Path("data") / "hf",
    show_default=True,
    help="Directory containing saved Hugging Face datasets (produced by the preprocessing script).",
)
@click.option(
    "--namespace",
    type=str,
    required=True,
    help="Hugging Face user or organization to publish under (e.g., 'your-username').",
)
def main(
    out_dir: Path,
    namespace: str,
) -> None:
    out_dir = out_dir.resolve()
    api = HfApi()

    ds_dirs = find_dataset_dirs(out_dir)
    if not ds_dirs:
        click.echo(f"No datasets found under: {out_dir}")
        raise SystemExit(0)

    for idx, ds_path in enumerate(ds_dirs, start=1):
        rel = ds_path.relative_to(out_dir)
        published_dataset_name = get_published_dataset_name(rel)
        repo_id = f"{namespace}/{published_dataset_name}"
        click.echo(f"[{idx}/{len(ds_dirs)}] Publishing {published_dataset_name}")

        api.create_repo(
            repo_id=repo_id, repo_type="dataset", private=False, exist_ok=True
        )
        ds = load_from_disk(str(ds_path))
        ds.push_to_hub(
            repo_id=repo_id, commit_message=f"Add {published_dataset_name} dataset"
        )


if __name__ == "__main__":
    main()
