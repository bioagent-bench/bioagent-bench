import json
import tarfile
from pathlib import Path
from urllib.parse import urlparse
import requests
from typing import List
import click


class BioAgentDataset:
    """Download and manage bioinformatics benchmark datasets."""
    
    def __init__(self, metadata_path: str = "src/task_metadata.json"):
        """Initialize with path to task metadata JSON file."""
        self.metadata_path = Path(metadata_path)
        self.tasks = self._load_metadata()
        self.base_dir = Path("tasks")
    
    def _load_metadata(self) -> List[dict]:
        """Load task metadata from JSON file."""
        with open(self.metadata_path, 'r') as f:
            return json.load(f)
    
    def _get_task(self, task_id: str) -> dict:
        """Get task metadata by task_id."""
        for task in self.tasks:
            if task["task_id"] == task_id:
                return task
        raise ValueError(f"Task '{task_id}' not found")
    
    def _download_file(self, url: str, output_path: Path) -> bool:
        """Download a single file from URL to output_path."""
        if not url or url.strip() == "":
            return False
            
        if output_path.exists():
            print(f"File already exists: {output_path}")
            return True
            
        print(f"Downloading {url} to {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            # Auto-extract tar.gz files
            if output_path.name.endswith('.tar.gz'):
                self._extract_tarfile(output_path)
                
            return True
            
        except Exception as e:
            print(f"Error downloading {url}: {e}")
            if output_path.exists():
                output_path.unlink()
            return False
    
    def _extract_tarfile(self, tar_path: Path) -> None:
        """Extract tar.gz file to the same directory."""
        extract_dir = tar_path.parent
        print(f"Extracting {tar_path} to {extract_dir}")
        
        try:
            with tarfile.open(tar_path, 'r:gz') as tar:
                tar.extractall(path=extract_dir)
            print(f"Successfully extracted {tar_path}")
        except Exception as e:
            print(f"Error extracting {tar_path}: {e}")
    
    def _download_urls(self, urls: List[str], target_dir: Path, category: str) -> bool:
        """Download multiple URLs to target directory."""
        if not urls:
            print(f"No {category} URLs to download")
            return True
            
        success = True
        for i, url in enumerate(urls):
            if not url or url.strip() == "":
                continue
                
            # Generate filename from URL or use index-based naming
            parsed_url = urlparse(url)
            filename = Path(parsed_url.path).name
            if not filename:
                filename = f"{category}_{i}.tar.gz"
                
            output_path = target_dir / filename
            if not self._download_file(url, output_path):
                success = False
                
        return success
    
    def download_data(self, task_id: str) -> bool:
        """Download data files for a specific task."""
        task = self._get_task(task_id)
        data_urls = task["download_urls"]["data"]
        target_dir = self.base_dir / task_id / "data"
        
        print(f"Downloading data for task: {task_id}")
        return self._download_urls(data_urls, target_dir, "data")
    
    def download_reference_data(self, task_id: str) -> bool:
        """Download reference data files for a specific task."""
        task = self._get_task(task_id)
        ref_urls = task["download_urls"]["reference_data"]
        target_dir = self.base_dir / task_id / "reference"
        
        print(f"Downloading reference data for task: {task_id}")
        return self._download_urls(ref_urls, target_dir, "reference")
    
    def download_results(self, task_id: str) -> bool:
        """Download result files for a specific task (for evaluation)."""
        task = self._get_task(task_id)
        result_urls = task["download_urls"]["results"]
        target_dir = self.base_dir / task_id / "results"
        
        print(f"Downloading results for task: {task_id}")
        return self._download_urls(result_urls, target_dir, "results")
    
    def download_task(self, task_id: str, include_reference: bool = False) -> bool:
        """Download all files for a task with optional reference data and results."""
        print(f"Downloading task: {task_id}")
        
        success = True
        
        if not self.download_data(task_id):
            success = False
            
        if include_reference:
            if not self.download_reference_data(task_id):
                success = False
                
        return success
    
    def download_all_tasks(self, include_reference: bool = False) -> bool:
        """Download all tasks with optional reference data and results."""
        success = True
        
        for task in self.tasks:
            task_id = task["task_id"]
            if not self.download_task(task_id, include_reference):
                success = False
        return success
    
    def list_tasks(self) -> List[str]:
        """List all available task IDs."""
        return [task["task_id"] for task in self.tasks]

    def download_all_results(self) -> bool:
        """Download results for all tasks."""
        success = True
        for task in self.tasks:
            task_id = task["task_id"]
            if not self.download_results(task_id):
                success = False
        return success


# -----------------------------
# CLI
# -----------------------------

@click.group()
@click.option(
    "--metadata",
    default="src/task_metadata.json",
    show_default=True,
    help="Path to task metadata JSON file.",
)
@click.pass_context
def cli(ctx: click.Context, metadata: str) -> None:
    """BioAgent dataset manager."""
    ctx.obj = BioAgentDataset(metadata_path=metadata)


@cli.command("list-tasks")
@click.pass_obj
def cli_list_tasks(dataset: BioAgentDataset) -> None:
    """List available task IDs."""
    for task_id in dataset.list_tasks():
        click.echo(task_id)


@cli.command("download")
@click.option("--task", "task_ids", multiple=True, help="Task ID(s) to download.")
@click.option("--all", "download_all", is_flag=True, help="Operate on all tasks.")
@click.option("--data/--no-data", default=True, show_default=True, help="Download data files.")
@click.option(
    "--reference/--no-reference",
    default=False,
    show_default=True,
    help="Include reference data.",
)
@click.option(
    "--results/--no-results",
    default=False,
    show_default=True,
    help="Include results files.",
)
@click.pass_obj
def cli_download(
    dataset: BioAgentDataset,
    task_ids: List[str],
    download_all: bool,
    data: bool,
    reference: bool,
    results: bool,
) -> None:
    """Download datasets for tasks (data, reference, results)."""
    if not download_all and not task_ids:
        click.echo("Specify at least one --task or use --all", err=True)
        raise SystemExit(1)

    # If only results are requested, optimize the flow
    if results and not data and not reference:
        if download_all:
            ok = dataset.download_all_results()
            raise SystemExit(0 if ok else 1)
        else:
            ok = True
            for tid in task_ids:
                ok = dataset.download_results(tid) and ok
            raise SystemExit(0 if ok else 1)

    target_ids = dataset.list_tasks() if download_all else list(task_ids)
    ok = True
    for tid in target_ids:
        if data:
            ok = dataset.download_task(tid, include_reference=reference) and ok
        elif reference:
            ok = dataset.download_reference_data(tid) and ok
        if results:
            ok = dataset.download_results(tid) and ok
    raise SystemExit(0 if ok else 1)


@cli.command("download-all-results")
@click.pass_obj
def cli_download_all_results(dataset: BioAgentDataset) -> None:
    """Download results for all tasks."""
    ok = dataset.download_all_results()
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    cli()