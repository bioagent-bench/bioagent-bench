from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union
import os
import urllib.request
from functools import partial
from urllib.parse import urlparse
import yaml
from tqdm import tqdm


class BioBenchTask:
    """
    A scikit-learn style dataset object for bioinformatics benchmarks.
    
    This class provides a standardized interface for downloading, processing,
    and evaluating bioinformatics datasets for LLM agent benchmarking.
    """
    
    def __init__(
        self,
        task_id: str,
        data_home: Optional[str] = None,
        download_if_missing: bool = True,
        force_redownload: bool = False,
        verbose: bool = False,
    ):
        """Initialize the dataset.
        
        Args:
            task_id: Identifier for the task/dataset
            data_home: Directory where data is stored. If None, uses ~/bioagent_data
            download_if_missing: If True, download data if not present
            force_redownload: If True, force redownload even if data exists
            verbose: If True, print progress information
        """
        self.task_id = task_id
        self.verbose = verbose
        
        if data_home is None:
            self.data_home = os.path.expanduser("~/bioagent_data")
        else:
            self.data_home = data_home
            
        self.task_path = os.path.join(self.data_home, task_id)
        self.data_path = os.path.join(self.task_path, "data")
        self.results_path = os.path.join(self.task_path, "results")
        
        # Load dataset metadata
        self.metadata = self._load_metadata()
        
        if download_if_missing or force_redownload:
            self.download_data(force=force_redownload)
            self.download_results(force=force_redownload)

    def _load_metadata(self) -> Dict:
        """Load dataset metadata from YAML file."""
        yaml_path = os.path.join(os.path.dirname(__file__), "datasets.yaml")
        with open(yaml_path, 'r') as f:
            metadata = yaml.safe_load(f)
        return metadata['datasets'].get(self.task_id, {})

    def download_data(self, force: bool = False) -> None:
        """Download the dataset files.
        
        Args:
            force: If True, force redownload even if files exist
        """
        if not os.path.exists(self.data_path) or force:
            os.makedirs(self.data_path, exist_ok=True)
            data_url = self.metadata.get('data_url')
            if data_url:
                if self.verbose:
                    print(f"Downloading data for {self.task_id}...")
                urllib.request.urlretrieve(data_url, 
                                        os.path.join(self.data_path, "data.zip"))
                # TODO: Add unzip functionality
            else:
                print(f"No data URL found for {self.task_id}")

    def download_results(self, force: bool = False) -> None:
        """Download the results files.
        
        Args:
            force: If True, force redownload even if files exist
        """
        if not os.path.exists(self.results_path) or force:
            os.makedirs(self.results_path, exist_ok=True)
            results_url = self.metadata.get('results_url')
            if results_url:
                if self.verbose:
                    print(f"Downloading results for {self.task_id}...")
                urllib.request.urlretrieve(results_url,
                                         os.path.join(self.results_path, "results.zip"))
                # TODO: Add unzip functionality
            else:
                print(f"No results URL found for {self.task_id}")

    def display_meta(self) -> None:
        """Display metadata for this dataset."""
        print(f"\nDataset: {self.task_id}")
        print("-" * 40)
        for key, value in self.metadata.items():
            print(f"{key}: {value}")

    @staticmethod
    def list_datasets() -> List[str]:
        """List all available datasets."""
        yaml_path = os.path.join(os.path.dirname(__file__), "datasets.yaml")
        with open(yaml_path, 'r') as f:
            metadata = yaml.safe_load(f)
        datasets = list(metadata['datasets'].keys())
        print("\nAvailable datasets:")
        print("-" * 40)
        for dataset in datasets:
            print(f"- {dataset}")
        return datasets
    

if __name__ == "__main__":
    print(dataset.list_datasets())