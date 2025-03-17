from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union
import json
import os
import hashlib
import urllib.request
from functools import partial
import shutil
import subprocess
from urllib.parse import urlparse


class TaskCategory(Enum):
    """Categories of bioinformatics tasks in the benchmark."""
    EXPERIMENTAL_EVOLUTION = "experimental_evolution"
    CYSTIC_FIBROSIS = "cystic_fibrosis"
    RNA_SEQ = "rna_seq"
    METAGENOMICS = "metagenomics"
    COMPARATIVE_GENOMICS = "comparative_genomics"


@dataclass
class BioTaskDescription:
    """Description of a bioinformatics task in the benchmark."""
    task_id: str
    name: str
    description: str
    difficulty: str 
    num_bio_tools: int
    num_packages: int


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