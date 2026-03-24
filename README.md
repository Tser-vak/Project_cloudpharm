# 🧬 Scalable Boltz-2 YAML Generator And Smile Cleaner

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.14-blue.svg)](#)

# - Boltz-2 -Yaml
A high-performance, multiprocessing pipeline designed to generate combinatorial YAML configuration files for Boltz-2. This tool takes a list of proteins (via FASTA and associated MSA files) and a list of ligands (via CSV/TSV) and orchestrates the creation of individual task files for downstream molecular modeling.

Designed to handle massive datasets (e.g., millions of combinations), it utilizes strict design patterns to ensure low memory footprint, process safety, and easy recoverability.

# - Smiles_cleaner
This utility is designed to clean and standardize SMILES strings within chemical CSV datasets. It is optimized for large-scale data, utilizing multiprocessing to ensure it doesn't bottleneck your pipeline.

** Note: **  This tool is highly specialized for personal workflows, specifically tailored for data that was carefully pre-filtered using SQLite on the ChEMBL database.

## Features

* **Multi-Core Processing:** Utilizes `ProcessPoolExecutor` to distribute row-level tasks across your CPU cores, drastically reducing processing time.
* **Modular Cleaning Strategies:** Decouples cleaning rules from the execution engine. Currently includes:
  * **Salt Remover:** Isolates the largest fragment in a SMILES string (e.g., stripping `[Na+]` from `CCO.[Na+]`).
  * **Whitespace Cleaner:** Removes invisible trailing spaces and tabs that can break chemical fingerprinting.
* **Deterministic Deduplication:** Tracks original indices to ensure "keep first" deduplication is strictly based on the original file order.
* **Phase-Based Sharding:** Optionally splits the cleaned dataset into separate CSVs based on clinical phase (`max_phase`) for downstream HPC pipelines.
* **Comprehensive Logging:** Features a silent `tqdm` console interface while maintaining deep, auditable logs for modifications, errors, and duplicates.

## ✨ Key Features

* **Massive Scalability:** Uses Python's `multiprocessing.Pool` and memory-efficient generators to orchestrate millions of tasks without overwhelming RAM.
* **Smart Sharding:** Automatically groups output YAML files into subdirectories based on the protein ID to prevent file system limits (e.g., inode exhaustion in single directories).
* **Idempotent Execution:** Safely interrupt and resume your pipeline. The script checks for existing outputs and automatically skips them, saving valuable compute time.
* **Resilient Architecture:** Implements the Repository and Strategy design patterns to separate data parsing from business logic. One failing combination will not crash the entire worker pool.
* **Real-time Progress:** Built-in integration with `tqdm` to monitor the pipeline's progress in real-time.

## 🚀 Installation

Clone the repository and install the dependencies. The pipeline relies on core libraries such as `pandas` for robust ligand data parsing and `tqdm` for tracking massive job pools. 

```bash
git clone <your-repo-url>
cd <your-repo-directory>

# Install the exact required dependencies 
pip install -r requirements.txt
