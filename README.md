# 🧬 Scalable Boltz-2 YAML Generator

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](#)

A high-performance, multiprocessing pipeline designed to generate combinatorial YAML configuration files for Boltz-2. This tool takes a list of proteins (via FASTA and associated MSA files) and a list of ligands (via CSV/TSV) and orchestrates the creation of individual task files for downstream molecular modeling.

Designed to handle massive datasets (e.g., millions of combinations), it utilizes strict design patterns to ensure low memory footprint, process safety, and easy recoverability.

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
