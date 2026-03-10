# Boltz-2 Scalable Pipeline Generator

A high-performance, concurrent Python pipeline designed to generate massive datasets of Boltz-2 compatible YAML configuration files. 

By combining protein sequences (FASTA + MSA) and ligands (CSV/TSV), this tool can orchestrate millions of protein-ligand pairs efficiently. It utilizes advanced design patterns including task sharding, generator-based memory management, and idempotent execution to handle massive workloads on local machines or HPC clusters.

## ✨ Key Features

* **Massive Scalability:** Uses Python's `multiprocessing.Pool` and memory-efficient generators to process millions of combinations (e.g., $1.4\times 10^6$ tasks) without exhausting RAM.
* **Smart Sharding:** Automatically groups output files into subdirectories based on the Protein ID. This prevents filesystem crashes caused by dumping hundreds of thousands of files into a single directory.
* **Idempotency (Safe Resumes):** Checks for existing files before processing. If a massive job is interrupted, simply rerun the command—it will skip completed tasks and resume exactly where it left off.
* **Fault Tolerant:** An error in generating one specific YAML file will be logged but will not crash the entire worker pool.
* **Dynamic Worker Allocation:** Automatically detects your CPU cores and optimizes the worker pool (defaulting to `n-2` cores to keep your system responsive).

## 📋 Prerequisites

Install the exact dependencies required for this pipeline using the provided `requirements.txt` file:

```bash
pip install -r requirements.txt
