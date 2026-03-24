# 🧬 Scalable Boltz-2 YAML Generator And Smile Cleaner

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.14-blue.svg)](#)

# - Boltz-2 -Yaml
A high-performance, multiprocessing pipeline designed to generate combinatorial YAML configuration files for Boltz-2. This tool takes a list of proteins (via FASTA and associated MSA files) and a list of ligands (via CSV/TSV) and orchestrates the creation of individual task files for downstream molecular modeling.

Designed to handle massive datasets (e.g., millions of combinations), it utilizes strict design patterns to ensure low memory footprint, process safety, and easy recoverability.

## ✨ Key Features

* **Massive Scalability:** Uses Python's `multiprocessing.Pool` and memory-efficient generators to orchestrate millions of tasks without overwhelming RAM.
* **Smart Sharding:** Automatically groups output YAML files into subdirectories based on the protein ID to prevent file system limits (e.g., inode exhaustion in single directories).
* **Idempotent Execution:** Safely interrupt and resume your pipeline. The script checks for existing outputs and automatically skips them, saving valuable compute time.
* **Resilient Architecture:** Implements the Repository and Strategy design patterns to separate data parsing from business logic. One failing combination will not crash the entire worker pool.
* **Real-time Progress:** Built-in integration with `tqdm` to monitor the pipeline's progress in real-time.

## 📂 Directory Structure

  ```bash
  Project Structure
  .
  ├── a3m_all/                # Source for Protein MSAs (.a3m)
  ├── gpcrs/                  # Input Protein data
  │   └── gpcrs.fasta         # Multi-FASTA file
  ├── Cleaned_data/           # Input Ligand data
  │   └── approved.csv        # CSV with 'chembl_id' and 'smiles'
  ├── logs/                   # Auto-generated runtime logs
  └── boltz_pipeline_v2.py    # Scalable YAML generator
  ``` 

  ### Advanced Configuration
You can customize paths and hardware utilization via command-line arguments:

| Argument | Default | Description |
| :--- | :--- | :--- |
| `--fasta` | `gpcrs/gpcrs.fasta` | Path to protein FASTA file. |
| `--msa_dir` | `a3m_all` | Directory where MSA (`.a3m`) files are stored. |
| `--ligand_path` | `.../approved.csv` | Path to ligand CSV/TSV. |
| `--out_dir` | `output_Approved` | Destination for generated YAMLs. |
| `--workers` | `CPU_COUNT - 2` | Number of parallel processes to spawn. |
 


# - Smiles_cleaner
This utility is designed to clean and standardize SMILES strings within chemical CSV datasets. It is optimized for large-scale data, utilizing multiprocessing to ensure it doesn't bottleneck your pipeline.

** Note: **  This tool is highly specialized for personal workflows, specifically tailored for data that was carefully pre-filtered using SQLite on the ChEMBL database.

##  Features 
* **Multi-Core Processing:** Utilizes `ProcessPoolExecutor` to distribute row-level tasks across your CPU cores, drastically reducing processing time.
* **Modular Cleaning Strategies:** Decouples cleaning rules from the execution engine. Currently includes:
  * **Salt Remover:** Isolates the largest fragment in a SMILES string (e.g., stripping `[Na+]` from `CCO.[Na+]`).
  * **Whitespace Cleaner:** Removes invisible trailing spaces and tabs that can break chemical fingerprinting.
* **Deterministic Deduplication:** Tracks original indices to ensure "keep first" deduplication is strictly based on the original file order.
* **Phase-Based Sharding:** Optionally splits the cleaned dataset into separate CSVs based on clinical phase (`max_phase`) for downstream HPC pipelines.
* **Comprehensive Logging:** Features a silent `tqdm` console interface while maintaining deep, auditable logs for modifications, errors, and duplicates.

## Usage

Run the script via the command line. At a minimum, you must provide an input CSV and an output destination.

### Basic Usage
Standardize SMILES and remove duplicates:

```bash
python "smiles_cleaner.py" --input "raw_molecules.csv" --output "cleaned_molecules.csv"
```

### Advanced Usage
Clean the data, force the script to use exactly 4 CPU cores, and shard the final output by clinical phase into a specific directory:

```bash
python "smiles_cleaner.py" --input "raw_molecules.csv" --output "cleaned_molecules.csv" --split_dir "phase_outputs/" --cores "4"
```
##Command Line Arguments 

| Argument | Type | Required | Default | Description |
| :--- | :--- | :--- | :--- | :--- |
| `--input` | String | **Yes** | None | Path to the raw input CSV file. |
| `--output` | String | **Yes** | None | Path where the final, cleaned CSV will be saved. |
| `--log` | String | No | `logs/smiles_cleaning.log` | Path for the general execution and audit log. |
| `--count_log` | String | No | `logs/data_count.log` | Path for the final dataset metrics and phase counts. |
| `--dupe_log` | String | No | `logs/duplicates_removed.csv` | Path to save the CSV containing removed duplicate rows. |
| `--split_dir` | String | No | None | Directory to save sharded phase-specific CSVs. |
| `--cores` | Integer | No | System CPU - 2 | Forces a specific number of CPU workers. |

## Output Example

``` bash 
name: PROTEIN_ID__LIGAND_ID
sequences:
  - protein:
      id: A
      sequence: MKWVTFISLLFLFSSAYS...
      msa: ./a3m_all/PROTEIN_ID.a3m
  - ligand:
      id: L
      smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O'
properties:
  - affinity:
      binder: L
```


## 🚀 Installation

Clone the repository and install the dependencies. The pipeline relies on core libraries such as `pandas` for robust ligand data parsing and `tqdm` for tracking massive job pools. 

```bash
git clone <your-repo-url>
cd <your-repo-directory>

# Install the exact required dependencies 
pip install -r requirements.txt
