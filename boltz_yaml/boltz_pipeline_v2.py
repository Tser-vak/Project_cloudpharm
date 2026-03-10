import os
import re
import logging
import argparse
from pathlib import Path
from multiprocessing import Pool
from dataclasses import dataclass
from typing import List, Generator, Tuple , Optional

# --- Constants & Logging ---
# Sharding prefix length: 2 characters (e.g., output/O0/...)
# This keeps ~4000-5000 files per folder instead of 1.4 million in one.
SHARD_PREFIX_LEN = 2

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger("BoltzPipeline")

# --- Models ---
@dataclass(frozen=True)
class Protein:
    id: str
    sequence: str
    msa_path: str

@dataclass(frozen=True)
class Ligand:
    id: str
    smiles: str

# --- Design Pattern: Repository ---
class ProteinRepository:
    """Handles parsing FASTA and linking MSA files."""
    def __init__(self, fasta_path: Path, msa_dir: Path):
        self.fasta_path = fasta_path
        self.msa_dir = msa_dir

    def get_all(self) -> List[Protein]:
        proteins = []
        if not self.fasta_path.exists():
            logger.error(f"FASTA file not found: {self.fasta_path}")
            return []

        logger.info(f"Parsing FASTA: {self.fasta_path}")
        with open(self.fasta_path, "r") as f:
            current_id = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        proteins.append(self._create_protein(current_id, "".join(current_seq)))
                    current_id = line[1:].split()[0] # Take first word as ID
                    current_seq = []
                else:
                    current_seq.append(line)
            
            if current_id: # Catch last one
                proteins.append(self._create_protein(current_id, "".join(current_seq)))
        
        # Filter out proteins without MSAs
        valid_proteins = [p for p in proteins if p is not None]
        logger.info(f"Loaded {len(valid_proteins)} proteins with valid MSAs.")
        return valid_proteins

    def _create_protein(self, p_id: str, seq: str) -> Optional[Protein]:
        msa_path = (self.msa_dir / f"{p_id}.a3m").resolve()
        if not msa_path.is_file():
            # In a production pipeline, we log this but don't stop the whole process
            logger.debug(f"Skipping {p_id}: MSA missing at {msa_path}")
            return None 
        # 2. The HPC path: formatted perfectly for the YAML file
        # This will output something like "./a3m_all/O94910.a3m"
        hpc_relative_path = f"./{self.msa_dir.name}/{p_id}.a3m"      
        return Protein(p_id, seq, hpc_relative_path)

class LigandRepository:
    """Handles parsing Ligand CSV/TSV data."""
    def __init__(self, file_path: Path):
        self.file_path = file_path

    def get_all(self) -> List[Ligand]:
        import pandas as pd
        all_ligands = []
        
        if not self.file_path.exists():
            logger.error(f"Ligand file not found: {self.file_path}")
            return []

        logger.info(f"Loading ligands from {self.file_path.name}")
        sep = "\t" if self.file_path.suffix == ".tsv" else ","
        df = pd.read_csv(self.file_path, sep=sep)
        
        for _, row in df.iterrows():
            all_ligands.append(Ligand(
                str(row['chembl_id']).strip(),
                str(row['smiles']).strip()
            ))
        
        logger.info(f"Loaded {len(all_ligands)} total ligands.")
        return all_ligands

# --- Design Pattern: Strategy (Task Worker) ---
def process_combination(task_args: Tuple[Protein, Ligand, Path]):
    """
    Worker function: Processes a single protein-ligand pair.
    Uses Sharding and Idempotency.
    """
    protein, ligand, output_root = task_args
    
    # 1. Generate safe name
    safe_p = re.sub(r"[^A-Za-z0-9._-]+", "_", protein.id)
    safe_l = re.sub(r"[^A-Za-z0-9._-]+", "_", ligand.id)
    job_name = f"{safe_p}__{safe_l}"
    
    # 2. Design Pattern: Sharding
    # We use the first 2 chars of protein ID as folder name if we input [:SHARD_PREFIX_LEN],without it it create for each gpcrs it own folder
    shard_dir = output_root / safe_p
    shard_dir.mkdir(exist_ok=True)
    
    out_path = shard_dir / f"{job_name}.yaml"
    
    # 3. Design Pattern: Idempotency (Resume capability)
    if out_path.exists():
        return False # Signal skipped

    # 4. Generate Content
    yaml_content = f"""version: 1
name: {job_name}

sequences:
  - protein:
      id: A
      sequence: {protein.sequence}
      msa: {protein.msa_path}

  - ligand:
      id: L
      smiles: "{ligand.smiles}"

properties:
  - affinity:
      binder: L
"""
    try:
        out_path.write_text(yaml_content, encoding="utf-8")
        return True # Signal success
    except Exception as e:
        # Don't let one failure crash the 1.4M pool
        return str(e)

# --- Orchestrator ---
def main():
    parser = argparse.ArgumentParser(description="Scalable Boltz-2 YAML Generator")
    parser.add_argument("--fasta", default="gpcrs/gpcrs.fasta")
    parser.add_argument("--msa_dir", default="a3m_all")
    parser.add_argument("--ligand_path", default="ligand_folder/Cleaned_phase_2.csv")
    parser.add_argument("--out_dir", default="output_phase_2")
    parser.add_argument("--workers", type=int, default=max(1 , (os.cpu_count() or 1) - 2), #check for cpu count and use that as default
                        help="Number of parallel processes (default: all cores)")
    args = parser.parse_args()

    # Paths
    root = Path.cwd()
    out_dir = root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Data Ingestion (Repository Pattern)
    prot_repo = ProteinRepository(root / args.fasta, root / args.msa_dir)
    lig_repo = LigandRepository(root / args.ligand_path)

    proteins = prot_repo.get_all()
    ligands = lig_repo.get_all()

    if not proteins or not ligands:
        logger.error("Insufficient data to proceed. Check your input paths.")
        return

    # 2. Producer Logic (Generating Task Arguments)
    # We use a generator to keep memory low (don't pre-calculate 1.4M tuples)
    def task_generator():
        for p in proteins:
            for l in ligands:
                yield (p, l, out_dir)

    total_tasks = len(proteins) * len(ligands)
    logger.info(f"Orchestrating {total_tasks} tasks across {args.workers} workers...")

    # 3. Execution (Consumer Pattern / Parallelism)
    count_new = 0
    count_skipped = 0
    count_error = 0

    # Optional: Use tqdm if you have it installed for a nice progress bar
    try:
        from tqdm import tqdm
        progress = tqdm(total=total_tasks, unit="job")
    except ImportError:
        progress = None

    with Pool(processes=args.workers) as pool:
        # imap_unordered is faster as it doesn't preserve order
        # chunksize=100 reduces IPC overhead for small tasks
        for result in pool.imap_unordered(process_combination, task_generator(), chunksize=100):
            if result is True:
                count_new += 1
            elif result is False:
                count_skipped += 1
            else:
                count_error += 1
                logger.error(f"Task failure: {result}")
            
            if progress:
                progress.update(1)

    if progress: progress.close()

    logger.info("--- Pipeline Summary ---")
    logger.info(f"Created: {count_new}")
    logger.info(f"Skipped (Existing): {count_skipped}")
    logger.info(f"Errors: {count_error}")
    logger.info(f"Results stored in: {out_dir}")

if __name__ == "__main__":
    main()
