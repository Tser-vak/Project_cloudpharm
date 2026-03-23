"""
smiles_cleaner.py

A high-performance, professional-grade utility for cleaning SMILES strings in CSV datasets.
Optimized with Multiprocessing and a silent tqdm interface for senior-level efficiency.

ARCHITECTURAL OVERVIEW:
1. Strategy Pattern: Decouples cleaning rules (salts, whitespace) from the execution engine.
2. Producer-Consumer Pattern: Uses a ProcessPoolExecutor to distribute row-level tasks across CPU cores.
3. Idempotency & Order: Tracks original indices to ensure "keep first" deduplication is deterministic.
4. Data Integrity: Performs phase-based sharding to prepare data for downstream HPC pipelines.
"""

import os
import logging
import argparse
import pandas as pd
import concurrent.futures
from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Any, Tuple
from tqdm import tqdm
from pathlib import Path

# =============================================================================
# DESIGN PATTERN: STRATEGY (SMILES CLEANERS)
# We use an Abstract Base Class (ABC) to define a strict interface for cleaning.
# This allows us to add new cleaning rules (e.g., chirality normalization)
# without ever touching the core ParallelProcessor logic.
# =============================================================================

class SmilesStrategy(ABC):
    """Abstract base class ensuring all cleaners implement the 'process' method."""
    @abstractmethod
    def process(self, smiles: str) -> Tuple[str, bool, str]:
        """
        Processes a SMILES string.
        Returns: (cleaned_smiles, was_modified_flag, audit_note)
        """
        pass

class SaltRemover(SmilesStrategy):
    """
    Rule: Chemical datasets often contain salts/ions separated by dots (e.g., 'CCO.[Na+]').
    Strategy: Isolate the largest fragment, assuming it is the primary organic molecule.
    """
    def process(self, smiles: str) -> Tuple[str, bool, str]:
        # Guard clause: If no dot, it's already a single fragment.
        if not smiles or not isinstance(smiles, str) or "." not in smiles:
            return smiles, False, ""

        # Split by the dot separator and select the fragment with the maximum length.
        fragments = smiles.split(".")
        longest_fragment = max(fragments, key=len)

        # If the longest fragment is different from the original, we've successfully stripped salts.
        if longest_fragment != smiles:
            return longest_fragment, True, f"Removed salts ({smiles} -> {longest_fragment})"
        return smiles, False, ""

class WhitespaceCleaner(SmilesStrategy):
    """
    Rule: Invisible characters (tabs, trailing spaces) can break chemical fingerprinting.
    Strategy: Perform a standard string strip to ensure clean indexing.
    """
    def process(self, smiles: str) -> Tuple[str, bool, str]:
        if not isinstance(smiles, str):
            return smiles, False, ""

        cleaned = smiles.strip()
        if cleaned != smiles:
            return cleaned, True, f"Trimmed whitespace ('{smiles}' -> '{cleaned}')"
        return smiles, False, ""

# =============================================================================
# WORKER FUNCTION (STATELESS & PICKLABLE)
# This function must reside at the top level so the 'multiprocessing' module
# can serialize (pickle) it and send it to child processes.
# =============================================================================

def process_row_worker(row_tuple: Tuple[int, pd.Series], strategies: List[SmilesStrategy]) -> Dict[str, Any]:
    """
    Independent worker function that processes a single CSV row.
    By making this stateless, we avoid race conditions and complex locking mechanisms.
    """
    index, row = row_tuple
    # Primary identifier lookup: Fallback to Row Index if chembl_id is missing.
    row_id = row.get('chembl_id') or row.get('ligand_id') or f"Row_{index}"
    original_smiles = str(row.get('smiles', ''))

    # Early exit for empty or corrupted data to prevent worker crashes.
    if not original_smiles or original_smiles.lower() == 'nan':
        return {"index": index, "id": row_id, "status": "ERROR", "msg": "Empty SMILES", "data": None}

    current_smiles = original_smiles
    any_modified = False
    audit_trail = []

    try:
        # Pipeline execution: Run the SMILES through every active strategy.
        for strategy in strategies:
            result_smiles, modified, note = strategy.process(current_smiles)
            if modified:
                any_modified = True
                audit_trail.append(note)
                current_smiles = result_smiles

        # Construct the new data row while preserving all original metadata.
        new_row = row.to_dict()
        new_row['smiles'] = current_smiles

        return {
            "index": index, # Essential for maintaining original file order during post-processing.
            "id": row_id,
            "status": "MODIFIED" if any_modified else "CLEAN",
            "msg": " | ".join(audit_trail) if any_modified else "Already clean",
            "data": new_row
        }
    except Exception as e:
        # Fail-safe: Capture the error and report it back to the orchestrator instead of crashing.
        return {"index": index, "id": row_id, "status": "ERROR", "msg": str(e), "data": None}

# =============================================================================
# CORE ORCHESTRATOR: ParallelProcessor
# This class manages the lifecycle of the data from ingestion to sharding.
# =============================================================================

class ParallelProcessor:
    """Manages high-speed, multi-core SMILES cleaning and post-processing."""

    def __init__(self, logger: logging.Logger, strategies: List[SmilesStrategy]):
        self.logger = logger
        self.strategies = strategies
        # Tracking metrics for the executive summary report.
        self.stats = {"total": 0, "modified": 0, "clean": 0, "errors": 0, "duplicates": 0}

    def run(self, input_path: Path, output_path: Path, dupe_log_path: Path, count_log_path: Path, split_dir: Optional[Path] = None, max_workers: Optional[int] = None):
        """Main execution loop for cleaning, deduplication, and phase-based splitting."""

        # 1. DATA INGESTION
        if not input_path.exists():
            self.logger.error(f"Input file missing: {input_path}")
            return

        self.logger.info(f"Loading data: {input_path}")
        try:
            df = pd.read_csv(input_path)
        except Exception as e:
            self.logger.critical(f"CSV Read Failure: {e}")
            return

        self.stats["total"] = len(df)
        results_with_metadata = []

        # 2. MULTIPROCESSING EXECUTION
        # Determine CPU count: Leave 2 cores free for OS/System stability.
        workers = max_workers or (os.cpu_count() - 2 if os.cpu_count() > 1 else 1)
        self.logger.info(f"Starting Multi-Core Clean (Workers: {workers})")

        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            # Map the dataframe rows into the worker pool.
            futures = [executor.submit(process_row_worker, (idx, row), self.strategies) for idx, row in df.iterrows()]

            # Progress Tracking: Use tqdm for a clean, professional CLI experience.
            with tqdm(total=len(df), desc="Cleaning (Multi-Core)", unit="row") as pbar:
                for future in concurrent.futures.as_completed(futures):
                    try:
                        result = future.result()
                        if result["status"] == "ERROR":
                            self.stats["errors"] += 1
                            self.logger.error(f"FAIL [{result['id']}]: {result['msg']}")
                        else:
                            # Tally success metrics.
                            if result["status"] == "MODIFIED":
                                self.stats["modified"] += 1
                                self.logger.info(f"MODIFIED [{result['id']}]: {result['msg']}")
                            else:
                                self.stats["clean"] += 1
                                self.logger.info(f"CLEAN [{result['id']}]: {result['msg']}")
                            results_with_metadata.append(result)
                    except Exception as e:
                        self.logger.critical(f"FATAL WORKER ERROR: {e}")
                    pbar.update(1)

        if not results_with_metadata:
            self.logger.warning("No valid results were generated.")
            return

        # 3. DETERMINISTIC RE-ORDERING
        # Because workers return asynchronously, we must sort by the original 'index'
        # to ensure the "Keep First" deduplication is truly based on the original file order.
        results_with_metadata.sort(key=lambda x: x["index"])
        final_data = [r["data"] for r in results_with_metadata]
        output_df = pd.DataFrame(final_data)

        # 4. DEDUPLICATION LOGIC
        # Identify rows where the cleaned SMILES string has already appeared.
        self.logger.info("Performing deduplication on 'smiles' column...")
        is_duplicate = output_df.duplicated(subset=['smiles'], keep='first')
        duplicates_df = output_df[is_duplicate]
        self.stats["duplicates"] = len(duplicates_df)

        if not duplicates_df.empty:
            dupe_log_path.parent.mkdir(parents=True, exist_ok=True)
            # Filter for requested audit columns to keep the log file focused.
            log_cols = ['chembl_id', 'pref_name', 'max_phase', 'smiles', 'inchi_key']
            actual_log_cols = [c for c in log_cols if c in duplicates_df.columns]
            duplicates_df[actual_log_cols].to_csv(dupe_log_path, index=False)
            self.logger.info(f"Logged {len(duplicates_df)} duplicates to: {dupe_log_path}")

        # Drop the duplicates from the primary dataset.
        output_df = output_df[~is_duplicate]

        # 5. FINAL SORTING & PERSISTENCE
        # Sort by chembl_id to provide a predictable output for the user.
        if 'chembl_id' in output_df.columns:
            output_df = output_df.sort_values(by='chembl_id')

        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_df.to_csv(output_path, index=False)
        self.logger.info(f"Saved primary cleaned data to: {output_path}")

        # 6. PHASE-BASED SHARDING (SPLITTING)
        # Group molecules by their clinical phase to facilitate targeted pipeline runs.
        phase_counts = {}
        if split_dir and 'max_phase' in output_df.columns:
            split_dir.mkdir(parents=True, exist_ok=True)
            for phase, group in output_df.groupby('max_phase'):
                # Handle cases where phase might be stored as a float/string.
                try:
                    phase_val = int(float(phase))
                except (ValueError, TypeError):
                    phase_val = phase
                
                count = len(group)
                phase_path = split_dir / f"phase_{phase_val}.csv"
                group.to_csv(phase_path, index=False)
                phase_counts[phase_val] = count
                self.logger.info(f"Generated: {phase_path} ({count} rows)")

        # 7. FINAL METRICS LOGGING
        # Write the final count and phase-specific counts to the dedicated log file.
        count_log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(count_log_path, 'w', encoding='utf-8') as f:
            f.write(f"Total Cleaned (No duplicates/salts): {len(output_df)}\n")
            if phase_counts:
                f.write("--- Phase Counts ---\n")
                for phase_val in sorted(phase_counts.keys()):
                    f.write(f"Phase {phase_val}: {phase_counts[phase_val]}\n")
        
        self.logger.info(f"Logged detailed data counts to: {count_log_path}")
        self._print_summary()

    def _print_summary(self):
        """Final executive summary for terminal output and logs."""
        summary = (
            "\n" + "=" * 40 + "\n"
            "         FINAL PERFORMANCE REPORT\n" +
            "=" * 40 + "\n"
            f"  Total Molecules:    {self.stats['total']}\n"
            f"  Perfect (Clean):    {self.stats['clean']}\n"
            f"  Improved (Fixed):   {self.stats['modified']}\n"
            f"  Duplicates Removed: {self.stats['duplicates']}\n"
            f"  Failures (Error):   {self.stats['errors']}\n" +
            "=" * 40
        )
        self.logger.info(summary)

# =============================================================================
# INFRASTRUCTURE: LOGGING & CLI
# =============================================================================

def setup_logging(log_path: Path) -> logging.Logger:
    """Configures high-fidelity logging with a silent console interface."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("SmilesCleaner")
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        # Detailed File Logging: Records EVERY operation for auditability.
        fh = logging.FileHandler(log_path, mode='w', encoding='utf-8')
        fh.setFormatter(logging.Formatter('%(asctime)s - [%(levelname)s] - %(message)s'))
        logger.addHandler(fh)

        # Silent Console: Only errors appear in the console to avoid breaking the tqdm bar.
        ch = logging.StreamHandler()
        ch.setLevel(logging.ERROR)
        ch.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(ch)

    return logger

def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(description="Multi-Core SMILES Processing Engine")
    parser.add_argument("--input", type=str, required=True, help="Path to input CSV")
    parser.add_argument("--output", type=str, required=True, help="Path to save results")
    parser.add_argument("--log", type=str, default="logs/smiles_cleaning.log", help="General log file")
    parser.add_argument("--count_log", type=str, default="logs/data_count.log", help="Log file for the final data count")
    parser.add_argument("--dupe_log", type=str, default="logs/duplicates_removed.csv", help="CSV to log removed duplicates")
    parser.add_argument("--split_dir", type=str, default=None, help="Directory to save split phase CSVs")
    parser.add_argument("--cores", type=int, default=None, help="Force CPU worker count")

    args = parser.parse_args()

    # Convert strings to Path objects for robust cross-platform path handling.
    input_path = Path(args.input)
    output_path = Path(args.output)
    log_path = Path(args.log)
    dupe_log_path = Path(args.dupe_log)
    split_dir = Path(args.split_dir) if args.split_dir else None

    # Initialize logging and cleaning strategies.
    logger = setup_logging(log_path)
    strategies = [WhitespaceCleaner(), SaltRemover()]

    # Execute the pipeline.
    count_log_path = Path(args.count_log)
    processor = ParallelProcessor(logger, strategies)
    processor.run(input_path, output_path, dupe_log_path, count_log_path, split_dir=split_dir, max_workers=args.cores)

if __name__ == "__main__":
    # Required for Windows multiprocessing compatibility.
    main()
