import os
import gemmi
import glob
import logging
import subprocess
import argparse
import sys
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from tqdm import tqdm


# ======= WAY TO RUN THE SCRIPT |||| python3 open_ba_split.py -i data/input -o data/output -w 16 -ph 7.0 ... ||||==================
# ====== Here we have A subprocess for linux proccess ||| Watch out the .venv is made for linux env , ubuntu  and to be able to eun the 
# script you need either to run it in a ,inux env or use the WSL terminal in the ID of you choosing 
# Writing to start the instance of env with source .venv/bin/activate  and then run the script with the command above ||| 
# You can also change the default input and output directories in the configuration section of the script

# =========================================================
# CONFIGURATION SECTION (Hardcode your defaults here)
# =========================================================
DEFAULT_INPUT_DIR = "data/no_2og"
LOG_BASE_DIR = "logs"

# =========================================================
# LOGGING SETUP
# =========================================================
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_dir = Path(LOG_BASE_DIR)
success_log_dir = log_dir / "success"
failure_log_dir = log_dir / "failure"

success_log_dir.mkdir(parents=True, exist_ok=True)
failure_log_dir.mkdir(parents=True, exist_ok=True)

success_log_path = success_log_dir / f"run_{timestamp}_success.log"
failure_log_path = failure_log_dir / f"run_{timestamp}_failure.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)

# Filter out logging to stream when using tqdm to avoid overlap
class TqdmLoggingHandler(logging.Handler):
    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)

# Update logging to use TqdmLoggingHandler for console
logger = logging.getLogger()
for handler in logger.handlers[:]:
    if isinstance(handler, logging.StreamHandler):
        logger.removeHandler(handler)
logger.addHandler(TqdmLoggingHandler())

def log_to_file(path, message):
    with open(path, "a") as f:
        f.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

def process_single_cif(cif_path, protein_dir, ligand_dir, success_log, failure_log, overwrite=False, ph="6.8", cofactor_logic=False, no_protein=False):
    cif_path = Path(cif_path)
    base_name = cif_path.stem
    
    prot_out = protein_dir / f"{base_name}_protein.pdb"
    lig_out = ligand_dir / f"{base_name}_ligand.sdf"
    
    # Skip if outputs already exist and overwrite is False
    if not overwrite and lig_out.exists() and (no_protein or prot_out.exists()):
        return f"Skipped: {base_name} (exists)"

    try:
        temp_lig_pdb = ligand_dir / f"{base_name}_temp_lig.pdb"

        doc = gemmi.cif.read_file(str(cif_path))
        st = gemmi.make_structure_from_block(doc.sole_block())
        
        # Extract Protein
        if not no_protein:
            st_protein = st.clone()
            st_protein.remove_ligands_and_waters()
            st_protein.write_pdb(str(prot_out))

        # Extract Ligand
        st_ligand = st.clone()
        st_ligand.remove_waters()
        
        # Determine target ligand based on directory
        folder_name = cif_path.parent.name
        target_res = None
        # !!!!!!Need to change the logic to let the user specify if it has cofactors or not |Now it works hardcoded|!!!!!!!!!
        if cofactor_logic.lower() == "y":
            target_res = "LIG2"
        elif cofactor_logic.lower() == "n":
            target_res = "LIG1"

        for model in st_ligand:
            for chain in model:
                for i in reversed(range(len(chain))):
                    res = chain[i]
                    if target_res:
                        # Keep only the target ligand (LIG1 or LIG2)
                        if res.name != target_res:
                            del chain[i]
                    else:
                        # Fallback to original logic if folder is unknown
                        res_info = gemmi.find_tabulated_residue(res.name)
                        if res_info.is_amino_acid() or res_info.is_nucleic_acid():
                            del chain[i]

        st_ligand.remove_empty_chains()
        
        if len(st_ligand) == 0 or len(st_ligand[0]) == 0:
            error_msg = "No ligand found in structure"
            log_to_file(failure_log, f"FAILED: {base_name} | {error_msg}")
            return f"No ligand: {base_name}"

        st_ligand.write_pdb(str(temp_lig_pdb))

        # Fix Chemistry with OpenBabel
        try:
            subprocess.run(
                ["obabel", str(temp_lig_pdb), "-O", str(lig_out), "-p", ph ],
                check=True, capture_output=True
            )
            log_to_file(success_log, f"SUCCESS: {base_name}")
        except subprocess.CalledProcessError as e:
            error_msg = f"OpenBabel error: {e.stderr.decode().strip()}"
            log_to_file(failure_log, f"FAILED: {base_name} | {error_msg}")
            return f"OB Error: {base_name}"
        finally:
            if temp_lig_pdb.exists():
                temp_lig_pdb.unlink()

        return f"Success: {base_name}"

    except Exception as e:
        log_to_file(failure_log, f"CRITICAL ERROR: {base_name} | {str(e)}")
        return f"Failed: {base_name}"

def main():
    # Set up Command Line Arguments
    parser = argparse.ArgumentParser(description="High-performance CIF protein-ligand splitter.")
    parser.add_argument("-i", "--input", default=DEFAULT_INPUT_DIR, help=f"Input directory containing .cif files (default: {DEFAULT_INPUT_DIR})")
    parser.add_argument("-l", "--logdir", default=LOG_BASE_DIR, help=f"Directory to store logs (default: {LOG_BASE_DIR})")
    parser.add_argument("-c" , "--cores" , type=int , default=max(1,(os.cpu_count()- 4)) , help="Number of CPU cores to use (default: all available)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files (default: False)")
    parser.add_argument("-ph", default="6.8", help="pH value for OpenBabel processing (default: 6.8)")
    parser.add_argument("-cf", default="n", choices=["y", "n"], help="Specify if the dataset contains cofactors (y/n, default: n)")
    parser.add_argument("--no-protein", action="store_true", help="Skip protein extraction (default: False)")
    args = parser.parse_args()

    input_dir = Path(args.input)
    
    # Logic to place output outside of the 'data' directory
    base_output_dir = input_dir.parent.parent if input_dir.parent.name == "data" else input_dir.parent
    
    protein_dir = base_output_dir / "proteins"
    ligand_dir = base_output_dir / "ligands"
    
    if not args.no_protein:
        protein_dir.mkdir(parents=True, exist_ok=True)
    ligand_dir.mkdir(parents=True, exist_ok=True)

    cif_files = list(input_dir.glob("*.cif"))
    
    if not cif_files:
        logging.error(f"No .cif files found in {input_dir}")
        return

    logging.info(f"Input: {input_dir}")
    if not args.no_protein:
        logging.info(f"Protein Output: {protein_dir}")
    logging.info(f"Ligand Output: {ligand_dir}")
    logging.info(f"Processing {len(cif_files)} files using {args.cores} cores...")

    executor = ProcessPoolExecutor(max_workers=args.cores)
    futures = {
        executor.submit(process_single_cif, f, protein_dir, ligand_dir, success_log_path, failure_log_path, args.overwrite, args.ph, args.cf, args.no_protein): f.name
        for f in cif_files
    }
    
    success_count = 0
    skipped_count = 0
    try:
        with tqdm(total=len(cif_files), desc="Processing CIFs", unit="file") as pbar:
            for future in as_completed(futures):
                result = future.result()
                if result.startswith("Success"):
                    success_count += 1
                elif result.startswith("Skipped"):
                    skipped_count += 1
                pbar.update(1)
                
    except KeyboardInterrupt:
        logging.warning("\nInterrupted by user (Ctrl+C). Cleaning up and exiting...")
        # Cancel all pending futures
        for future in futures:
            future.cancel()
        # Shutdown the executor immediately
        executor.shutdown(wait=False, cancel_futures=True)
        logging.info("Pending tasks canceled. Safe exit.")
        sys.exit(1)
    finally:
        executor.shutdown(wait=True)

    logging.info(f"Complete. Success: {success_count}/{len(cif_files)} (Skipped: {skipped_count})")
    logging.info(f"Logs stored in: {log_dir}/")

if __name__ == "__main__":
    main()
