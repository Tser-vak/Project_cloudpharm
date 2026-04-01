import os
import gemmi
import glob
import logging
import subprocess
import argparse
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path


# ======= WAY TO RUN THE SCRIPT |||| python3 open_ba_split.py -i data/input -o data/output -w 16 -p 7.0  ||||==================

# =========================================================
# CONFIGURATION SECTION (Hardcode your defaults here)
# =========================================================
DEFAULT_INPUT_DIR = "data/no_2og"
DEFAULT_OUTPUT_DIR = "data/processed"
LOG_BASE_DIR = "logs"

# =========================================================
# LOGGING SETUP
# =========================================================
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_dir = Path(LOG_BASE_DIR)
log_dir.mkdir(exist_ok=True)

success_log_path = log_dir / f"run_{timestamp}_success.log"
failure_log_path = log_dir / f"run_{timestamp}_failure.log"
main_log_path = log_dir / f"run_{timestamp}_main.log"

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(main_log_path),
        logging.StreamHandler()
    ]
)

def log_to_file(path, message):
    with open(path, "a") as f:
        f.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} - {message}\n")

def process_single_cif(cif_path, output_dir, success_log, failure_log):
    cif_path = Path(cif_path)
    base_name = cif_path.stem
    
    try:
        prot_out = output_dir / f"{base_name}_protein.pdb"
        lig_out = output_dir / f"{base_name}_ligand.sdf"
        temp_lig_pdb = output_dir / f"{base_name}_temp_lig.pdb"

        doc = gemmi.cif.read_file(str(cif_path))
        st = gemmi.make_structure_from_block(doc.sole_block())
        
        # Extract Protein
        st_protein = st.clone()
        st_protein.remove_ligands_and_waters()
        st_protein.write_pdb(str(prot_out))

        # Extract Ligand
        st_ligand = st.clone()
        st_ligand.remove_waters()
        
        for model in st_ligand:
            for chain in model:
                for i in reversed(range(len(chain))):
                    res = chain[i]
                    res_info = gemmi.find_tabulated_residue(res.name)
                    if res_info.is_amino_acid() or res_info.is_nucleic_acid():
                        del chain[i]

        st_ligand.remove_empty_chains()
        
        if len(st_ligand[0]) == 0:
            error_msg = "No ligand found in structure"
            log_to_file(failure_log, f"FAILED: {base_name} | {error_msg}")
            return f"No ligand: {base_name}"

        st_ligand.write_pdb(str(temp_lig_pdb))

        # Fix Chemistry with OpenBabel
        try:
            subprocess.run(
                ["obabel", str(temp_lig_pdb), "-O", str(lig_out), "-p", "7.4"],
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
    parser.add_argument("-o", "--output", default=DEFAULT_OUTPUT_DIR, help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})")
    parser.add_argument("-l", "--logdir", default=LOG_BASE_DIR, help=f"Directory to store logs (default: {LOG_BASE_DIR})")
    parser.add_argument("-c" , "--cores" , type=int , default=max(1,( os.cpu_count() - 2)) , help="Number of CPU cores to use (default: all available)")
    args = parser.parse_args()

    input_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    cif_files = list(input_dir.glob("*.cif"))
    
    if not cif_files:
        logging.error(f"No .cif files found in {input_dir}")
        return

    logging.info(f"Input: {input_dir} | Output: {output_dir}")
    logging.info(f"Processing {len(cif_files)} files using all available cores...")

    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(process_single_cif, f, output_dir, success_log_path, failure_log_path) 
            for f in cif_files
        ]
        
        success_count = 0
        for i, future in enumerate(futures):
            result = future.result()
            if result.startswith("Success"):
                success_count += 1
            
            if (i + 1) % max(1, len(cif_files) // 10) == 0:
                logging.info(f"Progress: {i+1}/{len(cif_files)} ({(i+1)/len(cif_files)*100:.1f}%)")

    logging.info(f"Complete. Success: {success_count}/{len(cif_files)}")
    logging.info(f"Logs stored in: {log_dir}/")

if __name__ == "__main__":
    main()

