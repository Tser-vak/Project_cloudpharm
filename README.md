# OpenBable: High-Performance CIF Protein-Ligand Splitter

A Python utility designed to efficiently extract and separate proteins and ligands from `.cif` (Crystallographic Information File) formats. It utilizes `gemmi` for structural parsing and `obabel` (OpenBabel) via subprocesses to correct ligand chemistry and assign protonation states at a specified pH. 

It supports multiprocessing for high-throughput tasks and includes robust error logging.

---

## ⚠️ Important Environment Requirements

**This project strictly requires a Linux environment to run correctly** because it relies on Linux-based subprocess calls to OpenBabel.

* **IDE Users:** If you are using an IDE (like VS Code, PyCharm, etc.) on a Windows machine, your integrated terminal **must** be set to a Linux or WSL (Windows Subsystem for Linux) profile (e.g., Bash/Zsh). It will fail if run in PowerShell or Command Prompt.
* **Virtual Environment:** The Python virtual environment (`.venv`) **must** be created entirely within the Linux/WSL environment. Do not create it using a Windows Python installation, or the paths and binaries will be fundamentally incompatible.

---

## 🛠 Prerequisites

Before running the script, ensure you have the following installed in your Linux/WSL environment:

1.  **Python 3.10+**
2.  **OpenBabel:** Must be installed at the system level so the `obabel` command is available in your system path.
    ```bash
    # Example for Ubuntu/Debian via WSL
    sudo apt-get update
    sudo apt-get install openbabel
    ```

## 🚀 Installation

1. **Clone the repository:**
   Open your Linux/WSL terminal and run:
   ```bash
   git clone [https://github.com/Tser-vak/openbable.git](https://github.com/Tser-vak/openbable.git)
   cd ..
   pip install requirements.txt
