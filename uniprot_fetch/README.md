# 🧬 UniProt Data Fetcher

A robust, GUI-based Python utility designed to query and download protein sequence data directly from the UniProt REST API. While pre-configured to fetch specific G-protein coupled receptors (GPCRs), the built-in query editor allows you to easily modify the search parameters for any custom UniProt data retrieval task.

---

## Key Features

* **Modern GUI:** A clean, high-DPI aware interface built with `tkinter` and `ttk`, featuring a scrollable query editor.
* **Real-Time Progress Tracking:** Features a visual progress bar displaying downloaded megabytes during large data streams.
* **Robust Networking:** Implements `requests.Session` with automatic retries for common HTTP errors (429, 500, 502, 503, 504) and stream processing to handle large datasets without consuming excessive RAM.
* **Atomic File Saving:** Writes data to a temporary file first and forces an OS-level hardware flush (`os.fsync`) before renaming. This ensures your data is never corrupted, even if the download is interrupted.
* **Automated Logging:** Automatically generates timestamped log files in a dedicated `logs/` directory for easy debugging and auditing.
* **Graceful Interruptions:** Safely handles `Ctrl+C` and window close events to prevent orphaned processes and clean up temporary files.

---

## Prerequisites

You will need **Python 3.10+** installed on your system. 

The script relies on the built-in `tkinter` module (which is included with most standard Python installations) and the third-party `requests` library.

## Installation

1. Clone this repository to your local machine.
2. Install the required external dependency using pip:

```bash
pip install requerements
```

---

## Usage

1. Launch the script from your terminal:

```bash
python fetch_GPCR_seq.py 
```

2. **Review or Edit the Query:** The application will open a window displaying the default UniProt query. You can edit this query string directly in the text box using standard UniProt search syntax.
3. **Run the Search:** Click the **Run Search Query** button.
4. **Choose Save Location:** A file dialog will prompt you to select a save location. By default, it suggests saving to a newly created `data/GPCR_Seq/` directory within the project root.
5. **Monitor Progress:** A dedicated progress window will appear, tracking the stream from the UniProt API to your local machine.

---

## The Default Query Explained

Out of the box, the tool is populated with a highly specific query targeting a subset of human GPCRs. 

**It explicitly filters for:**
* The Human Proteome (`UP000005640`)
* Reviewed entries (`reviewed:true`) with protein-level existence (`existence:1`)
* Receptors belonging to the G-protein coupled receptor or frizzled families.

**It explicitly excludes:**
* Olfactory, taste, opsin, vomeronasal, melanopsin, and mas-related proteins.
* Specific Trace amine-associated receptors (TAAR2-6, TAAR8, TAAR9).

---

## Project Structure

When you run the script, it will automatically organize your workspace by generating the following directories if they do not already exist:

* `logs/` - Contains timestamped `.log` files detailing network connections, successful downloads, and any encountered errors.
* `data/GPCR_Seq/` - The default recommended directory for saving your `.txt` or `.tsv` sequence results.
