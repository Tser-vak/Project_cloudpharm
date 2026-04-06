import sys
import os
import logging
import signal
import datetime
import tkinter as tk
from tkinter import messagebox, filedialog, scrolledtext, ttk

# --- High DPI Awareness (Fixes Blurriness on Windows) ---
try:
    from ctypes import windll
    windll.shcore.SetProcessDpiAwareness(1)
except Exception:
    pass

# --- Dynamic Path Resolution ---
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

# --- Instant Stop Handler ---
def handle_interrupt(sig, frame):
    """Forcefully exits the script."""
    print("\n[ABORT] User interrupted the process. Closing immediately...")
    os._exit(0)

signal.signal(signal.SIGINT, handle_interrupt)

# --- Enhanced Logging Configuration ---
def setup_logging():
    """Creates a 'logs' folder in the project root and a unique timestamped log file."""
    log_dir = os.path.join(PROJECT_ROOT, "logs")
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f"fetch_{timestamp}.log")
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, delay=True),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__), log_file

logger, current_log_path = setup_logging()

DEFAULT_QUERY = (
    "proteome:UP000005640 AND reviewed:true AND (existence:1) AND go:0004930 "
    "AND keyword:KW-0297 AND (family:G-protein coupled receptor OR family:frizzled) "
    "AND keyword:KW-0812 AND NOT (protein_name:olfactory "
    "OR protein_name:taste OR protein_name:opsin OR protein_name:vomeronasal "
    "OR protein_name:melanopsin OR protein_name:mas-related) "
    "AND NOT keyword:KW-0716 "
    "AND NOT (gene:TAAR2 OR gene:TAAR3 OR gene:TAAR4 OR gene:TAAR5 "
    "OR gene:TAAR6 OR gene:TAAR8 OR gene:TAAR9)"
)

class QueryDialog(tk.Tk):
    """Main Application Window with Modern Styling and High-DPI Support."""
    def __init__(self, title, initial_query):
        super().__init__()
        self.title(title)
        self.result = None
        
        self.configure(bg="#f5f5f7") 
        self.geometry("850x700")
        self.protocol("WM_DELETE_WINDOW", self.on_cancel)
        
        self.bind("<Control-c>", lambda e: handle_interrupt(None, None))
        
        self.style = ttk.Style(self)
        self.style.theme_use('clam')
        
        self.style.configure("TLabel", background="#f5f5f7", font=("Segoe UI", 11))
        self.style.configure("TButton", font=("Segoe UI", 10), padding=(15, 8))
        
        self.style.configure("Action.TButton", 
                            font=("Segoe UI", 10, "bold"), 
                            foreground="white", 
                            background="#007aff",
                            padding=(20, 10))
        self.style.map("Action.TButton", 
                      background=[('active', '#005bb5'), ('pressed', '#004494')])

        container = tk.Frame(self, bg="#f5f5f7", padx=40, pady=30)
        container.pack(expand=True, fill=tk.BOTH)

        # Added DNA emoji next to the title
        lbl = ttk.Label(container, text="🧬 UniProt Search Query", font=("Segoe UI", 18, "bold"))
        lbl.pack(pady=(0, 20), anchor="w")

        text_frame = tk.Frame(container, bg="#d1d1d6", padx=1, pady=1)
        text_frame.pack(expand=True, fill=tk.BOTH)
        
        self.text_area = scrolledtext.ScrolledText(
            text_frame, 
            wrap=tk.WORD, 
            width=80, 
            height=20, 
            font=("Consolas", 12),
            bg="white",
            fg="#1d1d1f",
            relief="flat",
            padx=15,
            pady=15,
            insertbackground="#007aff",
            undo=True,
            autoseparators=True
        )
        self.text_area.pack(expand=True, fill=tk.BOTH)
        self.text_area.insert(tk.END, initial_query)

        btn_frame = tk.Frame(container, bg="#f5f5f7")
        btn_frame.pack(pady=(30, 0), fill=tk.X)

        self.run_btn = ttk.Button(btn_frame, text="Run Search Query", style="Action.TButton", command=self.on_ok)
        self.run_btn.pack(side=tk.RIGHT)

    def on_ok(self):
        query_text = self.text_area.get("1.0", tk.END).strip()
        if not query_text:
            messagebox.showwarning("Empty Query", "Please enter a search query before running.")
            return
            
        self.result = query_text
        self.withdraw() 
        self.quit()     

    def on_cancel(self):
        self.destroy()
        os._exit(0)

class ProgressWindow:
    def __init__(self, parent):
        self.window = tk.Toplevel(parent)
        self.window.title("Downloading Data")
        self.window.geometry("550x220")
        self.window.configure(bg="#f5f5f7")
        self.window.resizable(False, False)
        self.window.attributes('-topmost', True)
        
        container = tk.Frame(self.window, bg="#f5f5f7", padx=30, pady=30)
        container.pack(expand=True, fill=tk.BOTH)

        self.label = ttk.Label(container, text="Connecting to UniProt API...", font=("Segoe UI", 12, "bold"))
        self.label.pack(pady=(0, 15), anchor="w")
        
        self.progress = ttk.Progressbar(container, orient=tk.HORIZONTAL, length=480, mode='determinate')
        self.progress.pack(pady=10)
        
        self.status_label = ttk.Label(container, text="Initializing stream...", font=("Segoe UI", 10), foreground="#6e6e73")
        self.status_label.pack(pady=(5, 0))
        
        self.window.protocol("WM_DELETE_WINDOW", lambda: handle_interrupt(None, None))
        self.window.update()

    def update(self, current, total):
        if total > 0:
            val = (current / total) * 100
            self.progress['value'] = val
            mb_total = total / (1024 * 1024)
            mb_current = current / (1024 * 1024)
            self.status_label.config(text=f"Downloaded {mb_current:.2f} MB of {mb_total:.2f} MB")
        else:
            self.progress['mode'] = 'indeterminate'
            self.progress.step(2)
            mb_current = current / (1024 * 1024)
            self.status_label.config(text=f"Downloaded {mb_current:.2f} MB (Total size unknown)")
        
        self.window.update()

    def close(self):
        if self.window.winfo_exists():
            self.window.destroy()

def setup_session():
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"]
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))
    return session

def fetch_prot(query, output_path, progress_ui):
    import requests
    
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "query": query,
        "fields": "accession,protein_name,gene_names,sequence",
        "format": "tsv"
    }
    
    timeout = (10, 60)
    session = setup_session()
    
    # Use a temporary file path for atomic writes
    temp_path = output_path + ".tmp"
    
    try:
        logger.info(f"Connecting to UniProt...")
        with session.get(url, params=params, stream=True, timeout=timeout) as response:
            if response.status_code != 200:
                logger.error(f"Error: UniProt returned status {response.status_code}")
                return False
            
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0
            
            # Atomic Write: Stream to temp file first
            with open(temp_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=512 * 1024):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        progress_ui.update(downloaded, total_size)
                
                # Production Grade: Force OS to flush buffer to physical hardware
                f.flush()
                os.fsync(f.fileno())
            
            # Atomic Swap: Replace existing file only if download finished successfully
            if os.path.exists(output_path):
                os.remove(output_path)
            os.rename(temp_path, output_path)

            logger.info(f"Successfully downloaded {downloaded} bytes to {output_path}")
            return True
            
    except Exception as e:
        logger.error(f"Network or File error: {e}")
        # Cleanup partial file on failure
        if os.path.exists(temp_path):
            try: os.remove(temp_path)
            except: pass
        return False
    finally:
        session.close()

def main():
    app = QueryDialog("UniProt Data Fetcher (Gemini 3 Enhanced)", DEFAULT_QUERY)
    
    while True:
        app.deiconify()
        app.mainloop() 
        
        query = app.result
        if not query:
            break

        default_dir = os.path.join(PROJECT_ROOT, "data", "GPCR_Seq")
        if not os.path.exists(default_dir):
            os.makedirs(default_dir)
            
        output_path = filedialog.asksaveasfilename(
            parent=app,
            title="Save Sequence Results",
            initialdir=default_dir,
            initialfile="GPCR_Sequence_results.txt",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("TSV files", "*.tsv")]
        )
        
        if output_path:
            progress_ui = ProgressWindow(app)
            success = fetch_prot(query, output_path, progress_ui)
            progress_ui.close()
            
            if success:
                messagebox.showinfo("Success", f"Data saved to:\n{output_path}")
            else:
                messagebox.showerror("Error", f"Fetch failed. Check the logs for details.")
        else:
            app.result = None

if __name__ == '__main__':
    main()
