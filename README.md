```text
    _____ _______ _____  _____ _   _ _______ 
  / ____|__   __|  __ \|_   _| \ | |__   __|
 | (___    | |  | |__) | | | |  \| |  | |   
  \___ \   | |  |  _  /  | | | . ` |  | |   
  ____) |  | |  | | \ \ _| |_| |\  |  | |   
 |_____/   |_|  |_|  \_\_____|_| \_|  |_|
```

# Strint
Strint: A single-cell preprocessing toolkit for full-length transcript identification, barcode demultiplexing, and UMI extraction.


- **Full-length transcript identification** 🧬  
- **Barcode identification** 🔖  
- **UMI identification** 🧩  

With Strint, researchers can streamline raw data processing before downstream single-cell RNA-seq analyses, ensuring reproducibility and accuracy.

Output Files

Strint generates the following output files:
  • BC_corrected.csv — Barcode-corrected information table
  • empty_bc_list.csv — Empty barcode list inferred from the barcode rank plot
  • knee_plot.png — Barcode rank plot for cell calling
  • matched_reads.fastq.gz — Reads assigned to valid barcodes
  • putative_bc.csv — Barcode information before correction
  • unmatched_reads.fastq.gz — Reads without valid barcodes
  • whitelist.csv — Final whitelist of selected barcodes

```
usage: main3.py [-h] --full_bc_whitelist FULL_BC_WHITELIST [--out_dir OUT_DIR] [--batch_size BATCH_SIZE] [--BC_fixed BC_FIXED] [--umi_fixed UMI_FIXED] [--putative_bc_out PUTATIVE_BC_OUT] [--out_whitelist_fn OUT_WHITELIST_FN]
                [--out_emptydrop_fn OUT_EMPTYDROP_FN] [--exp_cells EXP_CELLS] [--out_plot_fn OUT_PLOT_FN] [--DEFAULT_EMPTY_DROP_MIN_ED DEFAULT_EMPTY_DROP_MIN_ED] [--DEFAULT_EMPTY_DROP_NUM DEFAULT_EMPTY_DROP_NUM]
                [--fastq_out FASTQ_OUT] [--max_ed MAX_ED] [--minQ MINQ] [--threads THREADS]
                <input fastq filename/directory>

Strint: Single-cell preprocessing pipeline (full-length read identification, barcode splitting, UMI extraction)

positional arguments:
  <input fastq filename/directory>
                        Full-length sequencing fastq file (.fq or .fq.gz). If directory is given, all matching files are collected.

options:
  -h, --help            show this help message and exit
  --full_bc_whitelist FULL_BC_WHITELIST
                        Path to file containing all known barcodes (required).
  --out_dir OUT_DIR     Directory to save all output files (default: current directory).
  --batch_size BATCH_SIZE
                        Batch size for processing reads (default: 10000).
  --BC_fixed BC_FIXED   Fixed sequence between barcode parts (default: CCTTCC). Example: XXXXXXXXCCTTCCXXXXXXXX.
  --umi_fixed UMI_FIXED
                        Fixed sequence before UMI (default: CGATG). Example: CGATGXXXXXXXXXX.
  --putative_bc_out PUTATIVE_BC_OUT
                        Output filename for putative barcode table (default: putative_bc.csv).
  --out_whitelist_fn OUT_WHITELIST_FN
                        Output whitelist file containing selected barcodes (default: whitelist.csv).
  --out_emptydrop_fn OUT_EMPTYDROP_FN
                        Output file for empty droplet barcodes (default: empty_bc_list.csv).
  --exp_cells EXP_CELLS
                        Expected number of cells (default: 10000).
  --out_plot_fn OUT_PLOT_FN
                        Output filename for barcode rank plot (default: knee_plot.png).
  --DEFAULT_EMPTY_DROP_MIN_ED DEFAULT_EMPTY_DROP_MIN_ED
                        Minimum edit distance from empty drop BC to selected BC (default: 5).
  --DEFAULT_EMPTY_DROP_NUM DEFAULT_EMPTY_DROP_NUM
                        Maximum number of empty droplet barcodes to retain (default: 2000).
  --fastq_out FASTQ_OUT
                        Output fastq file containing reads assigned to whitelist barcodes (default: matched_reads.fastq.gz).
  --max_ed MAX_ED       Maximum edit distance threshold for barcode correction (default: 3).
  --minQ MINQ           Minimum quality score for barcode assignment (default: 10).
  --threads THREADS     Number of threads to used (default: 256).
```