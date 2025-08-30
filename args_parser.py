
import argparse
import os
import sys
from pathlib import Path
from utils import err_msg

def check_files_exist(file_list):
    if isinstance(file_list, str):
        file_list = [file_list]
    exit_code = 0
    for fn in file_list:
        if not os.path.exists(fn):
            exit_code = 1
            err_msg(f'Error: can not find {fn}', printit=True)
    if exit_code == 1:
        sys.exit(1)
    else:
        return True

# get file with a certian extensions
def get_files_by_suffix(search_dir, suffix, recursive=True):
    files = []
    if isinstance(suffix, str):
        suffix = [suffix]
    if recursive:
        for i in suffix:
            files.extend(Path(search_dir).rglob(i))
        return sorted(files)
    else:
        for i in suffix:
            files.extend(Path(search_dir).glob(i))
        return sorted(files)

def get_files_from_dir(fastq_dir):
    check_files_exist(fastq_dir)
    if os.path.isdir(fastq_dir):
        fastq_fns = get_files_by_suffix(
            search_dir=fastq_dir,
            suffix=['*.fastq', '*.fq', '*.fastq.gz', '*.fq.gz'],
            recursive=True
        )
        if not fastq_fns:
            err_msg(f"No FASTQ files found in directory: {fastq_dir}", printit=True)
            sys.exit(1)
        
    elif os.path.isfile(fastq_dir):
        fastq_fns = [fastq_dir]
    else:
        err_msg(f"File type of input {fastq_dir} is not supported.", printit=True)
        sys.exit(1)
    return list(map(str, fastq_fns))

def _to_outdir(out_dir, p):
    return p if os.path.isabs(p) else os.path.join(out_dir, os.path.basename(p))

def set_parser():
    parser = argparse.ArgumentParser(
        description="Strint: Single-cell preprocessing pipeline (full-length read identification, barcode splitting, UMI extraction)"
    )
     # 必选参数
    parser.add_argument(
        "fastq_fns", 
        metavar="<input fastq filename/directory>",
        help="Full-length sequencing fastq file (.fq or .fq.gz). If directory is given, all matching files are collected.",
        type=get_files_from_dir
    )

    

    parser.add_argument(
        "--full_bc_whitelist", 
        type=str, 
        required=True,
        help="Path to file containing all known barcodes (required)."
    )

    # 通用参数
    parser.add_argument(
        "--out_dir", 
        type=str, 
        default=".", 
        help="Directory to save all output files (default: current directory)."
    )

    parser.add_argument("--batch_size", type=int, default=10000,
        help="Batch size for processing reads (default: 10000).")
    parser.add_argument("--BC_fixed", type=str, default="CCTTCC",
        help="Fixed sequence between barcode parts (default: CCTTCC). Example: XXXXXXXXCCTTCCXXXXXXXX.")
    parser.add_argument("--umi_fixed", type=str, default="CGATG",
        help="Fixed sequence before UMI (default: CGATG). Example: CGATGXXXXXXXXXX.")
    parser.add_argument("--putative_bc_out", type=str, default="putative_bc.csv",
        help="Output filename for putative barcode table (default: putative_bc.csv).")
    parser.add_argument("--out_whitelist_fn", type=str, default="whitelist.csv",
        help="Output whitelist file containing selected barcodes (default: whitelist.csv).")
    parser.add_argument("--out_emptydrop_fn", type=str, default="empty_bc_list.csv",
        help="Output file for empty droplet barcodes (default: empty_bc_list.csv).")
    parser.add_argument("--exp_cells", type=int, default=10000,
        help="Expected number of cells (default: 10000).")
    parser.add_argument("--out_plot_fn", type=str, default="knee_plot.png",
        help="Output filename for barcode rank plot (default: knee_plot.png).")
    parser.add_argument("--DEFAULT_EMPTY_DROP_MIN_ED", type=int, default=5,
        help="Minimum edit distance from empty drop BC to selected BC (default: 5).")
    parser.add_argument("--DEFAULT_EMPTY_DROP_NUM", type=int, default=2000,
        help="Maximum number of empty droplet barcodes to retain (default: 2000).")
    parser.add_argument("--fastq_out", type=str, default="matched_reads.fastq.gz",
        help="Output fastq file containing reads assigned to whitelist barcodes (default: matched_reads.fastq.gz).")
    parser.add_argument("--max_ed", type=int, default=3,
        help="Maximum edit distance threshold for barcode correction (default: 3).")
    parser.add_argument("--minQ", type=int, default=10,
        help="Minimum quality score for barcode assignment (default: 10).")

    args = parser.parse_args()

    # 处理输出目录
    os.makedirs(args.out_dir, exist_ok=True)
    args.putative_bc_out  = _to_outdir(args.out_dir, args.putative_bc_out)
    args.out_whitelist_fn = _to_outdir(args.out_dir, args.out_whitelist_fn)
    args.out_emptydrop_fn = _to_outdir(args.out_dir, args.out_emptydrop_fn)
    args.out_plot_fn      = _to_outdir(args.out_dir, args.out_plot_fn)
    args.fastq_out        = _to_outdir(args.out_dir, args.fastq_out)

    print(f"[Strint] Collected {len(args.fastq_fns)} FASTQ files.")

    return args

# 自检入口
if __name__ == "__main__":
    args = set_parser()
    print("\n=== Parsed arguments ===")
    for k, v in vars(args).items():
        print(f"{k:20}: {v}")