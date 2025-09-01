from collections import namedtuple
import gzip
import pandas as pd
from collections import defaultdict, Counter
from tqdm import tqdm
import numpy as np 
from matplotlib import pyplot as plt
import zipfile
import io
from fast_edit_distance import edit_distance, sub_edit_distance
from io import StringIO
import multiprocessing as mp
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils import *
from args_parser import set_parser

def main():

    args = set_parser()
    fl_fq = args.fastq_fns
    out_dir = args.out_dir
    
    putative_bc_out = args.putative_bc_out
    putative_bc_csv = args.putative_bc_out
    out_whitelist_fn = args.out_whitelist_fn
    out_emptydrop_fn = args.out_emptydrop_fn
    out_plot_fn = args.out_plot_fn  
    fastq_out = args.fastq_out
    fastq_fns = args.fastq_fns
    full_bc_whitelist = args.full_bc_whitelist
    whitelsit_csv = args.full_bc_whitelist
    batch_size = args.batch_size
    batchsize = args.batch_size
    BC_fixed = args.BC_fixed
    umi_fixed = args.umi_fixed
    minQ = args.minQ
    exp_cells = args.exp_cells
    max_ed = args.max_ed
    
    DEFAULT_EMPTY_DROP_MIN_ED = args.DEFAULT_EMPTY_DROP_MIN_ED
    DEFAULT_EMPTY_DROP_NUM = args.DEFAULT_EMPTY_DROP_NUM
    n_process = args.threads
    
    ########################
    read_batchs = read_batch_generator(fl_fq, batch_size)
    
    read_ids = []
    putative_bcs = []
    putative_bc_min_qs = []
    bc_fixed_locs = [] #固定序列反向互补的5'位置 
    #raw_bc_pass = []
    umis = []
    umi_fixed_locs = []  #umi固定序列的5'端位置
    #trim_idxs = []
    post_umi_flankings = []
    polyA_starts = []
    BC_fixed = reverse_complement(BC_fixed) #reverse complement
    umi_fixed = reverse_complement(umi_fixed)
    
    for batch in read_batchs:
        for read_info in batch:
            part_id = read_info.id
            part_seq = read_info.seq[-30:]
            part_qv = read_info.q_letter[-30:]
            read_ids.append(part_id)
            putative_bc_min_q = None
            umi = None
            umi_fixed_loc= None
            post_umi_flanking = None
            polyA_start = None
            BC_fixed_loc = rfind_with_negative(part_seq, BC_fixed)
            #print(BC_fixed_loc)
            bc_fixed_locs.append(BC_fixed_loc)
            if BC_fixed_loc == -16:
                putative_bc = read_info.seq[-26:]
                putative_bcs.append(putative_bc)
                putative_bc_min_q = min([ord(x) for x in part_qv[-26:]]) -33
                putative_bc_min_qs.append(putative_bc_min_q)
                #locate umi
                find_umi_seq = read_info.seq[-36:-26] #barcode再往前10bp去找固定序列
                umi_fixed_loc_re = rfind_with_negative(find_umi_seq, umi_fixed) 
                if umi_fixed_loc_re != -1:
                    umi_fixed_loc = umi_fixed_loc_re - 26 #相对于read的位置
                    umi_fixed_locs.append(umi_fixed_loc)
                    umi = read_info.seq[umi_fixed_loc - 10 : umi_fixed_loc + 5]
                    umis.append(umi)
                    post_umi_flanking = read_info.seq[umi_fixed_loc - 10-5 : umi_fixed_loc - 10]
                    post_umi_flankings.append(post_umi_flanking)
                    #鉴定polyA的起始位置
                    seq_polyA = read_info.seq[umi_fixed_loc - 10 -100:umi_fixed_loc - 10]
                    last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
                    if last_polyA_idx: #可以检测到polyA
                        polyA_start =  last_polyA_idx - 10 + umi_fixed_loc #相对于整个read
                        polyA_starts.append(polyA_start)
                    else:
                        polyA_starts.append(polyA_start)
                else:
                    umis.append(umi)
                    umi_fixed_locs.append(umi_fixed_loc)
                    post_umi_flankings.append(post_umi_flanking)
                    polyA_starts.append(polyA_start)
                    
            elif BC_fixed_loc < -16:
                putative_bc = read_info.seq[BC_fixed_loc-10 : BC_fixed_loc+16]
                putative_bcs.append(putative_bc)
                putative_bc_min_q = min([ord(x) for x in read_info.q_letter[BC_fixed_loc-10:BC_fixed_loc+16]]) -33
                putative_bc_min_qs.append(putative_bc_min_q)
                #locate umi
                find_umi_seq = read_info.seq[BC_fixed_loc - 10 -10 :BC_fixed_loc - 10] #barcode再往前10bp去找固定序列
                umi_fixed_loc_re = rfind_with_negative(find_umi_seq, umi_fixed)
                
                    
                if umi_fixed_loc_re != -1:
                    umi_fixed_loc = umi_fixed_loc_re + BC_fixed_loc - 10 #相对于read的位置 修改了一个bug
                    umi_fixed_locs.append(umi_fixed_loc)
                    umi = read_info.seq[umi_fixed_loc - 10 : umi_fixed_loc + 5]
                    umis.append(umi)
                    post_umi_flanking = read_info.seq[umi_fixed_loc - 10-5 : umi_fixed_loc - 10]
                    post_umi_flankings.append(post_umi_flanking)
                    #if read_info.id == "250F302306011_11_158_8281_196090449_14000_1_14.58":
                    #    print(BC_fixed_loc)
                    #    print(umi_fixed_loc_re)
                    #    print(umi)
                    #    print(post_umi_flanking)
                    
                    #鉴定polyA的起始位置
                    seq_polyA = read_info.seq[umi_fixed_loc - 10 -100 :umi_fixed_loc - 10]
                    last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
                    if last_polyA_idx: #可以检测到
                        polyA_start =  last_polyA_idx - 10 + umi_fixed_loc #相对于整个read
                        polyA_starts.append(polyA_start)
                    else:
                        polyA_starts.append(polyA_start)
                    
                else:
                    umis.append(umi)
                    umi_fixed_locs.append(umi_fixed_loc)
                    post_umi_flankings.append(post_umi_flanking)
                    polyA_starts.append(polyA_start)
                
                
            elif BC_fixed_loc == -1: 
                #putative_bcs.append("No BC")
                putative_bcs.append(part_seq[-26:]) #直接输出后26bp
                putative_bc_min_qs.append(putative_bc_min_q)
                umi_fixed_locs.append(umi_fixed_loc)
                umis.append(umi)
                post_umi_flankings.append(post_umi_flanking)
                polyA_starts.append(polyA_start)
            else:#BC_fixed_loc>-16 右半段barcode不完全 #根据umi来判断barcode，即使是残缺的barcode
                #locate umi
                find_umi_seq = read_info.seq[BC_fixed_loc - 10 -10 :BC_fixed_loc - 10] #barcode再往前10bp去找固定序列
                umi_fixed_loc_re = rfind_with_negative(find_umi_seq, umi_fixed)
                if umi_fixed_loc_re != -1:
                    umi_fixed_loc = umi_fixed_loc_re + BC_fixed_loc - 10 #相对于read的位置 
                    umi_fixed_locs.append(umi_fixed_loc)
                    umi = read_info.seq[umi_fixed_loc - 10 : umi_fixed_loc + 5]
                    umis.append(umi)
                    post_umi_flanking = read_info.seq[umi_fixed_loc - 10-5 : umi_fixed_loc - 10]
                    post_umi_flankings.append(post_umi_flanking)
                    #根据umi来定位barcode 
                    putative_bc_start_loc = umi_fixed_loc + 5 
                    putative_bc = read_info.seq[putative_bc_start_loc:]
                    putative_bcs.append(putative_bc)
                    putative_bc_min_q = min([ord(x) for x in read_info.q_letter[putative_bc_start_loc:]]) -33
                    putative_bc_min_qs.append(putative_bc_min_q)
    
                    #if read_info.id == "250F302306011_13_7768_14189_229630721_8199_2_14.36":
                    #    print(BC_fixed_loc)
                    #    print(umi_fixed_loc_re)
                    #    print(umi)
                    #    print(post_umi_flanking)
                    #鉴定polyA的起始位置
                    seq_polyA = read_info.seq[umi_fixed_loc - 10 -100:umi_fixed_loc - 10]
                    last_polyA_idx = polyA_trimming_idx_neg(seq_polyA)
                    if last_polyA_idx: #可以检测到
                        polyA_start =  last_polyA_idx - 10 + umi_fixed_loc #相对于整个read
                        polyA_starts.append(polyA_start)
                    else:
                        polyA_starts.append(polyA_start)
                    
                else: #也没有找到umi
                    putative_bcs.append(part_seq[-26:])
                    putative_bc_min_qs.append(putative_bc_min_q)
                    umi_fixed_locs.append(umi_fixed_loc)
                    umis.append(umi)
                    post_umi_flankings.append(post_umi_flanking)
                    polyA_starts.append(polyA_start)
                    
            #break
    
        #break
    
    
    rst_df = pd.DataFrame(
            {'read_id': read_ids,
             'putative_bc': putative_bcs,
             'bc_fixed_locs':bc_fixed_locs,
             'putative_bc_min_qs': putative_bc_min_qs,
            'putative_umi': umis,
             'umi_fixed_locs':umi_fixed_locs,
             'post_umi_flankings':post_umi_flankings,
             'polyA_starts':polyA_starts
            }
            )
    
    rst_df.to_csv(putative_bc_out, index=False)
    
    dfs = pd.read_csv(putative_bc_out, chunksize=1_000_000)   #读取上一步输出的putative_bc.csv
    
    raw_bc_count = Counter()
    for df in tqdm(dfs, desc = 'Counting high-quality putative BC', unit='M reads'): #按块读取 CSV → 过滤掉低质量 putative_bc → 统计高质量 putative_bc 的出现次数 → 累加到总的 Counter。
            raw_bc_count += Counter(df[
            df.putative_bc_min_qs >=minQ].putative_bc.value_counts().to_dict()) #Counter({'GCGATAGTCACCTTCCGATGATTGCA': 166,'GACCTAATGGCCTTCCACGAGACTTG': 157,...}) 已经排序
    #print(raw_bc_count.values())
    
    #try:
    bc_whitelist, ept_bc = get_bc_whitelist(raw_bc_count=raw_bc_count, full_bc_whitelist=full_bc_whitelist, exp_cells=exp_cells, out_plot_fn = out_plot_fn, DEFAULT_EMPTY_DROP_MIN_ED=DEFAULT_EMPTY_DROP_MIN_ED, DEFAULT_EMPTY_DROP_NUM=DEFAULT_EMPTY_DROP_NUM)
    
    with open(out_whitelist_fn, 'w') as f:
        for k in bc_whitelist.keys():
            f.write(k+'\n')
    with open(out_emptydrop_fn, 'w') as f:
        for k in ept_bc:
            f.write(k+'\n')
    #except Exception as e: 
    #    print("Error: Failed to get whitelist. Please check the input files and settings.")
    #else:
    #    print(f'Whitelist saved as `{out_whitelist_fn}`!')
    
    
    demul_count_tot, count_tot, big_df = assign_read(fastq_fns=fastq_fns, fastq_out = fastq_out, putative_bc_csv=putative_bc_csv, 
                        whitelsit_csv=whitelsit_csv, max_ed=max_ed, n_process=n_process, batchsize=batchsize)
    
    #计算同时具有有效barcode和有效umi的read个数（也就是输出的fq的数量）
    mask = (big_df['BC_corrected'].notna()) & (big_df['BC_corrected'] != '') \
           & (big_df['putative_umi'].notna()) & (big_df['putative_umi'] != '')
    
    double_count = mask.sum()
    
    print(f"Number of full length reads: {count_tot}")
    print(f"Number of reads with valid barcode: {demul_count_tot}")
    print(f"Number of reads with valid barcode and valid umi: {double_count}")
    print(f"Proportion of reads with valid barcodes: {np.round(count_tot/demul_count_tot, 3)}")

if __name__ == "__main__":
    # Linux 通常默认 fork；若你换到 spawn/forkserver，请放在这里设置：
    # import multiprocessing as mp
    # try:
    #     mp.set_start_method("fork")
    # except RuntimeError:
    #     pass
    main()