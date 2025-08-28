from collections import namedtuple
import gzip
import pandas as pd
from collections import defaultdict, Counter
from tqdm import tqdm
import numpy as np 
from matplotlib import pyplot as plt
import zipfile
import io
from fast_edit_distance import edit_distance
from io import StringIO

# a light class for a read in fastq file
read_tuple = namedtuple('read_tuple', ['id', 'seq', 'q_letter'])
def fastq_parser(file_handle):
    while True:
        id = next(file_handle, None)
        if id is None:
            break
        seq = next(file_handle)
        next(file_handle) # skip  '+'
        q_letter = next(file_handle)
        yield read_tuple(id[1:].split()[0], seq.strip(), q_letter.strip()) #每次yield一条read的信息

# split any iterator in to batches  
def batch_iterator(iterator, batch_size):
    """generateor of batches of items in a iterator with batch_size.
    """
    batch = []
    i=0
    for entry in iterator:
        i += 1
        batch.append(entry)
        
        if i == batch_size:
            yield batch
            batch = []
            i = 0
    if len(batch):  #保证批次处理的时候，最后一批不满足batch_size的那些数据，也可以被yield
        yield batch

def read_batch_generator(fastq_fns, batch_size):   #输出batch size read info
    """Generator of barches of reads from list of fastq files

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 100.
    """
    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq = fastq_parser(handle)
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch
        else:
            with open(fn, "r") as handle:
                fastq = fastq_parser(handle)
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    yield batch

def reverse_complement(seq):
    '''
    Args: <str>
        queried seq
    Returns: <str>
        reverse_complement seq
    '''
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                    'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    letters = \
        [comp[base] if base in comp.keys() else base for base in seq]
    return ''.join(letters)[::-1]

def polyA_trimming_idx(seq, seed="AAAA", window=10, min_A=7, min_tail_len=8):
    """
    从 read 末端往前检测 polyA，返回 polyA 起始的绝对坐标（0-based）。
    若未检测到则返回 None。
    """
    s = seq.upper()
    anchor = s.rfind(seed)  # 最右侧 seed 的起点（绝对坐标）
    if anchor == -1:
        return None

    polyA_start = anchor  # 先把 polyA 起点放在 seed 起点
    i = anchor - 1        # 从 seed 之前的碱基开始，向左延伸

    while i >= 0:
        if s[i] == 'A':
            polyA_start = i
            i -= 1
            continue
        left = max(0, i - window + 1)
        if s[left:i+1].count('A') >= min_A:
            polyA_start = i
            i -= 1
            continue
        break

    # 确保 polyA 足够长
    if len(s) - polyA_start < min_tail_len:
        return None
    return polyA_start
    
def polyA_trimming_idx_neg(seq, **kwargs):
    idx_abs = polyA_trimming_idx(seq, **kwargs)  # 用上面的绝对坐标函数
    if idx_abs is None:
        return None
    return idx_abs - len(seq)  # 负数：从末尾往前的偏移

def get_bc_whitelist(raw_bc_count, full_bc_whitelist=None, exp_cells=None, out_plot_fn=out_plot_fn,empty_max_count = np.inf):
    percentile_count_thres = default_count_threshold_calculation
    whole_whitelist = []
    if full_bc_whitelist.endswith('.zip'):
        with zipfile.ZipFile(full_bc_whitelist) as zf:
            # check if there is only 1 file
            assert len(zf.namelist()) == 1

            with io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding="utf-8") as f:
                for line in f:
                    whole_whitelist.append(reverse_complement(line.strip()))
    else:
        with open(full_bc_whitelist, 'r') as f:
            for line in f:
                whole_whitelist.append(reverse_complement(line.strip()))
    
    whole_whitelist = set(whole_whitelist)
    raw_bc_count = {k:v for k,v in raw_bc_count.items() if k in whole_whitelist}
    #print(len(raw_bc_count))
    t = percentile_count_thres(list(raw_bc_count.values()), exp_cells) #t是
    knee_plot(list(raw_bc_count.values()), t, out_plot_fn)
    cells_bc = {k:v for k,v in raw_bc_count.items() if v > t}
    
    ept_bc = []
    ept_bc_max_count = min(cells_bc.values()) #空bc最大的read支持数
    ept_bc_max_count = min(ept_bc_max_count, empty_max_count)
    #print(ept_bc_max_count)

    ept_bc_candidate = [k for k,v in raw_bc_count.items() if v < ept_bc_max_count]
    #print(len(ept_bc_candidate))
    for k in ept_bc_candidate:
        if min([edit_distance(k, x, max_ed = DEFAULT_EMPTY_DROP_MIN_ED) for x in cells_bc.keys()]) >= DEFAULT_EMPTY_DROP_MIN_ED:
            ept_bc.append(k)
            #print(len(ept_bc))
        # we don't need too much BC in this list
        if len(ept_bc) >  DEFAULT_EMPTY_DROP_NUM:
            break
    return cells_bc, ept_bc

fl_fq = ["/home/liyy/2.project/ScRNA-Seq/Glycine_Sc/6011_4/6011.full-length.fq.gz"] #3,176,234 read count 
batch_size = 1000
read_batchs = read_batch_generator(fl_fq, batch_size)

def rfind_with_negative(s, sub):
    pos = s.rfind(sub)
    if pos == -1:
        return -1  # not found 
    return pos - len(s)

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
BC_fixed = "GGAAGG" #reverse complement
umi_fixed = "CATCG"

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

    break


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

rst_df.to_csv("./putative_bc.csv",index=False)

def default_count_threshold_calculation(count_array, exp_cells):
    top_count = np.sort(count_array)[::-1][:exp_cells]
    return np.quantile(top_count, 0.95)/20

def knee_plot(counts, threshold=None, out_fn = 'knee_plot.png'):
    """
    Plot knee plot using the high-confidence putative BC counts

    Args:
        counts (list): high-confidence putative BC counts
        threshold (int, optional): a line to show the count threshold. Defaults to None.
    """
    counts = sorted(counts)[::-1]
    plt.figure(figsize=(8, 8))
    plt.title(f'Barcode rank plot (from high-quality putative BC)')
    plt.loglog(counts,marker = 'o', linestyle="", alpha = 1, markersize=6)
    plt.xlabel('Barcodes')
    plt.ylabel('Read counts')
    plt.axhline(y=threshold, color='r', linestyle='--', label = 'cell calling threshold')
    plt.legend()
    plt.savefig(out_fn)


dfs = pd.read_csv("./putative_bc.csv", chunksize=1_000_000)   #读取上一步输出的putative_bc.csv
minQ = 10
out_whitelist_fn = "./whitelist.csv"
out_emptydrop_fn  = "./emtpy_bc_list.csv"
full_bc_whitelist = "/home/liyy/1.data/C4_Cyclone/C4_blaze/C4_26bp.txt"
exp_cells = 10000
out_plot_fn = "knee_plot.png"
DEFAULT_EMPTY_DROP_MIN_ED = 5
DEFAULT_EMPTY_DROP_NUM = 2000 


raw_bc_count = Counter()
for df in tqdm(dfs, desc = 'Counting high-quality putative BC', unit='M reads'): #按块读取 CSV → 过滤掉低质量 putative_bc → 统计高质量 putative_bc 的出现次数 → 累加到总的 Counter。
        raw_bc_count += Counter(df[
        df.putative_bc_min_qs >=minQ].putative_bc.value_counts().to_dict()) #Counter({'GCGATAGTCACCTTCCGATGATTGCA': 166,'GACCTAATGGCCTTCCACGAGACTTG': 157,...}) 已经排序
#print(raw_bc_count.values())

#try:
bc_whitelist, ept_bc = get_bc_whitelist(raw_bc_count=raw_bc_count, full_bc_whitelist=full_bc_whitelist, exp_cells=exp_cells)

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


