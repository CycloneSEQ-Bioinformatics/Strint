from collections import namedtuple, defaultdict, Counter
import gzip
import pandas as pd
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
import os

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

def rfind_with_negative(s, sub):
    pos = s.rfind(sub)
    if pos == -1:
        return -1  # not found 
    return pos - len(s)

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

def get_bc_whitelist(raw_bc_count, full_bc_whitelist=None, exp_cells=None, out_plot_fn = None,empty_max_count = np.inf, DEFAULT_EMPTY_DROP_MIN_ED=None, DEFAULT_EMPTY_DROP_NUM=None):
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
    #print(len(raw_bc_count))`
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

class read_fastq:
    """This class is for mimic the Bio.SeqIO fastq record. The SeqIO is not directly used because it's slow.
    """
    def __init__(self, title, sequence, qscore, quality_map = False):
        self.id = title.split()[0].strip("@")
        self.seq = sequence
        self.qscore = qscore

def _read_and_bc_batch_generator_with_idx(fastq_fns, putative_bc_csv, batch_size):
    """Generator of barches of reads from list of fastq files with the idx of the first read
    in each batch

    Args:
        fastq_fns (list): fastq filenames
        batch_size (int, optional):  Defaults to 1000.
    """
    read_idx = 0
    putative_bc_f = open(putative_bc_csv, 'r')
    putative_bc_header = next(putative_bc_f)

    for fn in fastq_fns:
        if str(fn).endswith('.gz'):
            with gzip.open(fn, "rt") as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in fastq_parser(handle))

                batch_iter = batch_iterator(fastq, batch_size=batch_size)
                
                for batch in batch_iter:
                    batch_len = len(batch)
                    batch_bc_df = pd.read_csv(
                        StringIO(
                            putative_bc_header + \
                            ''.join([next(putative_bc_f) for x in range(batch_len)])
                        ))
                    yield batch, read_idx, batch_bc_df
                    read_idx += batch_len
        else:
            with open(fn) as handle:
                fastq =\
                    (read_fastq(title, sequence, qscore) for title, sequence, qscore in fastq_parser(handle))
                read_batch = batch_iterator(fastq, batch_size=batch_size)
                for batch in read_batch:
                    batch_len = len(batch)
                    batch_bc_df = pd.read_csv(
                        StringIO( putative_bc_header + \
                        ''.join([next(putative_bc_f) for x in range(batch_len)]))
                    )

                    yield batch, read_idx, batch_bc_df #[read_id_seq_qv,....,] 0,1000,2000... df
                    read_idx += batch_len
    putative_bc_f.close()

def _match_bc_row(row, whitelist, max_ed, minQ):
    
    strand = '+'
    
    if minQ and row.putative_bc_qscore < minQ: #没有起到作用？
        return ['', '', '']

    if not row.putative_bc or row.putative_bc in whitelist: #若bc为空或者bc在白名单中，直接返回
        return [row.putative_bc, row.putative_umi, strand]
    else: #若这个barcode不为空，但是不在白名单中，需要尝试纠错
        bc = row.putative_bc
    
    best_ed = max_ed #2
    bc_hit = ''

    #对于不在white list中的barcode 这个是whitelist.csv
    #以下代码暂无问题 202508281418
    for i in whitelist:
        ed, end_idx = sub_edit_distance(i, bc, best_ed)  #直接和white list中的做比较
        if ed < best_ed:
            best_ed = ed
            bc_hit = i #最佳匹配到的一个white list中的barcode
        elif ed == best_ed:
            if not bc_hit:
                bc_hit = i
            else: 
                bc_hit = 'ambiguous'
                best_ed -= 1
                if best_ed < 0:
                    return ['', row.putative_umi, strand]
    
    if bc_hit == 'ambiguous' or bc_hit == '':
        return ['', row.putative_umi, strand] #如果矫正失败，那么将bc_corrected一列变为空，putative_umi不变
    else:
        pass
            
    out_umi = row.putative_umi
    return [bc_hit, out_umi, strand]

def assign_read_batches(r_batch, whitelist, max_ed, gz, minQ=0, emit_unmatched_fastq=True):
    read_batch, start_df_idx, df = r_batch
    df = df.fillna('')
    whitelist = set(whitelist)
    out_buffer = ''
    unmatched_fastq_buffer = ''
    
    new_cols = []
    for row in df.itertuples():
        new_cols.append(_match_bc_row(row, whitelist, max_ed, minQ))
    
    df[['BC_corrected','putative_umi', 'strand']] = new_cols
    demul_read_count = sum(df.BC_corrected!='') #最终所有的barcode（包括矫正回来的）

    for r, bc in zip(read_batch, df.itertuples()):
        try:
            assert bc.read_id == r.id            
        except AssertionError:
            err_msg("Different order in putative bc file and input fastq!", printit = True)
            sys.exit()
        if not bc.BC_corrected or not bc.putative_umi: #BC_corrected这一列有没有  并且putative_umi这一列有没有 目前是完全匹配，所以可能比例较低
            if emit_unmatched_fastq:
                putative_bc = getattr(bc, "putative_bc", "")
                if putative_bc:  # 确保不是空字符串
                    # 统计 A 的比例
                    a_count = putative_bc.count("A")
                    a_ratio = a_count / len(putative_bc)
        
                    if a_ratio <= 0.5:  # 只有 A 含量 ≤ 50% 才输出
                        seq = r.seq
                        qscore = r.qscore
                        cb_val = bc.BC_corrected if bc.BC_corrected else "NA"
                        umi_val = bc.putative_umi if bc.putative_umi else "NA"
                        header = (f"@{cb_val}_{umi_val}#{bc.read_id}_{getattr(bc, 'strand', '+')}"
                          f"\tCB:Z:{cb_val}\tUB:Z:{umi_val}")
                        unmatched_fastq_buffer += header + '\n'
                        unmatched_fastq_buffer += str(seq) + '\n'
                        unmatched_fastq_buffer += '+\n'
                        unmatched_fastq_buffer += qscore + '\n'
            continue
        #加上polyA的判断机制吧
        if bc.polyA_starts: #若polyA_starts不为空 目前为空的原因是umi固定序列左边的read太少，有可能是umi序列不完整
            seq = r.seq[:int(bc.polyA_starts)]
            qscore = r.qscore[:int(bc.polyA_starts)]
        else:
            if bc.umi_fixed_locs:
                seq = r.seq[:int(bc.umi_fixed_locs) -10 ] #如果没有找polyT,则根据umi位置进行裁剪
                qscore = r.qscore[:int(bc.umi_fixed_locs) - 10]
            else:
                continue
            
        # write to fastq
        out_buffer += f"@{bc.BC_corrected}_{bc.putative_umi}#{bc.read_id}_{bc.strand}\tCB:Z:{bc.BC_corrected}\tUB:Z:{bc.putative_umi}\n"
        out_buffer += str(seq) + '\n' 
        out_buffer += '+\n'
        out_buffer += qscore + '\n'
    
    if gz: #如果要求输出是gz文件
        b_out_buffer = gzip.compress(out_buffer.encode('utf-8'))
    else:
        b_out_buffer = out_buffer.encode('utf-8')

    if emit_unmatched_fastq:
        b_unmatched_fastq = (gzip.compress(unmatched_fastq_buffer.encode('utf-8'))
                             if gz else unmatched_fastq_buffer.encode('utf-8'))
    else:
        b_unmatched_fastq = None

    return df, b_out_buffer, demul_read_count, len(read_batch), b_unmatched_fastq


def assign_read(fastq_fns=None, fastq_out=None, putative_bc_csv=None, 
                    whitelsit_csv=None, max_ed=None,n_process=None, batchsize=None, minQ=0):
    
    gz = fastq_out.endswith('.gz') #判断输出文件是否为gz文件
    out_dir = os.path.dirname(fastq_out)       # 输出目录
    unmatched_out = os.path.join(out_dir, "unmatched_reads.fastq.gz")
        
    r_batches = \
        _read_and_bc_batch_generator_with_idx(fastq_fns, putative_bc_csv, batchsize)
    
    whitelist = [] 
    with open(whitelsit_csv, 'r') as f:
        for line in f:
            whitelist.append(line.split('-')[0].strip())

    if n_process == 1:
        demul_count_tot = 0
        count_tot = 0
        with open(fastq_out, 'wb') as output_handle, open(unmatched_out, 'wb') as unmatched_handle:
            pbar = tqdm(unit="Reads", desc='Processed')
            for r_batch in r_batches:
                _, b_fast_str, demul_count, read_count, b_unmatched_fastq = assign_read_batches(r_batch, whitelist, max_ed,  gz, minQ=0, emit_unmatched_fastq=True)
                demul_count_tot += demul_count
                count_tot += read_count
                output_handle.write(b_fast_str)
                if b_unmatched_fastq:
                    unmatched_handle.write(b_unmatched_fastq)
                pbar.update(read_count) #
        green_msg(f"Reads assignment completed. Demultiplexed read saved in {fastq_out}!")
        green_msg(f"Unmatched reads saved in: {unmatched_out}!")
        
    else:
        rst_futures = multiprocessing_submit(assign_read_batches, 
                            r_batches, 
                            n_process=n_process,
                            schduler = "process",
                            pbar_func=lambda x: len(x[0]),
                            whitelist = whitelist,
                            max_ed = max_ed,
                            gz = gz)

        demul_count_tot = 0
        count_tot = 0
        df_list = []
        with open(fastq_out, 'wb') as output_handle, open(unmatched_out, 'wb') as unmatched_handle:
            for f in rst_futures:
                df, b_fast_str, demul_count, read_count, b_unmatched_fastq = f.result()
                demul_count_tot += demul_count
                count_tot += read_count
                output_handle.write(b_fast_str)

                if b_unmatched_fastq:
                    unmatched_handle.write(b_unmatched_fastq)
                # 保存df
                df_list.append(df)
        big_df = pd.concat(df_list, ignore_index=True)
    
        green_msg(f"Reads assignment completed. Demultiplexed read saved in {fastq_out}!")
        green_msg(f"Unmatched reads saved in: {unmatched_out}!")
    
    return demul_count_tot, count_tot,big_df


def err_msg(msg, printit = False):
    CRED = '\033[91m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def warning_msg(msg, printit = False):
    CRED = '\033[93m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def green_msg(msg, printit = False):
    CRED = '\033[92m'
    CEND = '\033[0m'
    if printit:
        print(CRED + msg + CEND)
    else:
        return CRED + msg + CEND

def multiprocessing_submit(func, iterator, n_process=mp.cpu_count()-1 ,
                           pbar=True, pbar_unit='Read',pbar_func=len, 
                           schduler = 'process', *arg, **kwargs):
    """multiple processing or threading, 

    Args:
        func: function to be run parallely
        iterator: input to the function in each process/thread
        n_process (int, optional): number of cores or threads. Defaults to mp.cpu_count()-1.
        pbar (bool, optional): Whether or not to output a progres bar. Defaults to True.
        pbar_unit (str, optional): Unit shown on the progress bar. Defaults to 'Read'.
        pbar_func (function, optional): Function to calculate the total length of the progress bar. Defaults to len.
        schduler (str, optional): 'process' or 'thread'. Defaults to 'process'.

    Yields:
        return type of the func: the yield the result in the order of submit
    """
    class fake_future:
        # a fake future class to be used in single processing
        def __init__(self, rst):
            self.rst = rst
        def result(self):
            return self.rst

    if schduler == 'process':
        # make sure the number of process is not larger than the number of cores
        n_process = min(n_process-1, mp.cpu_count()-1)
        if n_process > 1:
            executor = concurrent.futures.ProcessPoolExecutor(n_process)
    elif schduler == 'thread':
        if n_process > 1:
            executor = concurrent.futures.ThreadPoolExecutor(n_process)
    else:
        green_msg('Error in multiprocessing_submit: schduler should be either process or thread', printit=True)
        sys.exit(1)

    if pbar:
        _pbar = tqdm(unit=pbar_unit, desc='Processed')
        
    # run in single process/thread if n_process < 1
    if n_process <= 1:
        for it in iterator:
            yield fake_future(func(it, *arg, **kwargs))
            if pbar:
                _pbar.update(pbar_func(it))
        return

    # A dictionary which will contain the future object
    max_queue = n_process
    futures = {}
    n_job_in_queue = 0
    
    # make sure the result is yield in the order of submit.
    job_idx = 0
    job_completed = {}

    # submit the first batch of jobs
    while n_job_in_queue < max_queue:
        i = next(iterator, None)
        if i is None:
            break
        futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
        job_idx += 1
        n_job_in_queue += 1
        job_to_yield = 0
    # yield the result in the order of submit and submit new jobs
    while True:
        # will wait until as least one job finished
        # batch size as value, release the cpu as soon as one job finished
        job = next(as_completed(futures), None)

        # yield the completed job in the order of submit  
        if job is not None:
            job_completed[futures[job][1]] = job, futures[job][0]
            del futures[job]

        # 
        if job is None and i is None and len(job_completed)==0:
            break

        # check order
        while job_to_yield in job_completed.keys():
            # update pregress bar based on batch size
            if pbar:
                _pbar.update(job_completed[job_to_yield][1])
            yield job_completed[job_to_yield][0]
            del job_completed[job_to_yield]
            
            # submit new job
            i = next(iterator, None)
            if i is not None:
                futures[executor.submit(func, i, *arg, **kwargs)] = (pbar_func(i),job_idx)
                job_idx += 1
                
            job_to_yield += 1


