import pandas as pd
import numpy as np
import pybedtools as pbt
import itertools
from multiprocessing import  Pool
from iteration_utilities import deepflatten


def get_name(s_file):
    """Return sample name from file path."""
    return s_file.split('/')[-1].replace('.bed', "").replace('.xl', "").replace('threshold_crosslinks_',"")


def parse_bed6_to_df(p_file):
    """Parse BED6 file to pandas.DataFrame."""
    return pd.read_csv(
        p_file,
        names=['chrom', 'start', 'end', 'name', 'score', 'strand'],
        sep='\t',
        header=None,
        dtype={'chrom': str, 'start': int, 'end': int, 'name': str, 'score': float, 'strand': str}
    )


def find_all(string,substring):
    """
    Function: Returning all the index of substring in a string
    Arguments: String and the search string
    Return:Returning a list
    """
    length = len(substring)
    c = 0
    indexes = []
    while c < len(string): # TODO: test with -len(substring)
        if string[c : c + length] == substring:
            indexes.append(c)
        c += 1
    return indexes



def find_motif(m, sequences, df, df_input, window, shift):
    m_pos_offshift = [find_all(x, m) for x in sequences]
    m_pos = [[x + shift - window for x in y] for y in m_pos_offshift]
    if df_input is not None:
        try:
            m_prtxn_str = df_input.loc[m, 'prtxn']
        except KeyError:
            m_prtxn_str = df_input.loc[m.replace('T', 'U'), 'prtxn']
        try:
            m_prtxn = [int(x) for x in m_prtxn_str.split(',')]
        except AttributeError:
            try:
                m_prtxn = [df_input.loc[m, 'mtxn']]
            except KeyError:
                m_prtxn = []        
    m_on_pos = []
    if df_input is not None:
        for pos in m_pos:
            pos_temp = []
            for p in pos:
                if p in m_prtxn:
                    pos_temp.append(p)
            m_on_pos.append(pos_temp)
        return m_on_pos
    else:
        return m_pos


def get_binding_sites(kmer_group, sequences, df, df_input, window, shift, k_length):
    for m in kmer_group:
        on_pos_list = find_motif(m, sequences, df, df_input, window, shift)
        df.reset_index(inplace=True, drop=True)
        motif_series = pd.Series(on_pos_list)
        df[m + '_on_positions'] = motif_series
        df_nan = df[df.isnull().any(1)]
        if len(df_nan):
            print(f' WARNING: {len(df_nan)} NaN in data')

    df['group_positions'] = df[[x for x in df.columns if '_on_positions' in x]].values.tolist()
    df.group_positions = df.group_positions.apply(itertools.chain.from_iterable)
    df.group_positions = df.group_positions.apply(list)

    appended_data = []
    for _, row in df.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        name = row['name']
        score = row['score']
        strand = row['strand']
        for p in row['group_positions']:
            if row['strand'] == '+':
                teple = (
                    chrom,
                    start + int(p) - shift,
                    end + int(p) + k_length - shift - 1,
                    name,
                    score,
                    strand
                )
            elif row['strand'] == '-':
                teple = (
                    chrom,
                    start - int(p) - k_length + shift + 1,
                    end - int(p) + shift,
                    name,
                    score,
                    strand
                )
            appended_data.append(teple)
    return appended_data


def merge_peaks(peaks, n):
    return peaks.sort().merge(d=n, s=True, c=[4, 5, 6], o=['distinct', 'sum', 'distinct'])


def parallelize(func, list_of_strings, sequences, df, df_input, myint1, myint2, myint3, n_cores):
    split_seqs = np.array_split(sequences, n_cores)
    split_df = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    iterable_args = [(list_of_strings, split_seq, split_df[i], df_input, myint1, myint2, myint3) for i, split_seq in enumerate(split_seqs)]
    results = pool.starmap(func, iterable_args)
    pool.close()
    pool.join()
    return results


def run_bs(kmer_group_T, xn_file, tsv_file, genome, genome_fai, window, k_length, max_merge_dist=30, n_cores=4, 
    chunk_size=10000, output_dir=None, consensus=None):
    print(xn_file)
    print('kmer group input ', kmer_group_T)
    kmer_group = [k.replace('U', 'T') for k in kmer_group_T]
    if consensus:
        file_name = consensus
    else:
        file_name = '_'.join(kmer_group[:3])
    print(file_name, consensus)
    print('kmer group with T ', kmer_group)
    df_in = parse_bed6_to_df(xn_file)
    print(f'{len(df_in)} crosslinks in input file')
    df_fai = pd.read_csv(genome_fai, sep='\t', header=None,
                         dtype={0 : str, 1 : int, 2 : int, 3 : int, 4 : int})
    chrom_names = set(df_fai[0].values)
    df = df_in[df_in.chrom.isin(chrom_names)].reset_index(drop=True)
    df_filtered = df_in[~df_in.chrom.isin(chrom_names)]
    print(f'{len(df_filtered)} crosslinks filtered and written to xl_filtered.bed')
    if len(df_filtered):
        df_filtered.to_csv('xl_filtered.bed', index=None, header=None, sep='\t')
    if tsv_file:
        df_input = pd.read_csv(tsv_file, sep='\t').set_index('Unnamed: 0')
    else:
        df_input = None
    shift = int((k_length + 1) / 2)
    n = len(df) // chunk_size + 1
    sites_chunks = np.array_split(df, n)

    print(f'There are {len(sites_chunks)} chunks')
    combined_results = []
    for chunk in sites_chunks:
        input_sites = pbt.BedTool.from_dataframe(chunk[['chrom', 'start', 'end', 'name', 'score', 'strand']])
        flank = input_sites.slop(l=window, r=window, s=True, g=genome_fai)
        seq_tab = flank.sequence(s=True, fi=genome, tab=True)
        sequences = [line.split("\t")[1].strip() for line in open(seq_tab.seqfn)]
        assert len(sequences) == len(chunk)
        results = parallelize(get_binding_sites, kmer_group, sequences, chunk, df_input, window, shift, k_length, n_cores)
        combined_results.append(results)
    appended_data_final = list(deepflatten(combined_results, types=list))
    print('Appended data length: ', len(appended_data_final))
    if not len(appended_data_final):
        print('Exiting, no binding sites found using specified inputs')
        return
    df_out = pd.DataFrame(appended_data_final, columns=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    
    input_name = get_name(xn_file)
    df_out = df_out[df_out.start > 1]
    print(f'df {input_name}_{file_name} done')
    print(f'df {input_name}_{file_name} saved')
    bed_out = pbt.BedTool.from_dataframe(df_out)
    bed_out_30 = merge_peaks(bed_out, max_merge_dist)
    if output_dir:
        bed_out_30.saveas(f'{output_dir}/merged30_{input_name}_{file_name}.bed')
    else:
        bed_out_30.saveas(f'merged30_{input_name}_{file_name}.bed')
    bed_out_0 = merge_peaks(bed_out, 0)
    if output_dir:
        bed_out_0.saveas(f'{output_dir}/merged0_{input_name}_{file_name}.bed')
    else:
        bed_out_0.saveas(f'merged0_{input_name}_{file_name}.bed')



if __name__ == "__main__":

    import sys
    
    motif_group = sys.argv[1].split(',')
    xl_in = sys.argv[2]
    if len(sys.argv) == 11:
        prtxn_file = None
        fasta = sys.argv[3]
        fai = sys.argv[4]
        window = int(sys.argv[5])
        kmer_len = int(sys.argv[6])
        num_cores = int(sys.argv[7])
        chunk_size = int(sys.argv[8])
        output_dir = sys.argv[9]
        consensus = sys.argv[11]
    if len(sys.argv) == 12:
        prtxn_file = sys.argv[3]
        fasta = sys.argv[4]
        fai = sys.argv[5]
        window = int(sys.argv[6])
        kmer_len = int(sys.argv[7])
        num_cores = int(sys.argv[8])
        chunk_size = int(sys.argv[9])
        output_dir = sys.argv[10]
        consensus = sys.argv[11]
    
    
    run_bs(
        motif_group,
        xl_in,
        prtxn_file,
        fasta,
        fai,
        window,
        kmer_len,
        n_cores=num_cores,
        chunk_size=chunk_size,
        output_dir=output_dir,
        consensus=consensus
    )

    
