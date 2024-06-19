#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 22:52:10 2021

@author: Peter Thomas, Laura McCoy's group, UCL 2019 - 2023.

Script is owned in entirety by UCL and is distributed freely according to fair use policy

Changes from V1 to V2:
    [X] - Pass a file to split into groups, rather than individual files
    [X] - Output gets written to a new folder in the original location (no argument to change)
    [X] - Do not write dendrograms to images, rather save square distance matrices (aim is to pass to R for plotting and stats)
    [X] - Duplicate CDRH3s removed, but repertoires are sampled with replacement, weighted by the CDRH3 frequency

"""

from pandas import read_csv, DataFrame, concat, merge
from numpy import repeat, array
import numpy as np
from collections import Counter
from datetime import date
from re import search
from os import mkdir
import argparse
from tqdm import tqdm

# %matplotlib inline


def transcribe(dna, dna_type = "cdna"):
    
    """
    
    Transcribes input DNA into RNA. Accepts dna or cdna arguments, with cdna as default.
    
    Arguments required:
        dna: A DNA string produced by read_fasta() with/without complimentary_dna()
        dna_type: cdna is defaulted, however function is compatible with dna if changed to "dna"
        
    Returns:
        rna: A string of RNA nucleotides
    
    """
    
    ### import necessary modules
    import re
    
    ### error handling, cannot accept dna types that are not dna or cdna. exits function
    if dna_type not in ["cdna", "dna"]:
        
        raise ValueError("unsupported dna type entered. accepted values are cdna or dna")
    
    ### value handling for incorrect data entry
    elif type(dna) is not str:
        
        raise ValueError("dna input should be of string type")
    
    else:
        
        ### if the dna_type is "cdna" (default)
        if dna_type == "cdna":
            
            ### define the comnplementarity dictionary
            dna2rna = {"A":"U",
                       "T":"A",
                       "G":"C",
                       "C":"G",
                       ".":"."}
            
            ### empty list for rna output
            rna = []
            
            ### loop through dna characters
            for nuc in dna:
                
                ### add rna values from dictionary to list
                rna.append(dna2rna[nuc])
            
            ### join into a single string and output
            return "".join(rna)
        
        else:
            
            ### if using direct from read_fasta, swap Ts with Us and return
            rna = re.sub("T", "U", dna)
            
    return rna
        
def translate(mrna, reading_frame = 1):
    
    """
    
    Translates an mRNA sequence in the specified reading frame, with frame 1 as the default.
    
    Arguments required:
        mrna: an mrna string produced by transcribe() function above
        reading_frame: a numerical value (1:3)
    
    Returns:
        A string of amino acids
    
    """
    
    ### value handling for incorrect data entry
    if type(mrna) is not str:
        
        raise ValueError("mrna input should be of string type")
        
    if type(reading_frame) is not int:
        
        raise ValueError("reading_frame input should be of integer type")
    
    ### create dictionary of codons
    codon2aa = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "CUU": "L",
                 "CUC": "L", "CUA": "L", "CUG": "L", "AUU": "I", "AUC": "I",
                 "AUA": "I", "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
                 "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S",
                 "AGC": "S", "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
                 "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCU": "A",
                 "GCC": "A", "GCA": "A", "GCG": "A", "UAU": "Y", "UAC": "Y",
                 "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAU": "N",
                 "AAC": "N", "AAA": "K", "AAG": "K", "GAU": "D", "GAC": "D",
                 "GAA": "E", "GAG": "E", "UGU": "C", "UGC": "C", "UGG": "W",
                 "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R",
                 "AGG": "R", "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
                 "AUG": "M", "UAA": "_", "UAG": "_", "UGA": "_"}
    
    ### use statements below to select start point of the mrna, and to exit if a non-reading frame value is supplied
    if reading_frame == 1:
        
        start_pos = 0
        
    elif reading_frame == 2:
        
        start_pos = 1
        
    elif reading_frame == 3:
        
        start_pos = 2
        
    else:
        
        raise ValueError("cannot translate outside frames 1, 2 or 3")
    
    ### initialise empty list for the protein output
    protein = []
    
    ### iterate in chunks of 3 over the mrna sequence input
    for codon in range(start_pos, len(mrna), 3):
        
        ### select the 3 residue section
        codon = mrna[codon:codon+3]
        
        ### ignore incomplete codons
        if len(codon) < 3:
            
            continue
        
        ### search the dictionary and add it to the protein list
        if codon not in codon2aa:
            
            protein.append('X')
            
        else:
            
            protein.append(codon2aa[codon])

    ### join the list and return
    protein = "".join(protein)
        
    return protein 

def findKmer(string, kmer):
    
    """
    Accessory function to kmerDecomp.
    Breaks the given string up into subsequences of the length given in kmer
    argument.
    
    Args:
        string- A character string of any length
        kmer- A numeric value to define the size of the kmer
        
    Returns:
        mer- A list of kmers from the input string
    
    """
    
    mer = []
    
    for position in range(0, len(string)):
        
        int_ = string[position:(position + kmer)]
        
        if len(int_) < kmer:
            
            int_ = int_ + '.' * (kmer - len(int_))
        
        mer.append(int_)
    
    return mer

def kmerDecomp(cdr3_list, kmer_size = 5):
    
    """
    Function to create a list of kmers per sequence.
    Contains findKmer function above.
    
    Function first finds the longest sequence in the list and then pads the
    remaining sequences with '.' to get a consistent length
    
    Args:
        cdr3_list- A list containing the strings to decompose into kmers
        kmer_size- Numeric value to define the size of the kmer
    
    Returns:
        kmer_list- A list the same length as cdr3_list, each element containing
        the relevant sequence in kmer form.
    
    """
    
    cdr3_lens = []
    
    for i in cdr3_list:
        
        cdr3_lens.append(len(i))
        
    longest_cdr3 = max(cdr3_lens)
    
    kmer_list = []
    
    for index, value in enumerate(cdr3_lens):
    
        pad_length = longest_cdr3 - value
        
        cdr3_list[index] = cdr3_list[index] + ''.join(list(repeat('.', pad_length)))
        
        kmer_list.append(findKmer(string = cdr3_list[index],
                                  kmer = kmer_size))
    
    return kmer_list

def countsPerPosition(kmer_list):
    
    """
    Function to count the kmer occurrences per position (largely deprecated,
    retained in case it's useful later on)
    
    Args:
        kmer_list- A list of kmers produced using kmerDecomp()
    
    Returns:
        kmer_data_frame- a pandas data frame containing the kmers per position
        in each sequence, and their frequency across the 'library'
    
    """
    
    n_kmer = len(kmer_list[0])
    
    kmers_by_position = []
    
    for position in range(0, n_kmer):
    
        kmer_position = []
        
        for index in range(0, len(kmer_list)):
            
            kmer_position.append(kmer_list[index][position])
    
        kmers_by_position.append(DataFrame.from_dict(Counter(kmer_position), orient='index').reset_index())
        kmers_by_position[position]['POSITION'] = position+1
        kmers_by_position[position] = kmers_by_position[position].rename(columns = {'index': 'KMER', 0: 'COUNT'})
        
        kmers_by_position[position]['FREQUENCY'] = kmers_by_position[position]['COUNT'] / sum(kmers_by_position[position]['COUNT'])
        

    kmer_data_frame = concat(kmers_by_position)
    
    return(kmer_data_frame)

def euclidean(l1, l2):
    
    """
    
    Function to calculate the euclidean distance between 2 kmer kidera factor scores.
    
    Args:
        l1- list of kidera factors for sequence 1
        l2- list of kidera factors for sequence 2
    Returns:
        The euclidean distance
    
    """
    
    out_ = [(i - j)**2 for i, j in zip(l1, l2)]
    return sum(out_)

def positionMatch(a, b, size):
    
    """
    Accessory function to getPositionalSimilarity. Scores kmers provided for
    identity per position. Matches get 1 point and mismatches get 0. If, because
    of sequence padding, two '.'s are being compared, the loop breaks to prevent comparison.
    
    Args:
        a- A list of kmers for an individual sequence
        b- A list of kmers for the sequence a is to be compared to
        size- The size of the kmers
        
    Returns:
        The similarity score based on kmer identity per position
    
    """
    
    score = 0
    
    valid_kmers = 0
    
    for position in range(0, len(a)):
        
        if a[position] == '.'*size and b[position] == '.'*size:
            break
        
        elif a[position] == b[position]:
            score += 1
        valid_kmers += 1
    if valid_kmers == 0:
        valid_kmers = 1
    return 1 - (score / valid_kmers)

def kideraDistance(a, b, size, kf_df):
    
    """
    
    Calculates the euclidean distances between the provided kmers, based on
    kidera factors (biochemical descriptors of amino acids) and summarises. 
    
    Currently in development.
    
    Args:
        a- A list of kmers for an individual sequence
        b- A list of kmers for the sequence a is to be compared to
        
    Returns:
        the sum of the euclidean distance measurements 
    
    """
    
    dist = []
    
    for position in range(0, len(a)):
    
        a_ = a[position]
        b_ = b[position]
        
        if a_ == '.'*size and b_ == '.'*size:
            break
        
        a_ = (kf_df.reindex(list(a_)).sum() / len(a_)).to_numpy()
        b_ = (kf_df.reindex(list(b_)).sum() / len(b_)).to_numpy()
        
        dist.append(euclidean(a_, b_))
    
    return sum(dist) / len(dist)

def kideraDistanceFast(kmer_list, kf_df, size):

    embeddings = []
    stop_at = []
    kmers = kmer_list
    
    for index_outer, seq in enumerate(kmers):
        int_ = []
        
        c = [len(seq)]
        
        for index_inner, k in enumerate(seq):
            
            if k == '.'*size:
                c.append(index_inner)
                
            int_.append((kf_df.reindex(list(k)).sum() / len(k)).to_numpy())
                
        stop_at.append(min(c)-1)
        embeddings.append(int_)
        
    euclidean_distances_ = []
    
    print('>Starting euclidean distance calculations...')
    
    for i in tqdm(range(0, len(embeddings))):
        
        arr_1 = array(embeddings[i])
        stop_1 = stop_at[i]
        for j in range(0, len(embeddings)):
            arr_2 = array(embeddings[j])
            stop_2 = stop_at[j]
            stop_ = max([stop_1, stop_2])
            euclidean_distances_.append(
                    np.mean(np.sqrt(np.sum((arr_1[0:stop_, :] - arr_2[0:stop_, :])**2, axis = 0)))
                    )
            
        
    distance_matrix = array(euclidean_distances_).reshape(len(embeddings), len(embeddings))
    
    return distance_matrix

def hammingDistance(a, b, size):
    
    """
    
    Calculates the hamming distances between the provided kmers and summarises. 
    
    Args:
        a- A list of kmers for an individual sequence
        b- A list of kmers for the sequence a is to be compared to
        
    Returns:
        the averaged hamming distance between sequences
    
    """
    
    percentage_matches = []
    
    for pos in range(0,len(a)):
        
        mismatch_int = 0
        
        kmer_a = a[pos]
        kmer_b = b[pos]
        
        if kmer_a == '.'*size and kmer_b == '.'*size:
            break
        
        else:
            for kmer_pos in range(0, len(kmer_a)):
            
                if kmer_a[kmer_pos] != kmer_b[kmer_pos]:
                    mismatch_int += 1
                    
                else:
                    mismatch_int += 0
    
        percentage_matches.append(mismatch_int / size)
    
    return sum(percentage_matches)/len(percentage_matches)

def similarityMeasure(data_input, method, names_, size):
    
    """
    Produces the similarity matrix used for hierarchical clustering
    
    kidera is currently in development
    
    Args:
        data_input- a list of kmers per sequence
        method- string value for similarity measure input. options are 'kidera',
        'positional' or 'jaccard'
        names_- a list of sequence ids
        size- a numeric variable for the size of the kmer to use
        
    Returns:
        A data frame of the pair-wise similarity scores per sequence
    
    """
    
    if method not in ['hamming', 'kidera', 'positional']:
        raise ValueError(f"\n>{method} is not a recognised similarity measure.\n>Options are 'hamming', 'kidera' or 'positional'")
    
    if method == 'kidera':
        kideras = read_csv('kidera.factors.tab',
                   sep = '\t', index_col = 0)
        
        distance_matrix = kideraDistanceFast(kmer_list = data_input,
                                             kf_df = kideras,
                                             size = size)
        return DataFrame(distance_matrix, columns = names_, index = names_)
    
    else:
        
        out_list = []
        
        for i in tqdm(range(0, len(data_input))):
            int_ = []
            for j in range(0, len(data_input)):
                if method == 'hamming':
                    int_.append(hammingDistance(data_input[i], data_input[j], size))
                elif method == 'positional':
                    int_.append(positionMatch(data_input[i], data_input[j], size))
            out_list.append(int_)
    
        return DataFrame(out_list, columns = names_, index = names_)
        
def runCode(data_path, similarity_,
            clone_select_loc = 'CLONE', region_select = 'CDR3_IMGT', translate_ = True,
            sequence_id = 'SEQUENCE_ID', split_column = 'MOUSE_ID',
            filter_column = 'SAMPLE_ID', filter_value = 'gp120 + mAb',
            kmer = 5, linkage_ = 'complete'):
    
    """
    Main code to run to decompose sequences into kmers and plot.
    """
    
    
    ### create folder in working directory
    
    mkdir(f'{date.today().strftime("%Y%m%d")}_{region_select}_{kmer}mer_{linkage_}')
    
    write_loc = f'{date.today().strftime("%Y%m%d")}_{region_select}_{kmer}mer_{linkage_}'
    
    if search('(\.tab|\.tsv)$', data_path):
        data_sep = '\t'
    elif search('\.csv$', data_path):
        data_sep = ','
    else:
        raise ValueError('>Accepted inputs are ".tsv", ".tab" and ".csv"')
    
    data = read_csv(data_path, data_sep)
    
    selection_cols = [clone_select_loc, region_select, sequence_id, split_column]
    
    for i in selection_cols:
        if i not in data.columns:
            raise ValueError(f'\n>Column {i} not in data frame\n')
    
    if filter_column is not None:
        data = data[data[filter_column] != filter_value]
    
    data.reset_index(drop = True, inplace = True)
    
    if translate_ == True:
        data[region_select] = [translate(transcribe(str(i), dna_type = 'dna')) for i in list(data[region_select])]
    
    region_group = data.groupby([split_column])[region_select].value_counts()
    region_group = DataFrame(region_group)
    region_group['COUNT'] = region_group[region_select]
    region_group[split_column] = [str(i).split(',')[0].split("'")[1] for i in region_group.index]
    region_group[region_select] = [str(i).split(',')[1].split("'")[1] for i in region_group.index]
    region_group.reset_index(drop = True, inplace = True)
    
    data = data.drop_duplicates([split_column, region_select]).reset_index(drop = True)
    data = data.sort_values([split_column, region_select])
    
    data = merge(data, region_group, how = 'outer')
    
    ### split into dictionary of dataframes
    split_values = data[split_column].unique()
    
    data_dict = {i : DataFrame() for i in split_values}
    
    for key in data_dict.keys():
        data_dict[key] = data[:][data[split_column] == key]
    
    ### get size of random sample
    sample_size = min([i.shape[0] for i in data_dict.values()])
    
    ### make random samples
    for key, value in data_dict.items():
        value = value.sample(n = sample_size, weights = value['COUNT'], replace = True, random_state = 42)        
        
        print(f'>Making {kmer}mers for {key}...')
        data_ = kmerDecomp(list(value[region_select]), kmer_size=kmer)
    
        print(f'>Making {similarity_} distance matrix...')
        similarity_table = similarityMeasure(data_,
                                             method = similarity_,
                                             names_ = list(value[region_select]),
                                             size = kmer)
        
        print(f'>Writing distance matrix to .tsv')                        
        similarity_table.to_csv(f'{write_loc}/{key}_{similarity_}-distance-matrix.tsv', sep = '\t')


#################################### DEFINE ARGUMENTS FROM CL ####################################

parser = argparse.ArgumentParser()
parser.add_argument("mode", help = "Define the similarity mode used for grouping.")
parser.add_argument("-f", "--file", help = "Path to ChangeO data file")
parser.add_argument("-c", "--clone", help = "Column containing the clonal family IDs")
parser.add_argument("-r", "--region", help = "Column containing the sequence region to analyse")
parser.add_argument("-t", "--translate", help = "Include to translate the dna sequences provided", action = 'store_true')
parser.add_argument("-s", "--sequence_id", help = "Column containing the unique sequence identifier")
parser.add_argument("-sc", "--splitcolumn", help = "Column containing the values to divide the data frame by")
parser.add_argument("-fc", "--filtercolumn", help = "Column for filtering out unnecessary data")
parser.add_argument("-fv", "--filtervalue", help = "Value to remove from filtercolumn", nargs = '+')
parser.add_argument("-k", "--kmer", help = "Numerical value for the size of the kmer")
parser.add_argument("-cm", "--cluster_method", help = "type of clustering for dendrogram. options are inputs to scipy dendrogram")

#################################### PARSE ARGUMENTS GIVEN AT CL ####################################

args = parser.parse_args()
sim_ = args.mode
file_ = args.file
clone_column = args.clone
region = args.region
tl = args.translate
seq_id = args.sequence_id
split_col = args.splitcolumn
filt_col = args.filtercolumn
filt_val = args.filtervalue
filt_val = ' '.join(filt_val)
kmer_ = args.kmer
cluster = args.cluster_method

"""

section here allows simpler definition of command line arguments for running in software like spyder

sim_ = 'kidera'
#sim_ = 'hamming'
#sim_ = 'positional'

#file_ = 'bulkDataForAnalysis.csv'
file_ = 'Sanger_Seq-Assay_results.csv'

clone_column = 'CLONE'
region = 'CDR3_IMGT'
tl = True
seq_id = 'SEQUENCE_ID'

#split_col = 'MOUSE_ID'
split_col = 'SAMPLE_ID'

#filt_col = 'SAMPLE_ID'
filt_col = 'EPITOPE_ADJUSTED'

#filt_val = 'gp120 + mAb'
filt_val = 'Non-specific'

kmer_ = 3
#kmer_ = 6

cluster = 'complete'

data_path = file_
similarity_ = sim_
clone_select_loc = clone_column
region_select = region
translate_ = tl
sequence_id = seq_id
split_column = split_col
filter_column = filt_col
filter_value = filt_val
kmer = kmer_
linkage_ = cluster

"""

#################################### MAKE AND SAVE KMERS ####################################
        
runCode(data_path = file_,
        similarity_ = sim_,
        clone_select_loc = clone_column,
        region_select = region,
        translate_ = tl,
        sequence_id = seq_id,
        split_column = split_col,
        filter_column = filt_col,
        filter_value = filt_val,
        kmer = kmer_,
        linkage_ = cluster
)

