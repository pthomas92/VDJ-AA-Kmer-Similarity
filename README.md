# VDJ-AA-Kmer-Similarity
Scripts to perform pairwise comparison of VDJ kmer sequences using various metrics, used during my PhD and then uploaded for tracking/future proofing. Inspired by Ostmeyer et al., (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1814-6).

_makeKmers.py:_
The script works by taking a repertoire of sequences (can define which section of the VDJ you're interested in (full length/CDR3 only)), decomposing them into overlapping kmers of a predefined length, and then performing all vs all pairwise comparisons to measure the similarity within repertoires. Multiple repertoires should be loaded in from the same source file, and are identified based on arguments supplied to the script. This allows down-sampling to the size of the smallest repertoire for more reliable comparisons.

Similarity measures are:
* Positional similarity- Kmer receives a score of 1 if the kmer is identical between comparisons, and 0 if not identical. Score is summed, and then divided by the number of kmers for analysis giving an overall percentage match.
* Hamming distance- Hamming distance is calculated between each kmer, and then averaged as above to give a percentage match per sequence comparison.
* Kidera distance- Kmers in sequences for comparison are arranged into n x 10 matrices, where n is each kmer, and 10 reflects the 10 kidera factor values assigned to each kmer. Numpy is then used to perform fast calculation of the euclidean distance between two sequences per position. This is then averaged to give an overall distance.

Input arguments:
- data_path = path to repertoire data file.
- similarity_ = similarity measure used (e.g., positional, kidera, hamming).
- clone_select_loc = column that has clonal definitions, to help deduplication.
- region_select = region of the VDJ that is analysed (e.g., CDR3_IMGT, if using CDR3s in changeo format).
- translate_ = boolean value, should the region be translated prior to analysis? Gaps are not allowed.
- sequence_id = sequence identifier used for final square matrix output.
- split_column = column containing grouping information (e.g., different groups for testing).
- filter_column = column containing values to remove data from (e.g., to only consider a antigen specific cells, I used the EPITOPE_ADJUSTED column here).
- filter_value = value to be removed from `filter_column` (e.g., to only consider antigen specific cells, I used the `Non-specific` value here).
- kmer = size of the kmer for analysis.
- linkage_ = Legacy code that is not removed. Used to write the output folder name.
  
_AnalyseKmerOutput.R:_
The script takes the distance matrix assembled by _makeKmers.py_ and performs hierarchical clustering according to user preferences (I use complete linkage). Cophenetic distances are then calculated as a measure of how diffuse or closely related the sequences within each repertoire are. It then creates heatmaps and boxplots summarising the similarity between the defined groups.

