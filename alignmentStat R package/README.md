The package is a wrapper for alignmentStat function. That function uses shuffling to pro-
duce a range of statistics and graphical outputs that help assess the statistical signfi-
cance of the alignment score of two protein or DNA sequences. The se-
quences are read from FASTA files. One of the sequences is iteratively shuf-
fled and scored against the other sequence to obtain a sample distribution of scores. The func-
tion estimates empirical p-value and e-value out of the obtained distibution. A Gumbel distribu-
tion is fitted on the sample and used to provide standardized Gumbel score and p-value.
