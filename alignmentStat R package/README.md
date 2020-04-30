The package is a wrapper for alignmentStat function. That function uses shuffling to produce a range of statistics and graphical outputs that help assess the statistical signficance of the alignment score of two protein or DNA sequences. The sequences are read from FASTA files. One of the sequences is iteratively shuffled and scored against the other sequence to obtain a sample distribution of scores. The function estimates empirical p-value and e-value out of the obtained distibution. A Gumbel distribution is fitted on the sample and used to provide standardized Gumbel score and p-value.

The zip and tarz folders "alignmentStat" contain the R package files:
	folder "man": help on package functions.
	folder "R": functions R scripts
	"DESCRIPTION": package description
	"NAMESPACE": dependencies

The PDF file "alignmentStat" is the reference manual of the package.

The Rmd file "Sample output of alignmenstat" is a parametrized R Markdown that lets user create dynamic reports with the package functions.

The R file "RenderScript Sample output of alignmentstat" is a script de to knit the Rmd "Sample output of alignmenstat"

The HTML files "example_DNA_morelessclose","example_proteins_close" y "example_proteins_distant" are examples of runs of the dynamic report.

Files "Q7IZ16.fasta", "P98052.fasta","Q96IY4.fasta", "B8D9R6.fasta", "gi32141095_N_1.fa", "gi32141095_N_0.fa" are FASTA files of example sequences.
