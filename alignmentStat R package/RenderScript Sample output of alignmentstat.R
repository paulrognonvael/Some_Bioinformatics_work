#Our package is not yet available on CRAN. To install it you need to unzip the file "alignmentStat.zip",
#install the package from the unzipped folder (path1) and load it, using the following code.

path1="~/Documents/MESIO/First Semester/Bioinformatica//Task 3/alignmentStat_doc_ex_english/alignmentStat"
install.packages(path1,repos = NULL, type="source")
library(alignmentStat)

rmarkdown::render("Sample output of alignmentstat.Rmd",
                  output_file = "example_DNA_morelessclose.html",
                  encoding="UTF-8",
                  params = list(  sequence1="gi32141095_N_1.fa",
                                  sequence2="gi32141095_N_0.fa",
                                  path="~/Documents/MESIO/First Semester/Bioinformatica/Task 3/alignmentStat_doc_ex_english",
                                  type="DNA",
                                  alignment="local",
                                  matrix="BLOSUM50",
                                  opening_penal=-3,
                                  ext_penal=-1,
                                  times_shuffling=10000,
                                  seq_shuffled=1)
                 )

rmarkdown::render("Sample output of alignmentstat.Rmd",
                  output_file = "example_proteins_close.html",
                  encoding="UTF-8",
                  params = list(  sequence1="Q7IZ16.fasta",
                                  sequence2="P98052.fasta",
                                  path="~/Documents/MESIO/First Semester/Bioinformatica/Task 3/alignmentStat_doc_ex_english",
                                  type="protein",
                                  alignment="local",
                                  matrix="BLOSUM80",
                                  opening_penal=-5,
                                  ext_penal=-1,
                                  times_shuffling=10000,
                                  seq_shuffled=1)
)


rmarkdown::render("Sample output of alignmentstat.Rmd",
                  output_file = "example_proteins_distant.html",
                  encoding="UTF-8",
                  params = list(  sequence1="B8D9R6.fasta",
                                  sequence2="Q96IY4.fasta",
                                  path="~/Documents/MESIO/First Semester/Bioinformatica/Task 3/alignmentStat_doc_ex_english",
                                  type="protein",
                                  alignment="global",
                                  matrix="PAM120",
                                  opening_penal=-3,
                                  ext_penal=-0.5,
                                  times_shuffling=10000,
                                  seq_shuffled=2)
)
