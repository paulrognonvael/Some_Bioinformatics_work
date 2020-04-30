setwd("~/MESIO/Bioinformatica/Trabajo final")

rmarkdown::render("Final Assignment - Paul Rognon.Rmd",
                  encoding="UTF-8",
                  output_file = "Final Assignment Paul Rognon.pdf",
                  params = list(path="~/MESIO/Bioinformatica/Trabajo final/",
                                wig="glycerol.wig",
                                gff="h37rv.gff3")
)
