params <-
list(path = "", wig = "", gff = "")

## ---- include=FALSE---------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, include = FALSE)


## ---- include=FALSE---------------------------------------------------------
# libraries

library(kableExtra)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(readxl)

# function to format printed numbers
formatgraph <- function(x) {
  formatC(x, format = "e", digits = 3)
}


# function to write matrix in markdown
write_matex2 <- function(x) {
  x <- as.matrix(x)
  begin <- "\\begin{bmatrix}"
  end <- "\\end{bmatrix}"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  paste(c(begin, X, end), collapse = "")
}


## ---------------------------------------------------------------------------
# Read the data
data <- read_table2(paste0(params$path, params$wig), col_names = FALSE, skip = 1)
colnames(data) <- c("variableStep", "count")


## ----include=TRUE-----------------------------------------------------------
# print sample data
set.seed(426)
data %>%
  sample_n(10) %>%
  kable("latex",
    booktabs = T,
    caption = "\\label{tab:tab1} Sample of data"
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ---------------------------------------------------------------------------
# States
states <- c("ES", "GD", "NE", "GA")
N <- length(states)


## ---------------------------------------------------------------------------
### Estimation of HMM model emission function and transition matrix

# Estimation of geometric distribution parameters for emission probabilities
reads <- data$count
reads_nz <- sort(reads[reads != 0 ])
size <- length(reads_nz)
mean_r <- mean(reads_nz[1:round(0.95 * size)])

mu <- c(1 / 0.99, 0.01 * mean_r + 2, mean_r, mean_r * 5.0)
L <- 1.0 / mu

fun1 <- function(i) {
  fun2 <- function(x) {
    dgeom(x, prob = L[i])
  }
  return(fun2)
}


B <- vector("list", 4)

for (i in 1:N) {
  B[[i]] <- fun1(i)
}


## ----include=TRUE-----------------------------------------------------------
# print sample data

df <- data.frame(states, formatgraph(L))
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab2} Emission distribution parameters",
  col.names = c("State", "$\\theta$"), escape = FALSE
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ---------------------------------------------------------------------------
# Estimation of the transition matrix
calculate_pins <- function(O) {
  non_ess_reads <- c()
  temp <- c()
  for (rd in O) {
    if (rd >= 1) {
      if (length(temp) < 10) {
        non_ess_reads <- c(non_ess_reads, temp)
      }
      non_ess_reads <- c(non_ess_reads, rd)
      temp <- c()
    } else {
      temp <- c(temp, rd)
    }
  }
  return(sum(non_ess_reads >= 1) / length(non_ess_reads))
}

pins <- calculate_pins(reads)


## ---------------------------------------------------------------------------
pnon <- 1.0 - pins

for (r in 0:100) {
  if ((pnon**r) < 0.01) {
    break
  }

  A <- matrix(nrow = N, ncol = N)
  a <- log1p(-(B[[3]](0)**r))
  # b <- log((1-exp(a))/3)
  b <- (r * log(B[[3]](0))) + log(1.0 / 3)
  for (i in 1:4) {
    A[i, ] <- rep(b, N)
    A[i, i] <- a
  }
}


## ---------------------------------------------------------------------------
# Initial probabilities

PI <- rep(0, N) # Initial state distribution
PI[1] <- 0.7
PI[2:N] <- 0.3 / (N - 1)


## ----include=TRUE-----------------------------------------------------------
# print initial probabilities
df <- data.frame(states, round(PI, 2))
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab3} Initial probabilities",
  col.names = c("State", "$\\pi_0$"), escape = FALSE
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ---------------------------------------------------------------------------
### Viterbi algorithm

# Notation: correspondance with the description of Viterbi algorithm by
# Cristianini et Hahn in "Introduction to Computational Genomics: A Case Studies Approach"

# B is the function giving E the emission matrix
# A is the transmission matrix T
# delta corresponds to the V matrix, the sequence of loglikelihood of the most probable hidden sequence
# Q is the pointer

O <- reads

N <- length(B)
T1 <- length(O)
delta <- matrix(ncol = T1, nrow = N)

# first run
b_o <- c()
for (i in 1:N) {
  b_o <- c(b_o, B[[i]](O[1]))
}

delta[, 1] <- log(PI) + log(b_o)

Q <- matrix(ncol = T1, nrow = N)


B <- vector("list", 4)

for (i in 1:N) {
  B[[i]] <- fun1(i)
}

for (t in 2:T1) {
  b_o <- c()
  for (i in 1:N) {
    b_o <- c(b_o, B[[i]](O[t]))
  }
  nus <- delta[, t - 1] + A
  delta[, t] <- apply(nus, 2, FUN = max) + log(b_o)
  Q[, t] <- apply(nus, 2, FUN = which.max)
}

# traceback
Q_opt <- c(which.max(delta[, T1]))

for (t in T1:2) {
  Q_opt <- c(Q[Q_opt[1], t], Q_opt)
}


## ----include=TRUE-----------------------------------------------------------
# print final site log likelihood
df <- data.frame(states, delta[, T1])
colnames(df) <- c("State", "Log-likelihood")
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab4} Final site log-likelihood"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ---------------------------------------------------------------------------
# print frequency of states in most probable path
df <- data.frame(round(100 * table(Q_opt) / length(Q_opt), 2))
df[, 1] <- states
colnames(df) <- c("State", "Total % of genome")
rownames(df) <- NULL
write.csv(df, paste0("state_freq_initprob", paste0(PI[1:2], collapse = "_"), ".csv"))
df2 <- read_csv("state_freq_initprob0.25_0.25.csv", col_types = cols(X1 = col_skip()))


## ----include=TRUE-----------------------------------------------------------
kable(list(df, df2), "latex",
  booktabs = T,
  caption = "\\label{tab:tab5} State frequency in TA sites with proposed (left) and equal (right) initial probabilities"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ---------------------------------------------------------------------------
### Results processing


data$state <- Q_opt
data$state <- factor(data$state, levels = c(1, 2, 3, 4), labels = states)

data <- data %>%
  arrange(variableStep)

data <- data %>%
  mutate(
    positive_count_flag = (count > 0),
    positive_count = case_when(count > 1 ~ count)
  )


# Summary TA x state

taxstate_summary <- data %>%
  group_by(state) %>%
  summarise(
    number_TA = n(),
    perc_TA = 100 * n() / nrow(data),
    total_count = sum(count),
    mean_count = mean(positive_count),
    insertion_dens = sum(positive_count_flag) / n()
  )


## ---------------------------------------------------------------------------
### Region boundaries

data <- data %>%
  mutate(region = rep(NA, nrow(data)))

# runs of equal values of state
rle_state <- rle(as.character(data$state))
rle_state <- data.frame(length = rle_state$lengths, state = rle_state$values)

data$region[1:rle_state$length[1]] <- paste0("Region", 1)

k <- 2
for (i in 2:nrow(rle_state)) {
  start <- sum(rle_state$length[1:(i - 1)]) + 1
  end <- start + rle_state$length[i] - 1
  data$region[start:end] <- paste0("Region", k)
  k <- k + 1
}


## ---------------------------------------------------------------------------
# Summary by region
region_summary <- data %>%
  group_by(region) %>%
  summarise(
    region.state = levels(data$state)[which.max(table(state))],
    number_TA = n(),
    total_count = sum(count),
    mean_count = max(mean(positive_count, na.rm = TRUE), 0, na.rm = TRUE),
    insertion_dens = sum(positive_count_flag) / n()
  )
# Summary region x state
regionxstate_summary <- region_summary %>%
  group_by(region.state) %>%
  summarise(
    mean_nb_TA = round(mean(number_TA), 1),
    mean_insertion_density = round(mean(insertion_dens), 3),
    mean_reads_count = round(mean(mean_count), 2)
  )


## ----include=TRUE-----------------------------------------------------------
# print statistics region:state

df <- regionxstate_summary %>% arrange(mean_insertion_density)

colnames(df) <- c("State", "Mean # TA sites", "Mean insertion density", "Mean read counts")
rownames(df) <- NULL
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab6} Statistics for state classification on regions"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ----include=TRUE,fig.cap="\\label{fig:fig1} Mean insertion density and read counts for regions", fig.pos = 'h!', fig.align="center"----
ggplot(region_summary) + geom_point(aes(
  x = insertion_dens,
  y = mean_count,
  col = region.state,
  shape = region.state
)) +
  xlab("Insertion frequency") +
  ylab("Mean read counts") +
  labs(col = "State", shape = "State")


## ---------------------------------------------------------------------------
#### Gene boundaries

data <- data %>%
  arrange(variableStep)

## Identified genes
gff3 <- readGFF(paste0(params$path, params$gff), version = 3)

hash_gene <- as_tibble(gff3) %>%
  filter(type == "gene") %>%
  arrange(start)

data <- data %>%
  mutate(
    gene.id = rep(NA, nrow(data)),
    gene.name = rep(NA, nrow(data)),
    gene.description = rep(NA, nrow(data))
  )

for (i in 1:nrow(hash_gene)) {
  start <- as.numeric(hash_gene[i, "start"])
  end <- as.numeric(hash_gene[i, "end"])
  data <- data %>%
    dplyr::filter(variableStep >= start & variableStep <= end) %>%
    mutate(
      gene.id = as.character(hash_gene[i, "gene_id"]),
      gene.name = as.character(hash_gene[i, "Name"]),
      gene.description = as.character(hash_gene[i, "description"])
    ) %>%
    rbind(data %>% dplyr::filter(!(variableStep >= start & variableStep <= end)))
}

genes_in_seq <- levels(factor(data$gene.id))
genes_in_gff <- hash_gene$gene_id
no_in_gff_notin_seq <- genes_in_gff[!genes_in_gff %in% genes_in_seq]


## ---------------------------------------------------------------------------
## Non-coding regions

data <- data %>%
  arrange(variableStep)

# runs of equal values of gene.id
rle_gene.id <- rle(is.na(data$gene.id))
rle_gene.id <- data.frame(length = rle_gene.id$lengths, isna = rle_gene.id$values)

even <- seq(2, nrow(rle_gene.id), 2)

k <- 1
for (i in even) {
  start <- sum(rle_gene.id$length[1:(i - 1)]) + 1
  end <- start + rle_gene.id$length[i] - 1
  data$gene.id[start:end] <- paste0("non_protein_coding", k)
  k <- k + 1
}


## ---------------------------------------------------------------------------
# Summary genes
gene_summary <- data %>%
  group_by(gene.id) %>%
  summarise(
    gene.state = levels(data$state)[which.max(table(state))],
    gene.name = max(gene.name),
    gene.description = max(gene.description),
    number_TA = n(),
    total_count = sum(count),
    mean_count = max(mean(positive_count, na.rm = TRUE), 0, na.rm = TRUE),
    insertion_dens = sum(positive_count_flag) / n()
  )
# Summary genes x state
genexstate_summary <- gene_summary %>%
  filter(!grepl("non_protein", gene.id)) %>%
  group_by(gene.state) %>%
  summarise(
    mean_nb_TA = round(mean(number_TA), 1),
    mean_insertion_density = round(mean(insertion_dens), 3),
    mean_reads_count = round(mean(mean_count), 2)
  )


## ----include=TRUE-----------------------------------------------------------
# print frequency of states in most probable path
df <- data.frame(round(100 * table(gene_summary$gene.state) / nrow(gene_summary), 2))
df[, 1] <- states
colnames(df) <- c("State", "Total % of genome")
rownames(df) <- NULL
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab7} State frequency in genes"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ----include=TRUE-----------------------------------------------------------
# print statistics region:state

df <- genexstate_summary %>% arrange(mean_insertion_density)

colnames(df) <- c("State", "Mean # TA sites", "Mean insertion density", "Mean read counts")
rownames(df) <- NULL
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab8} Statistics for state classification on genes"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ----include=TRUE,fig.cap="\\label{fig:fig2} Mean insertion density and read counts for regions", fig.pos = 'h!', fig.align="center"----
ggplot(gene_summary) + geom_point(aes(
  x = insertion_dens,
  y = mean_count,
  col = gene.state,
  shape = gene.state
)) +
  xlab("Insertion frequency") +
  ylab("Mean read counts") +
  labs(col = "State", shape = "State")


## ---------------------------------------------------------------------------
### Gumbel/EVD Variance and Expected run functions

# Variance
VarR <- function(n, p) {
  # VarR_n =  (pi^2)/(6*ln(1/p)^2) + 1/12 + r2(n) + E2(n) (Schilling, 1990)

  r2 <- getR2(n)
  E2 <- getE2(n)

  A <- pi**2 / (6 * log(1 / p)**2)
  V <- A + 1 / 12 + r2 + E2

  return(V)
}

# Expectation
ExpectedRuns <- function(n, p) {
  # ER_n =  log(1/p)(nq) + gamma/ln(1/p) -1/2 + r1(n) + E1(n) (Schilling, 1990)

  q <- 1 - p
  gamma <- getGamma()
  r1 <- getR1(n)
  E1 <- getE1(n)


  A <- log(n * q, base = 1.0 / p)
  B <- gamma / log(1 / p)
  ER <- A + B - 0.5 + r1 + E1

  return(ER)
}

getGamma <- function() {
  # Euler-Mascheroni constant ~ 0.577215664901
  return(0.5772156649015328606)
}

getR1 <- function(n) {
  # Small Correction term. Defaults to 0.000016
  return(0.000016)
}

getR2 <- function(n) {
  # Small Correction term. Defaults to 0.00006
  return(0.00006)
}

getE1 <- function(n) {
  # Small Correction term. Defaults to 0.01
  return(0.01)
}

getE2 <- function(n) {
  # Small Correction term. Defaults to 0.01
  return(0.01)
}


theta <- sum(data$count > 0) / nrow(data)


## ---------------------------------------------------------------------------
### Reassignment of essentiality based on EVD

# Summary genes x state
gene_summary2 <- data %>%
  group_by(gene.id) %>%
  summarise(
    gene.name = max(gene.name),
    gene.description = max(gene.description),
    number_TA = n(),
    total_count = sum(count),
    mean_count = max(mean(positive_count, na.rm = TRUE), 0, na.rm = TRUE),
    insertion_dens = sum(positive_count_flag) / n(),
    exp_max_len_es = ExpectedRuns(n(), 1 - theta),
    var_len = VarR(n(), 1 - theta),
    n_ES = table(state)["ES"],
    n_GD = table(state)["GD"],
    n_NE = table(state)["NE"],
    n_GA = table(state)["GA"],
    gene.state = levels(data$state)[which.max(table(state))]
  )

# Reassignment of state
gene_summary3 <- gene_summary2 %>%
  mutate(gene.state2 = case_when(
    n_ES == number_TA ~ "ES",
    n_ES >= exp_max_len_es + 3 * sqrt(var_len) ~ "ES",
    TRUE ~ gene.state
  ))


# Summary genes x state
genexstate_summary2 <- gene_summary3 %>%
  filter(!grepl("non_protein", gene.id)) %>%
  group_by(gene.state2) %>%
  summarise(
    mean_insertion_density = mean(insertion_dens),
    mean_reads_count = mean(mean_count)
  )

# differences before and after assignment by EVD
table(gene_summary3$gene.state, gene_summary3$gene.state2)


## ----include=TRUE-----------------------------------------------------------
# print differences before and after assignment by EVD
df <- as.data.frame.matrix(table(gene_summary3$gene.state, gene_summary3$gene.state2))
df <- df[, states]
df <- data.frame(State = c("ES", "GA", "GD", "NE"), df)
df <- df[c(1, 3, 4, 2), ]
colnames(df) <- c("", states)
rownames(df) <- NULL

kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab9} Contingency table of essentiality assignment"
) %>%
  add_header_above(c(" " = 1, "Refined assignment" = 4)) %>%
  pack_rows("Original", 1, 4) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ---------------------------------------------------------------------------
# Comparison to TraSH
non_essential_traSH <- read_excel(paste0(params$path, "mmi_3425_sm_tables3.xls"), skip = 1)
ne_genes_trash <- non_essential_traSH$`Rv designation`

trash_genes <- data.frame(gene.id = ne_genes_trash, gene.state.trash = "NE")

ne_genes_hmm <- gene_summary3 %>%
  filter(!grepl("non_protein", gene.id)) %>%
  filter(gene.state2 == "NE") %>%
  select(gene.id)
ne_genes_hmm <- ne_genes_hmm$gene.id
length(ne_genes_trash %in% ne_genes_hmm)

growth_defect_traSH <- read_excel(paste0(params$path, "mmi_3425_sm_tables2.xls"), skip = 1)
growth_defect_trash <- growth_defect_traSH$`Rv designation`
trash_genes <- rbind(trash_genes, data.frame(gene.id = growth_defect_trash, gene.state.trash = "GD"))

essential_traSH <- read_excel(paste0(params$path, "mmi_3425_sm_tables1.xls"), skip = 1)
essential_trash <- essential_traSH$`Rv designation`
trash_genes <- rbind(trash_genes, data.frame(gene.id = essential_trash, gene.state.trash = "ES"))

gene_summary4 <- gene_summary3 %>% full_join(trash_genes, by = "gene.id")

coding_genes <- gene_summary4 %>%
  filter(!grepl("non_protein", gene.id))

table(gene_summary4$gene.state.trash, gene_summary4$gene.state2)


## ----include=TRUE-----------------------------------------------------------
# pint differences with TraSH
df <- as.data.frame.matrix(table(gene_summary4$gene.state.trash, gene_summary4$gene.state2))
df <- df[, states]
df <- data.frame(State = c("NE", "GD", "ES"), df)
df <- df[c(3, 2, 1), ]
colnames(df) <- c("", states)
rownames(df) <- NULL
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab10} Contingency table of essentiality assignment with TraSH"
) %>%
  add_header_above(c(" " = 1, "HMM" = 4)) %>%
  pack_rows("TraSH", 1, 3) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )


## ----include=TRUE-----------------------------------------------------------
# Notable growth defect and growth advantage
notable_gd <- c("Rv0015c", "Rv0016c", "Rv0467", "Rv2379c", "Rv2380c", "Rv2381c", "Rv2382c", "Rv3841", "Rv0126", "Rv1097c", "Rv1098c", "Rv1099c")

df <- gene_summary3 %>%
  filter(gene.id %in% notable_gd) %>%
  select(gene.id, gene.state2, gene.name, insertion_dens, number_TA, mean_count)
colnames(df) <- c("Orf Ids", "State", "Included genes", "Insertion density", "Length", "Average reads")
df[, 4] <- round(df[, 4], 3)
df[, 6] <- round(df[, 6], 1)
rownames(df) <- NULL
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab11} Notable Growth-Defect genes"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )

notable_ga <- c("Rv3295", "Rv3296", "Rv2939", "Rv2940c", "Rv2941", "Rv2411c", "Rv0483", "Rv2930", "Rv2931", "Rv2932", "Rv2933", "Rv2934", "Rv2935", "Rv1843c", "Rv1844c", "Rv0554", "Rv0479c", "Rv0480c", "Rv0481c")

df <- gene_summary3 %>%
  filter(gene.id %in% notable_ga) %>%
  select(gene.id, gene.state2, gene.name, insertion_dens, number_TA, mean_count)
colnames(df) <- c("Orf Ids", "State", "Included genes", "Insertion density", "Length", "Average reads")
df[, 4] <- round(df[, 4], 3)
df[, 6] <- round(df[, 6], 1)
rownames(df) <- NULL
kable(df, "latex",
  booktabs = T,
  caption = "\\label{tab:tab12} Notable Growth-Advantage genes"
) %>%
  kable_styling(
    latex_options = c("hold_position"),
    position = "center"
  )

