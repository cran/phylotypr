## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# # devtools::install_github("mothur/phylotypr") #install development version
# install.packages("phylotypr") # install stable version from CRAN

## ----setup--------------------------------------------------------------------
library(phylotypr)

## -----------------------------------------------------------------------------
db <- build_kmer_database(rdp_trainset9_pds)

## -----------------------------------------------------------------------------
unknown <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTTAAGTCAGTGGTCAAATTGAGGGGCTCAACCCCTTCCCGCCATTGAAACTGGCGATCTTGAGTGGAAGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATGCCGGCTTCCTACTGACGCTCATGCACGAAAGTGTGGGTAACGAACAGG"

## ----echo = TRUE, eval = FALSE------------------------------------------------
# consensus <- classify_sequences(sequence = unknown, database = db, seed = 123)

## ----echo = FALSE-------------------------------------------------------------
consensus <- classify_sequences(sequence = unknown, database = db, seed = 123,
                                num_processors = 2)

## -----------------------------------------------------------------------------
consensus

## -----------------------------------------------------------------------------
filtered <- filter_taxonomy(consensus)

## -----------------------------------------------------------------------------
filtered

## -----------------------------------------------------------------------------
print_taxonomy(filtered)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# read_fasta(phylotypr_example("miseq_sop.fasta.gz")) |>
#   classify_sequences(database = db, seed = 123) |>
#   print_taxonomy() |>
#   head()

## ----echo = FALSE-------------------------------------------------------------
read_fasta(phylotypr_example("miseq_sop.fasta.gz")) |>
  classify_sequences(database = db, seed = 123, num_processors = 2) |>
  print_taxonomy() |>
  head()

## ----echo = TRUE, eval = FALSE------------------------------------------------
# strollur::new_dataset(dataset_name = "sop") |>
#   strollur::add(
#     table = phylotypr::read_fasta(phylotypr_example("miseq_sop.fasta.gz"))
#   ) |>
#   classify_sequences(db, seed = 123) |>
#   print_taxonomy() |>
#   head()

## ----echo = FALSE-------------------------------------------------------------
strollur::new_dataset(dataset_name = "sop") |>
  strollur::add(
    table = phylotypr::read_fasta(phylotypr_example("miseq_sop.fasta.gz"))
  ) |>
  classify_sequences(database = db, seed = 123, num_processors = 2) |>
  print_taxonomy() |>
  head()

## -----------------------------------------------------------------------------
parallelly::availableCores()

## -----------------------------------------------------------------------------
read_fasta(phylotypr_example("miseq_sop.fasta.gz")) |>
  classify_sequences(database = db, seed = 123, num_processors = 2) |>
  print_taxonomy() |>
  head()

## -----------------------------------------------------------------------------
classify_sequences(sequence = unknown, database = db, seed = 1) |>
  print_taxonomy()
classify_sequences(sequence = unknown, database = db, seed = 2) |>
  print_taxonomy()
classify_sequences(sequence = unknown, database = db, seed = 3) |>
  print_taxonomy()

## -----------------------------------------------------------------------------
classify_sequences(sequence = unknown, database = db, num_bootstraps = 1000,
                   seed = 123) |>
  print_taxonomy()

## ----eval = FALSE-------------------------------------------------------------
# devtools::install_github("mothur/phylotyprrefdata")
# library(phylotyprrefdata)

## ----eval = FALSE-------------------------------------------------------------
# data(package = "phylotyprrefdata")

