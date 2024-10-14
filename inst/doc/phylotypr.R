## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  # devtools::install_github("mothur/phylotypr") #install development version
#  install.packages("phylotypr") # install stable version from CRAN

## ----setup--------------------------------------------------------------------
library(phylotypr)

set.seed(19760620) # pat's birtday in YYYYMMDD format

## -----------------------------------------------------------------------------
db <- build_kmer_database(
  trainset9_pds$sequence,
  trainset9_pds$taxonomy
)

## -----------------------------------------------------------------------------
unknown <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTTAAGTCAGTGGTCAAATTGAGGGGCTCAACCCCTTCCCGCCATTGAAACTGGCGATCTTGAGTGGAAGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATGCCGGCTTCCTACTGACGCTCATGCACGAAAGTGTGGGTAACGAACAGG"

## -----------------------------------------------------------------------------
consensus <- classify_sequence(unknown = unknown, database = db)

## -----------------------------------------------------------------------------
consensus

## -----------------------------------------------------------------------------
filtered <- filter_taxonomy(consensus)

## -----------------------------------------------------------------------------
filtered

## -----------------------------------------------------------------------------
print_taxonomy(filtered)

## -----------------------------------------------------------------------------
library(dplyr)
library(purrr)

set.seed(19760620) # pat's birtday in YYYYMMDD format

miseq <- read_fasta(phylotypr_example("miseq_sop.fasta.gz"))

miseq |>
  dplyr::mutate(
    classification = purrr::map_chr(
      sequence,
      ~ classify_sequence(unknown = .x, database = db) |>
        filter_taxonomy() |>
        print_taxonomy(),
      .progress = TRUE
    )
  )

## ----eval = FALSE-------------------------------------------------------------
#  library(dplyr)
#  library(furrr)
#  
#  miseq <- read_fasta(phylotypr_example("miseq_sop.fasta.gz"))
#  
#  plan(strategy = multisession, workers = 4)
#  options(future.globals.maxSize = 10000000000)
#  
#  miseq |>
#    mutate(
#      classification = future_map_chr(
#        sequence,
#        ~ classify_sequence(unknown = .x, database = db) |>
#          filter_taxonomy() |>
#          print_taxonomy(),
#        .progress = TRUE,
#        .options = furrr_options(seed = 19760620)
#      )
#    )

## -----------------------------------------------------------------------------
map_chr(
  rep(unknown, 3),
  ~ classify_sequence(unknown = .x, database = db, num_bootstraps = 100) |>
    filter_taxonomy() |>
    print_taxonomy()
)

## -----------------------------------------------------------------------------
map_chr(
  rep(unknown, 3),
  ~ classify_sequence(unknown = .x, database = db, num_bootstraps = 100) |>
    filter_taxonomy() |>
    print_taxonomy()
)

## ----eval = FALSE-------------------------------------------------------------
#  devtools::install_github("mothur/phylotyprrefdata")
#  library(phylotyprrefdata)

## ----eval = FALSE-------------------------------------------------------------
#  data(package = "phylotyprrefdata")

