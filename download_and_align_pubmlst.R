# Download sequences from pubmlst, and align using mafft
library(tidyverse)
library(httr)
library(seqinr)

out_dir <- 'caged/pubmlst'

# Downloads alleles from pubmlst, sleeping for 20 seconds beteen each one
download_allele <- function(locus, output_dir, sleep_between = 20) {
  # base URL for pubmlst
  
  cat("Downloading locus ", locus, "\n")
  base_url <-
    paste0("https://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef/loci/", locus, "/alleles_fasta")

  # filename to write
  file_name <-  file.path(output_dir, paste0(locus, ".fas"))
  
  # download the file
  got <- httr::GET(url = base_url)
  
  content(got) |> write_file(file_name)

  # sleep to play nice
  Sys.sleep(sleep_between)
}

# create output folders
out_raw <- file.path(out_dir, 'raw')
out_aligned <- file.path(out_dir, 'aligned')
fs::dir_create(out_raw)
fs::dir_create(out_aligned)

# load in the genes
genes <- read_csv("cgmlst_genes.csv") %>%
  pull(gene)

genes <- genes[1:2]

map(genes, download_allele, sleep_between = 0, output_dir=out_raw)

# ok, now align the genes
fasta_files <- list.files(out_raw, "*.fas", full.names = TRUE)

align_gene <- function(input_fasta, output_dir, threads=1, direction='adjust') {
  extras <- ""
  if (direction == "adjust") {
    extras <- "--adjustdirection"
  } else if (direction == "accurate") {
    extras <- "--adjustdirectionaccurately"
  }
  output_fasta <- file.path(output_dir, fs::path_file(input_fasta))
  message("aligning gene ", input_fasta)
  system(paste0("mafft --thread ",
                threads,
                " --quiet ",
                extras,
                " ",
                input_fasta,
                " > ",
                output_fasta))
}

map(fasta_files, align_gene, output_dir=out_aligned, threads=8)
