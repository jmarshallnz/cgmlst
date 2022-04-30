# Download sequences from pubmlst, and align using mafft
library(tidyverse)
library(httr)
library(seqinr)

# Downloads alleles from pubmlst, sleeping for 20 seconds beteen each one
download_allele <- function(locus, output_dir, sleep_between = 20) {
  # base URL for pubmlst
  base_url <- "https://pubmlst.org/bigsdb?db=pubmlst_campylobacter_seqdef&page=downloadAlleles"
  
  # filename to write
  file_name <-  file.path(output_dir, paste0(locus, ".fas"))
  
  # parse the URL and add the locus to download
  parsed_url <- parse_url(base_url)
  parsed_url$query[['locus']] <- locus

  # download the file
  download.file(build_url(parsed_url), file_name)
  
  # sleep to play nice
  Sys.sleep(sleep_between)
}

# create output folders
fs::dir_create('raw')
fs::dir_create('aligned')

# load in the genes
genes <- read_csv("cgmlst_genes.csv") %>%
  pull(gene)

genes <- genes[1:2]

map(genes, download_allele, sleep_between = 0, output_dir="raw")

# ok, now align the genes
fasta_files <- list.files('raw', "*.fas", full.names = TRUE)

align_gene <- function(input_fasta, output_dir, threads=1, direction='none') {
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

map(fasta_files, align_gene, output_dir="aligned", threads=8)
