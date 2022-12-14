library(tidyverse)
library(seqinr)
library(DECIPHER)

# folder to read genes from
unaligned_folder <- "output_combined"

# read in the unaligned genes
fasta_files <- list.files(unaligned_folder, "*.fas.tab", full.names = TRUE)

# function for aligning our genes
align_gene <- function(input_fasta, output_fasta, threads=1, direction='none') {
  extras <- ""
  if (direction == "adjust") {
    extras <- "--adjustdirection"
  } else if (direction == "accurate") {
    extras <- "--adjustdirectionaccurately"
  }
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

align_table <- function(file_tab) {
  # read in the table
  cat("aligning ", file_tab, "\n")
  tab <- read.table(file_tab, header=TRUE, stringsAsFactors = FALSE)
  # create fasta input and output files
  fasta <- tempfile()
  aligned <- tempfile()
  write.fasta(sequences = as.list(tab$sequence),
              names = as.list(tab$file),
              file.out=fasta)
  align_gene(fasta, aligned, threads=8, direction="adjust")
  dna <- readDNAStringSet(aligned)
  tab$aligned_sequence = paste(dna)
  write.table(tab, file_tab, row.names=FALSE)
  return(0)
}

map(fasta_files, align_table)

#map(fasta_files, align_gene, output_dir=aligned_folder, threads=8)
