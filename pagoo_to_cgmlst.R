# code for blasting pagoo against cgmlst to find out which genes in the
# pagoo output align to which cgmlst genes.

# idea is then we can just use the pagoo output for everything.

library(tidyverse)
library(seqinr)

# read in the first isolate and it's cgmlst sequences. This will be
# used for the 'library' for blast.

labid <- 'SC0078'

input_path <- "output_combined"
sequence_files <- list.files(input_path, "*.fas.tab$")
get_cgmlst_sequence_for_isolate <- function(csv_file, isolate) {
  cat("processing file", csv_file, '\n')
  read.table(csv_file, header=TRUE) |>
    filter(file == paste0(isolate, ".fa")) |>
    pull(aligned_sequence) |>
    str_split(pattern='') |>
    unlist()
}

all_cgmlst_sequences <- 
  sequence_files |>
  set_names() |>
  map(~file.path(input_path, .)) |>
  map(~get_cgmlst_sequence_for_isolate(., isolate = labid))

cgmlst_file <- paste0("cgmlst_", labid, "_all.fas")
write.fasta(all_cgmlst_sequences, names=names(all_cgmlst_sequences),
            file.out = cgmlst_file)

# save this as a single fasta file

# read in the first isolate and it's pagoo sequences. This will be
# blasted against the library.

input_path <- "accessory/pagooOutput"

fasta_files <- list.files(input_path, "*.fas", recursive = TRUE)

# read in a pagoo sequence and pull out an isolate
get_pagoo_sequence_for_isolate <- function(fasta_file, isolate) {
  # read in the fasta file
  cat("processing file", fasta_file, '\n')
  gene <- read.fasta(fasta_file)
  gene_sequences <- tibble(isolate_id = names(gene),
                               sequence = gene) |>
    extract(isolate_id, into="isolate_clean", regex="(SC[0-9]*)", remove = FALSE)
  sequence <- gene_sequences |>
    filter(isolate_clean == isolate) |>
    pull(sequence)
  if (length(sequence)) {
    unlist(sequence) |>
      set_names(NULL)
  } else {
    return(NULL)
  }
}

all_pagoo_sequences <-
  fasta_files |>
  set_names() |>
  map(~file.path(input_path, .)) |>
  map(~ get_pagoo_sequence_for_isolate(., isolate = labid))

# save this as a single fasta file
pagoo_file <- paste0("pagoo_", labid, "_all.fas")
write.fasta(all_pagoo_sequences, names=names(all_pagoo_sequences),
            file.out = pagoo_file)

# create the blast database
system(paste0("makeblastdb -dbtype nucl -input_type fasta -in ", cgmlst_file," -out ","blast_db"))

# do the blasting
blast_file <- paste0("blast_out_", labid, "_pagoo.txt")
system(paste0("blastn -num_threads 4 -query ", pagoo_file, " -db ", "blast_db", " -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq' -out ",blast_file))

# read in the table
blast_csv_file <- paste0("blast_out_", labid, "_pagoo.csv")
blast_out <- read.table(blast_file, header=FALSE) |>
  set_names(c("pagoo", "cgmlst", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen","qseq", "sseq"))
blast_out |>
  write_csv(blast_csv_file)

blast_out |>
  dplyr::count(pagoo) |>
  pull(n) |> unique() # Hmm, some pagoo groups have multiple outputs

blast_out |>
  dplyr::count(cgmlst) |>
  pull(n) |> unique() # Hmm, some cgmlst groups have multiple outputs

blast_out |>
  group_by(cgmlst) |>
  summarise(n = n(), nd = n_distinct(pagoo)) |>
  filter(n > 1) |>
  print(n = 100)
