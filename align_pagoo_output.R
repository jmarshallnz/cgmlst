library(tidyverse)
source("align_fasta_files.R")

which <- "shell"
input_path <- file.path("accessory/pagooOutput", which)
output_path <- file.path("accessory/pagooOutput/aligned_adjust", which)

fasta_files <- list.files(input_path, "*.fas", full.names = TRUE)

fs::dir_create(output_path)

# about 6 seconds
#map(fasta_files, align_gene, output_dir=output_path, threads=8)

# about 9 seconds
system.time({
map(fasta_files, align_gene, output_dir=output_path, threads=8, direction='adjust')
})

# takes a long time...
#map(fasta_files, align_gene, output_dir=output_path, threads=8, direction='accurate')

fasta_files <- list.files(output_path, "*.fas", full.names = TRUE)

fasta_to_dataframe <- function(fasta_file) {
  raw <- read.fasta(fasta_file)
  gene_name <- fasta_file |> fs::path_file() |> fs::path_ext_remove()
  map(raw, paste, collapse='') |>
    unlist() |>
    enframe(name = 'labid', value = 'sequence') |>
    extract(labid, into='labid', regex='[_R]*([^_]*)__*') |>
    mutate(gene = gene_name) |>
    select(labid, gene, sequence)
}

core <- map(fasta_files, fasta_to_dataframe) |> list_rbind()

write_csv(core, file.path(output_path, "sequence_data.csv.gz"))
