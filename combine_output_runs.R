library(tidyverse)
# combine output runs from two folders

folder1 <- "output_original_run"
folder2 <- "output_ref"
outdir  <- "output_combined"

gene_files <- list.files(folder1, '*.fas.tab')

combine_tables <- function(gene_file, folder1, folder2, outdir) {
  tab1 <- read.table(file.path(folder1, gene_file), stringsAsFactors=FALSE, header=TRUE)
  tab2 <- read.table(file.path(folder2, gene_file), stringsAsFactors=FALSE, header=TRUE)
  comb_tab <- rbind(tab1, tab2)
  write.table(comb_tab, file.path(outdir, gene_file), row.names=FALSE)
  return(0)
}

fs::dir_create(outdir)

map(gene_files, combine_tables, folder1=folder1, folder2=folder2, outdir = outdir)


