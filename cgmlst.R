## Set up system

## sudo apt-get mafft
## sudo apt-get ncbi-tools

## Download the Campy cgMLST alleles from here: https://pubmlst.org/bigsdb?db=pubmlst_campylobacter_seqdef&page=downloadAlleles&tree=1
## I do this using the "batch link downloader" plugin for Chrome. You probably have a better way - maybe the API?

## Create this folder structure:

## Your_parent_folder
## Your_parent_folder/sequences     -     put the downloaded allele files in here, extensions should be ".fas"
## Your_parent_folder/db            -     empty folder, blast databases are written to here
## Your_parent_folder/output        -     empty folder, results go in here

# and for the data
## Folder_with_the_SACNZ_data     -     should have all the assembled genome fasta files in one folder. Extensions can be any of ".fa", ".fasta", ".fas" or ".fna".



## Load the functions in this file into the R environment
source('cgmlst_functions.R')

## Run this shiz
run_cgmlst(db_loc="pubmlst/sequences",
           qry_fld="data/SACNZ_genomes")

fix_alleles(results_folder = "pubmlst/output",
            threads=4)                      #### Add more threads if you have them

tab=summarise_alleles(results_folder = "pubmlst/output") 

