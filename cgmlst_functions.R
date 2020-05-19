library(dplyr)
library(DECIPHER) ## BIOCONDUCTOR
library(seqinr)

makedb <- function(ref=NULL){
  ##Makes the BLAST database for the reference allele data.
  print("Generating BLAST database")
  system(paste0("makeblastdb -dbtype nucl -input_type fasta -in ",ref," -out ",paste0(dirname(ref),"/../db/",basename(ref))))
  print("Done.")
}

cgmlst <- function(ref, isolate_dir, output_dir, temp_dir=NULL) {
  ## Calls BLAST against the reference allele database, and outputs the sequence data of EITHER
  ## 1) Exact matches - get "identity" = 1, and are assigned an allele number
  ## 2) If no exact match - best match, gets "identity" < 1, and allele number is NA.
  
  if (is.null(temp_dir)) temp_dir=tempdir()
  ref=normalizePath(ref)

  ## Create output folder if not already present (this should really be passed in, right?)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings=FALSE)
  }

  output_table <- file.path(output_dir,paste0(basename(ref),".tab"))
  if (!file.exists(output_table)) {

    temp_dir = normalizePath(temp_dir)
    lst = list.files(path=isolate_dir, pattern=".fa$|.fas$|.fasta|.fna", full.names = TRUE)
    nm = length(lst)

    message("The data folder contains ", nm, " files to be queried.")
    message("Reading reference fasta file.")

    tryCatch(expr={fas=Biostrings::readDNAStringSet(ref)},
             error=function(e){stop("ERROR: Please provide a valid fasta-format reference file.")})
    message("The reference file contains ", length(fas), " sequence(s).")

    db_file <- paste0(dirname(ref),"/../db/",basename(ref),".nhr")
    if (!file.exists(db_file)) {
      makedb(ref=ref)
    }

    print("BLASTing all query sequences, and extracting the best hits")
    res=list()
    for (x in 1:nm){
      print(paste0("BLASTing sequence ",x))
      pt2=file.path(output_dir, "tmp")
      system(paste0("blastn -num_threads 4 -query ",lst[x]," -db ", db_file, " -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq sseq' -out ",pt2))
      if(file.exists(pt2)){
        if(file.info(pt2)$size>0){
          tab=read.table(pt2,header=F,stringsAsFactors = F)
          colnames(tab)=c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen","qseq", "sseq")
          if(dim(subset(tab,tab$pident==100 & tab$length==tab$slen))[1]>0){
            tab=subset(tab,tab$pident==100 & tab$length==tab$slen)
            res[[x]]= data.frame(file=basename(lst[x]),closest_allele=tab$saccver,identity=1,allele=as.numeric(gsub(".*_","",tab$saccver)),sequence=tab$qseq)
            tab %>% dplyr::group_by(qaccver) %>% dplyr::mutate(best_hit=saccver[order(evalue,-rank(bitscore))[1]],qseq=qseq[order(evalue,-rank(bitscore))[1]],score=evalue[order(evalue,-rank(bitscore))[1]],length=length[order(evalue,-rank(bitscore))[1]],file=basename(lst[x]))
          } else {
            res[[x]]= tab  %>% dplyr::summarise(file=basename(lst[x]),closest_allele=saccver[order(evalue,-rank(bitscore))[1]],identity=pident[order(evalue,-rank(bitscore))[1]]/100*length[order(evalue,-rank(bitscore))[1]]/slen[order(evalue,-rank(bitscore))[1]],allele=NA,sequence=qseq[order(evalue,-rank(bitscore))[1]])
          }
        }else{
          res[[x]]= data.frame(file=basename(lst[x]),closest_allele=NA,identity=NA,allele=NA,sequence=NA)
        }}else{
          res[[x]]= data.frame(file=basename(lst[x]),closest_allele=NA,identity=NA,allele=NA,sequence=NA)
        }
    }
    file.remove(pt2)
    res=do.call(rbind,res)
    write.table(res, output_table, row.names = FALSE)
    return(res)
  }
}

run_cgmlst=function(db_dir, isolate_dir, output_dir, temp_dir=NULL){
  ## Runs cgmlst for all reference alleles (a loop)
  
  lst=list.files(path=db_dir, pattern="*.fas", full.names = TRUE)
  #####
  #####
  #####
  ## YOU CAN PARALLELISE THIS, IT WILL BE MUCH FASTER
  ## BUT DON'T SAVE THE OUTPUT VARIABLE TO MEMORY
  ## IT WILL GET BIG QUICKLY
  ## RESULTS ARE WRITTEN TO DISK
  #####
  #####
  #####
  for(x in 1:length(lst)){
    cgmlst(ref = lst[x], isolate_dir=isolate_dir, output_dir=output_dir, temp_dir=temp_dir)
  }
}

fix_alleles=function(results_folder=NULL,threads=4){
  ## Goes through the output data, and checks for non-exact matches (allele = NA)
  ## If there are any, it aligns all the alleles using MAFFT
  ## You have to do this because sometimes the sequences are reverse complements
  ## It then adds a "u_" unfront of a new allele number and overwrites the original
  ## output data with the newly assigned allele code (leaving the ones that were already there)
  ## and it writes all of the sequences as an alignment, so that they are directly comparable.
  
  lst=list.files(path = results_folder,pattern="*.tab",full.names = T)
  tmp=tempfile()
  tmp2=tempfile()
  for (x in 1:length(lst)){
    print(paste0("Checking file ",x," of ",length(lst)))
    tab=read.table(lst[x],header=T,stringsAsFactors = F)
    if(any(is.na(tab$allele))){
      print("Re-aligning...")
      write.fasta(sequences=as.list(tab$sequence),
                  names=as.list(tab$file),
                  file.out=tmp)
      system(paste0("mafft --thread ",threads," --quiet --adjustdirectionaccurately ",tmp," > ",tmp2))
      dna=readDNAStringSet(tmp2)
      alleles=match(paste(dna),unique(paste(dna)))
      alleles=paste0("u_",alleles)
      tab$allele[is.na(tab$allele)]=alleles[is.na(tab$allele)]
      tab$sequence=paste(dna)
      write.table(tab,lst[x],row.names = F)
    }
  }
}

summarise_alleles=function(results_folder=NULL){
  
  ## Tabulates the assigned allele codes for all result tables in the output folder.
  ## Then returns the result AND writes it to "summary.tab" in the output folder.
  
  lst=list.files(path = results_folder,pattern="*.tab",full.names = T)
  df=lapply(1:length(lst),function(x) {
    tab=read.table(lst[x],header=T,stringsAsFactors = F)
    return(data.frame(file=tab$file,allele=gsub(".fas.tab","",basename(lst[x])),value=tab$allele,stringsAsFactors = F))
  })
  df=do.call(rbind,df)
  tab=matrix("0",nrow = length(unique(df$file)),
             ncol = length(unique(df$allele)),
             dimnames = list(unique(df$file),unique(df$allele)))
  for(x in 1:dim(df)[1]){
    tab[df$file[x],df$allele[x]]=df$value[x]
  }
  write.table(tab,paste0(results_folder,"/summary.tab"))
  return(tab)
}
