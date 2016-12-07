#!/usr/bin/env Rscript 

#
# GenomeHash 
#  Reference implementation from ideas in 
#  Beyond MLST — Database-less Microbial Genotyping"
#
# 15/6/2016
# 

#helper function to test file existance and throw errors if fail
checkFileExist <-function(fname, var_name, dir.only=F) {
  x = ifelse(dir.only, "path", "file")
  cond = ifelse(dir.only, !file.exists(dirname(fname)), !file.exists(fname))
  fname = ifelse(dir.only, dirname(fname), fname)
  if(cond) {
      cat(sprintf("Given %s for %s:%s does not exist\n", x, var_name, fname))
      quit()
    }
}

#Get ref alleles from given allele database
extractRefAlleles <-function(allele_dbs, outfile, medoid=FALSE) {
  cat(sprintf("Determining reference alleles\n"))
  l<-pblapply(allele_dbs, function(allele_db) {
    checkFileExist(allele_db, "allele database")
  
    #Read in the allele_db
    alleles=readDNAStringSet(allele_db)
    
    #Get medoid allele or get first allele
    if(medoid) {
      cat(sprintf("Calculating distances between alleles (may take\n a while if number of alleles is large)\n"))
      #medoid= allele w/ min. ave. levenschtein dist to all other alleles in the db
      a=alleles[which.min(apply(as.matrix(stringdistmatrix(alleles)), 1, mean))]
    } else
      a=alleles[1]
    
    names(a)<-paste(file_path_sans_ext(basename(allele_db)), names(a), collapse="_", sep="_")
    return(a)
  })
  
  #Write out reference alleles
  checkFileExist(outfile, "output directory", dir.only = T)
  writeXStringSet(x=do.call('c', l), filepath=outfile)
}


#Read in custom BLAST blast output and make
#all returned alleles same direction
readBlastFile <- function(blast_file) {
  #Get best hits for each seed allele
  x<-read.table(blast_file, as.is=T)
  colnames(x)<-c("qseqid", "sseqid", "bitscore", "sseq", "sstrand", "length", "slen", "qstart", "qend", "qlen");
  #filter based on min lengthcontigs
  x <- x[x$length/x$slen > 0.5,]
  
  #flip all reverse strands
  I = which(x$sstrand=="minus")
  x$sseq[I] <- sapply(I, function(i) toString(reverseComplement(DNAString(x$sseq[i]))) )
  x$sstrand[I] <- "plus"  

  return(x)
}


# setup alignment function
align <- function(pattern, contig) {
  a1 = pairwiseAlignment(pattern, contig, type="global-local", gapOpening=10, gapExtension=5)
  a2 = pairwiseAlignment(pattern, reverseComplement(contig), type="global-local", gapOpening=10, gapExtension=5)
  if(a1@score>a2@score) 
    return(a1) 
  else 
    return(a2)
}


#Extraction
extractAllele <- function(assembly_file , ref_db, outfile) {
  cat(sprintf("Extracting alleles\n"))
  #Check that blast is in path
  if(any(sapply(c("makeblastdb", "blastn"), function(x) Sys.which(x)=="")))  {
    cat(sprintf("BLAST commands (makeblastdb, blastn) not found. Is it installed?\n"))
    cat(sprintf("Check for package 'ncbi-blast+' or go to ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST\n"))
    quit()
  }
  checkFileExist(ref_db, "reference alleles db")
  checkFileExist(assembly_file, "assemblies")
  checkFileExist(outfile, "allele output file", dir.only = T)

  #Create ref db blast database if it doesn't already exist 
  if(!file.exists(sprintf("%s.nhr", basename(ref_db)))) {
    cat(sprintf("Building BLAST database over given reference allele database\n"))
    cmd=sprintf("makeblastdb -in %s  -dbtype nucl -title %s -out %s", ref_db, basename(ref_db), basename(ref_db))
    system(cmd)
  }
  
  #Call blast over the given allele
  cat(sprintf("Use BLAST to detect alleles from reference alleles\n"))
  blast_file = paste(file_path_sans_ext(ref_db), "blast", sep=".")
  outfmt = "6 qseqid sseqid bitscore sseq sstrand length slen qstart qend qlen"
  system(sprintf("blastn -db %s -query %s -outfmt '%s' -word_size 15 -perc_identity 70  | sort -k 3 -g > %s", 
                 basename(ref_db), assembly_file,outfmt,  blast_file))

  checkFileExist(blast_file, "expected BLAST output")

  #read in the blast results and find the strongest matches (max bitscore)
  blast_res<-readBlastFile(blast_file)
  seeds <- sort(unique(blast_res$sseqid))
  best<-lapply(seeds, function(seed)  {
    subset(blast_res, sseqid==seed)[order(subset(blast_res, sseqid==seed)$bitscore, decreasing = T),][1,]
  })
  
  # Because we don't trust blast, conduct an exact alignment over whole contig ----
  cat(sprintf("Use Exact match to refine BLAST match\n"))
  contigs=readDNAStringSet(assembly_file)
  alleles<-lapply(best, function(b) {
    allele = toString(align(b$sseq[[1]], contigs[[b$qseqid]])@subject)
    return(sprintf(">%s %s\n%s\n", b$sseqid, paste(b[!(colnames(b) %in% c("sseqid", "sseq"))], collapse=" "), allele))
  })
  
  #Write alleles to a file
  write(unlist(alleles), file = outfile)
  return()
}


#Helper function for hash single allele
hashAlleles <- function(allele_db, outfile, chars=5, st=T, header=T) {
  cat(sprintf("Creating hashes of given alleles\n"))

  checkFileExist(allele_db, "Extracted alleles")
  checkFileExist(outfile, "allele output file", dir.only = T)

  alleles=readDNAStringSet(allele_db)
  hashes = unlist(lapply(alleles, hashAllele, chars=chars))
  header_txt= c("file",  sub(" .*", "", names(alleles)))

  if(st){
    st = hashAllele(paste(hashes, collapse=""), chars)
    hashes = c(hashes, st)
    header_txt= c(header_txt, "st")
  }

  f = file(outfile, 'w')
  if(header) {
      writeLines(paste(header_txt, collapse="\t"), con = f)
    }
  writeLines(paste(c(sub("\\..*", "", basename(allele_db)), hashes), collapse="\t"), con = f)
  close(f)
  return()
  
}


# Given an allele and an optional size of the resulting label (in characters), 
# produce as hash-based string.
hashAllele <- function(allele, chars=5) {
  if(class(allele)=="DNAString")
    allele=toString(allele)
    
  bitsPerChar=5
  bits=chars*bitsPerChar
  min_bytes=ceiling(bits/8)
  
  #create 32-bit alphabet. Remove characters that could be misread
  alph32 = setdiff(sapply(c(48:57,97:(97+25)), function(i) rawToChar(as.raw(i))), c("l","i","o","0"))
  
  packHashBits <- function(x) sum(2^(which(as.logical(x))-1))
  
  #Create a hash.
  #Shift so that we use the correct number of bits
  #Convert hash bits into alphabet 
  hash = substr(digest(allele, algo="md5"),1,min_bytes);
  pad = rawToBits(rep(as.raw(0),5))[1:(min_bytes*8-bits)]
  hash_bits = c(rawToBits(charToRaw(hash))[-c(1:(min_bytes*8-bits))], pad)
  hash_chars = alph32[sapply(seq(1,bits,bitsPerChar), 
                             function(i) packHashBits(hash_bits[i:(i+bitsPerChar-1)])
  )]
  return(paste(hash_chars, collapse=""))
}


#### Driver logic----
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.name <- ifelse(length(script.name), basename(script.name), "gh.R")

prereqs=c("Biostrings", "digest", "mclust", "docopt", "stringdist", "pbapply", "tools", "pacman", "stringr")

#Install packages if we run the "prereq" command
I = grep("--args", initial.options)
if(length(I) && initial.options[I+1]=="prereq"){
  if(length(grep("--install", initial.options)))
  {
      cat(sprintf("Install pre-requisite packages if required\n"))
      source("http://bioconductor.org/biocLite.R")
      if (!("pacman" %in% rownames(installed.packages()))) {cat("installing pacman"); biocLite("pacman"); library(pacman)}
      pacman::p_load_gh("hadley/stringr")
      pacman::p_load("docopt", "Biostrings", "digest", "mclust", "docopt", "stringdist", "pbapply", "tools")
      cat(sprintf("Finished install pre-requisite libraries\n"))
  } else {
      cat(sprintf("R packages required to run this script:\n%s\n", paste(prereqs, collapse="\n")))
      cat(sprintf("Run '%s prereq --install' to install any missing dependencies\n", script.name))
  }
  quit()
} 

#Otherwise assume libraries installed

#Setup command line parsing
doc<-sprintf("
GenomeHash - A minimal-database MLST implementation

Usage:
  %s prereq [--install]
  %s ref [--first|--medoid] (--output <ref_db>) <allele_db>...
  %s extract (--refdb <ref_db>) (--output <alleles>) <contigs>
  %s hash <alleles> [--nchar <char>] (--output <hashes>) [--st] [--header]
  %s -h | --help

Commands: 
  prereq    Install R packages required to run this script. NB: does not install BLAST.
  ref       Extract a set of 'reference alleles' from a given allele database
            Input 'allele_db' and output 'ref_db' are FASTA files.
  extract   Extract alleles by aligning to given seedds to a set of assemblies using BLAST. 
            Input 'contigs' and 'ref_db' and output 'alleles' are FASTA files.
  hash      Hash a given FASTA file of alleles.

Options:
  --help  -h          Show this screen.
  --install           List packages to be installed
  --output=<outfile>  Output filename
  --first             Use first allele as seed.
  --medoid            Use medoid allele as seed.
  --refdb=<ref_db>    Specify the database of reference alleles to use.
  --nchar=<char>      Number of characters in hash [default: 5].
  --header            Write out a header containing column (i.e. allele) names
  --st                Add a sequence type hash

Requirements: 
  NCBI blast+ 2.2.28-2 (available as ncbi-blast+ package in Ubuntu)

Examples: 
   Find the reference alleles for two loci for which we observed alleles using the 'first' method
      gh.R ref --first --output ./seeds_alleles.fa ./NEIS0001.txt ./NEIS0004.txt
   
   Extract alleles from assembled contigs using derived reference alleles
      gh.R extract --refdb ./seeds_alleles.fa --output ./10_6748.extract.fa ./10_6748.fas

   Create five character hashes for extracted alleles and write to file
      gh.R hash --nchar 5 --output ./10_6748_hashes.txt ./10_6748.extract.fa --header


", script.name, script.name, script.name, script.name, script.name)


if (!("pacman" %in% rownames(installed.packages()))){
  cat(doc)
  stop("Some prerequisite R packages not installed - run 'gh prereq'", call. = FALSE)
}

#Load just enough to parse command lines - faster if command line error
pacman::p_load_gh("hadley/stringr", install = F)
pacman::p_load("docopt")

#Try to parse command line args
if(!length(grep("--args", initial.options))) { docopt(doc, "--help" ); quit(); }

opts <- tryCatch(docopt(doc), error = function(e) { cat(sprintf("%sError parsing command line args\n", doc)); quit()})

if(opts$`--help`)
  quit()

#Load in rest of libraries after command args parse
#source("http://bioconductor.org/biocLite.R")
pacman::p_load("Biostrings", "digest", "mclust", "docopt", "stringdist", "pbapply", "tools", install=F)

#Create reference alleles
if(opts$ref) {
   zz<-extractRefAlleles(allele_dbs=opts$allele_db, outfile=opts$output, medoid=opts$medoid)
} else if(opts$extract) {
  zz<-extractAllele(assembly_file=opts$contigs, ref_db=opts$refdb, outfile=opts$output)
} else if(opts$hash) {
  zz<-hashAlleles(allele_db=opts$alleles, outfile=opts$output, chars=as.numeric(opts$nchar), st=opts$st, header=opts$header)
}
cat(sprintf("DONE\n"))

