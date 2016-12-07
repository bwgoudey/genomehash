GenomeHash - A minimal-database MLST implementation

Benjamin Goudey^{1,2,5}, Hannah Huckstep^1, Kelly Wyres^{3,4}, and Thomas Conway^{1,5}
1. IBM Research - Australia, Carlton, Australia
2. Centre For Epidemiology and Biostatistics, University of
Melbourne, Parkville, Australia
3. Centre for Systems Genomics, University of Melbourne,
Parkville, Australia
4. Department of Biochemistry and Molecular Biology, Bio21
Molecular Science and Biotechnology Institute, University of
Melbourne, Parkville, Australia
5. Department of Computing and Information Systems,
University of Melbourne, Parkville, Australia

Contact: Benjamin Goudey: bgoudey@au1.ibm.com
         Tom Conway: tconway@au1.ibm.com

Requirements: 
  NCBI blast+ 2.2.28-2 (available as ncbi-blast+ package in Ubuntu)
  Numerous third party R packages including:
   -Biostrings
   -digest 
   -mclust
   -docopt 
   -stringdist
   -pbapply
   -tools
   -pacman
   -stringr

Usage:
  gh.R prereq [--install]
  gh.R ref [--first|--medoid] (--output <ref_db>) <allele_db>...
  gh.R extract (--refdb <ref_db>) (--output <alleles>) <contigs>
  gh.R hash <alleles> [--nchar <char>] (--output <hashes>) [--st] [--header]
  gh.R -h | --help

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
    

Examples: 
   Find the reference alleles for two loci for which we observed alleles using the 'first' method
      gh.R ref --first --output ./seeds_alleles.fa ./NEIS0001.txt ./NEIS0004.txt
   
   Extract alleles from assembled contigs using derived reference alleles
      gh.R extract --refdb ./seeds_alleles.fa --output ./10_6748.extract.fa ./10_6748.fas

   Create five character hashes for extracted alleles and write to file
      gh.R hash --nchar 5 --output ./10_6748_hashes.txt ./10_6748.extract.fa --header 

License: 
    MIT License, see license.txt for details
