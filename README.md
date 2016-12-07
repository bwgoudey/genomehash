# GenomeHash - A minimal-database MLST implementation

* Benjamin Goudey
* Hannah Huckstep
* Kelly Wyres
* Thomas Conway

## About: 

GenomeHash is a simple method for using alleleic genotyping schemes (MLST, MLVA, Phage typing etc)
without the need for a reference database of alleles. 
This is achieved by replacing traditional numeric labels with hashes and by determining alleles for a given locus
using only a minimal allele database. 


## Contact: 

Benjamin Goudey: bgoudey@au1.ibm.com


## Requirements: 


 * NCBI blast+ 2.2.28-2 (available as ncbi-blast+ package in Ubuntu) 
*  R plus the following packages
   * Biostrings 
   * digest
   * mclust
   * docopt
   * stringdist
   * pbapply
   * tools
   * pacman
   * stringr

These packages will be installed by running the pre-req command

## Usage: 

* gh.R prereq [--install] gh.R ref --first|--medoid 
* gh.R extract (--refdb ) (--output ) 
* gh.R hash --nchar [--st] [--header] 
* gh.R -h | --help

**prereq:** Install R packages required to run this script. NB: does not install BLAST. 

**ref** Extract a set of 'reference alleles' from a given allele database Input 'allele_db' and output 'ref_db' are FASTA files. 

**extract** Extract alleles by aligning to given seedds to a set of assemblies using BLAST. Input 'contigs' and 'ref_db' and output 'alleles' are FASTA files. 

**hash** Hash a given FASTA file of alleles.

**Options:**

--help -h Show these instructions.

--install List packages to be installed

--output= Output filename 

--first Use first allele as seed.

--medoid Use medoid allele as seed.

--refdb= Specify the database of reference alleles to use.

--nchar= Number of characters in hash [default: 5].

--header Write out a header containing column (i.e. allele) names

--st Add a sequence type hash

## Examples: 

Find the reference alleles for two loci for which we observed alleles using the 'first' method

<code>gh.R ref --first --output ./seeds_alleles.fa ./NEIS0001.txt ./NEIS0004.txt</code>

Extract alleles from assembled contigs using derived reference alleles

<code>gh.R extract --refdb ./seeds_alleles.fa --output ./10_6748.extract.fa ./10_6748.fas</code>

Create five character hashes for extracted alleles and write to file 

<code>gh.R hash --nchar 5 --output ./10_6748_hashes.txt ./10_6748.extract.fa --header </code>

## License: 
MIT License, see license.txt for details
