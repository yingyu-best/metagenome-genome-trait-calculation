"
PROJECT: METAGENOMICS ANALYSIS ON PLASTISPHERE
CREATE DATE: 18 MARCH 2024
MODIFIED FROM JUSTIN LEE, CITY UNIVERSITY OF HONG KONG
DESCRIPTION: This is the first program script of entire data analysis, for assembly or MAG data preprocessing.
*Please cite our paper (will provide the info later) if you found the script usefull in your researc :)
"

# Step 0: Set working environment#
# (1) Import packages
##########################################
# source("https://bioconductor.org/biocLite.R"); biocLite("Biostrings")
# install.packages("seqinr")
# install.packages("reshape")
# install.packages("tidyverse")
# install.packages("scales")
# install.packages('data.table')
library("Biostrings")
library("seqinr")
library("tidyverse")
library("scales")
library('data.table')
##########################################

# (2) Functions definition
##########################################
A_in_Codon <- function (seq, forceToLower = TRUE, exact = FALSE, NA.A = NA, oldGC = FALSE) 
{ # FUNCTION:: a function to return base-A ration in a coding sequence. Copied from and modified from seqinr::GC().
  if (length(seq) == 1 && is.na(seq)) 
    return(NA)
  if (nchar(seq[1]) > 1) 
    stop("sequence is not a vector of chars")
  if (forceToLower) 
    seq <- tolower(seq)
  nc <- sum(seq == "c")
  ng <- sum(seq == "g")
  na <- sum(seq == "a")
  nt <- sum(seq == "t")
  if (oldGC) {
    warning("argument oldGC is deprecated")
    return((nc + ng)/length(seq))
  }
  if (!exact) {
    if (na + nc + ng + nt == 0) {
      result <- NA.A
    }
    else {
      result <- (na)/(na + nc + ng + nt)
    }
  }
  else {
    ngc <- ng + nc
    nat <- na + nt
    ngc <- ngc + sum(seq == "s")
    nat <- nat + sum(seq == "w")
    if (na + nc != 0) {
      nm <- sum(seq == "m")
      ngc <- ngc + nm * nc/(na + nc)
      nat <- nat + nm * na/(na + nc)
    }
    if (ng + nt != 0) {
      nk <- sum(seq == "k")
      ngc <- ngc + nk * ng/(ng + nt)
      nat <- nat + nk * nt/(ng + nt)
    }
    if (ng + na != 0) {
      nr <- sum(seq == "r")
      ngc <- ngc + nr * ng/(ng + na)
      nat <- nat + nr * na/(ng + na)
    }
    if (nc + nt != 0) {
      ny <- sum(seq == "y")
      ngc <- ngc + ny * nc/(nc + nt)
      nat <- nat + ny * nt/(nc + nt)
    }
    if (na + nc + ng != 0) {
      nv <- sum(seq == "v")
      ngc <- ngc + nv * (nc + ng)/(na + nc + ng)
      nat <- nat + nv * na/(na + nc + ng)
    }
    if (na + nc + nt != 0) {
      nh <- sum(seq == "h")
      ngc <- ngc + nh * nc/(na + nc + nt)
      nat <- nat + nh * (na + nt)/(na + nc + nt)
    }
    if (na + ng + nt != 0) {
      nd <- sum(seq == "d")
      ngc <- ngc + nd * ng/(na + ng + nt)
      nat <- nat + nd * (na + nt)/(na + ng + nt)
    }
    if (nc + ng + nt != 0) {
      nb <- sum(seq == "b")
      ngc <- ngc + nb * (nc + ng)/(nc + ng + nt)
      nat <- nat + nb * nt/(nc + ng + nt)
    }
    if (ngc + nat == 0) {
      result <- NA.A
    }
    else {
      result <- na/(ngc + nat)
    }
  }
  return(result)
}

A3 <- function (seq, frame = 0, ...) 
{# FUNCTION:: a function to return base-A ratio in 3rd-position of a codon
  if (nchar(seq[1]) > 1) {
    warning("sequence is not a vector of chars, I'm trying to cast it into one")
    seq <- s2c(seq[1])
  }
  if (frame != 0) 
    seq <- seq[(1 + frame):length(seq)]
  A_in_Codon(seq[seq(3, length(seq), by = 3)], ...)
}

get_nucleo_comp <- function(seq, nucleotide='a')
{# FUNCTION: to count te number of specific nucleotide in a sequence `seq`.
  if(!nucleotide %in% c("a", "t", "g", "c", "u")){
    print("`nucleotide` should be either on of a/t/g/c/u.")
    return()
  }
  seqinr::count(seq = seq, 1)[nucleotide]
}

get_atoms_comp <- function(seq_df, atom="O")
{# FUNCTION: to calculate codon-average number of atoms (O/H/C/N) of each CDS in a sequence dataframe
  switch (atom,
          O = ((seq_df$a * 0) + (seq_df$g * 1) + (seq_df$c * 1) + (seq_df$t * 2))/(seq_df$Length/3),
          N = ((seq_df$a * 5) + (seq_df$g * 5) + (seq_df$c * 3) + (seq_df$t * 2))/(seq_df$Length/3),
          H = ((seq_df$a * 5) + (seq_df$g * 5) + (seq_df$c * 5) + (seq_df$t * 4))/(seq_df$Length/3),
          C = ((seq_df$a * 5) + (seq_df$g * 5) + (seq_df$c * 4) + (seq_df$t * 5))/(seq_df$Length/3),
          'Please give correct atom(N/O/C/H)')
}

codon_counting_eff <- function(seq)
{# FUNCTION: to calculate effective number of codons in sequence `seq`.
  t(uco(seq, index = "eff"))
} 


isCDS <- function(seq)
{# FUNCTION: to check if the input sequence is a CDS (i.e. [length %% 3 == 0] AND [start with ATG] AND [stop with TAG/TAA/TGA])
  if (length(seq) %% 3 != 0){
    return(FALSE)
  }
  if (paste(seq[1:3], collapse="") != "atg") {
    return(FALSE)
  }
  if (!paste(seq[(length(seq)-2):length(seq)], collapse = "") %in% c("tag", "taa", "tga")){
    return(FALSE)
  }
  return(TRUE)
}
##########################################


# (3) Main Program
##########################################
# Step 0: Set working environment
setwd("/Users/baoyy/Desktop/test")
this.dir <- getwd()
INTERMETIDATE_dir <- "./intermediate_RDS/"
dir.create(INTERMETIDATE_dir, showWarnings = FALSE)
FFN <- "./ffn"   
FAA <- "./faa"   

# Step 1: Create a dataframe summerizing all .ffn.
print('STAGE 1:: Dataframe is being processed...')
seqs_df <- data.frame()
ffn_files <- list.files(FFN, pattern = 'ffn')#ffn
faa_files <- list.files(FAA, pattern = 'faa')
pb <- txtProgressBar(min = 0, max = length(ffn_files), char = "#", style = 3)
k <- 0               # progressbar counter
initializer <- TRUE  # dataframe-initializing indicator
for (ffn in ffn_files){
  k <- k + 1
  str.location <- gregexpr(pattern ='\\.', ffn)
  index <- substr(ffn, 1, str.location[[1]][2]-1)
  fasta_seqinr <- read.fasta(file = file.path(FFN, ffn))
  # extract data directly from .ffn file.
  new_df <- data.frame(seq_label = names(fasta_seqinr),
                       gc = sapply(X = fasta_seqinr, FUN = GC),
                       gc3 = sapply(X = fasta_seqinr, FUN = GC3),
                       a3 = sapply(X = fasta_seqinr, FUN = A3),
                       a = sapply(X = fasta_seqinr, FUN = get_nucleo_comp, nucleotide='a'),
                       g = sapply(X = fasta_seqinr, FUN = get_nucleo_comp, nucleotide='g'),
                       t = sapply(X = fasta_seqinr, FUN = get_nucleo_comp, nucleotide='t'),
                       c = sapply(X = fasta_seqinr, FUN = get_nucleo_comp, nucleotide='c'),
                       length = sapply(X = fasta_seqinr, FUN = length))
  new_df$Annotation <- sapply(getAnnot(fasta_seqinr), substr, 17, length(getAnnot(fasta_seqinr)))
  new_df$Carbon   <- ((new_df$a * 5) + (new_df$g * 5) + (new_df$c * 4) + (new_df$t * 5))/(new_df$length/3)   # U has 4 carbons and T has 5
  new_df$Hydrogen <- ((new_df$a * 5) + (new_df$g * 5) + (new_df$c * 5) + (new_df$t * 4))/(new_df$length/3)   # U has 4 Hydrogens and T has 6
  new_df$Nitrogen <- ((new_df$a * 5) + (new_df$g * 5) + (new_df$c * 3) + (new_df$t * 2))/(new_df$length/3)   # U and T have the same Nitrogen content = 2
  new_df$Oxygen   <- ((new_df$a * 0) + (new_df$g * 1) + (new_df$c * 1) + (new_df$t * 2))/(new_df$length/3)   # U and T have the same Oxygen content = 2
  new_df$Adenine  <- (new_df$a)/(new_df$length/3)
  new_df$Thymine  <- (new_df$t)/(new_df$length/3)
  new_df$Guanine  <- (new_df$g)/(new_df$length/3)
  new_df$Cytosine <- (new_df$c)/(new_df$length/3)
  new_df$gc3_to_gc <- new_df$gc3 / new_df$gc 
  new_df$sequence <- sapply(fasta_seqinr, getSequence)
  new_df$isCDS    <- sapply(new_df$sequence, isCDS) #ORFs must be started with ATG and ended with STOP-codon, as well as length%3!=0; CDSs is subset of ORFs.
  new_df$sample_id    <- index
  if (initializer){
    seqs_df <<- new_df
    initializer <- FALSE
  }else {
    seqs_df <- rbind(seqs_df, new_df)
  }
  setTxtProgressBar(pb, k)
}
seqs_df$Row.names <- NULL
print("======STAGE 1 completed!======")

# Step 2: Generate assembly-level or MAG-level dataframe (averaging by all ORFs from one sample's assembly or a MAG).
print("STAGE 2:: Computing assembly-level or MAG-level data.frames...")
seqs_df <- seqs_df %>% filter(isCDS == TRUE)
bins_df <- seqs_df %>% 
  select(sample_id, Oxygen, Nitrogen, Hydrogen, Carbon, 
         Adenine, Thymine, Guanine, Cytosine, length, gc, gc3, gc3_to_gc) %>% 
  group_by(sample_id) %>% 
  summarise(Sum_length=sum(length), gc_mean=mean(gc), gc3_mean=mean(gc3), GC_avg_ratio=mean(gc3_to_gc), cds_count=n(),
            Adenine=mean(Adenine), Thymine=mean(Thymine), Guanine=mean(Guanine), Cytosine=mean(Cytosine),
            Oxygen=mean(Oxygen), Nitrogen=mean(Nitrogen), Hydrogen=mean(Hydrogen), Carbon=mean(Carbon))

print("Saving all databases...")
saveRDS(seqs_df, file = paste(INTERMETIDATE_dir, "Database_CDSs.rds", sep=""))
saveRDS(bins_df, file = paste(INTERMETIDATE_dir, "Database_MAGs.rds", sep=""))
print("======STAGE 3 completed!======")
print('Data Preprocessing FINISHED!')
##########################################