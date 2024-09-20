"
PROJECT: METAGENOMICS ANALYSIS ON PLASTISPHERE
CREATE DATE: 18 MARCH 2024
MODIFIED FROM JUSTIN LEE, CITY UNIVERSITY OF HONG KONG
DESCRIPTION: This is a helper Rscript of the second program script of entire data analysis, for assembly or MAG data codon-analysis.
*Please cite our paper (will provide the info later) if you found the script usefull in your researc :)
"
# (1) Import packages
###########################
library('ape')
library('phytools')
library('dendextend')
#library('phylogram')
###########################


# (2) Define functions
###########################
# MATRIX: Element composition table of nucleotides (Justin)
nucleotides_element <- t(matrix(c(0, 1, 1, 2, 5, 5, 3, 2, 5, 5, 5, 4, 5, 5, 4, 5), nrow = 4, ncol = 4, byrow = TRUE, 
                                dimnames = list(c('Oxygen', 'Nitrogen', 'Hydrogen', 'Carbon'), c('a', 'g', 'c', 't'))))

# DATAFRAME: a handy template to calculate RSCU and CAI
template <- data.frame(codon = c('gca', 'gcc', 'gcg', 'gct',                # A [4]
                                 'tgc', 'tgt',                              # C [2]
                                 'gac', 'gat',                              # D [2]
                                 'gaa', 'gag',                              # E [2]
                                 'ttc', 'ttt',                              # F [2] 
                                 'gga', 'ggc', 'ggg', 'ggt',                # G [4]
                                 'cac', 'cat',                              # H [2] 
                                 'ata', 'atc', 'att',                       # I [3] 
                                 'aaa', 'aag',                              # K [2] 
                                 'cta', 'ctc', 'ctg', 'ctt', 'tta', 'ttg',  # L [6]
                                 'atg',                                     # M [1]
                                 'aac', 'aat',                              # N [2]
                                 'cca', 'ccc', 'ccg', 'cct',                # P [4] 
                                 'caa', 'cag',                              # Q [2] 
                                 'aga', 'agg', 'cga', 'cgc', 'cgg', 'cgt',  # R [6]
                                 'agc', 'agt', 'tca', 'tcc', 'tcg', 'tct',  # S [6]
                                 'aca', 'acc', 'acg', 'act',                # T [4]  
                                 'gta', 'gtc', 'gtg', 'gtt',                # V [4] 
                                 'tgg',                                     # W [1]
                                 'tac', 'tat',                              # Y [2] 
                                 'taa', 'tag', 'tga'),                      # STOP [3]
                       aminoacid = c('A', 'A', 'A', 'A', 
                                     'C', 'C', 
                                     'D', 'D', 
                                     'E', 'E', 
                                     'F', 'F', 
                                     'G', 'G', 'G', 'G', 
                                     'H', 'H', 
                                     'I', 'I', 'I', 
                                     'K', 'K', 
                                     'L', 'L', 'L', 'L', 'L', 'L', 
                                     'M', 
                                     'N', 'N', 
                                     'P', 'P', 'P', 'P', 
                                     'Q', 'Q', 
                                     'R', 'R', 'R', 'R', 'R', 'R', 
                                     'S', 'S', 'S', 'S', 'S', 'S', 
                                     'T', 'T', 'T', 'T', 
                                     'V', 'V', 'V', 'V', 
                                     'W', 
                                     'Y', 'Y', 
                                     '*', '*', '*'))

rscu_compute <- function(occurrences) 
{# FUNTION: to compute the Relative Synonymous Codon Usage - RSCU
  N <- length(occurrences)
  vaa <- c(4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 1, 2, 4, 2, 6, 6, 4, 4, 1, 2, 3)
  rscu <- vector()
  j <- 1
  k <- 1
  for (i in 1:N){
    rscu <- c(rscu, occurrences[i]/(sum(occurrences[j:(j+vaa[k]-1)])/vaa[k]))
    if (i/(j+vaa[k]-1) == 1){
      j <- i+1
      k <- k+1
    }
  } 
  return(rscu)
}

w_compute <- function (occurrences) 
{# FUNCTION: to compute the `w` value required to calculate the Codon Adaptation Index - CAI
  N <- length(occurrences)
  vaa <- c(4, 2, 2, 2, 2, 4, 2, 3, 2, 6, 1, 2, 4, 2, 6, 6, 4, 4, 1, 2, 3)
  w <- vector()
  j <- 1
  k <- 1
  for  (i in 1:N){
    w <- c(w, occurrences[i]/(max (occurrences[j:(j+vaa[k]-1)])))
    if (i/(j+vaa[k]-1)==1){
      j <- i+1
      k <- k+1
    }
  } 
  return(w)
}



dendroPlot_phyloRSCU <- function(nwk_file, rscu_values_txt, title_nwk, title_rscu)
{# FUNCTION: to plot a face-to-face phylo-RSCU dendrogram.
  rscu_values <- as.matrix(read.delim(file = rscu_values_txt, sep = '\t', header = T))
  clust_rscu <- t(rscu_values[c(-30,-59),]) # Removing M-ATG and  W-TGG codons
  d <- dist(clust_rscu)
  fit <- hclust(d, method = "ward.D2")
  RSCU <- as.dendrogram(fit)
  tree <- read.dendrogram(nwk_file)  #replace choosing root for conversion below.
  val = round(cor(cophenetic(RSCU), cophenetic(tree)), 6)    #covariance/correlation of two cophenetic distances.
  tanglegram(tree, RSCU, sort = FALSE, common_subtrees_color_lines = TRUE, 
             lwd = 2, dLeaf_left = 0, dLeaf_right = 0, margin_inner = 12.5, 
             highlight_distinct_edges  = FALSE, highlight_branches_lwd = FALSE, lab.cex = 1, 
             main_left = title_nwk, main_right = title_rscu, columns_width = c(4,1,4))
  mtext(paste("Cophenetic Correlation = ", val, sep=""), side=3, font=2, cex=2)
}
###########################

