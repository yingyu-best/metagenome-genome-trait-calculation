"
PROJECT: METAGENOMICS ANALYSIS ON PLASTISPHERE
CREATE DATE: 18 MARCH 2024
MODIFIED FROM JUSTIN LEE, CITY UNIVERSITY OF HONG KONG
DESCRIPTION: This is the second program script of entire data analysis, for assembly or MAG data codon-analysis.
*Please cite our paper (will provide the info later) if you found the script usefull in your researc :)
"

# (1) Import packages
##########################################
# install.packages('plotrix')
# install.packages('seqinr')
# install.packages('ggpubr')
#source("https://bioconductor.org/biocLite.R")
#biocLite('Biostrings')
# install.packages('data.table')
# install.packages('stringr')
# install.packages('ape')
# install.packages('tidyverse')
# install.packages('RColorBrewer')
# install.packages('gplots')
# install.packages("pheatmap)
library('plotrix')
library('seqinr')
library('ggpubr')
library('Biostrings')
library('data.table')
library('stringr')
library('ape')
library('tidyverse')
library('RColorBrewer')
library('gplots')
library("pheatmap")
##########################################


# (2) Functions definition
##########################################
get_HiExp_seqs <- function(my_seq = my_seq, annot_dir="emapper", strain_name = strain_name)
{# FUNCTION: to parse Ribosomal Protein CDSs from the genome based on the annotation file
  annot_dir <- sub(pattern = "/", replacement = "", annot_dir)  # Cut away the extra slash from input directory.
  strain_name <- gsub(pattern = ".fasta.ffn", replacement = "", my_seq)
  my_seq_fasta_temp <- readDNAStringSet(filepath = paste('ffn/', my_seq, sep = ''))
  emapper_ <- list.files(annot_dir)[grep(list.files(annot_dir), pattern = strain_name)]
  emapper_ <- paste(annot_dir, '/', emapper_, sep="")
  annotation_data <- read.delim(file = emapper_, sep = '\t', header = T, skip = 3 )
  query_names <- annotation_data$X.query_name[grep(pattern = 'ribosomal protein', x = annotation_data$eggNOG.free.text.desc.)]
  my_HiExp_fasta_temp <- DNAStringSet()
  if (str_detect(string = annotation_data[nrow(annotation_data), "X.query_name"], pattern = "Rate") == FALSE){
    print("FileError: Emapper annotation is not completed.")
    print(my_seq)
    return(my_HiExp_fasta_temp)
  }
  for (g in 1:length(query_names)){
    temp_seq <- my_seq_fasta_temp[grep(pattern = query_names[g], x = names(my_seq_fasta_temp))]
    my_HiExp_fasta_temp <- DNAStringSet(c(my_HiExp_fasta_temp, temp_seq))
  }
  writeXStringSet(my_HiExp_fasta_temp, filepath = paste(RESULT_dir, '/', strain_name, '_HiExp.fasta', sep = '')) 
  my_HiExp_fasta <- read.fasta(file = paste(RESULT_dir, '/', strain_name, '_HiExp.fasta', sep = ''), seqtype = 'DNA')
  return(my_HiExp_fasta)
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
##########################################


# (3) Environment Setting and Data Import
##########################################
# Envrionment settings
this.dir <- getwd()
date_ <- as.character(Sys.Date())
INTERMETIDATE_dir <- "./intermediate_RDS/"
dir.create(INTERMETIDATE_dir, showWarnings = FALSE)
RESULT_dir <- './ANALYSIS_RESULTS'
dir.create(RESULT_dir)

# Import useful preprocessed data.
source('FUNCTIONS.R')
bins_df <- readRDS(paste(INTERMETIDATE_dir, "Database_MAGs.rds", sep = ""))
seqs_df <- readRDS(paste(INTERMETIDATE_dir, "Database_CDSs.rds", sep = ""))
myassembly <- list.files(path = './ffn/') 
# Prepare empty RSCU matrices
rscu_genome_all_sample <- matrix(nrow = 64, ncol = length(myassembly))
rscu_HiExp_all_sample  <- matrix(nrow = 64, ncol = length(myassembly))
colnames(rscu_genome_all_sample) <- gsub(pattern = '\\.fasta.ffn', replacement = '', x = basename(myassembly))
colnames(rscu_HiExp_all_sample)  <- gsub(pattern = '\\.fasta.ffn', replacement = '', x = basename(myassembly))
##########################################


# (3) Main Program
##########################################
# Step 1: Codon Usage Analysis
completed_df_list <- list()  
stats_temp_list <- list()
codon_dataset_whole <- list()
n <- 0
pb <- txtProgressBar(min = 0, max = length(myassembly), char = "#", style = 3)
tic <- Sys.time()
print("STAGE 1:: Codon Usage Analysis ongoing...")
for(k in 1:length(myassembly)){
  my_seq <- myassembly[k]
  my_seq_fasta <- read.fasta(file = paste('./ffn/', my_seq, sep = ''), seqtype = 'DNA')
  strain_name <- gsub(pattern = '\\.ffn', replacement = '', x = my_seq)
  x_temp <- list.files(path = './emapper/', pattern = strain_name)
  my_HiExp_fasta <- if (length(x_temp) == 1) { get_HiExp_seqs(my_seq = my_seq, strain_name = strain_name) } else {my_seq_fasta}
  if (isEmpty(my_HiExp_fasta)){ # Exit program because of broken emapper file.
    exit()
  }
  # RSCU of whole genome
  eff_all <- sapply(my_seq_fasta, uco, index = "eff")
  dim(eff_all)[2] == length(my_seq_fasta) # Must be TRUE
  occurrences_genome_df <- data.frame(codon=rownames(eff_all), occurrences_genome = rowSums(eff_all))
  rownames(occurrences_genome_df) <- NULL
  codons_dataset <- merge(x = template, y = occurrences_genome_df, by = 'codon', sort = FALSE)
  codons_dataset$both_names <- paste(codons_dataset$aminoacid, '-', toupper(codons_dataset$codon),sep = '')
  codons_dataset$rscu_genome <- rscu_compute(codons_dataset$occurrences_genome)
  
  # RSCU of HiExp (Highly expressed proteins)
  eff_HiExp <- sapply(my_HiExp_fasta, uco, index = "eff")
  dim(eff_HiExp)[2] == length(my_HiExp_fasta) # Must be TRUE
  occurrences_HiExp_df <- data.frame(codon = rownames(eff_HiExp), occurrences_HiExp = rowSums(eff_HiExp))
  rownames(occurrences_HiExp_df) <- NULL
  codons_dataset$occurrences_HiExp <- occurrences_HiExp_df$occurrences_HiExp
  codons_dataset$rscu_HiExp <- rscu_compute(codons_dataset$occurrences_HiExp)
  codons_dataset$rscu_HiExp[is.nan(codons_dataset$rscu_HiExp)] <- 0
  
  # Plot out RSCU results
  low_RSCU <- rep(1.0, 64)
  postscript(file = paste(RESULT_dir, "/", strain_name, '_RSCU', '.eps', sep = ''), width = 15, height = 20)
  rscu_matrix_plot <- t(cbind(codons_dataset$rscu_genome, codons_dataset$rscu_HiExp))
  radial.plot(rscu_matrix_plot, labels = codons_dataset$both_names, label.prop=1.15,
              rp.type="s,p", 
              line.col=c(rgb(253,184,99, alpha = 255, maxColorValue = 255), rgb(128,115,172, alpha = 255, maxColorValue = 255)), 
              radial.lim=c(-1.5,max(rscu_matrix_plot)),
              show.grid.labels=T,
              lty=1, lwd=3, grid.col="grey90", 
              grid.bg="transparent", radlab=TRUE, 
              show.centroid=FALSE, mar = c(3, 3, 3, 3),
              main = '')
  radial.plot(low_RSCU, rp.type="l", line.col=c("red"), mar=c(3,3,3,3), 
              radial.lim=c(-1.5,max( c(codons_dataset$rscu_genome, codons_dataset$rscu_HiExp) )), 
              lwd=2, lty=2, add = TRUE)
  mtext(strain_name, side=3, adj = -0.05, line=1.2, cex=1.5, font=2)
  legend(x = 6, y = -5,  legend = c('Genome', 'HiExp'), 
         col = c( rgb(253,184,99, alpha = 255, maxColorValue = 255), rgb(128,115,172, alpha = 255, maxColorValue = 255) ), 
         lwd = 4, cex=0.9, bty = 'n')
  dev.off()
  
  # Attaching the RSCU columns of each sample to the "rscu_all_sample" matrix
  rscu_genome_all_sample[,k] <- codons_dataset$rscu_genome
  rscu_HiExp_all_sample[,k] <-  codons_dataset$rscu_HiExp
  
  # CAI genome
  # Use helper function `w_compute`
  codons_dataset$w_genome <- w_compute(codons_dataset$occurrences_genome)
  codons_dataset$w_HiExp <- w_compute(codons_dataset$occurrences_HiExp)
  codons_dataset$w_HiExp[is.nan(codons_dataset$w_HiExp)] <- 0
  write.table(x = codons_dataset, file = paste(RESULT_dir, '/', strain_name, '_RSCU_RESULTS.txt', sep = ''), 
              sep = '\t', row.names = FALSE)
  # reorder according to 'w' in seqinr
  data(caitab)
  w_order <- rownames(caitab)
  codons_dataset_ord <- codons_dataset[match(w_order, codons_dataset$codon),]
  # calculating CAI
  cai_all_genome <- sapply(my_seq_fasta, cai, w = codons_dataset_ord$w_genome)
  cai_all_HiExp <- sapply(my_seq_fasta, cai, w = codons_dataset_ord$w_HiExp)
  
  # Ploting CAI distribution (both genome and HiExp)
  pdf(file = paste(RESULT_dir, '/', strain_name, '_CAI', '.pdf', sep = ''), width = 10, height = 7)
  p1 <- hist(cai_all_genome, breaks = 10, plot = FALSE)
  p2 <- hist(cai_all_HiExp, breaks = 10, plot = FALSE)
  plot( p1, col = rgb(253,184,99, alpha = 150, maxColorValue = 255), xlim=c(0,1), 
        ylim=c(0, max( c(p1$counts, p2$counts) )), main = strain_name, xlab=" Codon Adaptation Index (CAI)")
  abline(v=median(cai_all_genome), col='red', lwd=3, lty=1)
  abline(v=mean(cai_all_genome), col='#4393c3', lwd=3, lty=1)
  abline(v=c(mean(cai_all_genome) + sd(cai_all_genome), mean(cai_all_genome) - sd(cai_all_genome)), col='#4393c3', lwd=2, lty=2)
  legend("topright", c(paste('Median: ', round(median(cai_all_genome), 3)), 
                       paste('Mean: ', round(mean(cai_all_genome), 3)),
                       paste('sd: ', round(sd(cai_all_genome), 3))), 
         col = c('red', '#4393c3', '#4393c3'),
         lwd = 1.5, lty = c(1,1,2), bty = "n")
  plot( p2, col= rgb(128,115,172, alpha = 150, maxColorValue = 255), xlim=c(0,1), ylim=c(0, max( c(p1$counts, p2$counts))), add=TRUE)
  abline(v=median(cai_all_HiExp), col='red', lwd=3, lty=1)
  abline(v=mean(cai_all_HiExp), col='#4393c3', lwd=3, lty=1)
  abline(v=c(mean(cai_all_HiExp) + sd(cai_all_HiExp), mean(cai_all_HiExp) - sd(cai_all_HiExp)), col='#4393c3', lwd=2, lty=2)
  legend("topleft", c(paste('Median: ', round(median(cai_all_HiExp), 3)), 
                      paste('Mean: ', round(mean(cai_all_HiExp), 3)),
                      paste('sd: ', round(sd(cai_all_HiExp), 3))), 
         col = c('red', '#4393c3', '#4393c3'),
         lwd = 1.5, lty = c(1,1,2), bty = "n")
  legend('left', legend = c('Genome', 'HiExp'), 
         col = c( rgb(253,184,99, alpha = 150, maxColorValue = 255), 
                  rgb(128,115,172, alpha = 150, maxColorValue = 255) ), lwd = 6, cex=1)
  dev.off()
  
  # GC vs GC3 & CAI
  gc   <- sapply(X = my_seq_fasta, FUN = GC)
  gc3  <- sapply(X = my_seq_fasta, FUN = GC3)
  L    <- sapply(X = my_seq_fasta, FUN = length)
  qGC  <- ecdf(gc)(gc) * 100
  qGC3 <- ecdf(gc3)(gc3) * 100
  qCAI_genome <- ecdf(cai_all_genome)(cai_all_genome) * 100
  qCAI_HiExp  <- ecdf(cai_all_HiExp)(cai_all_HiExp) * 100
  gc_gc3_cai_df <- cbind.data.frame(BIN = strain_name,
                                    CDS = names(my_seq_fasta),
                                    a = sapply(X = my_seq_fasta, FUN = get_nucleo_comp, nucleotide='a'),
                                    t = sapply(X = my_seq_fasta, FUN = get_nucleo_comp, nucleotide='t'),
                                    g = sapply(X = my_seq_fasta, FUN = get_nucleo_comp, nucleotide='g'),
                                    c = sapply(X = my_seq_fasta, FUN = get_nucleo_comp, nucleotide='c'),
                                    Length = L,
                                    GC = gc,
                                    qGC = qGC,
                                    GC3 = gc3, 
                                    qGC3 = qGC3,
                                    CAI_genome = cai_all_genome,
                                    qCAI_genome = qCAI_genome,
                                    CAI_HiExp = cai_all_HiExp,
                                    qCAI_HiExp = qCAI_HiExp)
  
  if (length(x_temp) == 1) {
    file_name <- paste('emapper/', list.files('emapper/')[grep(list.files('emapper/'), pattern = strain_name)], sep = '')
    LARGE_gc_gc3_cai_df <- merge(x = gc_gc3_cai_df, y = read.delim(file = file_name, sep = '\t', header = T, skip = 3 ),
                                 by.x= 'CDS', by.y='X.query_name', all.x = TRUE, sort = FALSE)
  } else {
    LARGE_gc_gc3_cai_df <- gc_gc3_cai_df
  }
  write.table(x = LARGE_gc_gc3_cai_df, file = paste(RESULT_dir, '/', strain_name, '_DATAFRAME.txt', sep = ''), 
              sep = '\t', row.names = FALSE)
  completed_df_list[[k]] <- gc_gc3_cai_df
  
  # Plot correlations
  postscript(file = paste(RESULT_dir, '/', strain_name, '_CORR.eps', sep = ''), width = 12, height = 20)
  par(mar=c(4, 5, 4, 4), mfrow=c(2,3), oma = c(1, 1, 1, 1))
  # Subplot 1: GC vs CAI_genome
  gcReg2 <- lm(CAI_genome ~ GC, data = gc_gc3_cai_df)
  cor2 <- cor.test(x = gc_gc3_cai_df$GC, y=gc_gc3_cai_df$CAI_genome, method = 'pearson')
  plot(x = gc_gc3_cai_df$GC, y=gc_gc3_cai_df$CAI_genome, main=strain_name, xlim = c(0,1), ylim = c(0,1), 
       pch=20, cex=0.6, xlab='GC', ylab='CAI_genome', col='#a6bddb', cex.lab=1.5, cex.axis=1.5)
  abline(reg = gcReg2, lty=1)
  abline(coef = c(0,1), lty=2)
  abline(v=c(mean(gc), mean(gc)+sd(gc), mean(gc)-sd(gc)), lty=c(3,3,3), lwd=1)
  abline(h=c(mean(gc_gc3_cai_df$CAI_genome), mean(gc_gc3_cai_df$CAI_genome)+sd(gc_gc3_cai_df$CAI_genome), 
             mean(gc_gc3_cai_df$CAI_genome)-sd(gc_gc3_cai_df$CAI_genome)), lty=c(3,3,3), lwd=1)
  mtext(text = paste(sep = '', 'PCC  ', round(cor2$estimate, 2), 
                     ',  ', ifelse(test = cor2$p.value<0.001, yes = "P < 0.001", no = paste("P = ", round(cor2$p.value, 2)))),
        side = 3, line=0, cex=0.8)
  box(lwd=1)
  # Subplot 2: GC3 vs CAI_genome
  gcReg3 <- lm(CAI_genome ~ GC3, data = gc_gc3_cai_df)
  cor3 <- cor.test(x = gc_gc3_cai_df$GC3, y=gc_gc3_cai_df$CAI_genome, method = 'pearson')
  plot(x = gc_gc3_cai_df$GC3, y=gc_gc3_cai_df$CAI_genome, main=strain_name, xlim = c(0,1), 
       ylim = c(0,1), pch=20, cex=0.6, xlab='GC3', ylab='CAI_genome', col='#a6bddb', cex.lab=1.5, cex.axis=1.5)
  abline(reg = gcReg3, lty=1)
  abline(coef = c(0,1), lty=2)
  abline(v=c(mean(gc3), mean(gc3)+sd(gc3), mean(gc3)-sd(gc3)), lty=c(3,3,3), lwd=1)
  abline(h=c(mean(gc_gc3_cai_df$CAI_genome), mean(gc_gc3_cai_df$CAI_genome)+sd(gc_gc3_cai_df$CAI_genome), 
             mean(gc_gc3_cai_df$CAI_genome)-sd(gc_gc3_cai_df$CAI_genome)), lty=c(3,3,3), lwd=1)
  mtext(text = paste(sep = '', 'PCC  ', round(cor3$estimate, 2), ',  ', 
                     ifelse(test = cor3$p.value<0.001, yes = "P < 0.001", no = paste("P = ", round(cor3$p.value, 2)))),
        side = 3, line=0, cex=0.8)
  box(lwd=1)
  # Subplot 3: GC vs GC3
  gcReg1 <- lm(GC3 ~ GC, data = gc_gc3_cai_df)
  cor1 <- cor.test(x = gc_gc3_cai_df$GC, y=gc_gc3_cai_df$GC3, method = 'pearson')
  plot(x = gc_gc3_cai_df$GC, y=gc_gc3_cai_df$GC3, main=strain_name, xlim = c(0,1), ylim = c(0,1), 
       pch=20, cex=0.6, xlab='GC', ylab='GC3', col='#c51b8a', cex.lab=1.5, cex.axis=1.5)
  abline(reg = gcReg1, lty=1)
  abline(coef = c(0,1), lty=2)
  abline(v=c(mean(gc), mean(gc)+sd(gc), mean(gc)-sd(gc)), lty=c(3,3,3), lwd=1)
  abline(h=c(mean(gc3), mean(gc3)+sd(gc3), mean(gc3)-sd(gc3)), lty=c(3,3,3), lwd=1)
  mtext(text = paste(sep = '', 'PCC  ', round(cor1$estimate, 2), ',  ', 
                     ifelse(test = cor1$p.value<0.001, yes = "P < 0.001", no = paste("P = ", round(cor1$p.value, 2)))),
        side = 3, line=0, cex=0.8)
  box(lwd=1)
  # Subplot 4: GC vs CAI_HiExp
  gcReg4 <- lm(CAI_HiExp ~ GC, data = gc_gc3_cai_df)
  cor4 <- cor.test(x = gc_gc3_cai_df$GC, y=gc_gc3_cai_df$CAI_HiExp, method = 'pearson')
  plot(x = gc_gc3_cai_df$GC, y=gc_gc3_cai_df$CAI_HiExp, main=strain_name, xlim = c(0,1), ylim = c(0,1), 
       pch=20, cex=0.6, xlab='GC', ylab='CAI_HiExp', col='#a1d99b', cex.lab=1.5, cex.axis=1.5)
  abline(reg = gcReg4, lty=1)
  abline(coef = c(0,1), lty=2)
  abline(v=c(mean(gc), mean(gc)+sd(gc), mean(gc)-sd(gc)), lty=c(3,3,3), lwd=1)
  abline(h=c(mean(gc_gc3_cai_df$CAI_HiExp), mean(gc_gc3_cai_df$CAI_HiExp)+sd(gc_gc3_cai_df$CAI_HiExp), 
             mean(gc_gc3_cai_df$CAI_HiExp)-sd(gc_gc3_cai_df$CAI_HiExp)), lty=c(3,3,3), lwd=1)
  mtext(text = paste(sep = '', 'PCC  ', round(cor4$estimate, 2), ',  ', 
                     ifelse(test = cor4$p.value<0.001, yes = "P < 0.001", no = paste("P = ", round(cor4$p.value, 2)))),
        side = 3, line=0, cex=0.8)
  box(lwd=1)
  # Subplot 5: GC3 vs CAI_HiExp
  gcReg5 <- lm(CAI_HiExp ~ GC3, data = gc_gc3_cai_df)
  cor5 <- cor.test(x = gc_gc3_cai_df$GC3, y=gc_gc3_cai_df$CAI_HiExp, method = 'pearson')
  plot(x = gc_gc3_cai_df$GC3, y=gc_gc3_cai_df$CAI_HiExp, main=strain_name, xlim = c(0,1), ylim = c(0,1), 
       pch=20, cex=0.6, xlab='GC3', ylab='CAI_HiExp', col='#a1d99b', cex.lab=1.5, cex.axis=1.5)
  abline(reg = gcReg5, lty=1)
  abline(coef = c(0,1), lty=2)
  abline(v=c(mean(gc3), mean(gc3)+sd(gc3), mean(gc3)-sd(gc3)), lty=c(3,3,3), lwd=1)
  abline(h=c(mean(gc_gc3_cai_df$CAI_HiExp), mean(gc_gc3_cai_df$CAI_HiExp)+sd(gc_gc3_cai_df$CAI_HiExp), 
             mean(gc_gc3_cai_df$CAI_HiExp)-sd(gc_gc3_cai_df$CAI_HiExp)), lty=c(3,3,3), lwd=1)
  mtext(text = paste(sep = '', 'PCC  ', round(cor5$estimate, 2), ',  ', ifelse(test = cor5$p.value<0.001, 
                                                                               yes = "P < 0.001", no = paste("P = ", round(cor5$p.value, 2)))), side = 3, line=0, cex=0.8)
  box(lwd=1)
  # Subplot 6: TO BE FIX - INSERT PLOT OF SLOPE AS 6th PLOT
  dev.off()
  
  # Saving STATS
  stats_temp <- data.frame(sample = strain_name, 
                           mean_GC         = mean(gc_gc3_cai_df$GC),         sd_GC         = sd(gc_gc3_cai_df$GC), 
                           mean_GC3        = mean(gc_gc3_cai_df$GC3),        sd_GC3        = sd(gc_gc3_cai_df$GC3), 
                           mean_CAI_genome = mean(gc_gc3_cai_df$CAI_genome), sd_CAI_genome = sd(gc_gc3_cai_df$CAI_genome),
                           mean_CAI_HiExp  = mean(gc_gc3_cai_df$CAI_HiExp),  sd_CAI_HiExp  = sd(gc_gc3_cai_df$CAI_HiExp),
                           
                           intercept_GC_GC3         = gcReg1$coefficients[1], slope_GC_GC3         = gcReg1$coefficients[2],
                           intercept_GC_CAI_genome  = gcReg2$coefficients[1], slope_GC_CAI_genome  = gcReg2$coefficients[2],
                           intercept_GC3_CAI_genome = gcReg3$coefficients[1], slope_GC3_CAI_genome = gcReg3$coefficients[2],
                           intercept_GC_CAI_HiExp   = gcReg4$coefficients[1], slope_GC_CAI_HiExp   = gcReg4$coefficients[2],
                           intercept_GC3_CAI_HiExp  = gcReg5$coefficients[1], slope_GC3_CAI_HiExp  = gcReg5$coefficients[2],
                           
                           GC_GC3_PCC         = cor1$estimate, GC_GC3_p_value         = cor1$p.value,
                           GC_CAI_genome_PCC  = cor2$estimate, GC_CAI_genome_p_value  = cor2$p.value,
                           GC3_CAI_genome_PCC = cor3$estimate, GC3_CAI_genome_p_value = cor3$p.value,
                           GC_CAI_HiExp_PCC   = cor4$estimate, GC_CAI_HiExp_p_value   = cor4$p.value,
                           GC3_CAI_HiExp_PCC  = cor5$estimate, GC3_CAI_HiExp_p_value  = cor5$p.value)
  
  stats_temp_list[[k]] <- stats_temp
  codons_dataset$sample <- strain_name
  codon_dataset_whole[[k]] <- codons_dataset
  setTxtProgressBar(pb, k) 
}
toc <- Sys.time()
print(paste("Time Elapsed: ", str(toc - tic), sep = ""))
print("======STAGE 1 completed!======")

# Step 2:Save resulting matrices
print("STAGE 2:: Dataframes saving...")
completed_df <- rbindlist(completed_df_list)
write.table(x = completed_df, file = paste(RESULT_dir, '/01_METRICS.txt', sep=""), sep = '\t', row.names = FALSE)
stats_df <- rbindlist(stats_temp_list) 
write.table(x = stats_df, file = paste(RESULT_dir, '/02_STATS.txt', sep=""), sep = '\t', row.names = FALSE)
codon_dataset_whole <- rbindlist(codon_dataset_whole)
write.table(x = codon_dataset_whole, file = paste(RESULT_dir, '/03_CODONS.txt', sep=""), sep = '\t', row.names = FALSE)
print("======STAGE 2 completed!======")

# Step 3: RSCU analysis visualization
print("STAGE 3:: Visualizations")
if (length(myassembly) > 2){
  rownames(bins_df) <- paste(bins_df$sample_id, "-", bins_df$Sample_type, sep = "")
  # RSCU genome matrix of all samples
  row.names(rscu_genome_all_sample) <- codons_dataset$both_names
  write.table(x = rscu_genome_all_sample, file = paste(RESULT_dir, '/rscu_genome_all_sample.txt', sep=""), quote = TRUE, sep = '\t', row.names = TRUE)
  rscu_genome_all_sample <- as.matrix(read.delim(file = paste(RESULT_dir, '/rscu_genome_all_sample.txt', sep=""), sep = '\t', header = T))
  # colnames(rscu_genome_all_sample) <- paste(colnames(rscu_genome_all_sample), "-", as.data.frame(taxa_df)[match(colnames(rscu_genome_all_sample), taxa_df$sample_id), "phyla"], sep = "")
  my_palette <- colorRampPalette(c("#ca0020", "#f7f7f7", "#0571b0"))(n = 29)
  col_breaks = c(seq(0,0.60,length=10),
                 seq(0.61,1.6,length=10),
                 seq(1.61,4,length=10))
  rbPal <- colorRampPalette(c('red','blue')) # Dark red=low-GC; dark blue=high-GC
  temp_codons <- str_sub(rownames(rscu_genome_all_sample), start = 2, end=5)
  codons_metrics <- data.frame(codons=rownames(rscu_genome_all_sample), 
                               oxy=str_count(temp_codons, "C")*2 + str_count(temp_codons, "G")*2 + str_count(temp_codons, "A") + str_count(temp_codons, "T") * 3)
  codons_metrics$color <- as.character(rbPal(10)[as.numeric(cut(as.numeric(codons_metrics$oxy),breaks = 10))])
  annotation_row = data.frame(OxygenInCodon = as.numeric(codons_metrics$oxy))
  rownames(annotation_row) <- codons_metrics$codons 
  annotation_col = data.frame(GC = bins_df$gc_mean, CodonAverage.Oxygen = bins_df$Oxygen, Type = bins_df$Sample_type) 
  rownames(annotation_col) = colnames(rscu_genome_all_sample)
  
  # Plot heatmap for genome data
  pdf(file = paste(RESULT_dir, '/03_HEATMAP_RSCU_genome.pdf', sep=""), width = 16, height = 8)
  pheatmap(rscu_genome_all_sample,
           fontsize_row = 6, 
           col= my_palette, 
           breaks = col_breaks, 
           annotation_row = annotation_row, 
           annotation_col = annotation_col, 
           main = expression("RSCU"["Genome"]),
           fontsize_col = 4)
  dev.off()
  
  # RSCU HiExp matrix of all samples
  rownames(rscu_HiExp_all_sample) <- codons_dataset$both_names
  write.table(x = rscu_HiExp_all_sample, file = paste(RESULT_dir, '/rscu_HiExp_all_sample.txt', sep=""), quote = TRUE, sep = '\t', row.names = TRUE)
  rscu_HiExp_all_sample <- as.matrix(read.delim(file = paste(RESULT_dir, '/rscu_HiExp_all_sample.txt', sep=""), sep = '\t', header = T))
  # colnames(rscu_HiExp_all_sample) <- paste(colnames(rscu_HiExp_all_sample), "-", as.data.frame(taxa_df)[match(colnames(rscu_HiExp_all_sample), taxa_df$sample_id), "phyla"], sep = "")
  
  # Plot heatmap for HiExp data
  pdf(file = paste(RESULT_dir, '/03_HEATMAP_RSCU_HiExp.pdf', sep=""), width = 16, height = 8)
  pheatmap(rscu_HiExp_all_sample,
           fontsize_row = 6, 
           col= my_palette, 
           breaks = col_breaks, 
           annotation_row = annotation_row, 
           annotation_col = annotation_col, 
           main = expression("RSCU"["HiExp"]),
           fontsize_col = 4)
  dev.off()
  
  # Clustering RSCU genome
  pdf(file = paste(RESULT_dir, '04_CLUSTERING_RSCU_GENOME.pdf', sep=""), width = 12, height = 12)
  clust_rscu_genome_all_sample <- t(rscu_genome_all_sample[c(-30,-59),]) # Removing M-ATG and  W-TGG codons
  d <- dist(clust_rscu_genome_all_sample)
  fit <- hclust(d, method = "ward.D2")
  saveRDS(fit, file = paste(INTERMETIDATE_dir, "04_CLUSTERING_RSCU_GENOME.rds", sep = ""))
  plot(fit, cex=0.3) # display dendogram
  dev.off()
  
  # Clustering RSCU HiExp
  pdf(file = paste(RESULT_dir, '/04_CLUSTERING_RSCU_HiExp.pdf', sep=""), width = 12, height = 12)
  clust_rscu_HiExp_all_sample <- t(rscu_HiExp_all_sample[c(-30,-59),]) # Removing M-ATG and  W-TGG codons
  d2 <- dist(clust_rscu_HiExp_all_sample)
  fit2 <- hclust(d2, method = "ward.D2")
  saveRDS(fit2, file = paste(INTERMETIDATE_dir, "04_CLUSTERING_RSCU_HiExp.rds", sep = ""))
  plot(fit2, cex=0.3) # display dendogram
  dev.off()
}
print("======STAGE 3 completed!======")
print('Analysis FINISHED!')
##########################################

