## --------------------------------------------
##
## Script name: .R
##
## Description: 
##
##
##
## Output: data for input in signature classification
##
## Date created: 04-2023
##
## Email: alcheng@nygenome.org
##
##
## Resources:
## --------------------------------------------

## ---------- Set working directory -----------



## ---------- Set environment -----------------

library(ggplot2); library(scales); library(pROC); #library(dplyr)
library(data.table); library(stringr); library(zoo); library(lsa);
library(deconstructSigs); library(BSgenome.Hsapiens.UCSC.hg38)
load("resources/signatures.genome.cosmic.v3.2.march2021.grch38.rda")
sig_ref = signatures.genome.cosmic.v3.2.march2021.grch38
bsg_ref = BSgenome.Hsapiens.UCSC.hg38
source("resources/theme_bw2.R")
## ----------- Functions ---------------------------------------------------------

get_denoised_duplex <- function(filename){
  x <- fread(filename)

  print(x)
  x$X_FLAGS= as.numeric(x$X_FLAGS)
  x$X_EDIST = as.numeric(x$X_EDIST)
  x = x[x$X_EDIST<=2, ]
  print(x)

  MI_top = x$MI[x$X_FLAGS == 0 | x$X_FLAGS==1024]
  MI_bot = x$MI[x$X_FLAGS == 16 | x$X_FLAGS==1040]
  print(MI_top)
  print(MI_bot)
  KEEP = MI_top[MI_top %in% MI_bot]
  print(KEEP)
  x = x[x$MI %in% KEEP, ]
  print(x)

  x <- x[!duplicated(x$MI_pos), ]
  df = x
  df <- df[df$duplex_bp2 !='N', ]
  df <- df[df$alt == df$duplex_bp2, ]
  
  df$pos <- as.numeric(df$pos)
  df$X_EDIST <- as.numeric(df$X_EDIST)
  df$X_FC1 <- as.numeric(df$X_FC1)
  df$X_FC2 <- as.numeric(df$X_FC2)
  df$X_FLAGS = as.numeric(df$X_FLAGS)
  df$X_FILTERED_COUNT <- as.numeric(df$X_FILTERED_COUNT)
  df$X_LENGTH <- as.numeric(df$X_LENGTH)
  df$X_MAPQ <- as.numeric(df$X_MAPQ)
  df$X_READ_COUNT<- as.numeric(df$X_READ_COUNT)
  df$X_SCORE <- as.numeric(df$X_SCORE)
  df$rq <- as.numeric(df$rq)
  df$rs <- as.numeric(df$rs) # Start position of the original read
  df$aD <- as.numeric(df$aD) # n molecules that create top-strand consensus
  df$aE <- as.numeric(df$aE) # error rate of top strand
  df$bD <- as.numeric(df$bD) # n molecules that creat bot-strand consensus
  df$bE <- as.numeric(df$bE) # error rate of bot strand
  df$cD <- as.numeric(df$cD) # depth of duplex consensus (~= aD + bD)
  df$cE <- as.numeric(df$cE) # error rate of duplex
  df$cU <- as.numeric(df$cU) # total number of reads that have the same duplex UMI
                                # this is always equal or greater than cD. It can be greater
                                # if some duplicates have different cigar strands. Reads with mismatched cigars
                                # are not used in duplex creation
  df$dS <- as.numeric(df$dS) # duplex start position
  df$FM_passed <- as.numeric(df$FM_passed) # number of duplicates that make it to the featuremap file
  df$num_N <- as.numeric(df$num_N) # number of Ns on the duplex read
  df$duplex_length <- as.numeric(df$duplex_length) #len of duplex
  df$duplex_PIR <- as.numeric(df$duplex_PIR) #position in read of duplex
  df$duplex_PIR2 <- as.numeric(df$duplex_PIR2) #position in read from 3' end (length - pir)
  df$DP <- as.numeric(df$DP) # depth (unfiltered) at that position
  df$AF <- as.numeric(df$AF) # unfiltered allele fraction
  df$SB <- as.numeric(df$SB) #unfiltered strand bioas
  df$ref_top <- as.numeric(df$ref_top) # number of top-mapping bases of REF
  df$ref_bot <- as.numeric(df$ref_bot) # number of bot-mapping bases of REF
  df$alt_top <- as.numeric(df$alt_top) # number of top-mapping bases of ALT
  df$alt_bot <- as.numeric(df$alt_bot) # number of bot-mapping bases of ALT
  df$pir <- as.numeric(df$pir)
  df$pir2<- as.numeric(df$pir2)
  
  df = df[df$alt_bot>=1 & df$alt_top >=1, ] # logically this should not filter
                                            # anything, but we do remove a few weirdly mapped reads this way
  df = df[df$FM_passed==df$cU, ] #only consider bases where all duplex-duplicates make it to the featuremap
  
  df <- df[df$pir >=25, ] # quality filter
  df <- df[df$pir2 >=25, ] # quality filter

  return(df)
}

get_triN_frequencies <- function(d){
  
  
  d <- d[, c("chrom", "pos","ref", "alt", 'sample_id')]
  
  d$sample_id <- 'sample_id'
  sigs.input <- mut.to.sigs.input(mut.ref = d, 
                                  sample.id = "sample_id", 
                                  chr = "chrom", 
                                  pos = "pos", 
                                  ref = "ref", 
                                  alt = "alt", bsg = BSgenome.Hsapiens.UCSC.hg38)
  v= whichSignatures(tumor.ref = sigs.input, 
                     sample.id = "sample_id",
                     contexts.needed = TRUE,
                     signature.cutoff = 0)
  return(v)
}

get_sig <- function(p, sample_id){
  df = get_denoised_duplex(p)
  df$sample_id <- sample_id
  
  df = df[df$AF<=0.2, ]  

  triN <- get_triN_frequencies(df)
  triN <- triN$tumor
  
  triN <- data.frame(t(triN))
  
  colnames(triN)<- c("frequency")
  triN$sample_id <- p
  triN$total_count <- nrow(df)
  triN$count <- round(triN$total_count)*triN$frequency
  
  triN$triN <- paste0(str_sub(rownames(triN), 1,1),
                      str_sub(rownames(triN), 3,3),
                      str_sub(rownames(triN), 7,7))
  
  triN$variant <- paste0(str_sub(rownames(triN), 3,5))
  triN <- triN[, c("sample_id", "triN", "variant", "frequency", "total_count", "count")]
  return(triN)
}
args = commandArgs(trailingOnly=TRUE)

p = args[1]
sample_id = args[2]
output=args[3]
output2 = args[4]
result <- get_sig(p, sample_id)

fwrite(x=result, file = output, sep = '\t')

result$triN_var = paste0(str_sub(result$triN,1,1), "[", result$variant, "]", str_sub(result$triN, 3,3))

triN_order = fread("resources/triN_order.R", header=FALSE)

result$triN_var = factor(result$triN_var, levels = triN_order$V1)

pdf(output2, height = 2, width = 5)
ggplot(data = result)+
  geom_col(aes(x=triN_var, y=frequency, fill = variant))+
    scale_fill_manual(values = c("#45BEEE", "#000000",
                                 "#E42D26", "#CDC9CA", 
                                 "#A2CD61", "#EDC6C5"))+
    theme_bw2()+
    xlab("Trinucleotide context")+
    ylab("Frequency")+
    theme(axis.text.x= element_blank())
dev.off()
