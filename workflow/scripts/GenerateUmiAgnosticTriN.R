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
library(data.table); library(parallel); library(stringr); library(zoo); library(lsa)
library(scales); library(deconstructSigs); library(BSgenome.Hsapiens.UCSC.hg38)
load("/gpfs/commons/home/aarora/UG_analysis/Itai_Error_Analysis/tri_nuc_scripts/signatures.genome.cosmic.v3.2.march2021.grch38.rda")
sig_ref = signatures.genome.cosmic.v3.2.march2021.grch38
bsg_ref = BSgenome.Hsapiens.UCSC.hg38
source('resources/theme_bw2.R')

## ----------- Functions ---------------------------------------------------------

get_denoised_umiagno <- function(filename){
  
  x <- fread(filename)
  df = x

  df$pos <- as.numeric(df$pos)
  df$X_EDIST <- as.numeric(df$X_EDIST)
  df$X_FC1 <- as.numeric(df$X_FC1)
  df$X_FC2 <- as.numeric(df$X_FC2)
  df$X_FLAGS<-as.numeric(df$X_FLAGS)
  
  df$X_FILTERED_COUNT <- as.numeric(df$X_FILTERED_COUNT)
  df$X_LENGTH <- as.numeric(df$X_LENGTH)
  df$X_MAPQ <- as.numeric(df$X_MAPQ)
  df$X_READ_COUNT<- as.numeric(df$X_READ_COUNT)
  df$X_SCORE <- as.numeric(df$X_SCORE)
  df$rq <- as.numeric(df$rq)
  df$rs <- as.numeric(df$rs)
  df$X_LENGTH <- as.numeric(df$X_LENGTH)
  df$PIR <- as.numeric(df$PIR)
  df$PIR2 <- as.numeric(df$PIR2)
  df$DP <- as.numeric(df$DP)
  df$AF <- as.numeric(df$AF)
  df$SB <- as.numeric(df$SB)
  df$ref_top <- as.numeric(df$ref_top)
  df$ref_bot <- as.numeric(df$ref_bot)
  df$alt_top <- as.numeric(df$alt_top)
  df$alt_bot <- as.numeric(df$alt_bot)
  
  df = df[(df$alt_bot>=1 & df$X_FLAGS == 16) | (df$alt_top>=1 & df$X_FLAGS == 0) | (df$alt_bot>=1 & df$X_FLAGS == 1040) | (df$alt_top>=1 & df$X_FLAGS == 1024), ] # logically this should not filter
                                            # anything, but we do remove a few weirdly mapped reads this way
  
  df <- df[df$PIR >=25, ] # quality filter
  df <- df[df$PIR2 >=25, ] # quality filter

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
  df = get_denoised_umiagno(p)
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


