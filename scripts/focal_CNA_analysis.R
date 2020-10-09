## Load abs copy number data and perform GISTIC CNA enrichment analysis
args = commandArgs(trailingOnly=TRUE)

# Load libraries
library(Biobase)
library(QDNAseqmod)
library(tidyverse)

# GISTIC parameters
sample_list <- as.character(readLines(snakemake@input[["list"]]))
seg_file_name <- snakemake@wildcards[["subset"]]
gistic_folder <- snakemake@config[["gistic"]] 
output_folder <- paste0(snakemake@config[["output_dir"]],seg_file_name,"_gistic_results/")
refgenefile <- paste0(gistic_folder,"refgenefiles/hg19.mat")

# Load abs data
rds <- snakemake@config[["rds"]]
abs_data <- readRDS(rds)

# Load meta data
meta <- read.table(file = snakemake@config[["meta"]],header = T,sep = "\t")

# Generate log2-like relative segments
# GISTIC defines the CN.seg value as (6) Seg.CN   =   (log2() -1 of copy number)
convert_to_rel <- function(segs = x,ploidy = n){
  segs[segs <= 0] <- 0.001
  segs <- log2(segs)-1
  # segs <- log2(segs)
  return(segs)
}

# Extract relative segment calls from QDNAseq segment data
extract_rel_segs <- function(abs_data = x){
  sample_n <- length(colnames(abs_data))
  to_use <- fData(abs_data)$use
  all_segs <- data.frame()
  for(i in 1:sample_n){
    cn_obj <- abs_data[to_use,i]
    sample_name <- colnames(cn_obj)
    segments <- assayDataElement(cn_obj,"segmented")
    segments <- convert_to_rel(segments,median(segments))
    rle_x <- rle(as.numeric(segments))
    end = cumsum(rle_x$lengths)
    start = c(1, lag(end)[-1] + 1)
    chr_pos_s <- fData(cn_obj)[start,c("chromosome")]
    chr_pos_e <- fData(cn_obj)[end,c("chromosome")]
    if(any(chr_pos_s != chr_pos_e)){
      stop("Inconsistent chromosome labels for segment")
    } else {
      chr_pos <- fData(cn_obj)[start,c("chromosome")]
    }
    start_pos <- fData(cn_obj)[start,c("start")]
    end_pos <- fData(cn_obj)[end,c("end")]
    seg_vals <- rle_x$values
    n_probes <- end - start +1
    segs <- data.frame(Sample=sample_name,
                       Chromosome=chr_pos,
                       Start=start_pos,
                       End=end_pos,
                       nProbes=n_probes,
                       Seg.CN=seg_vals)
    all_segs <- rbind(all_segs,segs)
  }
  return(all_segs)
}

# Subset metafile
meta <- meta[meta$SAMPLE_ID %in% sample_list,]

# Extract segments across all samples
abs_data <- abs_data[,colnames(abs_data) %in% meta$SAMPLE_ID]
seg_calls <- extract_rel_segs(abs_data = abs_data)

# Write segs to table
write.table(x = seg_calls,
            file = paste0(output_folder,seg_file_name,".txt"),
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

# call script that sets MCR environment and calls GISTIC executable - example gistic params
#gistic.def.cmd <- paste0(gistic_folder,"/gistic2 ",
#                        " -b ",output_folder,
#                        " -seg ",paste0(output_folder,seg_file_name,".txt"),
#                        " -refgene ",refgenefile,
#                        " -fname default",
#                        " -armpeel 1",
#                        " -brlen 0.5",
#                        " -savegene 1",
#                        " -genegistic 1",
#                        " -smallmem 1",
#                        " -broad 1",
#                        " -conf 0.90",
#                        " -gcm extreme",
#                        " -ta 0.1",
#                        " -cap 1.5",
#                        " -td 0.1",
#                        " -js 4",
#                        " -maxseg 2500",
#                        " -qvt 0.25",
#                        " -rx 1")

# call script that sets MCR environment and calls GISTIC executable - TCGA params

gistic.tcga.cmd <- paste0(gistic_folder,"gistic2",
                          " -b ",output_folder,
                          " -seg ",paste0(output_folder,seg_file_name,".txt"),
                          " -refgene ",refgenefile,
                          " -fname ",seg_file_name,
                          " -armpeel 1",
                          " -brlen 0.7",
                          " -savegene 1",
                          " -genegistic 1",
                          " -smallmem 1",
                          " -broad 1",
                          " -conf 0.99",
                          " -gcm extreme",
                          " -ta 0.1",
                          " -cap 1.5",
                          " -td 0.1",
                          " -js 4",
                          " -maxseg 2000",
                          " -qvt 0.25",
                          " -rx 0")

# Run gistic
#system(command = gistic.def.cmd)
system(command = gistic.tcga.cmd)
