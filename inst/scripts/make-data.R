### =========================================================================
### denyranges is an AnnotationHub package that stores genomic coordinates of
### problematic genomic regions as GRanges.
### -------------------------------------------------------------------------
###

# denyranges package currently contains 19 Rds objects. The original and 
# processed data are available at 
# https://drive.google.com/drive/folders/124DZtsU0YVWqkb7dgu8Nk6b3N8-ShVSC?usp=sharing

# The object names are structured as "<genome assembly>.<lab>.<original file name>", 
# e.g., "hg19.Birney.wgEncodeDacMapabilityConsensusExcludable".

# # ENCODE
# ## hg19
# 
# Mint_Blacklist_hg19.bed.gz
# https://www.encodeproject.org/files/ENCFF200UUD/, Bradley Bernstein, Broad
# wget https://www.encodeproject.org/files/ENCFF200UUD/@@download/ENCFF200UUD.bed.gz
# gzip --decompress --stdout ENCFF200UUD.bed.gz | bedtools sort -i - > hg19.Bernstein.Mint_Blacklist_hg19.bed
# 
# wgEncodeDacMapabilityConsensusExcludable.bed.gz, Ewan Birney, EBI
# https://www.encodeproject.org/files/ENCFF001TDO/, Ewan Birney, EBI
# wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz
# gzip --decompress --stdout ENCFF001TDO.bed.gz | bedtools sort -i - > hg19.Birney.wgEncodeDacMapabilityConsensusExcludable.bed
# 
# wgEncodeDukeMapabilityRegionsExcludable.bed.gz
# https://www.encodeproject.org/files/ENCFF001THR/, Gregory Crawford, Duke
# wget https://www.encodeproject.org/files/ENCFF001THR/@@download/ENCFF001THR.bed.gz
# gzip --decompress --stdout ENCFF001THR.bed.gz | bedtools sort -i - > hg19.Crawford.wgEncodeDukeMapabilityRegionsExcludable.bed
# 
# hg19mitoblack.bed.gz
# https://www.encodeproject.org/files/ENCFF055QTV/, Barbara Wold, Caltech
# wget https://www.encodeproject.org/files/ENCFF055QTV/@@download/ENCFF055QTV.bed.gz
# gzip --decompress --stdout ENCFF055QTV.bed.gz | bedtools sort -i - > hg19.Wold.hg19mitoblack.bed
# 
# eCLIP_blacklistregions.hg19.bed.gz
# https://www.encodeproject.org/files/ENCFF039QTN/, Gene Yeo, UCSD
# wget https://www.encodeproject.org/files/ENCFF039QTN/@@download/ENCFF039QTN.bed.gz
# gzip --decompress --stdout ENCFF039QTN.bed.gz | bedtools sort -i - > hg19.Yeo.eCLIP_blacklistregions.hg19.bed
# 
# ## hg38
# 
# Mint_Blacklist_GRCh38.bed.gz
# https://www.encodeproject.org/files/ENCFF023CZC/, Bradley Bernstein, Broad
# wget https://www.encodeproject.org/files/ENCFF023CZC/@@download/ENCFF023CZC.bed.gz
# gzip --decompress --stdout ENCFF023CZC.bed.gz | bedtools sort -i - > hg38.Bernstein.Mint_Blacklist_GRCh38.bed
# 
# GRCh38_unified_blacklist.bed.gz
# https://www.encodeproject.org/files/ENCFF356LFX/, Anshul Kundaje, Stanford
# wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
# gzip --decompress --stdout ENCFF356LFX.bed.gz | bedtools sort -i - > hg38.Kundaje.GRCh38_unified_blacklist.bed
# 
# GRCh38.blacklist.bed.gz
# https://www.encodeproject.org/files/ENCFF419RSJ/, Anshul Kundaje, Stanford
# wget https://www.encodeproject.org/files/ENCFF419RSJ/@@download/ENCFF419RSJ.bed.gz
# gzip --decompress --stdout ENCFF419RSJ.bed.gz | bedtools sort -i - > hg38.Kundaje.GRCh38.blacklist.bed
# 
# wgEncodeDacMapabilityConsensusExcludable.hg38.bed.gz
# https://www.encodeproject.org/files/ENCFF220FIN/, Tim Reddy, Duke
# wget https://www.encodeproject.org/files/ENCFF220FIN/@@download/ENCFF220FIN.bed.gz
# gzip --decompress --stdout ENCFF220FIN.bed.gz | bedtools sort -i - > hg38.Reddy.wgEncodeDacMapabilityConsensusExcludable.hg38.bed
# 
# hg38mitoblack.bed.gz
# https://www.encodeproject.org/files/ENCFF940NTE/, Barbara Wold, Caltech
# wget https://www.encodeproject.org/files/ENCFF940NTE/@@download/ENCFF940NTE.bed.gz
# gzip --decompress --stdout ENCFF940NTE.bed.gz | bedtools sort -i - > hg38.Wold.hg38mitoblack.bed
# 
# eCLIP_blacklistregions.hg38liftover.bed.fixed.bed.gz
# https://www.encodeproject.org/files/ENCFF269URO/, Gene Yeo, UCSD
# wget https://www.encodeproject.org/files/ENCFF269URO/@@download/ENCFF269URO.bed.gz
# gzip --decompress --stdout ENCFF269URO.bed.gz | bedtools sort -i - > hg38.Yeo.eCLIP_blacklistregions.hg38liftover.bed.fixed.bed
# 
# ## mm10
# 
# blacklist.full.bed.gz
# https://www.encodeproject.org/files/ENCFF790DJT/, Ross Hardison, PennState
# wget https://www.encodeproject.org/files/ENCFF790DJT/@@download/ENCFF790DJT.bed.gz
# gzip --decompress --stdout ENCFF790DJT.bed.gz | bedtools sort -i - > mm10.Hardison.blacklist.full.bed
# 
# psublacklist.mm10.bed.gz
# https://www.encodeproject.org/files/ENCFF226BDM/, Ross Hardison, PennState
# wget https://www.encodeproject.org/files/ENCFF226BDM/@@download/ENCFF226BDM.bed.gz
# gzip --decompress --stdout ENCFF226BDM.bed.gz | bedtools sort -i - > mm10.Hardison.psublacklist.mm10.bed
# 
# anshul.blacklist.mm10.bed.gz
# https://www.encodeproject.org/files/ENCFF999QPV/, Anshul Kundaje, Stanford
# wget https://www.encodeproject.org/files/ENCFF999QPV/@@download/ENCFF999QPV.bed.gz
# gzip --decompress --stdout ENCFF999QPV.bed.gz | bedtools sort -i - > mm10.Kundaje.anshul.blacklist.mm10.bed
# 
# mm10.blacklist.bed.gz
# https://www.encodeproject.org/files/ENCFF547MET/, Anshul Kundaje, Stanford
# wget https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz
# gzip --decompress --stdout ENCFF547MET.bed.gz | bedtools sort -i - > mm10.Kundaje.mm10.blacklist.bed
# 
# mm10mitoblack.bed.gz
# https://www.encodeproject.org/files/ENCFF759PJK/, Barbara Wold, Caltech
# wget https://www.encodeproject.org/files/ENCFF759PJK/@@download/ENCFF759PJK.bed.gz
# gzip --decompress --stdout ENCFF759PJK.bed.gz | bedtools sort -i - > mm10.Wold.mm10mitoblack.bed
# 
# ## mm9 
# 
# mm9mitoblack.bed.gz
# https://www.encodeproject.org/files/ENCFF299EZH/, Barbara Wold, Caltech
# wget https://www.encodeproject.org/files/ENCFF299EZH/@@download/ENCFF299EZH.bed.gz
# gzip --decompress --stdout ENCFF299EZH.bed.gz | bedtools sort -i - > mm9.Wold.mm9mitoblack.bed
# 
# # http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/
# 
# http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz
# wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/dm3-D.melanogaster/dm3-blacklist.bed.gz
# gzip --decompress --stdout dm3-blacklist.bed.gz | bedtools sort -i - > dm3.Kundaje.dm3-blacklist.bed
# 
# http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz
# wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/ce10-C.elegans/ce10-blacklist.bed.gz
# gzip --decompress --stdout ce10-blacklist.bed.gz | bedtools sort -i - > ce10.Kundaje.ce10-blacklist.bed

# The following example demonstrate how the coordinates of assembly-specific 
# problematic regions were converted into Rds objects

# Download results data from https://drive.google.com/drive/folders/124DZtsU0YVWqkb7dgu8Nk6b3N8-ShVSC?usp=sharing
library(GenomicRanges)
library(readr)
# Folder with results
dir_in <- "/Volumes/GoogleDrive/My Drive/denydata"
# All BED files
files <- list.files(path = dir_in, pattern = "bed$")
# In each subfolder
for (file in files) {
  # Read "fimo.bed" created by "fimo.qsub"
  denyBED <- read_tsv(file.path(dir_in, file), col_names = FALSE)
  # Assign column names depending on the number of columns
  if (ncol(denyBED) == 3) { # Only 3 columns
    colnames(denyBED) <- c("chr", "start", "stop")
  } else { # If more than 3 columns, consider the first 6
    colnames(denyBED) <- c("chr", "start", "stop", "name", "score", "strand")
  }
  # Convert to GRanges object
  denyGR <- GenomicRanges::makeGRangesFromDataFrame(denyBED, keep.extra.columns = TRUE)
  # Add seqinfo
  # Parse out genome ID from the file name, to get hg19, hg38, mm9, mm10, etc.
  genome_id <- strsplit(file, ".", fixed = TRUE)[[1]][1]
  # Get chromosome info and match it to the chromosome order in denyBED
  chrom_data <- GenomeInfoDb::getChromInfoFromUCSC(genome = genome_id)
  chrom_data <- chrom_data[chrom_data$chrom %in% seqlevels(denyGR), ]
  chrom_data <- chrom_data[match(seqlevels(denyGR), chrom_data$chrom), ]
  # Check if chromosome order is the same
  if (!all.equal(seqlevels(denyGR), chrom_data$chrom)) {
    print(paste("Chromosome order does not match for", genome_id, "genome."))
    break
  }
  # Assign seqinfo data
  seqlengths(denyGR) <- chrom_data$size
  isCircular(denyGR) <- chrom_data$circular
  genome(denyGR)     <- genome_id

  # Save as Rds object. file extension is changed to rds
  saveRDS(object = denyGR, file = file.path(dir_in, sub("bed$", "rds", file)))
  # denyGR <- readRDS(file = file.path(dir_in, sub("bed$", "rds", file)))
}