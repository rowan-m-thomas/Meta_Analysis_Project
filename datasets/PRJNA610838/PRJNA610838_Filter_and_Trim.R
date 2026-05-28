### Set your working directory to your dataset folder ###
setwd("/storage/home/hcoda1/4/rthomas376/rthomas376/PRJNA610838/")

### Check you're in the right place and list the files present ###
getwd()
list.files()

### Install and load necessary packages: ###
#install.packages("devtools")
library(devtools)
library(dada2)
packageVersion("dada2") #1.30.0

#BiocManager::install("Biostrings")
library(Biostrings)
packageVersion("Biostrings") #version 2.70.3

### Define the path variable ###
path <- getwd() #sets the path as the directory containing the fastq_files
list.files(path) #lists the files within here; verifies the path is correct

### Forward and reverse fastq filenames have format: SAMPLENAME_1.fastq and SAMPLENAME_2.fastq ###
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE)) 
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

### Check how many files found ###
cat("Found", length(fnFs), "forward files\n") #10 forward files
cat("Found", length(fnRs), "reverse files\n") #10 reverse files

### Extract sample names from the filenames ###
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

### Look at the quality profile of the forward and reverse reads ###
plotQualityProfile(fnFs[1:4]) #Great quality, starts to dip towards the end. Recommend trimming to 240bp.
plotQualityProfile(fnRs[1:4]) #Reverse reads aren't as high quality, recommend trimming to 180bp.

### Place filtered files in filtered/ subdirectory ###
#These lines create file paths for your filtered fastq files that will be generated after the filterandtrim step
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
#These lines assign sample names as labels to each file path
names(filtFs) <- sample.names
names(filtRs) <- sample.names

### Run filterAndTrim ###
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                            truncLen=c(240, 180),
                            maxN=0, 
                            maxEE=c(2, 3), 
                            truncQ=2, 
                            rm.phix=TRUE,
                            compress=TRUE, 
                            multithread=TRUE)

### ===== SAVE QUALITY PLOTS AND FILTERING REPORT ===== ###

### Create reports folder ###
reports_dir <- file.path(path, "reports")
dir.create(reports_dir, showWarnings=FALSE)

### Save quality profile plots ###
png(file.path(reports_dir, "quality_profile_forward.png"), width=800, height=600)
plotQualityProfile(fnFs[1:4])
dev.off()

png(file.path(reports_dir, "quality_profile_reverse.png"), width=800, height=600)
plotQualityProfile(fnRs[1:4])
dev.off()

### Create filtering report with plots ###
sink(file.path(reports_dir, "PRJNA610838_filtering_report.md"))

cat("# DADA2 Filtering Report\n\n")
cat("**Dataset:** PRJNA610838\n\n")
cat("**Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("## Quality Assessment\n\n")
cat("### Forward Reads Quality Profile\n\n")
cat("![Forward Quality](quality_profile_forward.png)\n\n")
cat("Forward reads show good quality throughout, with slight degradation toward the end.\n\n")

cat("### Reverse Reads Quality Profile\n\n")
cat("![Reverse Quality](quality_profile_reverse.png)\n\n")
cat("Reverse reads show lower quality compared to forward reads, particularly at the end.\n\n")

cat("## Filtering Parameters\n\n")
cat("Based on quality profiles above, the following parameters were used:\n\n")
cat("```\n")
cat("truncLen = c(240, 180)  # Truncate forward to 240bp, reverse to 180bp\n")
cat("maxN = 0                # No ambiguous bases allowed\n")
cat("maxEE = c(2, 3)         # Max expected errors (forward=2, reverse=3)\n")
cat("truncQ = 2              # Truncate at first Q score < 2\n")
cat("rm.phix = TRUE          # Remove PhiX contamination\n")
cat("```\n\n")

cat("## Summary Statistics\n\n")
cat("| Metric | Value |\n")
cat("|--------|-------|\n")
cat("| Total Samples | ", nrow(out), " |\n")
cat("| Mean Reads Input | ", round(mean(out[,1]), 0), " |\n")
cat("| Mean Reads Output | ", round(mean(out[,2]), 0), " |\n")
cat("| Mean Retention Rate | ", round(mean(out[,2]/out[,1]*100), 1), "% |\n")
cat("| Min Reads Output | ", min(out[,2]), " |\n")
cat("| Max Reads Output | ", max(out[,2]), " |\n\n")

sink()

### Save results as CSV backup ###
write.csv(filtering_results, file.path(reports_dir, "PRJNA610838_filtering_results.csv"), row.names=FALSE)

save.image(file = "PRJNA610838_filtered.RData")
