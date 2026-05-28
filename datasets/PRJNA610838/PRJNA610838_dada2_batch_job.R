#!/usr/bin/env Rscript
#DADA2 Batch Job - Paired-end reads
#Load libraries
library(dada2)
library(Biostrings)

# Set working directory
setwd("/storage/home/hcoda1/4/rthomas376/rthomas376/PRJNA610838")

# Define reports directory
reports_dir <- file.path(getwd(), "reports")
dir.create(reports_dir, showWarnings=FALSE)

# Load filtered files from filterAndTrim step
load(file.path(getwd(), "PRJNA610838_filtered.RData"))
print("Starting DADA2 pipeline...")

### ===== LEARN ERROR RATES ===== ###
print("Learning error rates...")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
print("Error rates learned")

# Save error rate plots
png(file.path(reports_dir, "error_rates_forward.png"), width=800, height=600)
plotErrors(errF, nominalQ=TRUE)
dev.off()

png(file.path(reports_dir, "error_rates_reverse.png"), width=800, height=600)
plotErrors(errR, nominalQ=TRUE)
dev.off()

### ===== SAMPLE INFERENCE (DADA) ===== ###
print("Running DADA algorithm...")
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
print("Denoising complete")

### ===== MERGE PAIRED READS ===== ###
print("Merging paired reads...")
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=FALSE)
print("Paired reads merged")

### ===== MAKE SEQUENCE TABLE ===== ###
print("Creating sequence table...")
seqtab <- makeSequenceTable(mergers)
print(paste("Sequence table dimensions:", nrow(seqtab), "samples,", ncol(seqtab), "ASVs"))

### ===== REMOVE CHIMERAS ===== ###
print("Removing chimeras...")
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
chimera_retained <- round(sum(seqtab.nochim)/sum(seqtab)*100, 2)
print(paste("Chimera removal: kept", chimera_retained, "% of sequences"))

### ===== TRACKING TABLE ===== ###
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "non-chim")
rownames(track) <- sample.names
print("Tracking table:")
print(head(track))

### ===== SAVE DATA ===== ###
# Save sequence table
seqtab_clean_PRJNA610838 <- seqtab.nochim
save(seqtab_clean_PRJNA610838, file = file.path(getwd(), "PRJNA610838_seqtab.RData"))
print("Sequence table saved!")

# Save entire workspace
save.image(file = file.path(getwd(), "PRJNA610838_ready_for_taxonomy.RData"))
print("Full workspace saved!")

### ===== UPDATE MARKDOWN REPORT ===== ###
sink(file.path(reports_dir, "PRJNA610838_filtering_report.md"), append=TRUE)

cat("\n## Error Rate Assessment\n\n")
cat("Error rates were learned separately for forward and reverse reads using the `learnErrors()` function.\n\n")

cat("### Forward Error Rates\n\n")
cat("![Forward Error Rates](error_rates_forward.png)\n\n")

cat("### Reverse Error Rates\n\n")
cat("![Reverse Error Rates](error_rates_reverse.png)\n\n")

cat("## DADA2 Denoising and ASV Inference\n\n")
cat("The DADA algorithm was applied to denoise sequences and infer Amplicon Sequence Variants (ASVs).\n\n")

cat("## Paired-End Read Merging\n\n")
cat("Forward and reverse reads were merged based on sequence overlap.\n\n")

cat("## Chimera Removal\n\n")
cat("Bimeric sequences (chimeras) were removed using the consensus method.\n\n")
cat("**Sequences retained after chimera removal:** ", chimera_retained, "%\n\n")

cat("## Sequence Table Summary\n\n")
cat("| Metric | Value |\n")
cat("|--------|-------|\n")
cat("| Total Samples | ", nrow(seqtab.nochim), " |\n")
cat("| Total ASVs | ", ncol(seqtab.nochim), " |\n")
cat("| Mean Reads per Sample | ", round(mean(rowSums(seqtab.nochim)), 0), " |\n")
cat("| Min Reads in Sample | ", min(rowSums(seqtab.nochim)), " |\n")
cat("| Max Reads in Sample | ", max(rowSums(seqtab.nochim)), " |\n\n")

cat("## Processing Pipeline Tracking\n\n")
cat("| Sample | Input | Filtered | Denoised | Merged | Non-Chim |\n")
cat("|--------|-------|----------|----------|--------|----------|\n")
for(i in 1:nrow(track)) {
  cat("| ", rownames(track)[i], " | ", 
      track[i,1], " | ", 
      track[i,2], " | ", 
      track[i,3], " | ", 
      track[i,4], " | ", 
      track[i,5], " |\n", sep="")
}

sink()

# Save tracking table as CSV
write.csv(track, file.path(reports_dir, "PRJNA610838_tracking_table.csv"))

print("Markdown report updated!")
print("DADA2 pipeline complete!")