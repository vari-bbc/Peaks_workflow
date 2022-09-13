#!/usr/bin/env Rscript

options(warn=1) # print all warnings as they happen

library(DiffBind)
library(readr)
library(dplyr)
library(stringr)
library(Rsamtools)

# variables
samplesheet <- snakemake@input[['samplesheet']] 
db_config <- list(cores=snakemake@threads)
out_dir <- snakemake@params[["outdir"]]
out_samplesheet <- snakemake@output[['samplesheet']]
DB_summits <- snakemake@params[['DB_summits']]
macs2_type <- snakemake@params[['macs2_type']]
subtract_controls <- as.logical(snakemake@params[['subtract_controls']])
curr_enriched_factor <- snakemake@params[['enriched_factor']]

# set up Diffbind samplesheet
samples <- read_tsv(samplesheet) %>%
        dplyr::mutate(SampleID=sample, 
                      bamReads=str_glue("analysis/filt_bams/{sample}_filt_alns.bam"), 
                      bamControl=ifelse(!is.na(control), str_glue("analysis/filt_bams/{control}_filt_alns.bam"), NA),
                      Peaks=str_glue("analysis/{macs2_type}/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak"),
                      PeakCaller="bed",
                      Factor=sample_group,
                      out_pref=ifelse(is.na(enriched_factor), "peaks", enriched_factor)) %>%
        dplyr::select(SampleID, bamReads, bamControl, Peaks, PeakCaller, Factor, out_pref) %>%
        unique() # some samples represented twice if they had more than 1 set of fastq files

if(any(!is.na(samples$bamControl))){
    samples <- samples %>% dplyr::filter(!is.na(bamControl)) # if project has controls, then we remove all the rows for the controls
} else{
    samples <- samples %>% dplyr::select(-bamControl)
}

# filter for just the samples related to the current enriched factor
stopifnot(curr_enriched_factor %in% samples$out_pref)
samples <- samples %>% dplyr::filter(out_pref == curr_enriched_factor)

# subtract controls?
if (!subtract_controls){
    samples$bamControl <- NA
}

write_tsv(samples, out_samplesheet)

print(samples)

message(str_glue("Making DBA for {curr_enriched_factor}"))

# make DBA object
diffbind <- dba(sampleSheet = samples %>% dplyr::select(-out_pref))

# set num cores and how many reads to store in memory
diffbind$config$cores <- db_config$cores
diffbind$config$yieldSize <- 200000000

# get counts
diffbind$config$scanbamparam <- ScanBamParam(flag = scanBamFlag(isDuplicate=FALSE, isSecondaryAlignment=FALSE, isSupplementaryAlignment=FALSE, isPaired=TRUE, isProperPair=TRUE, isNotPassingQualityControls=FALSE, isUnmappedQuery=FALSE, hasUnmappedMate=FALSE, isMinusStrand=NA, isMateMinusStrand=NA))

diffbind <- dba.count(diffbind, summits=DB_summits, bUseSummarizeOverlaps = TRUE, bParallel = TRUE)
print(diffbind)

# normalization
diffbind <- dba.normalize(diffbind, normalize = DBA_NORM_NATIVE, background = TRUE)

# print the normalization methods
dba.normalize(diffbind, bRetrieve=TRUE)

dir.create(out_dir, recursive=TRUE)
# write the DBA to an RDS file
saveRDS(diffbind, str_glue("{out_dir}/{curr_enriched_factor}.rds"))


# session info
sessionInfo()

