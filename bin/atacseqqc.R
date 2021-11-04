library(ATACseqQC)
library(stringr)
library(dplyr)
library(readr)
library(GenomicFeatures)

args <- commandArgs(trailingOnly=TRUE)

bamfile <- args[1]
out_pref <- args[2]
known_genes_lib <- args[3]
library(known_genes_lib, character.only=TRUE) # load the known genes for your organism

if (length(args)!=3) {
      stop("Three arguments must be supplied")
} 

if (length(args)==1 & !str_detect(bamfile, "bam$")) {
      stop("Arg 1 must be bam file")
} 


## bamfile tags to be read in
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                 "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                               "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                               "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                               "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                               "U2"))
library(Rsamtools)
bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags

gal <- readBamFile(bamfile, asMates=FALSE, tag=tags)#, bigFile=TRUE)
shiftedBamfile <- paste0(out_pref, "_shifted.bam")
#gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
gal1 <- gal
txs <- transcripts(eval(as.name(known_genes_lib)))

#gtfFile <- "/varidata/research/projects/bbc/versioned_references/2021-08-10_11.12.27_v6/data/mm10_gencode_plus_viruses_and_cfmedips/annotation/mm10_gencode_plus_viruses_and_cfmedips.gtf"
#chrominfo <- read_tsv("/varidata/research/projects/bbc/versioned_references/2021-08-10_11.12.27_v6/data/mm10_gencode_plus_viruses_and_cfmedips/sequence/mm10_gencode_plus_viruses_and_cfmedips.fa.fai", col_names=c("chrom", "length", "OFFSET", "LINEBASES", "LINEWIDTH")) %>% dplyr::select(chrom, length)

#txdb <- makeTxDbFromGFF(file=gtfFile,
             #chrominfo=chrominfo)

#txs <- transcripts(txdb)

# filter for standard chromosomes
txs <- keepStandardChromosomes(txs, pruning.mode = 'coarse')

tsse <- TSSEscore(gal1, txs)
write_rds(tsse, paste0(out_pref, "_tsse.rds"))
tsse$TSSEscore

pdf(paste0(out_pref, "_tsse.pdf"))
plot(100*(-9:10-.5), tsse$values, type="b", 
        xlab="distance to TSS",
        ylab="aggregate TSS score")
dev.off()


sessionInfo()
