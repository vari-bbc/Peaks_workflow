# for debugging
#save.image("/varidata/research/projects/bbc/research/JONR_20220120_cutnrun_VBCS-694/foo.Rdata")
#print(str(snakemake))

suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(Rsamtools))

# params from workflow
ref_fasta <- snakemake@params[["ref_fasta"]]#"/varidata/research/projects/bbc/versioned_references/2022-03-08_14.47.50_v9/data/mm10_gencode_plus_viruses_and_cfmedips/sequence/mm10_gencode_plus_viruses_and_cfmedips.fa"
outfile <- snakemake@output[[1]] #"output.txt"
print(outfile)
# make GenomicRanges from the genome
ref_gr <- as(seqinfo(FaFile(ref_fasta)), "GRanges")

# output the standard chromosome names as a space-separated text file
standardChromosomes(ref_gr) %>% paste(., collapse=" ") %>% writeLines(., outfile)

