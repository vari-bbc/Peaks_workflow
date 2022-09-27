import pandas as pd
import numpy as np
import os
import re
import itertools
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("6.15.0")

##### Editable variables #####

configfile: "bin/config.yaml"

bamCoverage_binsize = 10

macs2_species=config['macs2']['gsize']

ref_fa = config['ref']['sequence']
ref_fai = config['ref']['fai']
gtf = config['ref']['annotation']
bwa_index = config['ref']['index']
blacklist = config['ref']['blacklist']

mito_chrom=config['ref']['mito_chr']
fai_parsed = pd.read_table(ref_fai, names=['chr','len','offset','bases_per_line','bytes_per_line'])
if mito_chrom not in fai_parsed['chr'].values:
    raise Exception('{mito} not found in reference fai file.'.format(mito=mito_chrom))

chroms_no_mito = ' '.join(fai_parsed[fai_parsed['chr'] != mito_chrom]['chr'].values)

# When we filter alignments for peak calling, we remove contigs smaller than 'cutoff'. This is meant to be a way to retain only the nuclear chromosomes (not mitochondrial).
#chrom_min_bp = config['ref']['chrom_min_bp']
#chroms_gt_cutoff = ' '.join(fai_parsed[fai_parsed['len'] > chrom_min_bp]['chr'].values)

##### load config and sample sheets #####
samplesheet="bin/samples.tsv"
units = pd.read_table(samplesheet, dtype={"sample" : str, "sample_group" : str })
units['se_or_pe'] = ["SE" if x else "PE" for x in units['fq2'].isnull()]

samples = units[["sample","control","sample_group","enriched_factor","se_or_pe"]].drop_duplicates()
if not samples['sample'].is_unique:
    raise Exception('A sample has more than one combination of control, sample_group, enriched_factor, and/or se_or_pe.')

# Filter for sample rows that are not controls
controls_list = list(itertools.chain.from_iterable( [x.split(',') for x in samples['control'].values if not pd.isnull(x)] ))
samples_no_controls = samples[-samples['sample'].isin(controls_list)].copy()
samples_no_controls["enriched_factor"]=samples_no_controls["enriched_factor"].fillna("peaks")

# if SE data, make sure it is not ATACseq
if (config['atacseq'] and units['fq2'].isnull().any()):
    raise Exception('SE ATAC-seq not supported.')

snakemake_dir = os.getcwd() + "/"

# make a tmp directory for analyses
tmp_dir = os.path.join(snakemake_dir, "tmp")
if not os.path.exists(tmp_dir):
    os.mkdir(tmp_dir)

shared_snakemake_dir=config['shared_snakemake_repo'] + "/rules/"

include: os.path.join(shared_snakemake_dir, "post_alignment/flagstat")
include: os.path.join(shared_snakemake_dir, "post_alignment/CollectInsertSizeMetrics")
include: os.path.join(shared_snakemake_dir, "post_alignment/CollectAlignmentSummaryMetrics")

rule all:
    input:
        #expand("analysis/bwamem/{sample.sample}.bam", sample=samples.itertuples()),
        #expand("analysis/filt_bams/{sample.sample}_filt_alns.bam", sample=samples.itertuples()),
        #expand("analysis/deeptools_cov_rmdups/{sample.sample}_filt_alns_rmdups.bw", sample=samples.itertuples()),
        expand("analysis/merge_bigwigs/{group}.bw", group=pd.unique(samples['sample_group'])),
        "analysis/multiqc/multiqc_report.html",
        "analysis/deeptools_fingerprint/fingerprint.rawcts",
        #expand("analysis/peaks_rm_blacklist/{sample}_macs2_{size}_peaks.rm_blacklist.{size}Peak", sample=samples["sample"], size=["narrow","broad"]),
        #expand("analysis/hmmratac/{sample}_summits.bed", sample=samples["sample"]) if config['atacseq'] else [],
        expand("analysis/deeptools_plotenrichment/{sample}.pdf", sample=samples_no_controls["sample"]),
        expand("analysis/deeptools_bamcompare/{sample}_vs_control_log2ratio.bw", sample=samples_no_controls[pd.notnull(samples_no_controls['control'])]["sample"]),
        expand("analysis/filt_bams_nfr/CollectInsertSizeMetrics/{sample.sample}_filt_alns_nfr.insert_size_metrics.txt", sample=samples.itertuples()) if config['atacseq'] else [],
        "analysis/peaks_venn/report.html",
        "analysis/deeptools_heatmap_genes/genes.pdf",
        expand("analysis/atacseqc/{sample}{suffix}", sample=samples_no_controls["sample"], suffix=["_tsse.rds","_tsse.pdf"]) if config['atacseq'] and config['run_ATACseqQC'] else [],
        "analysis/deeptools_plotCorr/corr_ht.pdf",
        "analysis/deeptools_plotPCA/pca.pdf",
        expand("analysis/deeptools_cov_rmdups_{norm_type}/{sample}_filt_alns_rmdups.bw", norm_type=config['addtnl_bigwig_norms'], sample=samples_no_controls['sample']) if isinstance(config['addtnl_bigwig_norms'], list) else [],
        expand("analysis/homer_find_motifs/{sample}/homerMotifs.all.motifs", sample=samples_no_controls["sample"]) if config['homer']['run'] else [], #sample=samples[pd.notnull(samples['enriched_factor'])]['sample'])
        expand("analysis/diffbind_count/{factor}.rds", factor=pd.unique(samples_no_controls["enriched_factor"]))

def get_orig_fastq(wildcards):
    if wildcards.read == "R1":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq1"].values)
    elif wildcards.read == "R2":
            fastq = expand("raw_data/{fq}", fq = units[units["sample"] == wildcards.sample]["fq2"].values)
    return fastq 

rule rename_fastqs:
    """
    Rename fastqs by biologically meaningful name. Concatenate different runs of same library.
    """
    input:
        get_orig_fastq
    output:
        "analysis/renamed_data/{sample}_{read}.fastq.gz"
    log:
        stdout="logs/rename_fastqs/{sample}_{read}.o",
        stderr="logs/rename_fastqs/{sample}_{read}.e",
    benchmark:
        "benchmarks/rename_fastqs/{sample}_{read}.txt"
    params:
    threads: 1
    resources:
        mem_gb=8
    envmodules:
    shell:
        """
        if [ `printf '{input}' | wc -w` -gt 1 ]
        then
            cat {input} > {output}
        else
            ln -sr {input} {output}
        fi

        """

rule fastqc:
    """
    Run fastqc on raw_data/ files.
    """
    input:
        "analysis/renamed_data/{fq_pref}.fastq.gz"
    output:
        html="analysis/fastqc/{fq_pref}_fastqc.html",
        zip="analysis/fastqc/{fq_pref}_fastqc.zip"
    params:
        outdir="analysis/fastqc/"
    log:
        stdout="logs/fastqc/{fq_pref}.o",
        stderr="logs/fastqc/{fq_pref}.e"
    benchmark:
        "benchmarks/fastqc/{fq_pref}.txt"
    envmodules:
        config['modules']['fastqc']
    threads: 1
    resources:
        mem_gb = 32
    shell:
        """
        fastqc --outdir {params.outdir} {input}
        """

rule fastq_screen:
    """
    Run fastq_screen to detect any contamination from other species or excessive rRNA.
    """
    input:
        "analysis/renamed_data/{fq_pref}.fastq.gz"
    output:
        html = "analysis/fastq_screen/{fq_pref}_screen.html",
        txt = "analysis/fastq_screen/{fq_pref}_screen.txt",
    params:
    log:
        stdout="logs/fastq_screen/{fq_pref}.o",
        stderr="logs/fastq_screen/{fq_pref}.e"
    benchmark:
        "benchmarks/fastq_screen/{fq_pref}.txt"
    envmodules:
        config['modules']['fastq_screen']
    threads: 8
    resources:
        mem_gb = 32
    shell:
        """
        fastq_screen --threads {threads} --outdir analysis/fastq_screen/ {input}
        """

rule trim_galore_PE:
    """
    Run trim_galore on paired-end reads.
    """
    input:
        expand("analysis/renamed_data/{{sample}}_R{read}.fastq.gz", read=["1","2"])
    output:
        temp(expand("analysis/trim_galore/{{sample}}_R{ext}", ext=["1_val_1.fq.gz","2_val_2.fq.gz"])),
        expand("analysis/trim_galore/{{sample}}_R1{ext}", ext=[".fastq.gz_trimming_report.txt","_val_1_fastqc.html"]),
        expand("analysis/trim_galore/{{sample}}_R2{ext}", ext=[".fastq.gz_trimming_report.txt","_val_2_fastqc.html"]),
    params:
    log:
        stdout="logs/trim_galore/{sample}.o",
        stderr="logs/trim_galore/{sample}.e"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        mem_gb = 64
    shell:
        """
        trim_galore --paired {input} --output_dir analysis/trim_galore/ --cores {threads} --fastqc
        """

rule trim_galore_SE:
    """
    Run trim_galore on single-end reads.
    """
    input:
        "analysis/renamed_data/{sample}_R1.fastq.gz"
    output:
        temp("analysis/trim_galore/{sample}_R1_trimmed.fq.gz"),
        "analysis/trim_galore/{sample}_R1.fastq.gz_trimming_report.txt",
        expand("analysis/trim_galore/{{sample}}_R1_trimmed_fastqc{ext}", ext=['.html','.zip']),
    params:
    log:
        stdout="logs/trim_galore/{sample}.o",
        stderr="logs/trim_galore/{sample}.e"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        config['modules']['trim_galore']
    threads: 4
    resources:
        mem_gb = 64
    shell:
        """
        trim_galore {input} --output_dir analysis/trim_galore/ --cores {threads} --fastqc
        """

rule bwamem:
    input:
        lambda wildcards: expand("analysis/trim_galore/{sample}_R{read}_val_{read}.fq.gz", read=["1","2"], sample=wildcards.sample) if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "analysis/trim_galore/{sample}_R1_trimmed.fq.gz"
    output:
        outbam=temp("analysis/bwamem/{sample}.bam"),
        outbai="analysis/bwamem/{sample}.bam.bai",
        idxstat="analysis/bwamem/{sample}.bam.idxstat",
        samblaster_err="analysis/bwamem/{sample}.samblaster.e",
    log:
        stdout="logs/bwamem/{sample}.o",
        stderr="logs/bwamem/{sample}.e",
    benchmark:
        "benchmarks/bwamem/{sample}.txt"
    params:
        bwa_idx=bwa_index,
        samblaster_params=lambda wildcards: "--addMateTags" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "--ignoreUnmated"
    threads: 16
    envmodules:
        config['modules']['bwa'],
        config['modules']['samblaster'],
        config['modules']['samtools'],
    resources:
        mem_gb=180
    shell:
        """
        bwa mem \
        -t {threads} \
        {params.bwa_idx} \
        {input} | \
        samblaster {params.samblaster_params} 2>{output.samblaster_err} | \
        samtools sort \
        -m 6G \
        -@ {threads} \
        -O "BAM" \
        -o {output.outbam} \
        -

        echo "END bwamem"
        echo "END bwamem" 1>&2

        samtools index -@ {threads} {output.outbam}

        echo "END indexing"
        echo "END indexing" 1>&2
        
        samtools idxstats {output.outbam} > {output.idxstat}

        echo "END idxstats"
        echo "END idxstats" 1>&2
 
        """

rule atacseqc:
    input:
        bam="analysis/bwamem/{sample}.bam"
    output:
        expand("analysis/atacseqc/{{sample}}{suffix}", suffix=["_tsse.rds","_tsse.pdf"]),
    log:
        stdout="logs/atacseqc/{sample}.o",
        stderr="logs/atacseqc/{sample}.e",
    benchmark:
        "benchmarks/atacseqc/{sample}.txt"
    params:
        outpref="analysis/atacseqc/{sample}",
        knownGenesLib=config["atacseqc"]["known_genes_pkg"]
    threads: 4
    envmodules:
        "bbc/R/R-4.1.0-setR_LIBS_USER"
    resources:
        mem_gb=96
    shell:
        """
        Rscript --vanilla bin/atacseqqc.R {input.bam} {params.outpref} {params.knownGenesLib}
        """


rule filt_bams:
    """
    Filter for mapq score, properly paired, primary alignment, remove unmapped. Do not remove dup reads.
    """
    input:
        "analysis/bwamem/{sample}.bam"
    output:
        bam="analysis/filt_bams/{sample}_filt_alns.bam",
        bai="analysis/filt_bams/{sample}_filt_alns.bam.bai",
        idxstat="analysis/filt_bams/{sample}_filt_alns.bam.idxstat",
    log:
        stdout="logs/filt_bams/{sample}.o",
        stderr="logs/filt_bams/{sample}.e",
    benchmark:
        "benchmarks/filt_bams/{sample}.txt"
    params:
        mapq=30,
        flags_to_exclude="2828",
        flags_to_include=lambda wildcards: "-f 2" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "",
        keep_chroms = chroms_no_mito #chroms_gt_cutoff #chroms_no_mito if atacseq else ''
    threads: 8
    resources:
        mem_gb=80
    envmodules:
        config['modules']['samtools']
    shell:
        """
        samtools view \
              -@ {threads} \
              -q {params.mapq} \
              -F {params.flags_to_exclude} \
              {params.flags_to_include} \
              -b -o {output.bam} {input} {params.keep_chroms}

        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.idxstat}

        """

rule dedup_bams:
    """
    Dedup BAMs.
    """
    input:
        bam="analysis/filt_bams/{sample}_filt_alns.bam",
        bai="analysis/filt_bams/{sample}_filt_alns.bam.bai",
    output:
        dedup_bam=temp("analysis/dedup_bams/{sample}_filt_alns.dedup.bam"),
        dedup_bai="analysis/dedup_bams/{sample}_filt_alns.dedup.bam.bai",
    log:
        stdout="logs/dedup_bams/{sample}.o",
        stderr="logs/dedup_bams/{sample}.e",
    benchmark:
        "benchmarks/dedup_bams/{sample}.txt"
    params:
    threads: 8
    resources:
        mem_gb=80
    envmodules:
        config['modules']['samtools']
    shell:
        """
        samtools view -b -@ {threads} -F 1024 -o {output.dedup_bam} {input.bam}
        samtools index {output.dedup_bam}
        """

rule filt_bams_nfr:
    """
    Keep only paired alignments with <150 bp fragment size.
    """
    input:
        "analysis/dedup_bams/{sample}_filt_alns.dedup.bam"
    output:
        bam="analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam",
        bai="analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam.bai",
        idxstat="analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam.idxstat",
        #metrics="analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam.metrics"
    log:
        stdout="logs/filt_bams_nfr/{sample}.o",
        stderr="logs/filt_bams_nfr/{sample}.e",
    benchmark:
        "benchmarks/filt_bams_nfr/{sample}.txt"
    params:
        temp=os.path.join(snakemake_dir, "tmp"),
        maxFragmentSize = '150',
    threads: 8
    envmodules:
        config['modules']['bamtools'],
        #config['modules']['deeptools'],
        config['modules']['samtools']
    resources:
        mem_gb=80
    shell:
        """
        export TMPDIR={params.temp}
       
        bamtools filter  -insertSize "<={params.maxFragmentSize}" -in {input} -out {output.bam}
        
        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.idxstat}

        """


rule deeptools_cov_rmdups:
    input:
        bam="analysis/filt_bams/{sample}_filt_alns.bam"
    output:
        bigwig_rmdups="analysis/deeptools_cov_rmdups/{sample}_filt_alns_rmdups.bw"
    log:
        stdout="logs/deeptools_cov_rmdups/{sample}.o",
        stderr="logs/deeptools_cov_rmdups/{sample}.e"
    benchmark:
        "benchmarks/deeptools_cov_rmdups/{sample}.txt"
    params:
        blacklist=blacklist,
        binsize=bamCoverage_binsize,
        norm_method="CPM",
        sam_keep=lambda wildcards: "--samFlagInclude 64" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "", # count only first in pair if PE
        sam_exclude="1024",
        extend_reads=lambda wildcards: "--extendReads" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "", # extend to frag size if PE
        temp=os.path.join(snakemake_dir, "tmp")
    threads: 16
    envmodules:
        config['modules']['deeptools']
    resources:
        mem_gb=96
    shell:
        """
        export TMPDIR={params.temp}
        
        bamCoverage \
        -p {threads} \
        --binSize {params.binsize} \
        --bam {input.bam} \
        {params.extend_reads} \
        --blackListFileName {params.blacklist} \
        -o {output.bigwig_rmdups} \
        --normalizeUsing {params.norm_method} \
        --samFlagExclude {params.sam_exclude} \
        {params.sam_keep}

        """

rule get_std_chrom_names:
    input:
    output:
        "analysis/misc/std_chroms.txt"
    log:
        stdout="logs/get_std_chrom_names/out.o",
        stderr="logs/get_std_chrom_names/err.e",
    benchmark:
        "benchmarks/get_std_chrom_names/bench.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    threads: 1
    resources:
        mem_gb=20
    envmodules:
        config['modules']['R']
    script:
        "bin/scripts/get_standard_chrom_names.R"

rule csaw_norm_factors:
    input:
        bams=lambda wildcards: expand("analysis/filt_bams/{sample}_filt_alns.bam", sample=samples_no_controls[samples_no_controls['enriched_factor']==wildcards.enriched_factor]['sample']),
        std_chroms="analysis/misc/std_chroms.txt"
    output:
        bkgrd="analysis/bigwig_norm_factors/{enriched_factor}_csaw_bkgd.tsv", 
        hiAbund="analysis/bigwig_norm_factors/{enriched_factor}_csaw_hiAbund.tsv",
        rds=expand("analysis/bigwig_norm_factors/{{enriched_factor}}_{obj_nm}.rds", obj_nm=['binned','small_wins','filt_small_wins'])
    log:
        stdout="logs/csaw_norm_factors/{enriched_factor}.o",
        stderr="logs/csaw_norm_factors/{enriched_factor}.e",
    benchmark:
        "benchmarks/csaw_norm_factors/{enriched_factor}.txt"
    params:
        mito_chr=mito_chrom,
        blacklist=config['ref']['blacklist'],
        out_pref="analysis/bigwig_norm_factors/{enriched_factor}",
        bam_samp_names=lambda wildcards: samples_no_controls[samples_no_controls['enriched_factor']==wildcards.enriched_factor]['sample']
    threads: 8
    resources:
        mem_gb=120
    envmodules:
        config['modules']['R']
    script:
        "bin/scripts/calc_norm_factors.R"

rule cutnrun_ecoli_scale_factors:
    """
    Adapted from https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1
    """
    input:
        bam=lambda wildcards: expand("analysis/filt_bams/{sample}_filt_alns.bam", sample=samples['sample']),
    output:
        "analysis/bigwig_norm_factors/cutnrun_ecoli_scale_factors.tsv"
    log:
        stdout="logs/cutnrun_ecoli_scale_factors/out.o",
        stderr="logs/cutnrun_ecoli_scale_factors/err.e",
    benchmark:
        "benchmarks/cutnrun_ecoli_scale_factors/bench.txt"
    params:
        ecoli_chrom=config["ecoli_chrom"]
    threads: 8
    resources:
        mem_gb=120
    envmodules:
        config['modules']['samtools']
    shell:
        """
        printf "sample\\tecoli_frags\\tfinal.factors\\n" > {output}
        for bam in {input.bam}
        do
            sample=`basename $bam | perl -npe 'chomp; die unless /(_filt_alns.bam)/; s/$1//'`
            ecoli_frags=`samtools view -@{threads} -c -f 64 -F 1024 $bam {params.ecoli_chrom}` # count first read only, and ignore duplicates
            scale_factor=`echo "10000 / $ecoli_frags" | bc -l`
            printf "$sample\\t$ecoli_frags\\t$scale_factor\\n" >> {output}
        done
        """

def get_bigwig_norm_factors_file(wildcards):
    curr_enriched = samples[samples['sample']==wildcards.sample]['enriched_factor'].values[0]
    if (wildcards.norm_type == "csaw_bkgd"):
        return "analysis/bigwig_norm_factors/{enriched_factor}_csaw_bkgd.tsv".format(enriched_factor=curr_enriched)
    if (wildcards.norm_type == "csaw_hiAbund"):
        return "analysis/bigwig_norm_factors/{enriched_factor}_csaw_hiAbund.tsv".format(enriched_factor=curr_enriched)
    if (wildcards.norm_type == "ecoli"):
        return "analysis/bigwig_norm_factors/cutnrun_ecoli_scale_factors.tsv"

def get_bigwig_norm_factor(wildcards, input):
    df = pd.read_table(input.norm_factors)
    scalefactor = df[df['sample']==wildcards.sample]['final.factors'].values[0]
    return str(scalefactor)

rule alternate_norm_bigwigs:
    """
    Assume PE reads, so we set --extendReads by default and --samFlagInclude 64 (first read only)
    """
    input:
        bam="analysis/filt_bams/{sample}_filt_alns.bam",
        norm_factors=get_bigwig_norm_factors_file,
    output:
        bigwig_rmdups="analysis/deeptools_cov_rmdups_{norm_type}/{sample}_filt_alns_rmdups.bw"
    log:
        stdout="logs/deeptools_cov_rmdups_{norm_type}/{sample}.o",
        stderr="logs/deeptools_cov_rmdups_{norm_type}/{sample}.e",
    benchmark:
        "benchmarks/deeptools_cov_rmdups_{norm_type}/{sample}.txt"
    params:
        blacklist=blacklist,
        binsize=bamCoverage_binsize,
        scale_factor=get_bigwig_norm_factor,
        sam_exclude="1024",
        norm="None",
        temp="analysis/deeptools_cov_rmdups_{norm_type}/{sample}.tmp"
    threads: 8
    resources:
        mem_gb=120
    envmodules:
        config['modules']['deeptools']
    shell:
        """
        export TMPDIR={params.temp}
        mkdir -p {params.temp}

        bamCoverage \
        -p {threads} \
        --binSize {params.binsize} \
        -b {input.bam} \
        --extendReads \
        --blackListFileName {params.blacklist} \
        -o {output.bigwig_rmdups} \
        --normalizeUsing {params.norm} \
        --scaleFactor {params.scale_factor} \
        --samFlagExclude {params.sam_exclude} \
        --samFlagInclude 64

        rm -r {params.temp}
        """

rule deeptools_multiBWsummary:
    input:
        bws=expand("analysis/deeptools_cov_rmdups/{sample}_filt_alns_rmdups.bw", sample=samples['sample'].values),
        peaks="analysis/merge_all_peaks/all_merged_{peak_type}.bed".format(peak_type="macs2_nfr_broad" if config['atacseq'] else "macs2_narrow")
    output:
        "analysis/deeptools_multiBWsummary/results.npz"
    log:
        stdout="logs/deeptools_multiBWsummary/out.o",
        stderr="logs/deeptools_multiBWsummary/err.e"
    benchmark:
        "benchmarks/deeptools_multiBWsummary/bench.txt"
    params:
        blacklist=blacklist,
        temp=os.path.join(snakemake_dir, "tmp")
    threads: 8
    envmodules:
        config['modules']['deeptools']
    resources:
        mem_gb=96
    shell:
        """
        export TMPDIR={params.temp}
       
        multiBigwigSummary BED-file --BED {input.peaks} -p {threads} --blackListFileName {params.blacklist} -b {input.bws} -o {output}

        """

rule deeptools_plotCorr:
    input:
        "analysis/deeptools_multiBWsummary/results.npz"
    output:
        pdf="analysis/deeptools_plotCorr/corr_ht.pdf"
    log:
        stdout="logs/deeptools_plotCorr/out.o",
        stderr="logs/deeptools_plotCorr/err.e"
    benchmark:
        "benchmarks/deeptools_plotCorr/bench.txt"
    params:
        temp=os.path.join(snakemake_dir, "tmp")
    threads: 1
    envmodules:
        config['modules']['deeptools']
    resources:
        mem_gb=64
    shell:
        """
        export TMPDIR={params.temp}

        plotCorrelation --corData {input} -c spearman -p heatmap --skipZeros -o {output.pdf}

        """

rule deeptools_plotPCA:
    input:
        "analysis/deeptools_multiBWsummary/results.npz"
    output:
        pdf="analysis/deeptools_plotPCA/pca.pdf"
    log:
        stdout="logs/deeptools_plotPCA/out.o",
        stderr="logs/deeptools_plotPCA/err.e"
    benchmark:
        "benchmarks/deeptools_plotPCA/bench.txt"
    params:
        temp=os.path.join(snakemake_dir, "tmp")
    threads: 1
    envmodules:
        config['modules']['deeptools']
    resources:
        mem_gb=64
    shell:
        """
        export TMPDIR={params.temp}

        plotPCA --transpose -in {input} -o {output.pdf} --log2 --ntop 5000 
        """

rule deeptools_bamcompare:
    input:
        bam1="analysis/filt_bams/{sample}_filt_alns.bam",
        bam2=lambda wildcards: "analysis/filt_bams/{control}_filt_alns.bam".format(control=samples[samples['sample']==wildcards.sample]['control'].values[0]),
    output:
        bigwig="analysis/deeptools_bamcompare/{sample}_vs_control_log2ratio.bw"
    log:
        stdout="logs/deeptools_bamcompare/{sample}.o",
        stderr="logs/deeptools_bamcompare/{sample}.e"
    benchmark:
        "benchmarks/deeptools_bamcompare/{sample}.txt"
    params:
        blacklist=blacklist,
        binsize=bamCoverage_binsize,
        sam_keep=lambda wildcards: "--samFlagInclude 64" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "", # count only first in pair if PE
        sam_exclude="1024",
        extend_reads=lambda wildcards: "--extendReads" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "", # extend to frag size if PE
        temp=os.path.join(snakemake_dir, "tmp")
    threads: 16
    envmodules:
        config['modules']['deeptools']
    resources:
        mem_gb=96
    shell:
        """
        export TMPDIR={params.temp}
       
        bamCompare -b1 {input.bam1} -b2 {input.bam2} \
                -p {threads} \
                {params.extend_reads} \
                --binSize {params.binsize} \
                --blackListFileName {params.blacklist} \
                --samFlagExclude {params.sam_exclude} \
                {params.sam_keep} \
                -o {output.bigwig}

        """

rule prep_chromsizes_file:
    input:
        fai=config["ref"]["fai"]
    output:
        "analysis/prep_chromsizes_file/chrom_sizes.tsv"
    log:
        stdout="logs/prep_chromsizes_file/out.o",
        stderr="logs/prep_chromsizes_file/err.e",
    benchmark:
        "benchmarks/prep_chromsizes_file/bench.txt"
    params:
    threads: 1
    resources:
        mem_gb=8
    envmodules:
    shell:
        """
        cut -f 1,2 {input.fai} > {output}

        """

rule merge_bigwigs:
    input:
        bigwigs=lambda wildcards: expand("analysis/deeptools_cov_rmdups/{sample}_filt_alns_rmdups.bw", sample=samples[samples['sample_group']==wildcards.group]['sample']),
        chromsizes="analysis/prep_chromsizes_file/chrom_sizes.tsv"
    output:
        wig=temp("analysis/merge_bigwigs/{group}.wig"),
        bigwig="analysis/merge_bigwigs/{group}.bw"
    log:
        stdout="logs/merge_bigwigs/{group}.o",
        stderr="logs/merge_bigwigs/{group}.e"
    benchmark:
        "benchmarks/merge_bigwigs/{group}.txt"
    params:
        lambda wildcards, input, output:
        ("""
        # sum bigwigs at each interval
        wiggletools mean {in_bw} > {out_wig}

        # convert to bigwig
        wigToBigWig {out_wig} {in_chroms} {out_bw}

        """ if len(input.bigwigs) > 1 else """
        ln -rs {in_bw} {out_wig} # placeholder file. Will be removed
        ln -rs {in_bw} {out_bw}

        """).format(in_bw=input.bigwigs,
                   out_wig=output.wig,
                   out_bw=output.bigwig,
                   in_chroms=input.chromsizes)
    threads: 8
    envmodules:
        config['modules']['ucsc'],
        config['modules']['WiggleTools']
    resources:
        mem_gb=96
    shell:
        """
        {params}
        """

rule deeptools_heatmap_genes:
    input:
        bw=lambda wildcards: expand("analysis/merge_bigwigs/{sample_group}.bw", sample_group=pd.unique(samples["sample_group"])),
        genes=config['gtf_for_tss_heatmap'],
    output:
        compmat="analysis/deeptools_heatmap_genes/compmat.gz",
        compmatbed="analysis/deeptools_heatmap_genes/compmat.bed",
        heatmap="analysis/deeptools_heatmap_genes/genes.pdf"
    log:
        stdout="logs/deeptools_heatmap_genes/out.o",
        stderr="logs/deeptools_heatmap_genes/err.e"
    benchmark:
        "benchmarks/deeptools_heatmap_genes/bench.txt"
    envmodules:
        config['modules']['deeptools']
    params:
        after="2000",
        before="2000",
        binsize=bamCoverage_binsize,
        samp_labels=lambda wildcards, input: " ".join(os.path.basename(x).replace(".bw", "") for x in input.bw),
        temp="analysis/deeptools_heatmap_genes/ht_tmp",
        yaxislabel='"Mean CPMs"',
        blacklist=blacklist,
    threads: 16
    resources:
        mem_gb=100
    shell:
        """
        export TMPDIR={params.temp}

        computeMatrix \
        scale-regions \
        -p {threads} \
        -b {params.before} \
        -a {params.after} \
        --missingDataAsZero \
        --samplesLabel {params.samp_labels} \
        --binSize {params.binsize} \
        -R {input.genes} \
        -S {input.bw} \
        -o {output.compmat} \
        -bl {params.blacklist} \
        --transcriptID gene \
        --transcript_id_designator gene_id \
        --outFileSortedRegions {output.compmatbed}

        echo "END computeMatrix"
        echo "END computeMatrix" 1>&2

        plotHeatmap \
        --heatmapWidth 6 \
        --yAxisLabel {params.yaxislabel} \
        -m {output.compmat} \
        -out {output.heatmap} \

        echo "END plotHeatmap"
        echo "END plotHeatmap" 1>&2

        """

rule deeptools_fingerprint:
    input:
        expand("analysis/filt_bams/{sample.sample}_filt_alns.bam", sample=samples.itertuples())
    output:
        plot="analysis/deeptools_fingerprint/fingerprint.pdf",
        rawcts="analysis/deeptools_fingerprint/fingerprint.rawcts"
    log:
        stdout="logs/deeptools_fingerprint/combined.o",
        stderr="logs/deeptools_fingerprint/combined.e"
    benchmark:
        "benchmarks/deeptools_fingerprint/combined.txt"
    params:
        labels=lambda wildcards, input: ' '.join([os.path.basename(x).replace("_filt_alns.bam","") for x in input]),
        blacklist=blacklist,
        #sam_keep="64",
        extend_reads=lambda wildcards: "--extendReads" if all(x == "PE" for x in samples['se_or_pe'].values) else "", # extend to frag size if all samples are PE
        sam_exclude="1024",
        sam_keep=lambda wildcards: "--samFlagInclude 64" if all(x == "PE" for x in samples['se_or_pe'].values) else "", # count only first in pair if all samples are PE
        temp=tmp_dir
    threads: 16
    envmodules:
        config['modules']['deeptools']
    resources:
        mem_gb=160
    shell:
        """
        export TMPDIR={params.temp}

        plotFingerprint \
        -b {input} \
        -p {threads} \
        -o {output.plot} \
        --outRawCounts {output.rawcts} \
        {params.extend_reads} \
        --labels {params.labels} \
        --samFlagExclude {params.sam_exclude} \
        {params.sam_keep} \
        --blackListFileName {params.blacklist}

        """

def get_macs2_bams(wildcards):
    # pass deduped bams to macs2. Note that BAMs in 'filt_bams_nfr/' directory are deduped.
    if (config['atacseq'] and wildcards.macs2_type == "macs2_nfr"):
        macs2_bams = { 'trt': "analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam".format(sample=wildcards.sample) }
    elif (wildcards.macs2_type == "macs2"):
        macs2_bams = { 'trt': "analysis/dedup_bams/{sample}_filt_alns.dedup.bam".format(sample=wildcards.sample) }
    
    control = samples[samples['sample'] == wildcards.sample]['control'].values[0]
    
    if (not pd.isnull(control)):
        macs2_bams['control'] = expand("analysis/dedup_bams/{sample}_filt_alns.dedup.bam", sample = control.split(','))
    return macs2_bams

rule macs2:
    input:
        unpack(get_macs2_bams)
    output:
        multiext("analysis/{macs2_type}/{sample}_macs2_broad_peaks", ".xls", ".broadPeak"),
        multiext("analysis/{macs2_type}/{sample}_macs2_narrow_peaks", ".xls", ".narrowPeak"),
        "analysis/{macs2_type}/{sample}_macs2_narrow_summits.bed",
    log:
        stdout="logs/{macs2_type}/{sample}.o",
        stderr="logs/{macs2_type}/{sample}.e"
    benchmark:
        "benchmarks/{macs2_type}/{sample}.txt"
    params:
        control_param=lambda wildcards, input: "-c " + ' '.join(input.control) if 'control' in input.keys() else '',
        #lambda wildcards, input: print(len(input))
            #"-c {input.control}" if input.has_key('control') else '',
        species=macs2_species,
        q_cutoff="0.05",
        broad_name="{sample}_macs2_broad",
        narrow_name="{sample}_macs2_narrow",
        outdir="analysis/{macs2_type}/",
        macs2_format=lambda wildcards: "-f BAMPE" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "",
        temp_dir="tmp/"
    envmodules:
        config['modules']['macs2']
    threads: 1
    resources:
        mem_gb=100
    shell:
        """
        macs2 \
        callpeak \
        -t {input.trt} \
        {params.control_param} \
        {params.macs2_format} \
        --outdir {params.outdir} \
        -n {params.broad_name} \
        -g {params.species} \
        -q {params.q_cutoff} \
        --broad \
        --keep-dup 'all' \
        --tempdir {params.temp_dir} 

        echo "END broad peak calling" 1>&2
        echo "END broad peak calling"

        macs2 \
        callpeak \
        -t {input.trt} \
        {params.control_param} \
        {params.macs2_format} \
        --outdir {params.outdir} \
        -n {params.narrow_name} \
        -g {params.species} \
        -q {params.q_cutoff} \
        --keep-dup 'all' \
        --tempdir {params.temp_dir} 

        echo "END narrow peak calling" 1>&2
        echo "END narrow peak calling"
        
        """

rule bamtobed:
    """
    Dedup BAMs, name-sort and convert to BED for peak calling.
    """
    input:
        "analysis/dedup_bams/{sample}_filt_alns.dedup.bam"
    output:
        bam=temp("analysis/bamtobed/{sample}_filt_alns.dedup.nmsort.bam"),
        bedpe="analysis/bamtobed/{sample}.bedpe.gz",
        bed="analysis/bamtobed/{sample}.bed.gz"
    log:
        stdout="logs/bamtobed/{sample}.o",
        stderr="logs/bamtobed/{sample}.e",
    benchmark:
        "benchmarks/bamtobed/{sample}.txt"
    params:
    threads: 8
    envmodules:
        config['modules']['samtools'],
        config['modules']['bedtools']
    resources:
        mem_gb=80
    shell:
        """
        # name sort BAMs
        samtools sort -@ {threads} -n -o {output.bam} {input}

        # below adapted from https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_bam2ta.py (Latest commit 867cfe5)
        bedtools bamtobed -bedpe -mate1 -i {output.bam} | gzip -nc > {output.bedpe}

        zcat {output.bedpe} | awk 'BEGIN{{OFS="\\t"}}{{printf "%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n",$1,$2,$3,$9,$4,$5,$6,$10}}' | gzip -nc > {output.bed}


        """

rule macs2_ENCODE_atac:
    """
    For ATAC-seq, use the ENCODE method of running MACS2 via -f BED. Specific MACS2 params copied from PEPATAC. Both ENCODE and PEPATAC run only in narrow mode, but we run also in broad mode here for now to make it easier to fit in with the existing parts in this specific workflow.
    """
    input:
        trt="analysis/bamtobed/{sample}.bed.gz",
        chr_sizes="analysis/prep_chromsizes_file/chrom_sizes.tsv"
        #trt="analysis/filt_bams/{sample}_filt_alns.bam",
    output:
        broad=multiext("analysis/macs2_ENCODE_atac/{sample}_macs2_broad_peaks", ".xls", ".broadPeak", ".gappedPeak"),
        narrow=multiext("analysis/macs2_ENCODE_atac/{sample}_macs2_narrow_peaks", ".xls", ".narrowPeak"),
        narrow_summits="analysis/macs2_ENCODE_atac/{sample}_macs2_narrow_summits.bed",
        narrow_lambda_bedgraph=temp("analysis/macs2_ENCODE_atac/{sample}_macs2_narrow_control_lambda.bdg"),
        narrow_bedgraph=temp("analysis/macs2_ENCODE_atac/{sample}_macs2_narrow_treat_pileup.bdg"),
        narrow_bedgraph_clip=temp("analysis/macs2_ENCODE_atac/{sample}_macs2_narrow_treat_pileup.clip.bdg"),
        narrow_bigwig="analysis/macs2_ENCODE_atac/{sample}_macs2_narrow_treat_pileup.bw"
    log:
        stdout="logs/macs2_ENCODE_atac/{sample}.o",
        stderr="logs/macs2_ENCODE_atac/{sample}.e"
    benchmark:
        "benchmarks/macs2_ENCODE_atac/{sample}.txt"
    params:
        species=macs2_species,
        broad_name="{sample}_macs2_broad",
        narrow_name="{sample}_macs2_narrow",
        outdir="analysis/macs2_ENCODE_atac/",
        temp_dir="tmp/",
        q_cutoff="0.05"
    envmodules:
        config['modules']['macs2'],
        config['modules']['ucsc']
    threads: 1
    resources:
        mem_gb=100
    shell:
        """
        # ENCODE uses -p 0.01 but we use -q 0.05 here to be consistent with the rest of this workflow

        macs2 \
        callpeak \
        -t {input.trt} \
        -f 'BED' \
        --outdir {params.outdir} \
        -n {params.broad_name} \
        -g {params.species} \
        --shift -75 --extsize 150 \
        --nomodel \
        -q {params.q_cutoff} \
        --keep-dup all \
        --broad \
        --tempdir {params.temp_dir}

        echo "END broad peak calling" 1>&2
        echo "END broad peak calling"

        macs2 \
        callpeak \
        -t {input.trt} \
        -f 'BED' \
        --outdir {params.outdir} \
        -n {params.narrow_name} \
        -g {params.species} \
        --shift -75 --extsize 150 \
        --nomodel --call-summits \
        -q {params.q_cutoff} \
        --keep-dup all \
        -B \
        --SPMR \
        --tempdir {params.temp_dir}
        
        echo "END narrow peak calling" 1>&2
        echo "END narrow peak calling"

        # remove (-truncate leads to intervals where the end occurs before start) the intervals that extend past the end of a chromosome
        bedClip -verbose=2 {output.narrow_bedgraph} {input.chr_sizes} {output.narrow_bedgraph_clip}

        bedGraphToBigWig {output.narrow_bedgraph_clip} {input.chr_sizes} {output.narrow_bigwig} 

        echo "END narrow peak calling" 1>&2
        echo "END narrow peak calling"

        """

rule hmmratac:
    """
    Run HMMRATAC for peak calling.
    """
    input:
        "analysis/filt_bams/{sample}_filt_alns.bam",
    output:
        log="analysis/hmmratac/{sample}.log",
        genomeinfo="analysis/hmmratac/{sample}.genomeinfo",
        gappedpeak="analysis/hmmratac/{sample}_peaks.gappedPeak",
        gappedpeak_filt="analysis/hmmratac/{sample}_peaks.filteredPeaks.gappedPeak",
        gappedpeak_filt_open="analysis/hmmratac/{sample}_peaks.filteredPeaks.openOnly.bed",
        summits="analysis/hmmratac/{sample}_summits.bed",
        summits_filt="analysis/hmmratac/{sample}.filteredSummits.bed"
    log:
        stdout="logs/hmmratac/{sample}.o",
        stderr="logs/hmmratac/{sample}.e"
    benchmark:
        "benchmarks/hmmratac/{sample}.txt"
    envmodules:
        config['modules']['HMMRATAC'],
        config['modules']['samtools']
    params:
        outpref="analysis/hmmratac/{sample}",
        temp_dir="analysis/hmmratac/{sample}.tmp/",
        blacklist=config["ref"]["blacklist"],
        mito_chr=config["ref"]["mito_chr"]
    resources:
        mem_gb=120
    threads: 8
    shell:
        """
        # modified from https://github.com/LiuLabUB/HMMRATAC/blob/master/HMMRATAC_Guide.md
        # we remove the mitochondrial genome from the genominfo file in order to prevent it from being processed by HMMRATAC (see https://github.com/LiuLabUB/HMMRATAC/issues/50)
        samtools view -H {input} | perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){{print $1,"\\t",$2,"\\n"}}' | grep -Pv "^{params.mito_chr}\\t" > {output.genomeinfo}

        java -Xms80g -Xmx{resources.mem_gb}g -Djava.io.tmpdir={params.temp_dir} -jar $HMMRATAC \
        --output {params.outpref} \
        -b {input} \
        -i {input}.bai \
        -g {output.genomeinfo} \
        --blacklist {params.blacklist}

        # copied from https://github.com/LiuLabUB/HMMRATAC/blob/master/HMMRATAC_Guide.md
        awk -v OFS="\\t" '$13>=10 {{print}}' {params.outpref}_peaks.gappedPeak > {output.gappedpeak_filt}

        # custom one-liner to get only the open regions
        perl -F"\t" -lane 'print "$F[0]\\t$F[6]\\t$F[7]\\t$F[3]\\t$F[12]"' {output.gappedpeak_filt} > {output.gappedpeak_filt_open}

        # copied from https://github.com/LiuLabUB/HMMRATAC/blob/master/HMMRATAC_Guide.md
        awk -v OFS="\\t" '$5>=10 {{print}}' {params.outpref}_summits.bed > {output.summits_filt}
        """

def get_peaks_for_merging (wildcards):
    if wildcards.peak_type == "macs2_narrow":
        return(expand("analysis/macs2/{sample}_macs2_narrow_peaks.narrowPeak", sample=samples_no_controls['sample'].values))
    if wildcards.peak_type == "macs2_broad":
        return(expand("analysis/macs2/{sample}_macs2_broad_peaks.broadPeak", sample=samples_no_controls['sample'].values))
    if wildcards.peak_type == "macs2_nfr_narrow":
        return(expand("analysis/macs2_nfr/{sample}_macs2_narrow_peaks.narrowPeak", sample=samples_no_controls['sample'].values))
    if wildcards.peak_type == "macs2_nfr_broad":
        return(expand("analysis/macs2_nfr/{sample}_macs2_broad_peaks.broadPeak", sample=samples_no_controls['sample'].values))
    if wildcards.peak_type == "hmmratac_nfr":
        return(expand("analysis/hmmratac/{sample}_peaks.filteredPeaks.openOnly.bed", sample=samples_no_controls['sample'].values))
    if wildcards.peak_type == "macs2_ENCODE_atac_narrow":
        return(expand("analysis/macs2_ENCODE_atac/{sample}_macs2_narrow_peaks.narrowPeak", sample=samples_no_controls['sample'].values))
    if wildcards.peak_type == "macs2_ENCODE_atac_broad":
        return(expand("analysis/macs2_ENCODE_atac/{sample}_macs2_broad_peaks.broadPeak", sample=samples_no_controls['sample'].values))

rule merge_all_peaks:
    input:
        get_peaks_for_merging
    output:
        "analysis/merge_all_peaks/all_merged_{peak_type}.bed",
    log:
        stdout="logs/merge_all_peaks/{peak_type}.o",
        stderr="logs/merge_all_peaks/{peak_type}.e"
    benchmark:
        "benchmarks/merge_all_peaks/{peak_type}.txt"
    params:
    envmodules:
        config['modules']['bedops']
    threads: 4
    resources:
        mem_gb=100
    shell:
        """
        bedops --merge {input} 1> {output}
        """


rule rm_blacklist_peaks:
    """
    Remove blacklist regions from macs2 peaks using overlap of 1bp or more. In addition, remove nonstandard chromosomes.
    """
    input:
        broad="analysis/{macs2_type}/{sample}_macs2_broad_peaks.broadPeak",
        narrow="analysis/{macs2_type}/{sample}_macs2_narrow_peaks.narrowPeak",
        narrow_summits="analysis/{macs2_type}/{sample}_macs2_narrow_summits.bed",
        std_chroms="analysis/misc/std_chroms.txt"
    output:
        broad="analysis/{macs2_type}/rm_blacklist/{sample}_macs2_broad_peaks.rm_blacklist.broadPeak",
        narrow="analysis/{macs2_type}/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak",
        narrow_summits="analysis/{macs2_type}/rm_blacklist/{sample}_macs2_narrow_summits.rm_blacklist.bed"
    log:
        stdout="logs/{macs2_type}/rm_blacklist/{sample}.o",
        stderr="logs/{macs2_type}/rm_blacklist/{sample}.e"
    benchmark:
        "benchmarks/{macs2_type}/rm_blacklist/{sample}.txt"
    params:
        blacklist=blacklist,
    envmodules:
        config['modules']['bedtools'],
        config['modules']['R']
    threads: 4
    resources:
        mem_gb=100
    shell:
        """
        # below note the '||' condition for the grep statement that allows a return code of 0 even when no matches are found (as in the case of an empty bed file)
        bedtools intersect -v \
        -a {input.narrow} \
        -b {params.blacklist} | {{ grep -P "$(cat {input.std_chroms} | perl -lane 'print q:^:.join(q:\\t|^:, @F).q:\\t:')" || [[ $? == 1 ]]; }} > {output.narrow} 
    
        # subset the summits also based on the names of the retained narrow peaks
        Rscript -e 'library(magrittr); library(rtracklayer); keep_pks <- import("{output.narrow}")$name; gr <- import("{input.narrow_summits}"); gr[gr$name %in% keep_pks] %>% export(., "{output.narrow_summits}")'

        bedtools intersect -v \
        -a {input.broad} \
        -b {params.blacklist} | {{ grep -P "$(cat {input.std_chroms} | perl -lane 'print q:^:.join(q:\\t|^:, @F).q:\\t:')" || [[ $? == 1 ]]; }} > {output.broad}



        """

#def get_peaks_for_frip (wildcards):
#    if frip_peakset=="narrow":
#        return("analysis/macs2/{sample}_macs2_narrow_peaks.narrowPeak")
#    elif frip_peakset=="broad":
#        return("analysis/macs2/{sample}_macs2_broad_peaks.broadPeak")

rule diffbind_count:
    input:
        samplesheet=samplesheet,
        bams=expand("analysis/filt_bams/{sample}_filt_alns.bam", sample=samples['sample']),
        peaks=expand("analysis/macs2/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak", sample=samples_no_controls['sample']),

    output:
        outrds="analysis/diffbind_count/{enriched_factor}.rds",
        samplesheet="analysis/diffbind_count/{enriched_factor}_DB_samplesheet.tsv"
    log:
        stdout="logs/diffbind_count/{enriched_factor}.o",
        stderr="logs/diffbind_count/{enriched_factor}.e",
    benchmark:
        "benchmarks/diffbind_count/{enriched_factor}.txt"
    params:
        outdir="analysis/diffbind_count/",
        DB_summits=config['DiffBind']['summits'],
        enriched_factor="{enriched_factor}",
        subtract_controls=config['DiffBind']['subtract_controls'],
        macs2_type="macs2" if not config['atacseq'] else "macs2_ENCODE_atac"
    threads: 16
    resources:
        mem_gb=196
    envmodules:
        config['modules']['R']
    script:
        "bin/scripts/diffbind_count.R"


rule homer_find_motif:
    input:
        "analysis/macs2/rm_blacklist/{sample}_macs2_narrow_summits.rm_blacklist.bed"
    output:
        "analysis/homer_find_motifs/{sample}/homerMotifs.all.motifs"
    log:
        stdout="logs/homer_find_motif/{sample}.o",
        stderr="logs/homer_find_motif/{sample}.e"
    benchmark:
        "benchmarks/homer_find_motif/{sample}.txt"
    params:
        outdir="analysis/homer_find_motifs/{sample}/",
        genome=config['homer']['genome'],
        size=config['homer']['size']
    threads: 8
    envmodules:
        "bbc/HOMER/HOMER-4.11.1"
    resources:
        mem_gb=200
    shell:
        """
        findMotifsGenome.pl {input} {params.genome} {params.outdir} \
        -size {params.size} -p {threads} \
        -mask

        """

def get_peaks_for_venn (wildcards):
    peaks = {
            "broad": expand("analysis/macs2/rm_blacklist/{sample}_macs2_broad_peaks.rm_blacklist.broadPeak", sample=samples_no_controls['sample'].values),
            "narrow": expand("analysis/macs2/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak", sample=samples_no_controls['sample'].values)    
    }
    if (config['atacseq']):
        peaks['broad_nfr'] = expand("analysis/macs2_nfr/rm_blacklist/{sample}_macs2_broad_peaks.rm_blacklist.broadPeak", sample=samples_no_controls['sample'].values)
        peaks['narrow_nfr'] = expand("analysis/macs2_nfr/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak", sample=samples_no_controls['sample'].values)
        peaks['broad_ENCODE_atac'] = expand("analysis/macs2_ENCODE_atac/rm_blacklist/{sample}_macs2_broad_peaks.rm_blacklist.broadPeak", sample=samples_no_controls['sample'].values)
        peaks['narrow_ENCODE_atac'] = expand("analysis/macs2_ENCODE_atac/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak", sample=samples_no_controls['sample'].values)
        if (config['run_hmmratac']):
            peaks['hmmratac_nfr'] = expand("analysis/hmmratac/{sample}_peaks.filteredPeaks.openOnly.bed", sample=samples_no_controls['sample'].values)

    return peaks

rule peaks_venn:
    """
    Make Venn diagrams for peaksets of each peak type. Note that only first 5 samples are plotted.
    """
    input:
        unpack(get_peaks_for_venn)
    output:
        "analysis/peaks_venn/report.html"
    log:
        stdout="logs/peaks_venn/out.o",
        stderr="logs/peaks_venn/err.e"
    benchmark:
        "benchmarks/peaks_venn/benchmark.txt"
    params:
        out_dir="analysis/peaks_venn/peaks_venn_out_files/",
        sample_names=samples_no_controls['sample'].values,
    envmodules:
        config['modules']['R']
    threads: 1
    resources:
        mem_gb = 60
    script:
        "bin/peaks_venn.Rmd"

def get_merged_beds (wildcards):
    atac_peak_types = ['macs2_ENCODE_atac_broad','macs2_ENCODE_atac_narrow','macs2_nfr_broad','macs2_nfr_narrow','macs2_narrow','macs2_broad']
    if config['run_hmmratac']:
        atac_peak_types = atac_peak_types + ['hmmratac_nfr']# if config['run_hmmratac']
    merged_beds = expand("analysis/merge_all_peaks/all_merged_{peak_type}.bed", peak_type=atac_peak_types if config['atacseq'] else ['macs2_narrow','macs2_broad'])
    return (merged_beds)

rule deeptools_plotenrichment:
    input:
        bam="analysis/filt_bams/{sample}_filt_alns.bam",
        bed=expand("analysis/{macs2_type}/{{sample}}_macs2_{type}_peaks.{type}Peak", type=['narrow','broad'], macs2_type=['macs2','macs2_nfr','macs2_ENCODE_atac'] if config['atacseq'] else['macs2']),
        merged_beds=get_merged_beds,
        #expand("analysis/merge_all_peaks/all_merged_{peak_type}.bed", peak_type=['macs2_narrow','macs2_broad', 'macs2_nfr_narrow','macs2_nfr_broad','hmmratac_nfr'] if config['atacseq'] else ['macs2_narrow','macs2_broad'])],
    output:
        #merged_peaks="analysis/deeptools_plotenrichment/{sample}.narrowPeak",
        plot="analysis/deeptools_plotenrichment/{sample}.pdf",
        rawcts="analysis/deeptools_plotenrichment/{sample}.rawcts"
    log:
        stdout="logs/deeptools_plotenrichment/{sample}.o",
        stderr="logs/deeptools_plotenrichment/{sample}.e"
    benchmark:
        "benchmarks/deeptools_plotenrichment/{sample}.txt"
    params:
        blacklist=blacklist,
        temp=tmp_dir,
        extend_reads=lambda wildcards: "--extendReads" if all(x == "PE" for x in samples['se_or_pe'].values) else "", # extend to frag size if all samples are PE
        sam_keep=lambda wildcards: "--samFlagInclude 64" if all(x == "PE" for x in samples['se_or_pe'].values) else "", # count only first in pair if all samples are PE
        sam_exclude="1024",    
        #sam_keep="64",
    threads: 16
    envmodules:
        config['modules']['bedops'],
        config['modules']['deeptools']
    resources:
        mem_gb=100
    shell:
        """
        export TMPDIR={params.temp}

        plotEnrichment \
        --bamfiles {input.bam} \
        --BED {input.bed} {input.merged_beds} \
        --smartLabels \
        --variableScales \
        --outRawCounts {output.rawcts} \
        --blackListFileName {params.blacklist} \
        {params.extend_reads} \
        --samFlagExclude {params.sam_exclude} \
        {params.sam_keep} \
        -p {threads} \
        -o {output.plot} 

        """

rule rm_supplementary_alns:
    """
    Run samtools view to remove supplementary alignments (no bit set in SAM flag 2304 / -F 2304).
    """
    input:
        "analysis/bwamem/{bam_name}.bam"
    output:
        temp("analysis/rm_supplementary_alns/{bam_name}.F2304.bam")
    params:
    log:
        stdout="logs/rm_supplementary_alns/{bam_name}.o",
        stderr="logs/rm_supplementary_alns/{bam_name}.e"
    benchmark:
        "benchmarks/rm_supplementary_alns/{bam_name}.txt"
    envmodules:
        config['modules']['samtools']
    threads: 8
    resources:
        mem_gb = 32
    shell:
        """
        samtools view -@ {threads} -F 2304 -b -o {output} {input}
        """

rule preseq_complexity:
    """
    Run preseq c_curve and lc_extrap on the BAMs after removal of supplementary alignments.
    """
    input:
        "analysis/rm_supplementary_alns/{sample}.F2304.bam",
    output:
        ccurve="analysis/preseq_complexity/{sample}.c_curve.txt",
        lcextrap="analysis/preseq_complexity/{sample}.lc_extrap.txt"
    log:
        stdout="logs/preseq_complexity/{sample}.o",
        stderr="logs/preseq_complexity/{sample}.e"
    benchmark:
        "benchmarks/preseq_complexity/{sample}.txt"
    envmodules:
        config['modules']['preseq']
    params:
        paired=lambda wildcards: "-P" if samples[samples['sample']==wildcards.sample]['se_or_pe'].values=="PE" else "",
    resources:
        mem_gb=100
    threads: 8
    shell:
        """
        # preseq doesn't process supplmentary alignments properly. Remove these using samtools.

        preseq c_curve \
        -v \
        {params.paired} \
        -bam \
        -o {output.ccurve} \
        {input}

        echo "Finished c_curve." >&1
        echo "Finished c_curve." >&2

        preseq lc_extrap \
        -v \
        {params.paired} \
        -bam \
        -o {output.lcextrap} \
        {input}

        echo "Finished lc_extrap." >&1
        echo "Finished lc_extrap." >&2

        """

rule qualimap:
    """
    Run qualimap on unfiltered alignments.
    """
    input:
        "analysis/bwamem/{sample}.bam",
    output:
        touch("analysis/qualimap/{sample}/done")
    log:
        stdout="logs/qualimap/{sample}.o",
        stderr="logs/qualimap/{sample}.e"
    benchmark:
        "benchmarks/qualimap/{sample}.txt"
    envmodules:
        config['modules']['qualimap']
    params:
    resources:
        mem_gb=100
    threads: 8
    shell:
        """
        qualimap bamqc -bam {input} --java-mem-size={resources.mem_gb}G --paint-chromosome-limits -outdir analysis/qualimap/{wildcards.sample} -nt {threads}

        """

rule multiqc:
    input:
        expand("analysis/fastqc/{sample.sample}_R1_fastqc.html", sample=samples.itertuples()),
        expand("analysis/fastqc/{sample.sample}_R2_fastqc.html", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/trim_galore/{sample.sample}_R{read}_val_{read}_fastqc.html", sample=samples[samples['se_or_pe']=="PE"].itertuples(), read=["1","2"]),
        expand("analysis/trim_galore/{sample.sample}_R1_trimmed_fastqc.html", sample=samples[samples['se_or_pe']=="SE"].itertuples()),
        expand("analysis/preseq_complexity/{sample.sample}.c_curve.txt", sample=samples.itertuples()),
        expand("analysis/preseq_complexity/{sample.sample}.lc_extrap.txt", sample=samples.itertuples()),
        expand("analysis/bwamem/flagstat/{sample.sample}.flagstat", sample=samples.itertuples()),
        expand("analysis/bwamem/CollectInsertSizeMetrics/{sample.sample}.insert_size_metrics.txt", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/bwamem/{sample.sample}.bam.idxstat", sample=samples.itertuples()),
        expand("analysis/filt_bams/{sample.sample}_filt_alns.bam.idxstat", sample=samples.itertuples()),
        expand("analysis/filt_bams/flagstat/{sample.sample}_filt_alns.flagstat", sample=samples.itertuples()),
        expand("analysis/filt_bams/CollectInsertSizeMetrics/{sample.sample}_filt_alns.insert_size_metrics.txt", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/dedup_bams/flagstat/{sample.sample}_filt_alns.dedup.flagstat", sample=samples.itertuples()),
        expand("analysis/bwamem/CollectAlignmentSummaryMetrics/{sample.sample}.aln_metrics.txt", sample=samples.itertuples()),
        expand("analysis/bwamem/{sample.sample}.samblaster.e", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R1_screen.html", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R2_screen.html", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/deeptools_plotenrichment/{sample.sample}.rawcts", sample=samples_no_controls.itertuples()),
        expand("analysis/qualimap/{sample.sample}/done", sample=samples.itertuples()),
    output:
        "analysis/multiqc/multiqc_report.html",
    log:
        stdout="logs/multiqc/multiqc.o",
        stderr="logs/multiqc/multiqc.e"
    benchmark:
        "benchmarks/multiqc/multiqc.txt"
    params:
        workdir="analysis/multiqc/",
        dirs=["analysis/trim_galore/",
        "analysis/fastqc/",
        "analysis/qualimap/",
        "analysis/preseq_complexity/",
        #"analysis/fastp/",
        "analysis/bwamem/flagstat/",
        "analysis/bwamem/",
        "analysis/fastq_screen/",
        "analysis/bwamem/CollectAlignmentSummaryMetrics/",
        "analysis/filt_bams/",
        "analysis/filt_bams/flagstat/",
        "analysis/dedup_bams/flagstat/",
        "analysis/deeptools_plotenrichment/"],
        PE_dirs=["analysis/filt_bams/CollectInsertSizeMetrics/"] if not all(x == "SE" for x in samples['se_or_pe'].values)  else [],
        outfile="multiqc_report"
    envmodules:
        config['modules']['multiqc']
    threads: 4
    resources:
        mem_gb=100
    shell:
        """
        multiqc \
        --force \
        --config bin/multiqc_config.yaml \
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs} {params.PE_dirs} analysis/macs2/*narrow*

        """

