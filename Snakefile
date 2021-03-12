import pandas as pd
import numpy as np
import os
import re
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.28.0")

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
#chroms_no_mito = ' '.join(fai_parsed[fai_parsed['chr'] != mito_chrom]['chr'].values)
chroms_gt_5Mb = ' '.join(fai_parsed[fai_parsed['len'] > 5000000]['chr'].values)

##### load config and sample sheets #####

samples = pd.read_table("bin/samples.tsv")

# Filter for sample rows that are not controls
samples_no_controls = samples[-samples['sample'].isin(samples['control'].values)]

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
        #expand("analysis/filt_bams/{sample.sample}_filt_alns.bam", sample=samples.itertuples()),
        expand("analysis/deeptools_cov_keepdups/{sample.sample}_filt_alns_keepdups.bw", sample=samples.itertuples()),
        "analysis/multiqc/multiqc_report.html",
        "analysis/deeptools_fingerprint/fingerprint.rawcts",
        #expand("analysis/peaks_rm_blacklist/{sample}_macs2_{size}_peaks.rm_blacklist.{size}Peak", sample=samples["sample"], size=["narrow","broad"]),
        #expand("analysis/hmmratac/{sample}_summits.bed", sample=samples["sample"]) if config['atacseq'] else [],
        expand("analysis/deeptools_plotenrichment/{sample}.pdf", sample=samples["sample"]),
        expand("analysis/filt_bams_nfr/CollectInsertSizeMetrics/{sample.sample}_filt_alns_nfr.insert_size_metrics.txt", sample=samples.itertuples()) if config['atacseq'] else [],
        "analysis/peaks_venn/report.html",

def get_orig_fastq(wildcards):
    if wildcards.read == "R1":
            fastq = expand("raw_data/{fq}", fq = samples[samples["sample"] == wildcards.sample]["fq1"].values)
    elif wildcards.read == "R2":
            fastq = expand("raw_data/{fq}", fq = samples[samples["sample"] == wildcards.sample]["fq2"].values)
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
        "bbc/fastqc/fastqc-0.11.9"
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
        "bbc/fastq_screen/fastq_screen-0.14.0"
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
        expand("analysis/trim_galore/{{sample}}_R1{ext}", ext=["_val_1.fq.gz",".fastq.gz_trimming_report.txt","_val_1_fastqc.html"]),
        expand("analysis/trim_galore/{{sample}}_R2{ext}", ext=["_val_2.fq.gz",".fastq.gz_trimming_report.txt","_val_2_fastqc.html"]),
    params:
    log:
        stdout="logs/trim_galore/{sample}.o",
        stderr="logs/trim_galore/{sample}.e"
    benchmark:
        "benchmarks/trim_galore/{sample}.txt"
    envmodules:
        "bbc/trim_galore/trim_galore-0.6.0"
    threads: 4
    resources:
        mem_gb = 64
    shell:
        """
        trim_galore --paired {input} --output_dir analysis/trim_galore/ --cores {threads} --fastqc
        """

rule bwamem:
    input:
        expand("analysis/trim_galore/{{sample}}_R{read}_val_{read}.fq.gz", read=["1","2"])
    output:
        outbam="analysis/bwamem/{sample}.bam",
        outbai="analysis/bwamem/{sample}.bam.bai",
        idxstat="analysis/bwamem/{sample}.bam.idxstat"
    log:
        stdout="logs/bwamem/{sample}.o",
        stderr="logs/bwamem/{sample}.e",
    benchmark:
        "benchmarks/bwamem/{sample}.txt"
    params:
        bwa_idx=bwa_index,
    threads: 16
    envmodules:
        "bbc/bwa/bwa-0.7.17",
        "bbc/samblaster/samblaster-0.1.24",
        "bbc/samtools/samtools-1.9",
    resources:
        mem_gb=180
    shell:
        """
        bwa mem \
        -t {threads} \
        {params.bwa_idx} \
        {input} | \
        samblaster | \
        samtools sort \
        -m 6G \
        -@ {threads} \
        -O "BAM" \
        -o {output.outbam} \
        -

        echo "END bwamem" >> {log.stdout}
        echo "END bwamem" >> {log.stderr}

        samtools index -@ {threads} {output.outbam}

        echo "END indexing" >> {log.stdout}
        echo "END indexing" >> {log.stderr}
        
        samtools idxstats {output.outbam} > {output.idxstat}

        echo "END idxstats" >> {log.stdout}
        echo "END idxstats" >> {log.stderr}
 
        """

rule filt_bams:
    input:
        "analysis/bwamem/{sample}.bam"
    output:
        bam="analysis/filt_bams/{sample}_filt_alns.bam",
        bai="analysis/filt_bams/{sample}_filt_alns.bam.bai",
        idxstat="analysis/filt_bams/{sample}_filt_alns.bam.idxstat"
    log:
        stdout="logs/filt_bams/{sample}.o",
        stderr="logs/filt_bams/{sample}.e",
    benchmark:
        "benchmarks/filt_bams/{sample}.txt"
    params:
        mapq=20,
        flags_to_exclude="2308",
        keep_chroms = chroms_gt_5Mb #chroms_no_mito if atacseq else ''
    threads: 8
    resources:
        mem_gb=80
    envmodules:
        "bbc/samtools/samtools-1.9"
    shell:
        """
        samtools view \
              -@ {threads} \
              -q {params.mapq} \
              -F {params.flags_to_exclude} \
              -b -o {output.bam} {input} {params.keep_chroms}

        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.idxstat}

        """

rule filt_bams_nfr:
    """
    Keep only paired alignments with <150 bp fragment size.
    """
    input:
        "analysis/filt_bams/{sample}_filt_alns.bam"
    output:
        bam="analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam",
        bai="analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam.bai",
        idxstat="analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam.idxstat",
        metrics="analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam.metrics"
    log:
        stdout="logs/filt_bams_nfr/{sample}.o",
        stderr="logs/filt_bams_nfr/{sample}.e",
    benchmark:
        "benchmarks/filt_bams_nfr/{sample}.txt"
    params:
        temp=os.path.join(snakemake_dir, "tmp"),
        maxFragmentSize = 150,
    threads: 8
    envmodules:
        "bbc/deeptools/deeptools-3.4.3",
        "bbc/samtools/samtools-1.9"
    resources:
        mem_gb=80
    shell:
        """
        export TMPDIR={params.temp}
        
        alignmentSieve --bam {input} \
        --outFile {output.bam} \
        -p {threads} \
        --filterMetrics {output.metrics} \
        --maxFragmentLength {params.maxFragmentSize} \
        2> {log.stderr} 
        
        samtools index {output.bam}
        samtools idxstats {output.bam} > {output.idxstat}

        """


rule deeptools_cov_keepdups:
    input:
        bam="analysis/filt_bams/{sample}_filt_alns.bam"
    output:
        bigwig_keepdups="analysis/deeptools_cov_keepdups/{sample}_filt_alns_keepdups.bw"
    log:
        stdout="logs/deeptools_cov_keepdups/{sample}.o",
        stderr="logs/deeptools_cov_keepdups/{sample}.e"
    benchmark:
        "benchmarks/deeptools_cov_keepdups/{sample}.txt"
    params:
        blacklist=blacklist,
        binsize=bamCoverage_binsize,
        norm_method="CPM",
        sam_keep="64",
        temp=os.path.join(snakemake_dir, "tmp")
    threads: 16
    envmodules:
        "bbc/deeptools/deeptools-3.4.3"
    resources:
        mem_gb=96
    shell:
        """
        export TMPDIR={params.temp}
        
        # calculate the coverage
        ## Since this is PE data, we can use '--extendReads' to extend the reads to connect the two mates.
        bamCoverage \
        -p {threads} \
        --binSize {params.binsize} \
        --bam {input.bam} \
        --extendReads \
        --blackListFileName {params.blacklist} \
        -o {output.bigwig_keepdups} \
        --binSize {params.binsize} \
        --normalizeUsing {params.norm_method} \
        --samFlagInclude {params.sam_keep}

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
        labels=' '.join(samples["sample"].values),
        blacklist=blacklist,
        sam_keep="64",
        temp=tmp_dir
    threads: 16
    envmodules:
        "bbc/deeptools/deeptools-3.4.3"
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
        --extendReads \
        --minMappingQuality 20 \
        --labels {params.labels} \
        --samFlagInclude {params.sam_keep} \
        --blackListFileName {params.blacklist}

        """

def get_macs2_bams(wildcards):
    if (config['atacseq'] and wildcards.macs2_type == "macs2_nfr"):
        macs2_bams = { 'trt': "analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam".format(sample=wildcards.sample) }
    elif (wildcards.macs2_type == "macs2"):
        macs2_bams = { 'trt': "analysis/filt_bams/{sample}_filt_alns.bam".format(sample=wildcards.sample) }
    
    control = samples[samples['sample'] == wildcards.sample]['control'].values[0]
    
    if (not pd.isnull(control)):
        macs2_bams['control'] = "analysis/filt_bams/{sample}_filt_alns.bam".format(sample=control)
    return macs2_bams

rule macs2:
    input:
        unpack(get_macs2_bams)
        #trt="analysis/filt_bams/{sample}_filt_alns.bam",
        #control=get_macs2_control,
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
        control_param=lambda wildcards, input: "-c "+input.control if len(input) > 1 else '',
        #lambda wildcards, input: print(len(input))
            #"-c {input.control}" if input.has_key('control') else '',
        species=macs2_species,
        q_cutoff="0.05",
        broad_name="{sample}_macs2_broad",
        narrow_name="{sample}_macs2_narrow",
        outdir="analysis/{macs2_type}/",
        temp_dir="tmp/"
    envmodules:
        "bbc/macs2/macs2-2.2.7.1"
    threads: 1
    resources:
        mem_gb=100
    shell:
        """
        macs2 \
        callpeak \
        -t {input.trt} \
        {params.control_param} \
        -f BAMPE \
        --outdir {params.outdir} \
        -n {params.broad_name} \
        -g {params.species} \
        -q {params.q_cutoff} \
        --broad \
        --keep-dup 1 \
        --tempdir {params.temp_dir} 

        echo "END broad peak calling" >> {log.stdout}
        echo "END broad peak calling" >> {log.stderr}

        macs2 \
        callpeak \
        -t {input.trt} \
        {params.control_param} \
        -f BAMPE \
        --outdir {params.outdir} \
        -n {params.narrow_name} \
        -g {params.species} \
        -q {params.q_cutoff} \
        --keep-dup 1 \
        --tempdir {params.temp_dir} 

        echo "END narrow peak calling" >> {log.stdout}
        echo "END narrow peak calling" >> {log.stderr}
        
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
        "bbc/HMMRATAC/HMMRATAC-1.2.10",
        "bbc/samtools/samtools-1.9"
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
        perl -lane 'print "$F[0]\\t$F[6]\\t$F[7]\\t$F[3]\\t$F[12]"' {output.gappedpeak_filt} > {output.gappedpeak_filt_open}

        # copied from https://github.com/LiuLabUB/HMMRATAC/blob/master/HMMRATAC_Guide.md
        awk -v OFS="\\t" '$5>=10 {{print}}' {params.outpref}_summits.bed > {output.summits_filt}
        """

def get_peaks_for_merging (wildcards):
    if wildcards.peak_type == "macs2_narrow":
        return(expand("analysis/macs2/{sample}_macs2_narrow_peaks.narrowPeak", sample=samples['sample'].values))
    if wildcards.peak_type == "macs2_broad":
        return(expand("analysis/macs2/{sample}_macs2_broad_peaks.broadPeak", sample=samples['sample'].values))
    if wildcards.peak_type == "macs2_nfr_narrow":
        return(expand("analysis/macs2_nfr/{sample}_macs2_narrow_peaks.narrowPeak", sample=samples['sample'].values))
    if wildcards.peak_type == "macs2_nfr_broad":
        return(expand("analysis/macs2_nfr/{sample}_macs2_broad_peaks.broadPeak", sample=samples['sample'].values))
    if wildcards.peak_type == "hmmratac_nfr":
        return(expand("analysis/hmmratac/{sample}_peaks.filteredPeaks.openOnly.bed", sample=samples['sample'].values))


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
        "bbc/bedops/bedops-2.4.37"
    threads: 4
    resources:
        mem_gb=100
    shell:
        """
        bedops --merge {input} 1> {output} 2> {log.stderr}
        """


rule rm_blacklist_peaks:
    """
    Remove blacklist regions from macs2 peaks using overlap of 1bp or more
    """
    input:
        broad="analysis/{macs2_type}/{sample}_macs2_broad_peaks.broadPeak",
        narrow="analysis/{macs2_type}/{sample}_macs2_narrow_peaks.narrowPeak"
    output:
        broad="analysis/{macs2_type}/rm_blacklist/{sample}_macs2_broad_peaks.rm_blacklist.broadPeak",
        narrow="analysis/{macs2_type}/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak"
    log:
        stdout="logs/{macs2_type}/rm_blacklist/{sample}.o",
        stderr="logs/{macs2_type}/rm_blacklist/{sample}.e"
    benchmark:
        "benchmarks/{macs2_type}/rm_blacklist/{sample}.txt"
    params:
        blacklist=blacklist,
    envmodules:
        "bbc/bedtools/bedtools-2.29.2"
    threads: 4
    resources:
        mem_gb=100
    shell:
        """
        cat {input.narrow} | \
        bedtools intersect -v \
        -a "stdin" \
        -b {params.blacklist} > {output.narrow} 

        cat {input.broad} | \
        bedtools intersect -v \
        -a "stdin" \
        -b {params.blacklist} > {output.broad}

        """

#def get_peaks_for_frip (wildcards):
#    if frip_peakset=="narrow":
#        return("analysis/macs2/{sample}_macs2_narrow_peaks.narrowPeak")
#    elif frip_peakset=="broad":
#        return("analysis/macs2/{sample}_macs2_broad_peaks.broadPeak")

def get_peaks_for_venn (wildcards):
    peaks = {
            "broad": expand("analysis/macs2/rm_blacklist/{sample}_macs2_broad_peaks.rm_blacklist.broadPeak", sample=samples_no_controls['sample'].values),
            "narrow": expand("analysis/macs2/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak", sample=samples_no_controls['sample'].values)    
    }
    if (config['atacseq']):
        peaks['broad_nfr'] = expand("analysis/macs2_nfr/rm_blacklist/{sample}_macs2_broad_peaks.rm_blacklist.broadPeak", sample=samples_no_controls['sample'].values)
        peaks['narrow_nfr'] = expand("analysis/macs2_nfr/rm_blacklist/{sample}_macs2_narrow_peaks.rm_blacklist.narrowPeak", sample=samples_no_controls['sample'].values)
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
        "bbc/R/R-4.0.2-setR_LIBS_USER"
    threads: 1
    resources:
        mem_gb = 60
    script:
        "bin/peaks_venn.Rmd"

rule deeptools_plotenrichment:
    input:
        bam="analysis/filt_bams/{sample}_filt_alns.bam",
        bed=[expand("analysis/{macs2_type}/{{sample}}_macs2_{type}_peaks.{type}Peak", type=['narrow','broad'], macs2_type=['macs2','macs2_nfr'] if config['atacseq'] else['macs2']),
        expand("analysis/merge_all_peaks/all_merged_{peak_type}.bed", peak_type=['macs2_narrow','macs2_broad', 'macs2_nfr_narrow','macs2_nfr_broad','hmmratac_nfr'] if config['atacseq'] else ['macs2'])],
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
        sam_keep="64",
        #sam_exclude="1024",    
    threads: 16
    envmodules:
        "bbc/bedops/bedops-2.4.37",
        "bbc/deeptools/deeptools-3.4.3"
    resources:
        mem_gb=100
    shell:
        """
        export TMPDIR={params.temp}

        plotEnrichment \
        --bamfiles {input.bam} \
        --BED {input.bed} \
        --smartLabels \
        --variableScales \
        --outRawCounts {output.rawcts} \
        --samFlagInclude {params.sam_keep} \
        --blackListFileName {params.blacklist} \
        -p {threads} \
        -o {output.plot} 

        """

rule rm_supplementary_alns:
    """
    Run samtools view to remove supplementary alignments (no bit set in SAM flag 2048 / -F 2048).
    """
    input:
        "analysis/bwamem/{bam_name}.bam"
    output:
        temp("analysis/rm_supplementary_alns/{bam_name}.F2048.bam")
    params:
    log:
        stdout="logs/rm_supplementary_alns/{bam_name}.o",
        stderr="logs/rm_supplementary_alns/{bam_name}.e"
    benchmark:
        "benchmarks/rm_supplementary_alns/{bam_name}.txt"
    envmodules:
        "bbc/samtools/samtools-1.9"
    threads: 8
    resources:
        mem_gb = 32
    shell:
        """
        samtools view -@ {threads} -F 2048 -b -o {output} {input}
        """

rule preseq_complexity:
    """
    Run preseq c_curve and lc_extrap on the BAMs after removal of supplementary alignments.
    """
    input:
        "analysis/rm_supplementary_alns/{sample}.F2048.bam",
    output:
        ccurve="analysis/preseq_complexity/{sample}.c_curve.txt",
        lcextrap="analysis/preseq_complexity/{sample}.lc_extrap.txt"
    log:
        stdout="logs/preseq_complexity/{sample}.o",
        stderr="logs/preseq_complexity/{sample}.e"
    benchmark:
        "benchmarks/preseq_complexity/{sample}.txt"
    envmodules:
        "bbc/preseq/preseq-2.0.3"
    params:
    resources:
        mem_gb=100
    threads: 8
    shell:
        """
        # preseq doesn't process supplmentary alignments properly. Remove these using samtools.

        preseq c_curve \
        -v \
        -P \
        -bam \
        -o {output.ccurve} \
        {input}

        echo "Finished c_curve." >&1
        echo "Finished c_curve." >&2

        preseq lc_extrap \
        -v \
        -P \
        -bam \
        -o {output.lcextrap} \
        {input}

        echo "Finished lc_extrap." >&1
        echo "Finished lc_extrap." >&2

        """


rule multiqc:
    input:
        expand("analysis/fastqc/{sample.sample}_R{read}_fastqc.html", sample=samples.itertuples(), read=["1","2"]),
        expand("analysis/trim_galore/{sample.sample}_R{read}_val_{read}_fastqc.html", sample=samples.itertuples(), read=["1","2"]),
        expand("analysis/preseq_complexity/{sample.sample}.c_curve.txt", sample=samples.itertuples()),
        expand("analysis/preseq_complexity/{sample.sample}.lc_extrap.txt", sample=samples.itertuples()),
        #expand("analysis/fastp/{sample.sample}_fastp.html", sample=samples.itertuples()),
        expand("analysis/bwamem/flagstat/{sample.sample}.flagstat", sample=samples.itertuples()),
        expand("analysis/filt_bams/{sample.sample}_filt_alns.bam.idxstat", sample=samples.itertuples()),
        expand("analysis/filt_bams/CollectInsertSizeMetrics/{sample.sample}_filt_alns.insert_size_metrics.txt", sample=samples.itertuples()),
        expand("analysis/bwamem/CollectAlignmentSummaryMetrics/{sample.sample}.aln_metrics.txt", sample=samples.itertuples()),
        #expand("analysis/align/{sample.sample}_vs_mm10.bam.flagstat", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R{read}_screen.html", sample=samples.itertuples(), read=["1","2"]),
        expand("analysis/deeptools_plotenrichment/{sample.sample}.rawcts", sample=samples.itertuples()),
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
        "analysis/preseq_complexity/",
        #"analysis/fastp/",
        "analysis/bwamem/flagstat/",
        "analysis/fastq_screen/",
        "analysis/bwamem/CollectAlignmentSummaryMetrics/",
        "analysis/filt_bams/",
        "analysis/filt_bams/CollectInsertSizeMetrics/",
        "analysis/deeptools_plotenrichment/"],
        outfile="multiqc_report"
    envmodules:
        "bbc/multiqc/multiqc-1.9"
    threads: 4
    resources:
        mem_gb=100
    shell:
        """
        multiqc \
        --force \
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs}

        """

