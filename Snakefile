import pandas as pd
import numpy as np
import os
import re
import itertools
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
if mito_chrom not in fai_parsed['chr'].values:
    raise Exception('{mito} not found in reference fai file.'.format(mito=mito_chrom))

chroms_no_mito = ' '.join(fai_parsed[fai_parsed['chr'] != mito_chrom]['chr'].values)

# When we filter alignments for peak calling, we remove contigs smaller than 'cutoff'. This is meant to be a way to retain only the nuclear chromosomes (not mitochondrial).
#chrom_min_bp = config['ref']['chrom_min_bp']
#chroms_gt_cutoff = ' '.join(fai_parsed[fai_parsed['len'] > chrom_min_bp]['chr'].values)

##### load config and sample sheets #####

units = pd.read_table("bin/samples.tsv")
units['se_or_pe'] = ["SE" if x else "PE" for x in units['fq2'].isnull()]

samples = units[["sample","control","sample_group","se_or_pe"]].drop_duplicates()
if not samples['sample'].is_unique:
    raise Exception('A sample has more than one combination of control, smaple_group and/or se_or_pe.')

# Filter for sample rows that are not controls
controls_list = list(itertools.chain.from_iterable( [x.split(',') for x in samples['control'].values if not pd.isnull(x)] ))
samples_no_controls = samples[-samples['sample'].isin(controls_list)]

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
        outbam="analysis/bwamem/{sample}.bam",
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
    """
    Filter for mapq score, properly paired, primary alignment, remove unmapped. Do not remove dup reads.
    """
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
        mapq=30,
        flags_to_exclude="780",
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
        --binSize {params.binsize} \
        --normalizeUsing {params.norm_method} \
        --samFlagExclude {params.sam_exclude} \
        {params.sam_keep}

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
        bedgraph=temp("analysis/merge_bigwigs/{group}.bedGraph"),
        sorted_bg=temp("analysis/merge_bigwigs/{group}.sort.bedGraph"),
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
        bigWigMerge {in_bw} {out_bg}

        sort -k1,1 -k2,2n {out_bg} > {out_sorted_bg}

        # convert to bigwig
        bedGraphToBigWig {out_sorted_bg} {in_chroms} {out_bw}

        """ if len(input.bigwigs) > 1 else """
        ln -rs {in_bw} {out_bg} # placeholder file. Will be removed
        ln -rs {in_bw} {out_sorted_bg} # placeolder file. Will be removed
        ln -rs {in_bw} {out_bw}

        """).format(in_bw=input.bigwigs,
                   out_bg=output.bedgraph,
                   out_sorted_bg=output.sorted_bg,
                   out_bw=output.bigwig,
                   in_chroms=input.chromsizes)
    threads: 8
    envmodules:
        config['modules']['ucsc']
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
        "bbc/deeptools/deeptools-3.5.1"
    params:
        after="2000",
        before="2000",
        binsize=bamCoverage_binsize,
        samp_labels=lambda wildcards, input: " ".join(os.path.basename(x).replace(".bw", "") for x in input.bw),
        temp="analysis/deeptools_heatmap_genes/ht_tmp",
        yaxislabel='"Sum CPMs"',
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
        {params.sam_keep} \
        --blackListFileName {params.blacklist}

        """

def get_macs2_bams(wildcards):
    if (config['atacseq'] and wildcards.macs2_type == "macs2_nfr"):
        macs2_bams = { 'trt': "analysis/filt_bams_nfr/{sample}_filt_alns_nfr.bam".format(sample=wildcards.sample) }
    elif (wildcards.macs2_type == "macs2"):
        macs2_bams = { 'trt': "analysis/filt_bams/{sample}_filt_alns.bam".format(sample=wildcards.sample) }
    
    control = samples[samples['sample'] == wildcards.sample]['control'].values[0]
    
    if (not pd.isnull(control)):
        macs2_bams['control'] = expand("analysis/filt_bams/{sample}_filt_alns.bam", sample = control.split(','))
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
        --keep-dup 1 \
        --tempdir {params.temp_dir} 

        echo "END broad peak calling" >> {log.stdout}
        echo "END broad peak calling" >> {log.stderr}

        macs2 \
        callpeak \
        -t {input.trt} \
        {params.control_param} \
        {params.macs2_format} \
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
        config['modules']['bedtools']
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
    atac_peak_types = ['macs2_narrow','macs2_broad', 'macs2_nfr_narrow','macs2_nfr_broad']
    if config['run_hmmratac']:
        atac_peak_types = atac_peak_types + ['hmmratac_nfr']# if config['run_hmmratac']
    merged_beds = expand("analysis/merge_all_peaks/all_merged_{peak_type}.bed", peak_type=atac_peak_types if config['atacseq'] else ['macs2_narrow','macs2_broad'])
    return (merged_beds)

rule deeptools_plotenrichment:
    input:
        bam="analysis/filt_bams/{sample}_filt_alns.bam",
        bed=expand("analysis/{macs2_type}/{{sample}}_macs2_{type}_peaks.{type}Peak", type=['narrow','broad'], macs2_type=['macs2','macs2_nfr'] if config['atacseq'] else['macs2']),
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
        #sam_keep="64",
        #sam_exclude="1024",    
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


rule multiqc:
    input:
        expand("analysis/fastqc/{sample.sample}_R1_fastqc.html", sample=samples.itertuples()),
        expand("analysis/fastqc/{sample.sample}_R2_fastqc.html", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/trim_galore/{sample.sample}_R{read}_val_{read}_fastqc.html", sample=samples[samples['se_or_pe']=="PE"].itertuples(), read=["1","2"]),
        expand("analysis/trim_galore/{sample.sample}_R1_trimmed_fastqc.html", sample=samples[samples['se_or_pe']=="SE"].itertuples()),
        expand("analysis/preseq_complexity/{sample.sample}.c_curve.txt", sample=samples.itertuples()),
        expand("analysis/preseq_complexity/{sample.sample}.lc_extrap.txt", sample=samples.itertuples()),
        expand("analysis/bwamem/flagstat/{sample.sample}.flagstat", sample=samples.itertuples()),
        expand("analysis/filt_bams/{sample.sample}_filt_alns.bam.idxstat", sample=samples.itertuples()),
        expand("analysis/filt_bams/CollectInsertSizeMetrics/{sample.sample}_filt_alns.insert_size_metrics.txt", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/bwamem/CollectAlignmentSummaryMetrics/{sample.sample}.aln_metrics.txt", sample=samples.itertuples()),
        expand("analysis/bwamem/{sample.sample}.samblaster.e", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R1_screen.html", sample=samples.itertuples()),
        expand("analysis/fastq_screen/{sample.sample}_R2_screen.html", sample=samples[samples['se_or_pe']=="PE"].itertuples()),
        expand("analysis/deeptools_plotenrichment/{sample.sample}.rawcts", sample=samples_no_controls.itertuples()),
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
        "analysis/bwamem/*.samblaster.e",
        "analysis/fastq_screen/",
        "analysis/bwamem/CollectAlignmentSummaryMetrics/",
        "analysis/filt_bams/",
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
        --outdir {params.workdir} \
        --filename {params.outfile} \
        {params.dirs} {params.PE_dirs}

        """

