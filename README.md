# Peaks Workflow

This workflow is designed to align, perform basic QC and call peaks for peak-based methods such as ChIP-seq, CUT&RUN and ATAC-seq.

# Usage
1. `git clone <repo> <new_directory_name>`

2. Put the 'fastq.gz' files in the `raw_data/` subdirectory. You may use symlinks.

3. Look over the config file to ensure that the correct indexes and other species-specific settings are used.

4. Set up the sample sheet (samples.tsv). You may find it helpful to run `./make_samples_template.sh` from inside the `bin/` subdirectory to get a template file based on the files in `raw_data/`. The columns are:

   i.`sample` -- Name of the sample. Fastqs will be renamed to this. You may use the same 'sample' name in multiple rows in this file to represent fastqs that should be `cat` together.
   
   ii. `control` -- 'sample' name for the control to be used, for example, as control for MACS2. If not applicable, use 'NA'. Leaving it blank should also be fine.
   
   iii. `fq1` -- R1 file
    
    iv. `fq2` -- R2 file; fill in with NA if SE data.

   v. `sample_group` -- Grouping variable for replicates.

5. From the root directory of the project, run `qsub bin/run_snakemake.sh`.
