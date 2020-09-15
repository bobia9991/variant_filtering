# variant_filtering

##This directory contains scripts to process and annotate variant call files

There are two sub directories:
1) scripts - contains all code

2) test - contains test data

##To run the test data on the Tufts HPC:

1. `git clone https://github.com/rbatorsky/variant_filtering.git`
2. `cd variant_filtering/scripts`
3. The script `sbatch_variant_filtering_test.sh` runs the script `run_pipeline.sh` with the test files. To submit to slurm in batch mode: `sbatch sbatch_variant_filtering_test.sh` To run in interactive mode, first get an interactive session, e.g.: `srun --pty --mem=100Gb --cpus=4 bash`.
Then run `sh sbatch_variant_filtering_test.sh`

VCF file:
In order to run this on a new vcf, edit the file `sbatch_variant_filtering_test.sh` put the full path to the gzipped vcf on the line:
`-v ../test/test_raw.vcf.gz \`

Genelist:
Genelist argument is optional.
To change the genelist, specify the full path to the genelist in this line:
`-g ../test/genelist.txt \`

Genome build:
Build GRCh38 and b37 (broad's version of GRCh37) are supported.
Geneome build can be specified this way:
`-b GRCh38 \`

The output directory can be specified in the last line:
`-o out`

##Output files

When run on the test data with the example run script, the output directory contains the following files:
-test_raw.split.hmtnote.pickgene-gencode.formatcsq.unfiltered.tsv
This contains all variants that were in the initial vcf, one line per variant per transcript.
That is, variants will be listed multiple times, once for each transcript for which they have a VEP consequence.

-test_raw.split.hmtnote.pickgene-gencode.formatcsq.genelist.tsv
If a genelist is provided, this variant list is restricted to the genelist

-test_raw.split.hmtnote.pickgene-gencode.formatcsq.genelist.removecols.stringent-filter.tsv
Stringent filtered variant list that contains variants that meet the following criteria

        Variant has only Pathogenic or Likely_pathogenic clinvar reports with no conflicts (which means no Benign or Likely_benign reports)
        &
        There is assertion criteria in the Clinvar review status for the variant
        &
        MAX_AF is < 0.01 or missing
        &
        VEP IMPACT is High or Moderate

-test_raw.split.hmtnote.pickgene-gencode.formatcsq.genelist.removecols.stringent-filter.biobank.tsv
Same as above but with biobank format

-test_raw.split.hmtnote.pickgene-gencode.formatcsq.genelist.removecols.relaxed-filter.tsv
Relaxed filtered variant list that contains variants that meet the following criteria

        Variant has only Pathogenic or Likely_pathogenic clinvar reports with no conflicts (which means no Benign or Likely_benign reports)
        OR
        variant has NO Benign or Likely_benign reports  & MAX_AF is < 0.01 or missing  & VEP IMPACT is High or Moderate

-test_raw.split.hmtnote.pickgene-gencode.formatcsq.genelist.removecols.relaxed-filter.biobank.tsv
Same as above but with biobank format

-tmp
Folder with intermediate outputs, can be ignored but is useful for debugging.
