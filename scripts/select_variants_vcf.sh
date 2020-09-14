#!/bin/bash

## This script takes as argument the name of a vcf file with extension vcf.gz
## The result is an index file and another vcf where all multi-allelic sites are split

input_vcf=$1
output_dir=$2
ref=$3

output_base=$( basename $input_vcf )
output_vcf=${output_base%vcf.gz}
output_vcf=${output_vcf%vcf}split.vcf

## generate vcf index with tabix

/cluster/tufts/bio/tools/bcbio/1.1.5/bin/tabix -p vcf $input_vcf

/cluster/tufts/bio/tools/GATK/gatk-4.1.2.0/gatk LeftAlignAndTrimVariants \
 -R $ref \
 -V $input_vcf \
 -O ${output_dir}/${output_vcf} \
 --max-indel-length 535 \
 --split-multi-allelics

