#!/bin/bash

#module load anaconda/3
#source activate /cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/98/

if [ $# -lt 5 ]
  then
    echo "Error: no argument supplied. Usage: sh run_vep_pickgene_genecode.sh input.vcf ref cache_dir cadd_files clinvar_file"
    exit 1
fi

input_vcf=$1
REF=$2
cache_dir=$3
cadd_files=$4
clinvar=$5

annotated_vcf=${input_vcf%vcf}pickgene-gencode.vcf
species=homo_sapiens


#cache_dir=/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache/
#cadd_files=/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache/CADD/grch38/whole_genome_SNVs.tsv.gz,/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache/CADD/grch38/InDels.tsv.gz
#clinvar=/cluster/tufts/bio/data/clinvar/15oct19/clinvar.vcf.gz

## Annotate vcf, output vcf, filter common variants
vep -i $input_vcf \
--force_overwrite \
--fork 4 \
--vcf \
--per_gene \
--no_intergenic \
--canonical \
--cache \
--fasta $REF \
--max_af \
--af_1kg \
--af_esp \
--af_gnomad \
--sift b \
--polyphen b \
--plugin CADD,$cadd_files \
--offline \
--dir_cache $cache_dir \
--species $species \
--custom $clinvar,ClinVar,vcf,exact,0,CLNHGVS,GENEINFO,CLNSIG,CLNREVSTAT,CLNDN, \
-o $annotated_vcf

echo "Done."
