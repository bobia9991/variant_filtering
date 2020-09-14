#!/bin/bash

## Run example
## ./run_pipeline.sh -v ../test/test_raw.vcf.gz -g ../test/genelist.txt -b GRCh38 -o ../outdir

while getopts ":v:g:o:b:" opt; do
  case $opt in
    v) vcf="$OPTARG"
    ;;
    g) genelist="$OPTARG"
    ;;
    b) build="$OPTARG"
    ;;
    o) outdir="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ -z $vcf ]
then
  exit "Error: You must provide a vcf with option -v. Example usage: sh run_pipeline.sh --vcf input.vcf --genelist genelist.txt --build genome-build --outdir outdir"
else
  printf "Argument vcf is %s\n" "$vcf"
fi

if [ -z $build ]
then
  exit "Error: You must provide a genome build with option -b. Example usage: sh run_pipeline.sh --vcf input.vcf --genelist genelist.txt --build genome-build --outdir outdir"
else
  printf "Argument build is %s\n" "$build"
fi

if [ -z $outdir ]
then
  exit "Error: You must provide an output dir with option -o. Example usage: sh run_pipeline.sh --vcf input.vcf --genelist genelist.txt --build genome-build --outdir outdir"
else
  printf "Argument outdir is %s\n" "$outdir"
fi

printf "Argument genelist is %s\n" "$genelist"

## get script dir
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## starting analysis
module load anaconda/3
module load java/1.8.0_60

date=`date +%m-%d-%Y`

if [ $build = 'GRCh38' ]
then
  ref=/cluster/tufts/bio/data/genomes/HomoSapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa
  cache_dir=/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache/
  cadd_files=/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache/CADD/grch38/whole_genome_SNVs.tsv.gz,/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache/CADD/grch38/InDels.tsv.gz
  clinvar=/cluster/tufts/bio/data/clinvar/15oct19/clinvar.vcf.gz
elif [ $build = 'b37' ]
then
  ref=/cluster/tufts/bio/data/genomes/HomoSapiens/Broad/b37/Sequence/WholeGenomeFasta/hs37d5.fa.gz
  cache_dir=/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache-grch37/
  cadd_files=/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache/CADD/grch37/whole_genome_SNVs.tsv.gz,/cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/cache/CADD/grch37/InDels.tsv.gz
  clinvar=/cluster/tufts/bio/data/clinvar/grch37/3Aug20/clinvar.vcf.gz
else
  exit "Supported builds are b37 and GRCh38. You entered build $build which is not supported"
fi

base_out=$(basename $vcf)
base_out_remove=${base_out%.vcf.gz}
base_out_remove=${base_out_remove%.vcf}
out_prefix=${outdir}/${base_out_remove}

mkdir -p ${outdir}

echo "----- Starting split multi-allelic variants $vcf ----"

sh $DIR/select_variants_vcf.sh $vcf $outdir $ref

echo "--- Starting Hmtnote annotation ----"
source activate /cluster/tufts/bio/tools/conda_envs/ensembl-vep-versions/98/

hmtnote annotate ${out_prefix}.split.vcf ${out_prefix}.split.hmtnote.vcf --variab --offline

echo "--- Starting VEP annotation ----"
sh $DIR/run_vep_pickgene_gencode.sh ${out_prefix}.split.hmtnote.vcf $ref $cache_dir $cadd_files $clinvar

echo "--- Starting conversion to TSV ----"
sh $DIR/variants_to_table.sh ${out_prefix}.split.hmtnote.pickgene-gencode.vcf $ref

echo "---- Starting parsing and filtering with Python  -----"

python $DIR/formatcsq.py -tsv ${out_prefix}.split.hmtnote.pickgene-gencode.tsv -vcf ${out_prefix}.split.hmtnote.pickgene-gencode.vcf

if [ -z $genelist ]
then
  python $DIR/filter.py -tsv ${out_prefix}.split.hmtnote.pickgene-gencode.formatcsq.tsv -build $build
else
  python $DIR/filter.py -tsv ${out_prefix}.split.hmtnote.pickgene-gencode.formatcsq.tsv -genelist $genelist -build $build
fi

# move intermediate files to tmp location
mkdir -p ${outdir}/tmp
mv ${out_prefix}.split.hmtnote.pickgene-gencode.tsv ${outdir}/tmp
mv ${out_prefix}.split.hmtnote.pickgene-gencode.vcf ${outdir}/tmp 
mv ${out_prefix}.split.hmtnote.pickgene-gencode.vcf_summary.html ${outdir}/tmp
mv ${out_prefix}.split.hmtnote.pickgene-gencode.formatcsq.tsv ${outdir}/tmp
mv ${out_prefix}.split.vcf ${outdir}/tmp
mv ${out_prefix}.split.vcf.idx ${outdir}/tmp
mv ${out_prefix}.split.hmtnote.vcf ${outdir}/tmp
