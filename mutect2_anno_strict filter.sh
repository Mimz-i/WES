#!/bin/bash

# mutect2_anno.sh - Modified to filter AD >= 5 and AF >= 0.05

# Maintainer: E Chambers

############################################################
# Help                                                     #
############################################################
Help()
{
   echo
   echo "*** mutect2_anno.sh ***"
   echo "Script to filter Mutect2 VCF output based on PASS, AD >= 5, AF >= 0.05, and annotate with ANNOVAR."
   echo "Syntax: scriptTemplate [-a|s|h]"
   echo "options:"
   echo "a     Set annovar directory - default = annovar/"
   echo "s     Set sampleSheet.csv"
   echo "o     Set output directory - default = mutect2_anno/"
   echo "i     Set input directory - default results/"
   echo "h     Print this Help."
}

# Default paths
annovar_dir=$PWD/annovar
output_dir=$PWD/mutect2_anno2
input_dir=$PWD/results

# Get the options
while getopts ":hs:a:o:i:" option; do
   case $option in
      h) Help; exit;;
      a) annovar_dir=${OPTARG};;
      s) sample_sheet=${OPTARG};;
      o) output_dir=${OPTARG};;
      i) input_dir=${OPTARG};;
      ?) echo "Error: Invalid option"; Help; exit;;
   esac
done

echo "Running mutect2_anno with additional filtering for AD >= 5 and AF >= 0.05"

# Pre-flight checks
if [ -z ${sample_sheet} ]; then
   echo "Please specify input CSV sample sheet using the -s option."
   exit 1
elif [ ! -f "${sample_sheet}" ]; then
   echo "Error: Cannot find sample sheet at ${sample_sheet}."
   exit 1
fi

if ! command -v bcftools &> /dev/null; then
   echo 'Error: bcftools not found.' >&2
   exit 1
fi

if [ ! -d "${annovar_dir}" ]; then
   echo "Error: Cannot find annovar at ${annovar_dir}. Use -a to set annovar dir."
   exit 1
fi

if [ ! -d "${input_dir}" ]; then
   echo "Error: Cannot find results at ${input_dir}. Use -i to set nextflow results dir."
   exit 1
fi

mkdir -p $output_dir
humandb=${annovar_dir}/humandb

# Loop through samples
for sample in $(cut -d',' -f4 ${sample_sheet} | tail -n +2)
do
    echo "Processing ${sample}..."
    mkdir -p $output_dir/$sample
    
    # Define file paths
    input_vcf=${input_dir}/variant_calling/mutect2/${sample}/${sample}.mutect2.filtered.vcf.gz
    pass_vcf=${output_dir}/${sample}/${sample}.mutect2.filtered.variants.PASS.vcf.gz
    filtered_vcf=${output_dir}/${sample}/${sample}.mutect2.filtered.AD_AF.vcf.gz
    
    # Step 1: Filter by PASS
    bcftools view -Oz -f PASS ${input_vcf} > ${pass_vcf}
    bcftools index ${pass_vcf}
    
    # Step 2: Extract AD and AF from FORMAT column and apply filtering
    bcftools view -i 'FMT/AD[*:1]>=5 & FMT/AF[*]>=0.05' ${pass_vcf} -Oz -o ${filtered_vcf}
    bcftools index ${filtered_vcf}
    
    echo "Filtering complete. Annotating..."
    
    # Step 3: Convert to ANNOVAR format and annotate
    ${annovar_dir}/convert2annovar.pl -format vcf4 $filtered_vcf > ${filtered_vcf}.avinput
    ${annovar_dir}/annotate_variation.pl -filter -infoasscore -otherinfo -dbtype gnomad_exome -buildver hg38 ${filtered_vcf}.avinput $humandb
    ${annovar_dir}/table_annovar.pl ${filtered_vcf}.avinput $humandb -buildver hg38 -remove -protocol refGene,cytoBand,dbnsfp30a,clinvar_20221231,nci60,cosmic70,gnomad_exome,exac03 -operation gx,r,f,f,f,f,f,f -nastring NA
done
