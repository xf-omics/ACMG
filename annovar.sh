#! /bin/bash


# source $HOME/.bash_profile
if [ "$#" -lt 1 ]; then
	echo 'Please provide input file without the file extension .vcf'
    exit 1
else
	INPUTFILE=$1
    OUTPUT_DIR=$2
fi

ANNOVAR_DIR=/share/data/resources/hg19/ANNOVAR_humandb
annovar=/home/xf2193/software/annovar

cd $OUTPUT_DIR

$annovar/table_annovar.pl $INPUTFILE.vcf $ANNOVAR_DIR -buildver hg19 -out $INPUTFILE -remove -protocol refGene,gnomad_exome,gnomad_genome,avsnp147,intervar_20180118,clinvar_20190305,revel,cadd14coding,spliceAI,dbnsfp31a_interpro -operation g,f,f,f,f,f,f,f,f,f -nastring . -vcfinput

# Rscript ACMG.R $INPUTFILE

