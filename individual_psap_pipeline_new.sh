#!/bin/bash
# Requires vcf file as input
# $1 = vcf file, $2 = output file, $3 = ped file, $4 = genome build version

PSAP_PATH=/scratch/dclab/psap4/psap/
ANNOVAR_PATH=/scratch/dclab/annovar/
FORCE=0
#FORCE="-F"
# forces to output p-values for all possible coding variants instead of the top ones for each gene

echo $#

echo $PWD
module load R
if [ $# -gt 0 ] && [ $1 == "-h" ]
then
	echo "arg1 =  VCF file"
	echo "arg2 = output file name" 
	echo "arg3 = family pedigree file"
	echo "arg4 = human genome build version (eg. hg19, hg38)"
	echo "arg5 = indels annotated with CADD"
	echo "arg6 = CNV calls"
	echo "Example: popScore_analysis_pipeline.sh INDV001.vcf INDV.OUTFILE INDV001.ped HG_BUILD"
	exit
fi

        OUTFILE=$2 # Name of output file (no directory should be included here)
        PED_FILE=$3 # Name of pedigree file (directory should be included here)
	HG_BUILD=$4
	INDEL_FILE=$5 # Name of indel file (directory should be included here)
	CNV_FILE=$6 # Name of CNV file (directory should be included here)


if [ $# -ge 4 ]
then
# Establish Genocde version
	if [[ $HG_BUILD == "hg38" ]]
        then
                GENCODE_VER=20
        else
                GENCODE_VER=19
        fi
# Check that all required ANNOVAR annotation files are in the humandb directory
        MISSING=0
        for FILE in "${HG_BUILD}_ALL.sites.2015_08.txt" "${HG_BUILD}_cadd13.txt" "${HG_BUILD}_esp6500siv2_all.txt" "${HG_BUILD}_wgEncodeGencodeBasicV${GENCODE_VER}Mrna.fa" "${HG_BUILD}_wgEncodeGencodeBasicV${GENCODE_VER}.txt" "${HG_BUILD}_exac03.txt" "${HG_BUILD}_gnomad_exome.txt" "${HG_BUILD}_clinvar_20180603.txt"
        do
                if [ ! -f ${ANNOVAR_PATH}humandb/$FILE ]
                then
                        MISSING=$(( $MISSING+1 ))
                fi
        done
# If any of the required annotation files are missing, exit and print error  message
        if [ $MISSING -gt 0 ]
        then
                echo "ERROR: Missing required ANNOVAR annotations.  Please run get_annovar_annos.sh prior to running this script."
                exit
        fi

# Extract and move to VCF file directory
        FILE_LOC=${1%/*.vcf} # Extract location of VCF file
        cd $FILE_LOC # Use location of  VCF file as working directory, this is where all output will be written
        echo $PWD
        VCF=${1##/*/} # Extract VCF file name
	
# Convert vcf file to annovar file
        echo "PROGRESS: Converting VCF file to annovar input"
        cut -f1-8 $VCF > tmp #make small file just of sites, annotate this instead of full VCF
       perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old tmp -outfile ${OUTFILE}.avinput -includeinfo
       rm tmp
# Write column names from VCF file to header file (will be used later)
        grep '#' $VCF | tail -n 1 > ${OUTFILE}.avinput.header # Extract all lines of the VCF header.  The last line of the VCF header contains coumn names - write columna names to .avinput.header file

# If there is no annotated directory create annotated directory
        if [ $(ls -d $PWD/*/ | grep -c -w "annotated") == 0 ]
        then
                echo "Creating directory annotated/ to store annotation and analysis output"
                mkdir annotated
        fi

# Annotate with ANNOVAR
	echo "PROGRESS: Annotating data with ANNOVAR"
#	perl ${ANNOVAR_PATH}table_annovar.pl ${OUTFILE}.avinput -remove -outfile annotated/${OUTFILE}.avinput ${ANNOVAR_PATH}humandb/ -buildver $HG_BUILD -protocol wgEncodeGencodeBasicV${GENCODE_VER},gnomad_exome,exac03,esp6500siv2_all,1000g2015aug_all,cadd13,clinvar_20180603 -operation g,f,f,f,f,f,f -nastring NA -otherinfo -argument -separate,,,,,,


# Split fusion genes and duplicate the variants (work for >=2 gene fused together)
        echo "PROGRESS: Split fusion genes and duplicate the variants"
        python3 ${PSAP_PATH}pythonscript/fusion_gene_split.py annotated/${OUTFILE}.avinput.${HG_BUILD}_multianno.txt


# EXTRACT INDIVIDUAL IDS
	echo "PROGRESS: Extracting individual IDs"
	IDS=($(awk '{print $2}' $PED_FILE))
	IDX=1

# RUN PSAP_individual.R for each individual
	echo "PROGRESS: Starting PSAP annotation" 
	
NSAMP=${#IDS[@]}
echo "#!/bin/bash " >psap.sbatch
echo "#SBATCH --mem=12g " >>psap.sbatch
echo "#SBATCH --array=1-$NSAMP" >> psap.sbatch
echo "SAMPLE=\$( sed -n \${SLURM_ARRAY_TASK_ID}p $PED_FILE | cut -f2 -d\" \")" >> psap.sbatch
echo "\${SAMPLE}" >> psap.sbatch
echo "ml vcftools">> psap.sbatch
echo "ml R" >> psap.sbatch
echo "vcftools --indv \${SAMPLE} --vcf $VCF --recode --out \${SAMPLE} ">> psap.sbatch
echo "perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old \${SAMPLE}.recode.vcf -outfile \${SAMPLE}.avinput -includeinfo" >> psap.sbatch
echo "Rscript ${PSAP_PATH}RScripts/PSAP_individual.R ${OUTFILE} \${SAMPLE} $PSAP_PATH $PED_FILE $HG_BUILD $FORCE $INDEL_FILE" >> psap.sbatch
echo "rm \${SAMPLE}.avinput \${SAMPLE}.recode.vcf" >> psap.sbatch

       
sbatch psap.sbatch


# Generate report file - will look for variants present in two or more affected with PSAP < 1e-3 and not observed in unaffected
	#echo "PROGRESS: Generating report file for all individuals"
	#Rscript ${PSAP_PATH}RScripts/unrelated_candidate_analysis.R ${OUTFILE}.avinput $PED_FILE $PSAP_PATH $HG_BUILD

else
	echo "ERROR: Incorrect number of arguments." $# "arguments provided"
	echo "Please provide a VCF file, output file name, pedigree file and genome build version.  Please use the -h argument for help."
fi
