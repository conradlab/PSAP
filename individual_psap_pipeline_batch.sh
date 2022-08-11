#!/bin/bash
# Requires vcf file as input
# $1 = vcf file, $2 = output file, $3 = ped file, $4 = genome build version
bcftools=/opt/installed/bcftools-1.6/bin/bcftools
PSAP_PATH=/home/groups/ConradLab/daniel/gemini_phase_2_analysis/PSAP/
ANNOVAR_PATH=/home/groups/ConradLab/daniel/annovar/
FORCE=0
#FORCE="-F"
# forces to output p-values for all possible coding variants instead of the top ones for each gene

echo $#

echo $PWD
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
                GENCODE_VER=35
        else
                GENCODE_VER=19
        fi
# Check that all required ANNOVAR annotation files are in the humandb directory
        MISSING=0
        for FILE in "${HG_BUILD}_wgEncodeGencodeBasicV${GENCODE_VER}Mrna.fa" "${HG_BUILD}_wgEncodeGencodeBasicV${GENCODE_VER}.txt" "${HG_BUILD}_gnomad211_exome.txt" "${HG_BUILD}_clinvar_20210501.txt"
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


# PROJET MANAGEMENT
	if [ ! -d avinput ]; then 
        mkdir avinput
	fi

	if [ ! -d annotated ]; then 
	mkdir annotated
	fi

# Extract and move to VCF file directory
        FILE_LOC=${1%/*.vcf.gz} # Extract location of VCF file
#        cd $FILE_LOC # Use location of  VCF file as working directory, this is where all output will be written
#        echo $PWD
        VCF=${1##/*/} # Extract VCF file name
	VCF=${1}

	
# Convert vcf file to annovar file
        echo "PROGRESS: Converting VCF file to annovar input"
        $bcftools view $VCF |  cut -f1-8  > tmp #make small file just of sites, annotate this instead of full VCF
        perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old tmp -outfile ${OUTFILE}.avinput -includeinfo
        rm tmp
# Write column names from VCF file to header file (will be used later)
        $bcftools view $VCF | grep '#' | tail -n 1 > ${OUTFILE}.avinput.header # Extract all lines of the VCF header.  The last line of the VCF header contains coumn names - write columna names to .avinput.header file

# If there is no annotated directory create annotated directory
        if [ $(ls -d $PWD/*/ | grep -c -w "annotated") == 0 ]
        then
                echo "Creating directory annotated/ to store annotation and analysis output"
                mkdir annotated
        fi

# Annotate with ANNOVAR
	echo "PROGRESS: Annotating data with ANNOVAR"
	perl ${ANNOVAR_PATH}table_annovar.pl ${OUTFILE}.avinput -remove -outfile annotated/${OUTFILE}.avinput ${ANNOVAR_PATH}humandb/ -buildver $HG_BUILD -protocol wgEncodeGencodeBasicV${GENCODE_VER},gnomad211_exome,clinvar_20210501 -operation g,f,f -nastring NA -otherinfo -argument -separate,,

# Annotate with CADD
        echo "PROGRESS: Annotate CADD scores"
	# insert your conda env here with necessary packages
        conda activate daniel_CADD
        python3 /home/groups/ConradLab/daniel/gemini_phase_2_analysis/annotate_CADD/extract_scores.py -p /home/groups/ConradLab/daniel/gemini_phase_2_analysis/CADD1.6/CADD1.6_hg38_whole_genome_SNVs.tsv.gz -i annotated/${OUTFILE}.avinput.${HG_BUILD}_multianno.txt --found_out annotated/${OUTFILE}.avinput.${HG_BUILD}_multianno_CADD.txt

# Split fusion genes and duplicate the variants (work for >=2 gene fused together)
        echo "PROGRESS: Split fusion genes and duplicate the variants"
        python3 ${PSAP_PATH}pythonscript/fusion_gene_split.py annotated/${OUTFILE}.avinput.${HG_BUILD}_multianno_CADD.txt


# EXTRACT INDIVIDUAL IDS
	echo "PROGRESS: Extracting individual IDs"
	IDS=($(awk '{print $2}' $PED_FILE))
	IDX=1




# RUN PSAP_individual.R for each individual
	echo "PROGRESS: Starting PSAP annotation" 
NSAMP=${#IDS[@] + 1}
echo "#!/bin/bash " >psap1.sbatch
echo "#SBATCH --mem=100g " >>psap1.sbatch
echo "#SBATCH --array=2-1000%500" >> psap1.sbatch
echo "SAMPLE=\$( sed -n \${SLURM_ARRAY_TASK_ID}p $PED_FILE | cut -f2 )" >> psap1.sbatch
echo "\${SAMPLE}" >> psap1.sbatch
echo "if [ ! -f avinput/\${SAMPLE}.avinput ]" >> psap1.sbatch
echo "then" >> psap1.sbatch
echo "$bcftools view $VCF -s \${SAMPLE} -c 1 | $bcftools annotate -x ^INFO/AC >  \${SAMPLE}.tmp.vcf ">> psap1.sbatch
echo "perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old \${SAMPLE}.tmp.vcf -outfile avinput/\${SAMPLE}.avinput -includeinfo" >> psap1.sbatch
echo "fi" >> psap1.sbatch
echo "Rscript ${PSAP_PATH}RScripts/PSAP_individual.R ${OUTFILE} \${SAMPLE} $PSAP_PATH $PED_FILE $HG_BUILD $FORCE $INDEL_FILE" >> psap1.sbatch
echo "rm \${SAMPLE}.avinput \${SAMPLE}.tmp.vcf" >> psap1.sbatch

echo "#!/bin/bash " >psap2.sbatch
echo "#SBATCH --mem=100g " >>psap2.sbatch
echo "#SBATCH --array=1001-2000%500" >> psap2.sbatch
echo "SAMPLE=\$( sed -n \${SLURM_ARRAY_TASK_ID}p $PED_FILE | cut -f2 )" >> psap2.sbatch
echo "\${SAMPLE}" >> psap2.sbatch
echo "if [ ! -f avinput/\${SAMPLE}.avinput ]" >> psap2.sbatch
echo "then" >> psap2.sbatch
echo "$bcftools view $VCF -s \${SAMPLE} -c 1 | $bcftools annotate -x ^INFO/AC >  \${SAMPLE}.tmp.vcf ">> psap2.sbatch
echo "perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old \${SAMPLE}.tmp.vcf -outfile avinput/\${SAMPLE}.avinput -includeinfo" >> psap2.sbatch
echo "fi" >> psap2.sbatch
echo "Rscript ${PSAP_PATH}RScripts/PSAP_individual.R ${OUTFILE} \${SAMPLE} $PSAP_PATH $PED_FILE $HG_BUILD $FORCE $INDEL_FILE" >> psap2.sbatch
echo "rm \${SAMPLE}.avinput \${SAMPLE}.tmp.vcf" >> psap2.sbatch

echo "#!/bin/bash " >psap3.sbatch
echo "#SBATCH --mem=100g " >>psap3.sbatch
echo "#SBATCH --array=2001-$NSAMP%500" >> psap3.sbatch
echo "SAMPLE=\$( sed -n \${SLURM_ARRAY_TASK_ID}p $PED_FILE | cut -f2 )" >> psap3.sbatch
echo "\${SAMPLE}" >> psap3.sbatch
echo "if [ ! -f avinput/\${SAMPLE}.avinput ]" >> psap3.sbatch
echo "then" >> psap3.sbatch
echo "$bcftools view $VCF -s \${SAMPLE} -c 1 | $bcftools annotate -x ^INFO/AC >  \${SAMPLE}.tmp.vcf ">> psap3.sbatch
echo "perl ${ANNOVAR_PATH}convert2annovar.pl -format vcf4old \${SAMPLE}.tmp.vcf -outfile avinput/\${SAMPLE}.avinput -includeinfo" >> psap3.sbatch
echo "fi" >> psap3.sbatch
echo "Rscript ${PSAP_PATH}RScripts/PSAP_individual.R ${OUTFILE} \${SAMPLE} $PSAP_PATH $PED_FILE $HG_BUILD $FORCE $INDEL_FILE" >> psap3.sbatch
echo "rm \${SAMPLE}.avinput \${SAMPLE}.tmp.vcf" >> psap3.sbatch

       
#sbatch psap1.sbatch
echo "submitted first batch"
#sbatch psap2.sbatch
echo "done with second batch"
#sbatch psap3.sbatch
echo "done with third batch"

# Generate report file - will look for variants present in two or more affected with PSAP < 1e-3 and not observed in unaffected
	#echo "PROGRESS: Generating report file for all individuals"
	#Rscript ${PSAP_PATH}RScripts/unrelated_candidate_analysis.R ${OUTFILE}.avinput $PED_FILE $PSAP_PATH $HG_BUILD

else
	echo "ERROR: Incorrect number of arguments." $# "arguments provided"
	echo "Please provide a VCF file, output file name, pedigree file and genome build version.  Please use the -h argument for help."
fi

