# PSAP
Pipeline code for calculating population sampling probabilities

## CONTENTS
THE PSAP PACKAGE CONTAINS BASH SCRIPTS, AND ALL DEPENDANT PYTHON AND R SCRIPTS AND LOOKUP TABLES. EACH SCRIPT IS DESCRIBED BELOW.

1) ```individual_analysis_pipeline.sh```: Calls ANNOVAR to annotate data, calls a python script to appropriately annotate variants that overlap transcripts for multiple genes, calls an Rscript that performs some basic cleaning steps (mendelian inheritance filter - allows de novos, PSAP calibration filter, missing data filter, allele frequency discrepancy filter) and annotates all individuals with PSAP, and calls an R script that will report out candiate variants (inheritance pattern consistent with disease model)

## REQUIRED SOFTWARE
This pipeline uses the R statistical software and ANNOVAR.  Please ensure R (http://r-project.org) and ANNOVAR (http://annovar.openbioinformatics.org) are installed.  Paths to all other accessory softwares/scripts are hard coded to the directories within the PSAP directory. Note that ANNOVAR output formats have changed over time, and thus it is recommended to use the verion of ANNOVAR used in current PSAP development (Version: Thu, 24 Oct 2019). 


### PEDIGREE FILE FORMAT (SPACE SEPARATED, NO HEADER):

```
FAMILY ID
INDIVIDUAL ID
PATERNAL ID (0 IF NO FATHER)
MATERNAL ID (0 IF NO MOTHER)
GENDER (1 FOR MALE, 2 FOR FEMALE)
CASE-CONTROL STATUS (1 FOR UNAFFECTED, 2 FOR AFFECTED)
ETHNICITY (1 for generic model ["ALL"]; 2 for African ancestry; 3 for non-Finnish European model)
```
A new feature in this implementation of PSAP is ethnicity-specific models. If the ethnicity of your subject is unknown, then 
it is best to use coding "1"; if your sample has substantial (> 50%) ancestry from sub-Saharan Africa, we recommend coding "2"; coding "3" will be appropriate for most Caucasian samples (samples from population with strong isolation/drift/inbreeding may fare better with "1". In short, when in doubt, use model "1". 

##### NOTE: The individual ID must match the ID used in the VCF header for that individual


##### CONTENTS: 

PSAP is intended to be used as a pipeline for generating single-sample reports from a batch VCF file, taking advantage of a 
cluster environment for parallelization. The pipeline code on github consists of:
1. an example shell script individual_psap_pipeline_new.sh, which could be modified as needed for local needs, or it could 
simply serve as inspiration. This script produces a batch file for each individual exome to analyzed; each batch file is 
submitted using the slurm management system. If you don't have a cluster environment, or use a different one, this section
will need to be modified. 

2. R Scripts that, together, implement PSAP for individual samples (PSAP_individual.R and apply_PSAP.R). In theory these
scripts, and the PSAP lookup tables are all that are needed to run PSAP, assuming that the input data are complete and
formatted properly).

3. helper scripts and PSAP lookup tables. 


##### USAGE: 

Run "individual_psap_pipeline.sh -h" for an details of how to use the PSAP pipeline example script for a single sample VCF

"individual_psap_pipeline_batch.sh" is an example script is tailored to the slurm cluster environment; changes may be necessary. 

The full PSAP pipeline provides a method by which to intergrate analysis of SNVs, indels, and CNVs in the same statistical 
framework. In the current implmementation, we require SNVs and indels to be annotated with CADD1.3 phred scores. This means 
that indels must be provided in a separate file from SNVs (in most cases this indel file will be obtained directly from the
CADD web server: https://cadd.gs.washington.edu/score). The PSAP pipeline will work directly with the default, gzipped, file
format produced by the web server; no reformatting is necessary. 

Currently all of our analyses are based on CADD1.3, so it is important to use that version when annotating indels manual 
with the CADD webserver. The file format for CADD1.4 is different (as of fall of 2019) and will cause the PSAP pipeline to
break.

CNVs: documentation on how to incorporate CNV calls will be added here.


## WHEN USING THIS SCRIPT PLEASE CITE
Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data. Nucleic Acids Research, 38:e164, 2010

Wilfert AB, Chao K, Kaushal M, Jain S, Zöllner S, Adams DR and Conrad DF.  Genome-wide significance testing of variation from single case exomes. Nature Genetics. doi:10.1038/ng.3697. 2016
