# PSAP
Pipeline code for calculating population sampling probabilities

## REQUIRED SOFTWARE
This pipeline uses the R statistical software and ANNOVAR.  Please ensure R (http://r-project.org) and ANNOVAR (http://annovar.openbioinformatics.org) are installed.  Paths to all other accessory softwares/scripts are hard coded to the directories within the PSAP directory.


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


## WHEN USING THIS SCRIPT PLEASE CITE
Wang K, Li M, Hakonarson H. ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data. Nucleic Acids Research, 38:e164, 2010

Wilfert AB, Chao K, Kaushal M, Jain S, Zöllner S, Adams DR and Conrad DF.  Genome-wide significance testing of variation from single case exomes. Nature Genetics. doi:10.1038/ng.3697. 2016
