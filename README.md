# PSAP
Pipeline code for calculating population sampling probabilities



PEDIGREE FILE FORMAT (SPACE SEPARATED, NO HEADER):
FAMILY ID
INDIVIDUAL ID
PATERNAL ID (0 IF NO FATHER)
MATERNAL ID (0 IF NO MOTHER)
GENDER (1 FOR MALE, 2 FOR FEMALE)
CASE-CONTROL STATUS (1 FOR UNAFFECTED, 2 FOR AFFECTED)
NOTE: The individual ID must match the ID used in the VCF header for that individual
