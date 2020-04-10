#!/bin/Rscript

#Rscript PSAP_individual.R ${OUTFILE}.avinput $i $PSAP_PATH $PED_FILE $HG_BUILD $Force INDEL_FILE CNV_DEL_FILE CNV_DUP_FILE
print(getwd())
## NAME OF FILE FOR ANALYSIS, PROVIDED AS AN ARGUMENT WHEN CALLING THE RSCRIPT
arg <- commandArgs(trailingOnly=T) ## arg[1] = family - must be annovar annotated (avinput.hg19_multianno.txt) and have a separate header file for the vcf colums(header); arg[2] = individual ID; arg[3] = path to psap directory; arg[4] = ped file; arg[5] = genome build version; arg[6] whether report all the variants instead output the worst in each gene; arg[7] = optional input file of indel calls annotated by CADD online; arg[8] = optional input file of CNV-DEL calls; arg[9] = optional input file of CNV-DUP calls
print(arg[1])
print(arg[2])
print(arg[3])
print(arg[4])
print(arg[5])
print(arg[6])
print(arg[7])
print(arg[8])
print(arg[9])

CNVCALLS<-0 ## VARIABLE WILL BE SET TO 1 IF CNV CALLS ARE AVAILABLE FOR THIS sAMPLE


#ARGUMENTS FOR TESTING VIA EMACS
#arg<-c("GEMINI4","H_VZ-P16-P16","/Users/conradon/psap/PSAP/","/Users/conradon/gemini/analysis/dec/GEMINI_190304_all_genotyped_samples_clean.ped","hg19","FALSE","/Users/conradon/gemini/analysis/dec/GEMINI_190505_variants_FINAL_biallelicINDELs.hg19_cut_CADD_CHR.tsv","/Users/conradon/gemini/analysis/dec/Gemini4-sq60-DEL-NOA-no-exclusions-no-outliers-maf-1pc-h19-liftover-final.hg19_multianno_dup.txt","/Users/conradon/gemini/analysis/dec/Gemini4-NOA-sq60-DUP-LoF-no-exclusions-no-outliers-maf-1pc-hg19-liftover.hg19_multianno_dup.txt")



fam.id <- strsplit(arg[1],".avinput",fixed=T)
indv.id <- arg[2] # Individual ID - ASSUMES only one individual is being analyzed/annotated
dir <- arg[3]
ped <- read.table(arg[4],stringsAsFactors=F)
build_ver <-arg[5]
Force <- "FALSE"
indel.file <- arg[7]
cnv.del.file <- arg[8]
cnv.dup.file <- arg[9]

#find ethnicity information
if(ped$V7[which(ped$V2 == indv.id)] == 2) {
	eth <- "_AFR"
} else if (ped$V7[which(ped$V2 == indv.id)] == 3) {
	eth <- "_NFE"
} else {
	eth <- ""
}
## gene version and correspondingn columns
if(build_ver == "hg38"){
    gencode_ver = "wgEncodeGencodeBasicV20"
}else{
    gencode_ver = "wgEncodeGencodeBasicV19"
}
gene_col <- paste("Gene.",gencode_ver,sep="")
func_col <- paste("Func.",gencode_ver,sep="")
loc_col <- paste("ExonicFunc.",gencode_ver,sep="")
aa_col <- paste("AAChange.",gencode_ver,sep="")
score <- "CADD13_PHRED"
scale <- seq(0,1,7.14285e-4)
# read in lookup gene list and XY gene list
lookup.genes <- scan(paste(dir,"lookups/lookup_genes_041219_rm4.txt",sep=""),what = "character") #Need to make this (this is the whole gene list with XY genes)
XY_genelist <- scan(paste(dir,"lookups/Gender_genelist_041219_rm4.txt",sep=""),what = "character") #Need to put this





cnv<-NULL

if (is.na(cnv.del.file)==FALSE){
## ADDITIONAL EXPERIMENTAL INTEGRATION OF CNVS
#SHOULD BE PROCESSED WITH PYTHON SPLIT GENE SCRIPT BEFOREHAND
cnv<-read.table(cnv.del.file,header=F,sep="\t",quote="",comment.char="",skip=1,stringsAsFactors=F)
header<-read.table(cnv.del.file,header=F,sep="\t",nrow=1,stringsAsFactors=F)
header[length(header)]<-"ID"
header<-c(header,"TYPE","REGION","INFO","FILTER","QUAL","Alt_Read","Total_Read")

cnv.gt<-rep("het",dim(cnv)[1])
#QUICK AND DIRTY CLASSIFICATION OF HOMS
cnv.gt[which(cnv[,dim(cnv)[2]]==0)]<-"hom"

if (dim(cnv)[2]>length(header)){cnv<-cnv[,c(1:length(header))];}
colnames(cnv)<-as.character(header)
cnv$gt<-cnv.gt
cnv<-cnv[,-which(colnames(cnv)=="REGION")]

colnames(cnv)[dim(cnv)[2]]<-indv.id
cnv$ExonicFunc.wgEncodeGencodeBasicV19<-"frameshift deletion"
}

if (is.na(cnv.dup.file)==FALSE){
### NOW DUPlICATIONS
dup<-read.table(cnv.dup.file,header=F,sep="\t",quote="",comment.char="",skip=1,stringsAsFactors=F)
dup.header<-read.table(cnv.dup.file,header=F,sep="\t",nrow=1,stringsAsFactors=F)
dup.header[length(dup.header)]<-"ID"
dup.header<-c(dup.header,"TYPE","INFO","FILTER","QUAL","Alt_Read","Total_Read")

colnames(dup)<-dup.header
#stopifnot(colnames(dup)==colnames(cnv))

dup$ExonicFunc.wgEncodeGencodeBasicV19<-"frameshift insertion"
dup$gt<-"het"
colnames(dup)[dim(dup)[2]]<-indv.id

if (cnv.del.file!=0){
cnv<-rbind(cnv,dup) } else {cnv<-dup}

}


if (length(cnv)>0){
### CHECK THAT THIS PERSON HAS ANY CNV CALLS
if (indv.id %in% cnv$ID){
CNVCALLS<-1  
if(ped$V5[which(ped$V2 == indv.id)] == 1){ #male
	gender = "male"
	cnv[which(cnv[,gene_col] %in% XY_genelist),indv.id] = "hem" #hem in non-PAR region and Y

} else if (ped$V5[which(ped$V2 == indv.id)] == 2){#female
	gender = "female" 
        cnv[which(cnv[,"Chr"] == "Y" | cnv[,"Chr"] == "chrY"),indv.id] = "ref" #!!! add to include the chrY
	cnv[which(cnv[,indv.id]=="het" & cnv[,gene_col] %in% XY_genelist),indv.id] = "het_female" #het in non-PAR region
	cnv[which(cnv[,indv.id]=="hom" & cnv[,gene_col] %in% XY_genelist),indv.id] = "hom_female" #hom in non-PAR region	



      }

cnv<-cnv[which(cnv$ID==indv.id),]
Key <- paste(cnv$Chr,":",cnv$Start,"_",cnv$Ref,">",cnv$Alt, sep = "")
cnv<-cbind(Key,cnv)
cnv$FILTER<-"PASS"
}
}

########### END CNV PREP WORK 

# Read in exomes
print("read in exomes")



#READ IN ANNOTATIONS
annotations <-read.table(paste("annotated/",fam.id,".avinput.",build_ver,"_multianno_dup.txt",sep=""),sep="\t",stringsAsFactors=F,skip=1,quote = "")

header <-read.table(paste("annotated/",fam.id,".avinput.",build_ver,"_multianno_dup.txt",sep=""),sep="\t",stringsAsFactors=F,nrow=1,quote = "")
vcf.header <- read.table(paste(fam.id,".avinput.header",sep=""), sep="\t",stringsAsFactors=F,comment.char="@")
n.annos=ncol(header)

#CHECK THAT SPECIFIED INDIVIDUAL IS IN THE DATA
stopifnot(indv.id %in% vcf.header) 
names(annotations)<-header
#names(annotations)=c(header[-n.annos],vcf.header[1:8]) #take out the last column name of header (which is Otherinfo)	
start.anno<-n.annos-7
names(annotations[start.anno:n.annos])<-vcf.header[1:8]

#CLEAN UP OBJECT
annotations<-annotations[,-c(37,38,40:44)]

# For count for variants, comparison, finding grids later on
annotations$Key <- paste(annotations$Chr,":",annotations$Start,"_",annotations$Ref,">",annotations$Alt, sep = "")



#####READ IN INDIVIDUAL EXOME DATA FILE
exome.raw<-read.table(file=paste("avinput/",indv.id,".avinput",sep=""),stringsAsFactors=F,sep="\t",col.names=c("Chr","Start","End","Ref","Alt","c2","s2","ID","r2","a2","QUAL","FILTER","BLAH","INFO","GT"))

exome.raw$ID <- paste(exome.raw$Chr,":",exome.raw$Start,"_",exome.raw$Ref,">",exome.raw$Alt, sep = "")


exome.raw<-exome.raw[,c("ID","QUAL","FILTER","INFO","GT")]



## MERGING TO ACCOMMODATE SPLIT GENE (INSTEAD OF CBIND)
## ANNOTATIONS HAVE SPLIT GENE, RAW EXOME DOES NOT
#ONLY KEEP ANNOTATIONS FOR SITES IN THE EXOME BEING ANALYZED
annotations<-annotations[which(annotations$Key %in% exome.raw$ID),]
tmp<-merge(annotations,exome.raw,by.x=dim(annotations)[2],by.y=1,all.x=T)
exome.raw<-tmp

# !!! added to remove the duplication issue in CHET modes
exome.raw <- exome.raw[! duplicated(exome.raw[,c("Key",gene_col)]),]


colnames(exome.raw)[dim(exome.raw)[2]]<-indv.id

exome.raw$ID <- indv.id




## This pieces of codes are used to fix some problem from different VCF file sources
print("split ALT and REF")
# Add the columns for Alt and Total reads (additional codes)
exome.raw$Alt_Read <- sapply(exome.raw[,indv.id], function(x) unlist(strsplit(unlist(strsplit(x, split=":"))[2], split=","))[2])
exome.raw$Total_Read <- sapply(exome.raw[,indv.id], function(x) unlist(strsplit(x, split=":"))[3])	
  
## Extracts genotype info for specified individual
print("Extract genotypes")
a1 <- substr(exome.raw[,indv.id],1,1)
a2 <- substr(exome.raw[,indv.id],3,3)
exome.raw[indv.id] <- NA
exome.raw[which(a1 == "." | a2 == "."),indv.id] <- "miss" #add to remove the missing base
if(length(which(a1 != a2)) > 0){
	exome.raw[which(a1 != a2 & a1 != "." & a2 != "."),indv.id] = "het" # add condition to remove the missing base
}
if(length(which(!a1 %in% c(0,".") & !a2 %in% c(0,".") & a1 == a2)) > 0){
	exome.raw[which(!a1 %in% c(0,".") & !a2 %in% c(0,".") & a1 == a2),indv.id] = "hom"
}
if(length(which(a1 == 0 & a2 == 0)) > 0){
	exome.raw[which(a1 == 0 & a2 == 0),indv.id] = "ref"
}



if(ped$V5[which(ped$V2 == indv.id)] == 1){ #male
	gender = "male"
	exome.raw[which(exome.raw[,indv.id] %in% c("hom") & exome.raw[,gene_col] %in% XY_genelist),indv.id] = "hem" #hem in non-PAR region and Y
	exome.raw[which(exome.raw[,indv.id] %in% c("het") & exome.raw[,gene_col] %in% XY_genelist),indv.id] = "ref" #need to add this for male
} else if (ped$V5[which(ped$V2 == indv.id)] == 2){#female
	gender = "female" 
         exome.raw[which(exome.raw[,"Chr"] == "Y" | exome.raw[,"Chr"] == "chrY"),indv.id] = "ref" #!!! add to include the chrY
	exome.raw[which(exome.raw[,indv.id]=="het" & exome.raw[,gene_col] %in% XY_genelist),indv.id] = "het_female" #het in non-PAR region
	exome.raw[which(exome.raw[,indv.id]=="hom" & exome.raw[,gene_col] %in% XY_genelist),indv.id] = "hom_female" #hom in non-PAR region	

      }

######### MERGE CNV AND EXOME BEFORE CLEANING
if (CNVCALLS){
cnv<-cnv[,which(colnames(cnv) %in% colnames(exome.raw))]
exome.raw<-rbind(exome.raw,cnv)
}



### CLEAN DATA: 1) REMOVE BLACKLIST GENES, 2) REMOVE VARIANTS WITH AF DISCREPANCIES, 3) REMOVE GENES NOT INCLUDED IN LOOKUP TABLES, 4) REMOVE LINES WHERE ALL AFs ARE MISSING, 5) KEEP THE CODING AND SPLICING SITES, 5) REMOVE VARIANTS THAT DO NOT PASS TRANCHE FILTER 6) REMOVE VARIANTS THAT DO NOT PASS QUALITY FILTER OR HAVE MISSING pCADD SCORES
print("clean data process")
# Make score and gnomAD frequency information as numeric class
class(exome.raw[,score]) = "numeric"
class(exome.raw$gnomAD_exome_ALL) = "numeric"
# 1) REMOVE BLACKLISTED GENES; how about gnomAD low coverage list?
# bl <- scan(paste(dir,"psap/lookups/blacklist_081818.txt",sep=""),what="character") #Need to check file name
# bl.remove <- unique(c(which(exome.raw[,gene_col] %in% bl),grep("^HLA", exome.raw[,gene_col]), grep("^MUC", exome.raw[,gene_col]), grep("^KRT", exome.raw[,gene_col]), grep("^TRBV", exome.raw[,gene_col])))
bl.remove <- unique(c(grep("^HLA", exome.raw[,gene_col]), grep("^MUC", exome.raw[,gene_col]), grep("^KRT", exome.raw[,gene_col]), grep("^TRBV", exome.raw[,gene_col])))
# 2) REMOVE AF DISCREPANCIES (ANYTHING THAT IS MISSING IN ExAC BUT PRESENT IN 1000GP OR ESP AT GREATER THAN 5% FREQUENCY
af.remove <- which(is.na(exome.raw$gnomAD_exome_ALL) == T & exome.raw[,"1000g2015aug_all"] > 0.05 | is.na(exome.raw$gnomAD_exome_ALL) == T & exome.raw$esp6500siv2_all > 0.05)
# 3) REMOVE GENES NOT IN LOOKUP TABLES
lookup.remove <- which(! exome.raw[,gene_col] %in% lookup.genes)
# 4) REMOVE LINES WHERE ALL AFs ARE MISSING
af <- ped$V2[which(ped$V6>0)]
af<-af[which(af %in% names(exome.raw))]
missing.remove<-NULL
if (length(af)>0){
  missing.remove <- which(apply(exome.raw[af],1,function(row) return(sum(is.na(row)))) == length(af))}


# Remove these variants
tmp.exome <- exome.raw[-unique(c(bl.remove,af.remove,lookup.remove,missing.remove)),]
# 5) JUST RETAINING PROTEIN CODING SITES AND SPLICE SITES
keep <- unique(c(grep("splic",tmp.exome[,func_col]),which(is.na(tmp.exome[,loc_col])==FALSE)))
exome <- tmp.exome[keep,]

# 5b) SCORE INDELS (provide indel scores before the last filter step)
#lookup.lof = read.table(file=paste(dir,"psap/lookups/full.lof.CADD13phred.gencodeV19.allsites_031318.txt.gz",sep=""),stringsAsFactors=F)
#lookup.lof = read.table(file=paste(dir,"lookups/LOF_lookup.CADD13phred_031318.txt.gz",sep=""),stringsAsFactors=F)
#indels = grep("^frameshift",exome[,loc_col])
#gene.index = as.integer(factor(exome[,gene_col][indels],levels=lookup.lof[,1]))
#exome[,score][indels] = lookup.lof[gene.index,2]

if (is.na(indel.file)==FALSE){
#### LIINA INDEL CADD SCORE FIX
# 5b LN) UPDATE INDEL SCORING FOR INDELS. CADD scores manually retrieved for all INDELs via CADD website (v1.3).
indel.scores = read.table(indel.file, sep='\t', quote='', header=T, comment.char='',skip=1)
if (any(! names(indel.scores)[1:3]==c("X.CHROM","POS","REF"))){

print("Error with CADD indel input file format\n Please see PSAP documentation for details\n Expecting second line to
	      contain #CHROM POS REF ALT RawScore PHRED\n")

}


#Try to detect whether the contig naming system for INDELs is same as VCF. 
ichr<-grep("chr",indel.scores$X.CHROM)
vchr<-grep("chr",exome$Chr)

#ADD "CHR" TO INDEL CONTIG NAMES SINCE THIS CONVENTION IS USED IN THE VCF
if (length(vchr)>0 & length(ichr)==0){
  indel.scores$X.CHROM<-paste("chr",indel.scores$X.CHROM,sep="")
}


# Try to detect whether INDEL allele encoding is same as VCF. annotations from 0-based to 1-based to match the PSAP output
adash<-grep("-",indel.scores$ALT)
rdash<-grep("-",indel.scores$REF)

if (length(adash)==0 & length(rdash)==0){

# Remove the first character in Ref and Alt to match the PSAP output. Add +1 nt to Pos of deletions, insertions unchanged
indel.scores$REF <- as.character(indel.scores$REF)
indel.scores$ALT <- as.character(indel.scores$ALT)
indel.scores$POS <- as.integer(apply(indel.scores, 1, function(x) ifelse(nchar(x[3])>1,as.numeric(x[2])+1,x[2])))
indel.scores$REF <- apply(indel.scores, 1, function(x) ifelse(nchar(x[3])>1,substring(x[3],2),"-"))
indel.scores$ALT <- apply(indel.scores, 1, function(x) ifelse(nchar(x[4])>1,substring(x[4],2),"-"))
indel.scores <- unique(indel.scores[,c("X.CHROM","POS","REF","ALT","RawScore","PHRED")])
}

drops <- c("CADD13_RawScore","CADD13_PHRED")
indels = exome[grep("^frameshift",exome[,loc_col]), !names(exome) %in% drops]
indels = merge(indels, indel.scores, by.x=c("Chr","Start","Ref","Alt"), by.y=c("X.CHROM","POS","REF","ALT"), all.x=T, all.y=F)
names(indels)[names(indels)=="RawScore"] <- "CADD13_RawScore"
names(indels)[names(indels)=="PHRED"] <- score
indels <- indels[names(exome)]
exome <- rbind(exome[-grep("^frameshift",exome[,loc_col]),], indels)
}


# Normalize scores for grids later on
MAX <- 62
MIN <- 0
exome[,score][which(exome[,score] > MAX)] <- MAX
exome[,score][which(exome[,score] < MIN)] <- MIN
exome$score_norm <- (exome[,score]-MIN)/(MAX-MIN)


# 6) REMOVE VARIANTS THAT DO NOT PASS QUALITY FILTER OR HAVE MISSING pCADD SCORES
info <- exome[which(exome$FILTER=="PASS" & is.na(exome[,score]) == F | exome$FILTER=="." & is.na(exome[,score]) == F),
c("ID","Key","Chr","Start","Ref","Alt","Total_Read","Alt_Read",paste(c("Gene.","Func.","ExonicFunc.","AAChange.","GeneDetail."),gencode_ver,sep=""),"gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_NFE","ExAC_ALL","1000g2015aug_all","esp6500siv2_all","score_norm",score,"CLNSIG","CLNDN",indv.id)]



#TRACK THE DATA/DATA NOT INCLUDED IN ANY OF THE ABOVE ANALYSES AND OUTPUT IN THE END
id.raw <- exome.raw$Key
id.final <- info$Key
not_include <- unique(exome.raw[which(! id.raw %in% id.final),c("ID","Key","Chr","Start","Ref","Alt","Total_Read","Alt_Read",paste(c("Gene.","Func.","ExonicFunc.","AAChange.","GeneDetail."),gencode_ver,sep=""),"gnomAD_exome_ALL","gnomAD_exome_AFR","gnomAD_exome_NFE","ExAC_ALL","1000g2015aug_all","esp6500siv2_all",score,"CLNSIG","CLNDN",indv.id)])

# Remove information that will not be further used
rm(list=c("keep","exome","tmp.exome","exome.raw","af.remove","lookup.remove","bl.remove","lookup.genes"))
# Replace the column  
tmp <- info[-which(names(info) == indv.id)]
tmp["Geno"] <- info[,which(names(info) == indv.id)]
print(dim(tmp))

### START PSAP ANNOTATION
print("data cleaned, beginning PSAP annotation")
out <- data.frame()
## SOURCES CODE THAT WILL FORMAT AND ANALYZE THE DATA FOR EACH MODE OF INHERITANCE (AD, AR, CHET or HEM) MODEL
source(paste(dir,"RScripts/apply_PSAP.R",sep=""))
print(dim(out))  
# Remove non-necessary columns
out <- out[which(!names(out) %in% c("i","j","score_norm"))]
# reorder
out <- out[order(out$Dz.Model, out[,gene_col]),]


## WRITE OUTPUT
print("writing file")
print(dim(out))
# Output worst variants in each gene with PSAP annotation
if (Force == "FALSE") {
	write.table(out,paste("annotated/",fam.id,"_",indv.id,"_PSAP.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
} else {
	write.table(out,paste("annotated/",fam.id,"_",indv.id,"_PSAP_Force.txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
}
# Output the variants being filtered out in the cleaning process
write.table(not_include,file=paste("annotated/",fam.id,"_",indv.id,"_missing_data.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=FALSE)
