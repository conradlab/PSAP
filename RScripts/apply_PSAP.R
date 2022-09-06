## IDENTIFIES CANDIDATE VARIANTS FOR AD, AR, CHET, HEM
#tmp$f <- NA # for assign frequency for HW calculation 
tmp$Dz.Model <- NA #for assign Disease models
##For AD (Het)
##subset the exome
exome.AD <- subset(tmp, Geno == "het"|Geno == "het_female")
if (nrow(exome.AD)>0) {
	## add disease mode information #modify for gender pipeline
	exome.AD$Dz.Model[which(exome.AD$Geno == "het")] <- "DOM-het"
	exome.AD$Dz.Model[which(exome.AD$Geno == "het_female")] <- "X-linked-het"
	##ANNOTATE DATASETS WITH PSAP
	# Read in lookup table
	lookup <-read.table(paste(dir,"lookups_CADD16/full.het.CADD16phred.gencodeV35.allsites",eth,"_030122.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE)
	# GETS COLUMN (SCORE INTERVAL) INFO FOR LOOKUP TABLE
	exome.AD$j<-findInterval(exome.AD[,"score_norm"],scale)+1
	# find row
	exome.AD$i<-as.integer(factor(exome.AD[,gene_col],levels=lookup[,1]))
	# PSAP p(C > c)
	exome.AD$PSAP<-unlist(apply(exome.AD,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
	if (Force == "FALSE") {
		###pick the worst one
		tmp.AD <- do.call(rbind,by(exome.AD,exome.AD[,gene_col],function(x) x[which.max(as.numeric(x[,"score_norm"])),]))
		out <- rbind(out, tmp.AD)
	} else {
		out <- rbind(out, exome.AD)
	}
	## remove information that will not be used anymore
	rm(list = c("exome.AD"))
}
print("AD model complete")
##For AR (Hom)
##subset the exome
exome.AR <- subset(tmp, Geno == "hom"|Geno == "hom_female")
if (nrow(exome.AR)>0) {
	## Add disease mode information
	exome.AR$Dz.Model[which(exome.AR$Geno == "hom")] <- "REC-hom"
	exome.AR$Dz.Model[which(exome.AR$Geno == "hom_female")] <- "X-linked-hom"
	##ANNOTATE DATASETS WITH PSAP
	# Read in lookup table
	lookup <-read.table(paste(dir,"lookups_CADD16/full.hom.CADD16phred.gencodeV35.allsites",eth,"_030122.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE)
	# GETS COLUMN (SCORE INTERVAL) INFO FOR LOOKUP TABLE
	exome.AR$j<-findInterval(exome.AR[,"score_norm"],scale)+1
	# find row
	exome.AR$i<-as.integer(factor(exome.AR[,gene_col],levels=lookup[,1]))
	# PSAP p(C > c)
	exome.AR$PSAP<-unlist(apply(exome.AR,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
	# whether output all variants or not
	if (Force == "FALSE") {
		#pick the worst one
		tmp.AR <- do.call(rbind,by(exome.AR,exome.AR[,gene_col],function(x) x[which.max(as.numeric(x[,"score_norm"])),]))
		out <- rbind(out, tmp.AR)
	} else {
		out <- rbind(out, exome.AR)
	}
	## remove information that will not be used anymore
	rm(list = c("exome.AR"))
}
print("REC-hom model complete")
##For CHET (chet)
##subset the exome
exome.CHET <- subset(tmp, Geno == "het"|Geno == "het_female") # For CHET: just select genes with more than one variants
dups <- names(table(exome.CHET[,gene_col])[which(table(exome.CHET[,gene_col]) >= 2)])
exome.CHET <- exome.CHET[exome.CHET[,gene_col] %in% dups,]
if (nrow(exome.CHET) > 0) {
	#change genotype from het to chet
	exome.CHET$Geno[which(exome.CHET$Geno == "het")] <- "chet"
	exome.CHET$Geno[which(exome.CHET$Geno == "het_female")] <- "chet_female"
	## Add disease mode information
	exome.CHET$Dz.Model[which(exome.CHET$Geno == "chet")] = "REC-chet"
	exome.CHET$Dz.Model[which(exome.CHET$Geno == "chet_female")] = "X-linked-chet"
	##ANNOTATE DATASETS WITH PSAP
	# Read in lookup table
	lookup <-read.table(paste(dir,"lookups_CADD16/full.chet.CADD16phred.gencodeV35.allsites",eth,"_030122.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE)
	# GETS COLUMN (SCORE INTERVAL) INFO FOR LOOKUP TABLE
	exome.CHET$j<-findInterval(exome.CHET[,"score_norm"],scale)+1
	# find row
	exome.CHET$i<-as.integer(factor(exome.CHET[,gene_col],levels=lookup[,1]))
	# PSAP p(C > c)
	exome.CHET$PSAP<-unlist(apply(exome.CHET,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
	# whether output all variants or not
	if (Force == "FALSE") {
		###pick the worst pair #Need to test
		tmp.CHET <- do.call(rbind,by(exome.CHET,exome.CHET[,gene_col],function(x,score) { x[,"score_norm"]<-as.numeric(x[,"score_norm"]); x<-x[order(x[,"score_norm"],decreasing=T,na.last=NA),]; return(x[1:2,]) },score="score_norm")) 
		# Assign the second worst PSAP values to the CHET pair
		genes <- unique(tmp.CHET[,gene_col])
		for (gene in genes){
			tmp.CHET$PSAP[which(tmp.CHET[,gene_col] == gene)] <- tmp.CHET$PSAP[which(tmp.CHET[,gene_col] == gene)][2]	
		}			
		out <- rbind(out, tmp.CHET)
	} else {
		out <- rbind(out, exome.CHET)
	}
	## remove information that will not be used anymore
	rm(list = c("exome.CHET"))
}
print("REC-chet model complete")
## For male: additional mode HEM
if (gender == "male"){
	##subset the control exome
	exome.HEM <- subset(tmp, Geno == "hem") # for HEM mode
	if (nrow(exome.HEM) > 0) {
		## add disease mode information
		exome.HEM$Dz.Model[which(exome.HEM$Geno == "hem" & exome.HEM$Chr == "X")] <- "X-linked-hem"
		exome.HEM$Dz.Model[which(exome.HEM$Geno == "hem" & exome.HEM$Chr == "Y")] <- "Y-linked-hem"	
		 exome.HEM$Dz.Model[which(exome.HEM$Geno == "hem" & exome.HEM$Chr == "chrX")] <- "X-linked-hem" #!!! Wu-Lin add this line April 2019
                exome.HEM$Dz.Model[which(exome.HEM$Geno == "hem" & exome.HEM$Chr == "chrY")] <- "Y-linked-hem" #!!! Wu-Lin add this line April 2019
	

		##ANNOTATE DATASETS WITH PSAP
		# Read in lookup table
		lookup <-read.table(paste(dir,"lookups_CADD16/full.hem.CADD16phred.gencodeV35.allsites",eth,"_030122.txt.gz",sep=""),stringsAsFactors=FALSE,header=FALSE)
		# GETS COLUMN (SCORE INTERVAL) INFO FOR LOOKUP TABLE
		exome.HEM$j<-findInterval(exome.HEM[,"score_norm"],scale)+1
		# find row
		exome.HEM$i<-as.integer(factor(exome.HEM[,gene_col],levels=lookup[,1]))
		# PSAP p(C > c)
		exome.HEM$PSAP<-unlist(apply(exome.HEM,1,function(x,tab) {i<-as.numeric(x["i"]); j<-as.numeric(x["j"]); y<-tab[i,j]; return(y)},lookup))
		# whether output all variants or not
		if (Force == "FALSE") {
			#pick the worst one
			tmp.HEM <- do.call(rbind,by(exome.HEM,exome.HEM[,gene_col],function(x) x[which.max(as.numeric(x[,"score_norm"])),]))
			out <- rbind(out, tmp.HEM)
		} else {
			out <- rbind(out, exome.HEM)
		}
	}
print("HEM model complete")
}










