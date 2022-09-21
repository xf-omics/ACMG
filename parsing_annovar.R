# ACMG annotations based on InterVar

#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Please provide input file name without the file extension")
} else {
	filename = args[1]
}

options(warn=-1)

# ANNOVAR
annovar<-read.csv(paste(filename,"hg19_multianno.txt",sep="."),header=TRUE,quote="",sep="\t",stringsAsFactors=FALSE)

# Function
annovar["Func"] = 0
annovar$Func[grep("synonymous",annovar$ExonicFunc.refGene)] = 1
annovar$Func[annovar$ExonicFunc.refGene=="unknown"] = 2.5
annovar$Func[annovar$ExonicFunc.refGene=="nonsynonymous SNV"] = 3
annovar$Func[grep("splicing", annovar$Func.refGene)] = 4
annovar$Func[grep("stop|start", annovar$ExonicFunc.refGene)] = 5
annovar$Func[grep("frameshift", annovar$ExonicFunc.refGene)] = 6
annovar$Func[grep("nonframeshift", annovar$ExonicFunc.refGene)] = 2

# function filtering
annovar$MAX_SPLICE[annovar$MAX_SPLICE=='.']=0
annovar$MAX_SPLICE = as.numeric(annovar$MAX_SPLICE)

# ACMG 28 pieces of evidence
for (i in 40:67) {
  annovar[annovar[,i]=='.',i]=0
  annovar[,i] <- as.numeric(annovar[,i])
}


# null variant in a gene where LOF is a known mechanism of disease
annovar$PVS1 = 0
annovar$PVS1[annovar$Func>3] = 1

# Same amino acid change as a previously established pathogenic variant
annovar$PS1 = 0
annovar$PP5 = 0
annovar$PS1[grep("athogenic", annovar$CLNSIG)] = 1
annovar$PS1[grep("Conflicting",annovar$CLNSIG)] = 0
a = intersect(grep("single_|no_assertion_", annovar$CLNREVSTAT), which(annovar$PS1==1))
annovar$PP5[a] = 1
annovar$PS1[a] = 0
annovar$PP5[annovar$PS1==1] = 0

# Hot spot or functional domain
a = grep("[a-z]", annovar$Interpro_domain)
annovar$PM1[a] = 1

# frequency
for (i in c(11:18,23)) {
  annovar[annovar[,i]=='.',i]=0
  annovar[,i] <- as.numeric(annovar[,i])
  annovar[annovar[,i]>0.5,i] = 1 - annovar[annovar[,i]>0.5,i]
}
annovar["gnomAD_max"] = apply(annovar[,c(11:18,23)], 1, max)

annovar$PM2[annovar$gnomAD_max<0.0001] = 1
annovar$PM2[annovar$gnomAD_max>=0.0001] = 0

# Multiple lines of computational evidence support a deleterious effect
annovar$CADD14_PHRED[annovar$CADD14_PHRED=='.'&annovar$Func>3]=1
annovar$CADD14_PHRED[annovar$CADD14_PHRED=='.'&annovar$Func<=3]=0
annovar$CADD14_PHRED = as.numeric(annovar$CADD14_PHRED)
annovar$REVEL[annovar$REVEL=='.'&annovar$Func>3]=1
annovar$REVEL[annovar$REVEL=='.'&annovar$Func<=3]=0
annovar$REVEL = as.numeric(annovar$REVEL)
annovar$PP3[annovar$CADD14_PHRED>=20] = annovar$PP3[annovar$CADD14_PHRED>=20]+1
annovar$PP3[annovar$REVEL>=0.5] = annovar$PP3[annovar$REVEL>=0.5]+1
annovar$PP3[annovar$MAX_SPLICE>0.5] = annovar$PP3[annovar$MAX_SPLICE>0.5] + 1
annovar$PP3[annovar$PP3==1] = 0
annovar$PP3[annovar$PP3>1] = 1

# BA1: Allele frequency is >5% in Exome Sequencing Project, 1000 Genomes Project, or Exome Aggregation Consortium
annovar$BA1 = 0
annovar$BA1[annovar$gnomAD_max>=0.05] = 1

# BS1: Allele frequency is greater than expected for disorder
annovar$BS1 = 0
annovar$BS1[annovar$gnomAD_max>=0.01&annovar$BA1==0] = 1

# Well-established in vitro or in vivo functional studies how no damaging effect
annovar$BS3[grep("enign", annovar$CLNSIG)] = 1

# Multiple lines of computational evidence suggest no impact
annovar$BP4[annovar$CADD14_PHRED>=20] = 0
annovar$BP4[annovar$REVEL>=0.5] = 0
annovar$BP4[annovar$MAX_SPLICE>=0.5] = 0
annovar$BP4[annovar$REVEL<0.3&annovar$CADD14_PHRED<10&annovar$MAX_SPLICE<0.5&annovar$Func<4] = 1
# a synonymous variant has no impact to splicing
annovar$BP7[annovar$Func==1&annovar$MAX_SPLICE<0.5] = 1

# aggregate evidence
i=41
annovar["PS"] = rowSums(annovar[,i:(i+3)])
annovar["PM"] = rowSums(annovar[,(i+4):(i+9)])
annovar["PP"] = rowSums(annovar[,(i+10):(i+14)])

annovar["BS"] = rowSums(annovar[,(i+16):(i+19)])
annovar["BP"] = rowSums(annovar[,(i+20):(i+26)])

# pext
reference = "/share/terra/xf2193/genome/pext/chr/pext_chr"
filename=paste(filename,"var",sep=".")
colnames(annovar)[1:2] = c("CHROM","POS")
annovar = annovar[,c(1,2,4:10,38:67,71:73,75,88:96)]

for (i in unique(annovar$CHROM)) {
	file = paste0(reference,i)
	data <- read.csv(paste(file,"txt.avinput.bz2",sep='.'),header=FALSE,sep="\t",stringsAsFactors=FALSE)[,c(1,2,8)]
	colnames(data) = c("CHROM","POS","pext")
	data = data[data$pext!="NaN",]
	data$pext = as.numeric(data$pext)
	data = data[!is.na(data$pext),]
	data = merge(annovar[annovar$CHROM==i,],data,all.x=TRUE)
	data$pext[is.na(data$pext)] = 0.5
	data = unique(data)
# these commands below will change the order of columns
#        data1 = aggregate(data$pext,by=list(data$CHROM,data$POS),FUN="mean")
#        colnames(data1) = c("CHROM","POS","pext")
#	if (nrow(data)>nrow(data1)) {
#		data = merge(data,data1)
#	}

	data$PVS1[data$pext<0.1] = 0
	# data$BA1[data$pext<0.1] = data$BA1[data$pext<0.1] + 1
	data$BS[data$pext<0.2] = data$BS[data$pext<0.2] + 1
	data$BP[data$pext>=0.2&data$pext<0.5] = data$BP[data$pext>=0.2&data$pext<0.5] + 1

	# ACMG
	data$ACMG = "VUS"
	data$ACMG[data$BA1==1|data$BS>=2] = "benign"
	data$ACMG[(data$BS+data$BP)>=2&data$ACMG=="VUS"] = "likely benign"

	data$ACMG[data$PVS1==1&data$PS>=1&data$ACMG!="VUS"] = paste(data$ACMG[data$PVS1==1&data$PS>=1&data$ACMG!="VUS"],"pathogenic",sep='-')
	data$ACMG[data$PVS1==1&data$PS>=1&data$ACMG=="VUS"] = "pathogenic"

	data$ACMG[data$PVS1==1&(data$PM+data$PP)>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PVS1==1&(data$PM+data$PP)>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")],"pathogenic",sep='-')
	data$ACMG[data$PVS1==1&(data$PM+data$PP)>=2&data$ACMG=="VUS"] = "pathogenic"

	data$ACMG[data$PS>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PS>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")],"pathogenic",sep='-')
	data$ACMG[data$PS>=2&data$ACMG=="VUS"] = "pathogenic"

	data$ACMG[data$PS==1&data$PM>=3&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PS==1&data$PM>=3&(data$ACMG=="benign"|data$ACMG=="likely benign")],"pathogenic",sep='-')
        data$ACMG[data$PS==1&data$PM>=3&data$ACMG=="VUS"] = "pathogenic"

	data$ACMG[data$PS==1&data$PM>=2&data$PP>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PS==1&data$PM>=2&data$PP>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")],"pathogenic",seq='-')
        data$ACMG[data$PS==1&data$PM>=2&data$PP>=2&data$ACMG=="VUS"] = "pathogenic"

	data$ACMG[data$PS==1&data$PM>=1&data$PP>=4&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PS==1&data$PM>=1&data$PP>=4&(data$ACMG=="benign"|data$ACMG=="likely benign")],"pathogenic",sep='-')
	data$ACMG[data$PS==1&data$PM>=1&data$PP>=4&data$ACMG=="VUS"] = "pathogenic"

	data$ACMG[data$PVS1==1&data$PM==1&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PVS1==1&data$PM==1&(data$ACMG=="benign"|data$ACMG=="likely benign")],"likely pathogenic",sep='-')
	data$ACMG[data$PVS1==1&data$PM==1&data$ACMG=="VUS"] = "likely pathogenic"

	data$ACMG[data$PS==1&data$PM>=1&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PS==1&data$PM>=1&(data$ACMG=="benign"|data$ACMG=="likely benign")],"likely pathogenic",sep='-')
        data$ACMG[data$PS==1&data$PM>=1&data$ACMG=="VUS"] = "likely pathogenic"

        data$ACMG[data$PS==1&data$PP>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PS==1&data$PP>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")],"likely pathogenic",sep='-')
	data$ACMG[data$PS==1&data$PP>=2&data$ACMG=="VUS"] = "likely pathogenic"

	data$ACMG[data$PM>=3&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PM>=3&(data$ACMG=="benign"|data$ACMG=="likely benign")],"likely pathogenic",sep='-')
        data$ACMG[data$PM>=3&data$ACMG=="VUS"] = "likely pathogenic"

        data$ACMG[data$PM==2&data$PP>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PM==2&data$PP>=2&(data$ACMG=="benign"|data$ACMG=="likely benign")],"likely pathogenic",sep='-')
	data$ACMG[data$PM==2&data$PP>=2&data$ACMG=="VUS"] = "likely pathogenic"

	data$ACMG[data$PM==1&data$PP>=4&(data$ACMG=="benign"|data$ACMG=="likely benign")] = paste(data$ACMG[data$PM==1&data$PP>=4&(data$ACMG=="benign"|data$ACMG=="likely benign")],"likely pathogenic",sep='-')
	data$ACMG[data$PM==1&data$PP>=4&data$ACMG=="VUS"] = "likely pathogenic"

	data$ACMG[data$PS1==1&data$ACMG!="pathogenic"&data$ACMG!="likely pathogenic"] = paste(data$ACMG[data$PS1==1&data$ACMG!="pathogenic"&data$ACMG!="likely pathogenic"],"surrender",sep='-')

	OR = 4.4^data$PVS1 * 4.31^data$PS * 2.09^data$PM * 1.45^data$PP * 0.05^data$BA1 * 0.223^data$BS * 0.316^data$BP # OR
	data["Prob"] = OR / (OR + 1)
	
	if (i==1) {
		write.table(data,filename,row.names=FALSE,sep="\t")	
	} else {
		write.table(data,filename,row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
	}
}


