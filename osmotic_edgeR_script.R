library(limma)
library(edgeR)
library(magrittr)
#library(NMF)
library(stringr)

##read in counts, alter/fix column names. 
fcc <- read.table("/home/nreid/rnaseq/featurecounts/osmotic.counts.meta.fc",stringsAsFactors = FALSE,header = TRUE)
	colnames(fcc) <- gsub("X.home.nreid.rnaseq.alignments.merged.","",colnames(fcc))
	colnames(fcc) <- gsub(".bam","",colnames(fcc))

##separate counts from metadata
fcc2 <- as.matrix(fcc[,7:136])

##table of factors. process them. 
tags <- read.table("/home/nreid/rnaseq/bwa.bams.list")
	tags <- gsub("/home/nreid/rnaseq/alignments/merged/","",tags[,1])
	tags <- gsub(".bam","",tags)
	tags2 <- str_extract(tags,".$")
	tags3 <- str_extract(tags,regex("[^_]+(?=_.$)"))
	tags4 <- str_extract(tags,regex(".*(?=_[^_]+_.$)"))
	tags5 <- c("m","f","b","b","m","m","m","f","f","f","b","f","f","m","m","f","m")
		names(tags5) <- unique(tags4)

	tags <- cbind(tags4,tags3,tags5[tags4],tags2)
	colnames(tags) <- c("species","treatment","physiology","replicate")

	facs <- paste(tags[,1],tags[,2],sep = ".")


##image counts
image(log(fcc2[order(rowSums(fcc2),decreasing = T),]))

##remove misidentified fish:
	#F. diaphanus FW 1 is actually a heteroclitus.
	#F. sciadicus BW 2 is actually a chrysotus.
	toss <- which(colnames(fcc2)%in%c("F_diaphanus_FW_1","F_sciadicus_BW_2"))
	fcc <- fcc[,-(toss+6)]
	fcc2 <- fcc2[,-toss]
	tags <- tags[-toss,]
	facs <- facs[-toss]


##create annotation matrix
annot <- fcc[,1:6]
annot[,2] <- gsub(";.*","",annot[,2])
annot <- annot[,-(3:5)]
	
#create design matrix
TS <- factor(paste(tags[,1],tags[,2],sep = "."))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)

###start analysis. turn counts into DGEList object. 
fcc3 <- DGEList(counts = fcc2, genes = annot, group = TS)

##toss out some low coverage genes
keep <- rowSums(cpm(fcc3))>20
fcc3 <- fcc3[keep,,keep.lib.sizes=FALSE]
dim(fcc3)

##TMM normalization. creates sample weights as object in fcc3 list. 
fcc3 <- calcNormFactors(fcc3)
image(log(fcc3$counts[order(rowSums(fcc3$counts),decreasing = T),]))

##MDS plot of samples for 20 most variable genes
plotMDS(fcc3,labels = colnames(fcc2),top = 10000, col = as.numeric(factor(tags[,3])), gene.selection = "common",prior.count = 5)

##estimate GLM dispersions
fcc3 <- estimateGLMCommonDisp(fcc3, design, verbose=TRUE)
fcc3 <- estimateGLMTrendedDisp(fcc3, design)
fcc3 <- estimateGLMTagwiseDisp(fcc3, design)

##fit glm model
fit <- glmFit(fcc3,design)
fitQL <- glmQLFit(fcc3, design, robust=TRUE)

##contrasts to be evaluated.
##YOU FORGOT ZEBRINUS!
cont.matrix <- makeContrasts(
 	MDPL_FvB = (F_heteroclitus_MDPL.FW-F_heteroclitus_MDPL.BW),
 	MDPP_FvB = (F_heteroclitus_MDPP.FW-F_heteroclitus_MDPP.BW),
 	PLvPP_FvB =(F_heteroclitus_MDPL.FW-F_heteroclitus_MDPL.BW)-(F_heteroclitus_MDPP.FW-F_heteroclitus_MDPP.BW),
 	M_FvB = (A_xenica.FW-A_xenica.BW)+(F_grandis.FW-F_grandis.BW)+(F_heteroclitus_MDPL.FW-F_heteroclitus_MDPL.BW)+(F_heteroclitus_MDPP.FW-F_heteroclitus_MDPP.BW)+(F_similis.FW-F_similis.BW)+(L_parva.FW-L_parva.BW),
 	M_FvT = (A_xenica.FW-A_xenica.transfer)+(F_grandis.FW-F_grandis.transfer)+(F_heteroclitus_MDPL.FW-F_heteroclitus_MDPL.transfer)+(F_heteroclitus_MDPP.FW-F_heteroclitus_MDPP.transfer)+(F_similis.FW-F_similis.transfer)+(L_parva.FW-L_parva.transfer),
 	M_treatment = (A_xenica.FW-A_xenica.BW)+(F_grandis.FW-F_grandis.BW)+(F_heteroclitus_MDPL.FW-F_heteroclitus_MDPL.BW)+(F_heteroclitus_MDPP.FW-F_heteroclitus_MDPP.BW)+(F_similis.FW-F_similis.BW)+(L_parva.FW-L_parva.BW)+(A_xenica.FW-A_xenica.transfer)+(F_grandis.FW-F_grandis.transfer)+(F_heteroclitus_MDPL.FW-F_heteroclitus_MDPL.transfer)+(F_heteroclitus_MDPP.FW-F_heteroclitus_MDPP.transfer)+(F_similis.FW-F_similis.transfer)+(L_parva.FW-L_parva.transfer),
 	F_FvB = (F_catanatus.FW-F_catanatus.BW) + (F_notatus.FW-F_notatus.BW) + (F_olivaceous.FW-F_olivaceous.BW) + (F_rathbuni.FW-F_rathbuni.BW) + (F_sciadicus.FW-F_sciadicus.BW) + (L_goodei.FW-L_goodei.BW),
 	F_FvT = (F_catanatus.FW-F_catanatus.transfer) + (F_notatus.FW-F_notatus.transfer) + (F_olivaceous.FW-F_olivaceous.transfer) + (F_rathbuni.FW-F_rathbuni.transfer) + (F_sciadicus.FW-F_sciadicus.transfer) + (L_goodei.FW-L_goodei.transfer),
 	FvB = ((A_xenica.FW-A_xenica.BW)+(F_grandis.FW-F_grandis.BW)+(F_heteroclitus_MDPL.FW-F_heteroclitus_MDPL.BW)+(F_heteroclitus_MDPP.FW-F_heteroclitus_MDPP.BW)+(F_similis.FW-F_similis.BW)+(L_parva.FW-L_parva.BW) + (F_catanatus.FW-F_catanatus.BW) + (F_notatus.FW-F_notatus.BW) + (F_olivaceous.FW-F_olivaceous.BW) + (F_rathbuni.FW-F_rathbuni.BW) + (F_sciadicus.FW-F_sciadicus.BW) + (L_goodei.FW-L_goodei.BW))/12,
 	MvB_FvB = ((A_xenica.FW-A_xenica.BW)+(F_grandis.FW-F_grandis.BW)+(F_heteroclitus_MDPL.FW-F_heteroclitus_MDPL.BW)+(F_heteroclitus_MDPP.FW-F_heteroclitus_MDPP.BW)+(F_similis.FW-F_similis.BW)+(L_parva.FW-L_parva.BW))/6 - ((F_catanatus.FW-F_catanatus.BW) + (F_notatus.FW-F_notatus.BW) + (F_olivaceous.FW-F_olivaceous.BW) + (F_rathbuni.FW-F_rathbuni.BW) + (F_sciadicus.FW-F_sciadicus.BW) + (L_goodei.FW-L_goodei.BW))/6,
 	levels = design
 	)


##evaluate all contrasts. store results. FOR LRT
ncont <- dim(cont.matrix)[2]
lrt <- list()
for(i in 1:ncont){
	lrt[[i]] <- glmLRT(fit, contrast=cont.matrix[,i])
	cat(i," ")
	}

outvals <- list()
outvals[["FDR"]] <- numeric(0)
outvals[["logFC"]] <- numeric(0)
outvals[["logCPM"]] <- numeric(0)

for(i in 1:ncont){
	tmp <- topTags(lrt[[i]],sort.by="none",n=dim(fcc3)[1])
	outvals[["FDR"]] <- cbind(outvals[["FDR"]],tmp$table$FDR)
	outvals[["logFC"]] <- cbind(outvals[["logFC"]],tmp$table$logFC)
	outvals[["logCPM"]] <- cbind(outvals[["logCPM"]],tmp$table$logCPM)
	}

outvals[["sig"]] <- rowSums(outvals[["FDR"]]<0.05)>0

resum <- table(apply(outvals[["FDR"]]<0.05,MAR=1,FUN=function(x){paste(as.numeric(x),collapse=".")}))
resum <- cbind(do.call(rbind,strsplit(split="\\.",names(resum))),as.numeric(resum))
class(resum)<-"numeric"
resum<-resum[order(resum[,ncont+1],decreasing=TRUE),]
cbind(colnames(cont.matrix),colSums(outvals[["FDR"]]<0.05))


###evaluate all contrasts. store results. FOR QL F test

qlf <- list()
for(i in 1:ncont){
	qlf[[i]] <- glmQLFTest(fitQL, contrast=cont.matrix[,i])
	cat(i," ")
	}

outvalsQL <- list()
outvalsQL[["FDR"]] <- numeric(0)
outvalsQL[["logFC"]] <- numeric(0)
outvalsQL[["logCPM"]] <- numeric(0)

for(i in 1:ncont){
	tmp <- topTags(qlf[[i]],sort.by="none",n=dim(fcc3)[1])
	outvalsQL[["FDR"]] <- cbind(outvalsQL[["FDR"]],tmp$table$FDR)
	outvalsQL[["logFC"]] <- cbind(outvalsQL[["logFC"]],tmp$table$logFC)
	outvalsQL[["logCPM"]] <- cbind(outvalsQL[["logCPM"]],tmp$table$logCPM)
	}

outvalsQL[["sig"]] <- rowSums(outvalsQL[["FDR"]]<0.05)>0

resumQL <- table(apply(outvalsQL[["FDR"]]<0.05,MAR=1,FUN=function(x){paste(as.numeric(x),collapse=".")}))
resumQL <- cbind(do.call(rbind,strsplit(split="\\.",names(resumQL))),as.numeric(resumQL))
class(resumQL)<-"numeric"
resumQL<-resumQL[order(resumQL[,ncont+1],decreasing=TRUE),]
cbind(colnames(cont.matrix),colSums(outvalsQL[["FDR"]]<0.05),colSums(outvals[["FDR"]]<0.05))




##create table of treatment mean CPM

	tmeans <- c()
	for(i in unique(tags[,1])){
		for(j in unique(tags[,2])){
			tmeans <- cbind(
				tmeans,
				rowSums(cpm(fcc3[,tags[,1]==i&tags[,2]==j]))/sum(tags[,1]==i&tags[,2]==j)
				)
			}	
		}
	
	tmeans <- tmeans[,colSums(is.nan(tmeans))<1]
		colnames(tmeans) <- unique(paste(tags[,1],tags[,2],sep="_"))


#save results
save.image(file="/home/nreid/rnaseq/osmoticALL.Rdata")
save(fit,tmeans,lrt,outvals,fcc3,tags,file="/home/nreid/rnaseq/edgeR_genes.Robj")


