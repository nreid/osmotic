library(limma)
library(edgeR)
library(magrittr)
library(NMF)
library(mscr)

cbind(colnames(cont.matrix),colSums(outvalsQL[["FDR"]]<0.03),colSums(outvals[["FDR"]]<0.03))

##plot some results. 
treat <- cbind(c(1,cumsum(table(paste(tags[,1], tags[,2])))[-48]+1),cumsum(table(paste(tags[,1], tags[,2])))) %>%
	apply(., MAR=1, FUN=function(x){x[1]:x[2]})
names(treat) <- NULL

#lowcov <- log(colSums(fcc3$counts),base=10)<5.25

subCm <- c()
for(i in colnames(design)){
		subi <- paste(tags[,1],tags[,2],sep=".")==i
		subCm <- cbind(
			subCm,
			rowSums(cpm(fcc3[,subi]))/sum(subi)
			)
		}	
	

colnames(subCm) <- colnames(design)

subC <- cpm(fcc3,log=TRUE,prior.count=1)

subi <- (!grepl("999",colnames(fcc3$counts)))
subi <- (!grepl("999",colnames(fcc3$counts)))&(!tags[,2]=="999")
subg <- which(outvalsQL[["FDR"]][,5]<0.03)
subg <- which(outvalsQL[["FDR"]][,5]<0.03&!outvalsQL[["FDR"]][,4]<0.05&outvalsQL[["logFC"]][,5]>0)
#subg <- which(rowSums(outvalsQL$FDR[fcc3x$prior.df<8,]<0.05)!=0)
#subg <- which(outvalsQL[["FDR"]][,3]>0.1&outvals[["FDR"]][,3]<0.05)
#subg <- which((outvals[["FDR"]][,1]<0.05)&(outvals[["FDR"]][,3]>0.05))
mCPM <- log(rowSums(cpm(fcc3[,subi],log=FALSE,prior.count=1))/sum(subi),base=2)


#fold change from CPM average per individual
subC2 <- subC[subg,subi]-mCPM[subg]
subC2[subC2>3] <- 3
subC2[subC2<(-3)] <- -3

ann <- as.data.frame(tags[subi,c(1,2,3)],stringsAsFactors=FALSE)
aheatmap(subC2,annCol=ann,distfun="pearson")

#fold change from CPM average per species, per individual

subC2 <- cpm(fcc3)[subg,subi]
tmp <- cpm(fcc3)
mCPM <- apply(tags[subi,],MAR=1,FUN=
	function(x){
		rowSums(tmp[subg,grep(x[1],colnames(fcc3))])/sum(grepl(x[1],colnames(subC2)))
		}
	)
subC2 <- log((subC2+5)/(mCPM+5),base=2)
subC2[subC2>3] <- 3
subC2[subC2<(-3)] <- -3

ann <- as.data.frame(tags[subi,c(3,2)][order(tags[subi,3],tags[subi,1],tags[subi,2]),],stringsAsFactors=FALSE)
aheatmap(subC2[,order(tags[subi,3],tags[subi,1],tags[subi,2])],annCol=ann,distfun="pearson",Colv=NA)


plotMDS(subC2,
labels = colnames(fcc2)[subi],top = 500, 
#pch=as.numeric(factor(tags[,2][subi]))+19,
col = as.numeric(factor(tags[,3][subi])), 
gene.selection = "common",cex=.75)


#fold change from CPM average per treatment
subi2 <- grepl("hetero",colnames(subCm))
mCPM <- log(rowSums(subCm[,subi2])/sum(subi2),base=2)
subC3 <- log(subCm[subg,subi2],base=2)-mCPM[subg]
subC3[subC3>2] <- 2
subC3[subC3<(-2)] <- -2

ann <- as.data.frame(tags[subi2,c(1,2,3)])
aheatmap(subC3,distfun="pearson",annCol=ann)


plotMDS(fcc3$counts[subg,subi],
labels = colnames(fcc2)[subi],top = 200, 
col = as.numeric(factor(tags[,3][subi])), 
gene.selection = "common")

aheatmap(subC[subg,subi],distfun="pearson")


vennDiagram(outvalsQL[["FDR"]][,c(4,5)]<0.03)


plot(cpm(fcc3)["30217",],col=factor(tags[,3]),pch=as.numeric(factor(tags[,2]))+19)
plot(cpm(fcc3,log=TRUE,prior.count=1)["30217",],col=factor(tags[,3]),pch=as.numeric(factor(tags[,2]))+19)

plotMDS(fcc3,labels = colnames(fcc2),top = 100, col = as.numeric(factor(tags[,4])), gene.selection = "common",prior.count = 5)


