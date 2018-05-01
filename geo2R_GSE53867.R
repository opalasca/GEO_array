# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Mon Apr 30 04:11:11 EDT 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

acc="GSE53867"
gset <- getGEO("GSE53867", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL14149", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "20101"
sml <- c()
labels <- c("CD","UC","Control")
conditions <- c("CD","UC","C")
for (i in 1:nchar(gsms)) { sml[i] <- conditions[as.integer(substr(gsms,i,i))+1] }

fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
#cont.matrix <- makeContrasts(G2-G0, G1-G0, G1-G2, levels=design)
cont.matrix <- makeContrasts(UC-CD, C-UC, C-CD, levels=design)
cont_names<-c("UC_vs_CD", "C_vs_UC", "C_vs_CD")
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

for (i in 1:length(cont_names)) {
  tname=paste("T",cont_names[i],sep="_")
  t=data.frame()
  t=topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=5000)
  #t$Species.Scientific.Name[t$Species.Scientific.Name == "Homo sapiens"] <- "AA Homo sapiens"
  #t<-t[order(t$Sequence,t$Species.Scientific.Name),]
  #t<-t[!duplicated(t$Sequence),]
  #t<-t[t$adj.P.Val<0.9,]
  #t<-t[order(-t$B),]
  assign(tname, t)
  write.table(t, file=paste(acc,"_",tname,".txt",sep=""), row.names=F, sep="\t")
}

summ <- summary(decideTests(fit2, adjust.method="BH", p.value=0.05))
summ

#dim(exprs(gset))
#head(exprs(gset))
#dim(pData(gset))
#t<-pData(gset)
