# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sun Apr 15 08:55:48 EDT 2018


# http://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day1/bioc-intro.pdf

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
acc="GSE48957"
gset <- getGEO("GSE48957", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL14613", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
#gsms <- "222222222200000000001111111"
gsms <- "000000000011111111112222222"
#           G0          G1      G2
sml <- c()
labels <- c("Control","UC_active","UC_inactive")
conditions <- c("C","UCa","UCi")
for (i in 1:nchar(gsms)) { sml[i] <- conditions[as.integer(substr(gsms,i,i))+1] }

#for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
#sml <- paste("G", sml, sep="")    # set group names

fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
#cont.matrix <- makeContrasts(G2-G0, G1-G0, G1-G2, levels=design)
cont.matrix <- makeContrasts(UCa-C, UCa-UCi, UCi-C, levels=design)
cont_names<-c("UCa_vs_C", "UCa_vs_UCi", "UCi_vs_C")
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

for (i in 1:length(cont_names)) {
  tname=paste("T",cont_names[i],sep="_")
  t=data.frame()
  t=topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=5000)
  t$Species.Scientific.Name[t$Species.Scientific.Name == "Homo sapiens"] <- "AA Homo sapiens"
  t<-t[order(t$Sequence,t$Species.Scientific.Name),]
  t<-t[!duplicated(t$Sequence),]
  t<-t[t$adj.P.Val<0.1,]
  t<-t[order(-t$B),]
  assign(tname, t)
  write.table(t, file=paste(acc,"_",tname,".txt",sep=""), row.names=F, sep="\t")
}
#tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","miRNA_ID_LIST","SPOT_ID"))

summ <- summary(decideTests(fit2, adjust.method="BH", p.value=0.05))
summ


#dim(exprs(gset))
#head(exprs(gset))
#dim(pData(gset))
#t<-pData(gset)

probe_boxplot <- function(myrow){
  mygene <- exprs(gset)[myrow,]
  groups <- pData(gset)$characteristics_ch1.2
  df <- data.frame(values = mygene,
                   vars = groups)
  print(groups)
  par(cex.axis=0.8) 
  boxplot(values ~ vars, data = df)
}
#probe_boxplot("characteristics_ch1.1","SNORD123_st")
probe_boxplot("SNORD123_st")
probe_boxplot("hsa-miR-378c_st")



