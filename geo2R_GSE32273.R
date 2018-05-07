# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Apr 12 05:31:16 EDT 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

setwd("/Users/jmx139/Desktop/IBD/geo_microarray_analysis")
#load external functions
source(functions.R)

# load series and platform data from GEO
acc="GSE32273"
gset <- getGEO("GSE32273", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL8786", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("00000000000000000000XX11111111111111111111XX222222",
               "22222222222222XX33333333333333333333XX444444444444",
               "44444444XX55555555555555555555XX")

sml <- c()
labels <- c("UC_platelets","control_platelets", "UC_microvesicles","control_microvesicles", "UC_PBMC", "control_PBMC")
conditions <- c("UCpl","Cpl","UCm","Cm","UCP","CP")
for (i in 1:nchar(gsms)) { sml[i] <- conditions[as.integer(substr(gsms,i,i))+1] }
#for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

annot <- fData(gset)
annot <- annot[annot$SEQUENCE!="",]
annot <- subset(annot, select=c("ID","Species.Scientific.Name","SEQUENCE"))
write.table(annot, file=paste(acc,"_probe_annot.txt",sep=""), row.names=F, col.names=F,sep="\t")


# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# set up the data and proceed with analysis
#sml <- paste("G", sml, sep="")    # set group names
#gname <-("")
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
#cont.matrix <- makeContrasts(G1-G0, G3-G2, G5-G4, levels=design)
cont.matrix <- makeContrasts(UCpl-Cpl, UCm-Cm, UCP-CP, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
cont_names <-c("UCpl_vs_Cpl", "UCm_vs_Cm", "UCP_vs_CP")

#write_tables(cont_names, fit2, acc)
for (i in 1:length(cont_names)) {
  tname=paste("T",cont_names[i],sep="_")
  t=data.frame()
  t=topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=5000)
  t$Species.Scientific.Name[t$Species.Scientific.Name == "Homo sapiens"] <- "AA Homo sapiens"
  t<-t[order(t$SEQUENCE,t$Species.Scientific.Name),]
  t<-t[!duplicated(t$SEQUENCE),]
  t<-t[t$adj.P.Val<0.1,]
  t<-t[order(-t$B),]
  assign(tname, t)
  write.table(t, file=paste(acc,"_",tname,".txt",sep=""), row.names=F, sep="\t")
}

probe_boxplot <- function(myrow){
  mygene <- exprs(gset)[myrow,]
  #groups <- pData(gset)[,40]
  groups <- paste (pData(gset)[ ,40], pData(gset)[ ,41], sep = "" )
  df <- data.frame(values = mygene, vars = groups)
  print(groups)
  par(cex.axis=0.8) 
  boxplot(values ~ vars, data = df)
}
#probe_boxplot("characteristics_ch1.1","SNORD123_st")
probe_boxplot("hsa-miR-501-5p_st")
probe_boxplot("hsa-miR-27a-star_st")


summ <- summary(decideTests(fit2, adjust.method="BH", p.value=0.05))
summ

