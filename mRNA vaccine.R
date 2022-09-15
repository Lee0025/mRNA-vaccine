####This project was used to fliter BRCA mRNA vaccine candidates https://doi.org/10.1186/s12943-021-01310-0
####The DEGs and its location landscape, muatation, fraction altered was drawed in online kits 
###GEPIA, http://gepia2.cancer-pku.cn & http://www.cbioportal.org. We also dowloaded mutation and DEGs data.
#1. Now we search the mRNA candidates through intersecting over-expressed degs and mutations.
table_degenes <- read.delim("~/project2022/Figure1/table_degenes.txt")
overexp <- table_degenes[table_degenes$Log2.Fold.Change. > 0,]
Mutated_Genes <- read.delim("~/project2022/Figure1/Mutated_Genes.txt")
Mutated_Genes <- Mutated_Genes[!duplicated(Mutated_Genes$Gene),]
ovMut <- intersect(overexp$Gene.Symbol,Mutated_Genes$Gene)
#--------------------------------------------------------------------------------------------------------------------
#2.Searching the target genes through intersecting RFS & OS & ovMut. 
table_survival.OS <- read.delim("~/project2022/Figure1/table_survival-OS.txt")
table_survival.RFS <- read.delim("~/project2022/Figure1/table_survival-RFS.txt")
ovMutRFS <- intersect(ovMut,table_survival.RFS$Gene.Symbol)# 20 genes were maped and imported into GEPIA wepsite to 
# map the OS-related genes. 4 target genes were filtered. (CENPW poor prognosis, ZMYND10, RSPH1, TMEM229B)
ovMutOS <- intersect(ovMut,table_survival.OS$Gene.Symbol) #27 genes were maped and imported into GEPIA wepsite to 
# map the RFS-related genes. Finally, 6 target genes were filtered. (SLC7A5, CHPF, CCNE1 poor prognosis, SKAP1,CXCL9, SERPINA1)
# Those 47mRNAs were inputed into http://gepia2.cancer-pku.cn to search to the target genes.
#3. Relationships between immune cellular concentration & tumor puity score
#----------------------------------------------------------------------------------------------------
load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/TCGA/ICI/BRCA_estimate_socre.Rdata")
load("/Users/llls2012163.com/TCGA/TCGA_BRCA_FPKM/pcExpdata.Rdata")
load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/TCGA/ssgsva.Rdata")
exp <- pcExpdata$Exp
exp <- exp[c('CENPW','SLC7A5','CHPF','CCNE1'),]
exp <- t(exp)
exp <- as.data.frame(exp)
exp <- apply(exp,2,as.numeric)
exp <- log2(exp+1)
barcode <- colnames(pcExpdata$Exp)
rownames(exp) <- barcode
ssgsva <- t(ssgsva)
a <- ssgsva[,c(13,21,24)]
scores <- as.data.frame(scores)
scores$purity <- cos(0.6049872018+0.0001467884*scores$ESTIMATEScore)
preData <- cbind(exp,a,purity=scores$purity)
preData <- as.data.frame(preData)
setwd("/Users/llls2012163.com/project2022/Figure3")
library(ggExtra)
library(ggplot2)
library(ggprism)
library(ggpubr)
p <- ggplot(preData,aes(CCNE1,`purity`))+geom_point(size=0.8,color='red')+ theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),panel.border = element_rect(fill = NA,colour = 'black',size = 2))+
  geom_smooth(method = "lm")+
  stat_cor(data=preData,aes(CCNE1,`purity`),method = 'pearson')+
  theme(axis.text = element_text(size = 11))
p# the plot were drew and saved one by one
save(preData,file = 'preData.Rdata')
###-------------------------------------------------------------------------------------------------------------------------
## here we merged the rna seq data from TCGA and Metabric cohort and removed the batch effects through sva::combat function.
load("/Users/llls2012163.com/TCGA/TCGA_BRCA_FPKM/pcExpdata.Rdata")
load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/ICIScore validata/metabric/metabricRNAseq.Rdata")
exp <- pcExpdata$Exp[,113:ncol(pcExpdata$Exp)]
id <- rownames(exp)
exp <- apply(exp,2,as.numeric)
exp <- log2(exp+1)
exp <- as.data.frame(exp)
rownames(exp) <- id
id1 <- rownames(MetaBric_Corhort)
id2 <- intersect(id1,id)
Metacohort <- cbind(exp[id2,],MetaBric_Corhort[id2,])
library(sva)
batch <- paste0('batch',rep(c(1,2),c(1110,1904)))
s <- ComBat(dat=Metacohort, batch=batch, mod=NULL, par.prior=TRUE,  prior.plots=FALSE)
save(s,file = 'noneBatch.Rdata')
###pca plot to see batch effects
Metacohort <- na.omit(Metacohort)
pcab <- prcomp(t(Metacohort[1:100,]),center = T, scale. = T)
pcab.data <- as.data.frame(pcab$x)
platform <- c(rep('TCGA',1110),rep('METABRIC',1904))
pcab.data <- cbind(platform,pcab.data)
library(ggplot2)
library(ggprism)
library(patchwork)
p1 <- ggplot(pcab.data,aes(PC1,PC2,color=platform))+geom_point(size=0.3)+theme_prism()
s <- na.omit(s)
pcaa <- prcomp(t(s[1:100,]),center = T,scale. = T)
pcaa.data <- as.data.frame(pcaa$x)
pcaa.data <- cbind(platform,pcaa.data)
p2 <- ggplot(pcaa.data,aes(PC1,PC2,color=platform))+geom_point(size=0.3)+theme_prism()
p1+p2+plot_layout(ncol = 2,nrow = 2)
#--------------------------------------------------------------------------------------------------------------------------
#load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/TCGA/ICI/BRCA_estimate_socre.Rdata")########Now the old data were discarded
#load("/Users/llls2012163.com/TCGA/TCGA_BRCA_FPKM/pcExpdata.Rdata")
#load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/TCGA/ssgsva.Rdata")#########Now the old data were discarded. We used the Metacohort
# rna-seq data to estimate the estimatescore and immune cell abundance and the figure 2 were redrew.
load("~/project2022/noneBatch.Rdata")
load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/TCGA/ImmuneGenesetList.Rdata")# This immune gene List were download immune AI
library(GSVA)
ssgsva <- gsva(expr = s,gset.idx.list = z,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
save(ssgsva,file = 'MetaImmunecellabun.Rdata')
# estimate immune and stromal scores.
library(estimate)
estimate <- function(dat,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(dat,file = input.f,sep = '\t',quote = F)#####重点
  library(estimate)
  filterCommonGenes(input.f=input.f,
    output.f=output.f ,
    id="GeneSymbol")
  estimateScore(input.ds = output.f,
    output.ds=output.ds,
    platform="illumina") ## 注意这个platform参数
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='Metacohort'
scores=estimate(s,pro)
# Tumor purity=cos(0.6049872018+0.0001467884*ESTIMATEScore)
scores <- as.data.frame(scores)
scores$Purity <- cos(0.6049872018+0.0001467884*scores$ESTIMATEScore)
save(scores,file = 'Metaestimate.Rdata')# this data were saved and used in the Figure8.
#------------------------------------------------------------------------------------------------------
# immune cell abudance were obtained from my own research, which performed through ssgsva and the geneset was 
# assessed from ru bei bei et al. more details see doi: 10.3389/fonc.2022.844082.
load("~/project2022/immcellAbudance.Rdata")
clinical <- read.delim("~/TCGA/TCGA_BRCA/clinical.cart.2020-11-27/clinical.tsv")
clinical <- clinical[-c(1:24),]
TractClinicaldata <- function(x){
  m <- c("case_submitter_id","age_at_index","gender","vital_status","days_to_last_follow_up",
    "days_to_death","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_t","figo_stage")
  x <- x[,m]
  x <- x[!duplicated(x$case_submitter_id),]
  Day <- ifelse(x$days_to_death >0,x$days_to_death,x$days_to_last_follow_up)
  x <- cbind(x[,1:6],Day,x[,7:ncol(x)])
  x <- x[x$Day >0,]
  colnames(x) <-c('ID','Age','Gender','statu','Follow','days','Day','M','N',"T",'Stage')
  return(x)
}
clinical <- TractClinicaldata(clinical)
save(clinical,file = 'clinicTCGA.Rdata')
ssgsva <- ssgsva[,113:ncol(ssgsva)]
# 4.Now we performed consensus cluster perfile to identify the immune clusters
# 4.1immune gene List were downloaded from immune port: https://www.immport.org/shared/genelists
Cytokines <- read.delim2("~/Downloads/Cytokines.txt")
Interleukins <- read.delim2("~/Downloads/Interleukins.txt")
Interferons <- read.delim2("~/Downloads/Interferons.txt")
TNF <- read.delim2("~/Downloads/TNF.txt")
TGF <- read.delim2("~/Downloads/TGF.txt")
Chemokines <- read.delim2("~/Downloads/Chemokines.txt")
Cytokine_Receptors <-  read.delim2("~/Downloads/Cytokine_Receptors.txt")
Interleukins_Receptors <- read.delim2("~/Downloads/Interleukins_Receptors.txt")
Interferons_Receptors <- read.delim2("~/Downloads/Interferons_Receptors.txt")
TNF_Family_Members_Receptors <- read.delim2("~/Downloads/TNF_Family_Members_Receptors.txt")
TGF_b_Family_Members_Receptors <- read.delim2("~/Downloads/TGF_b_Family_Members_Receptors.txt")
Chemokine_Receptors <- read.delim2("~/Downloads/Chemokine_Receptors.txt")
TCR_Signaling_Pathway <- read.delim2("~/Downloads/TCR_Signaling_Pathway.txt")
BCR_Signaling_Pathway <- read.delim2("~/Downloads/BCR_Signaling_Pathway.txt")
Natural_Killer_Cell <-  read.delim2("~/Downloads/Natural_Killer_Cell.txt")
Antigen_Processing_and_Presentation <- read.delim2("~/Downloads/Antigen_Processing_and_Presentation.txt")
Antimicrobials <- read.delim2("~/Downloads/Antimicrobials.txt")
tarctIgene <- function(x){
  library(limma)
  a <- x[,2]
  b <- strsplit2(x[,3],',')
  c <- c(a,b)
  return(c)
}
Chemokines <- tarctIgene(Chemokines)
Interleukins <- tarctIgene(Interleukins)
Interferons <- tarctIgene(Interferons)
TNF<- tarctIgene(TNF)
TGF<- tarctIgene(TGF)
Cytokines<- tarctIgene(Cytokines)
Cytokine_Receptors<- tarctIgene(Cytokine_Receptors)
Interleukins_Receptors<- tarctIgene(Interleukins_Receptors)
TNF_Family_Members_Receptors<- tarctIgene(TNF_Family_Members_Receptors)
TGF_b_Family_Members_Receptors<- tarctIgene(TGF_b_Family_Members_Receptors)
Chemokine_Receptors <- tarctIgene(Chemokine_Receptors)
TCR_Signaling_Pathway <- tarctIgene(TCR_Signaling_Pathway)
BCR_Signaling_Pathway <- tarctIgene(BCR_Signaling_Pathway)
Natural_Killer_Cell<- tarctIgene(Natural_Killer_Cell)
Antigen_Processing_and_Presentation<- tarctIgene(Antigen_Processing_and_Presentation)
Antimicrobials<- tarctIgene(Antimicrobials)
Interferons_Receptors <- tarctIgene(Interferons_Receptors)
gene <- c(Antigen_Processing_and_Presentation,Antimicrobials,BCR_Signaling_Pathway,Chemokine_Receptors,
  Chemokines,Cytokine_Receptors,Cytokines,Interferons,Interferons_Receptors,Interleukins,
  Interleukins_Receptors,Natural_Killer_Cell,TCR_Signaling_Pathway,TGF,TGF_b_Family_Members_Receptors,
  TNF,TNF_Family_Members_Receptors)
gene <- gene[!duplicated(gene)]
save(gene,file = 'ImmuenGeneList.Rdata')
rm(list = ls())
# 4.2 we merge the TCGA & Metabric cohort, and remove the batch effects through sva::combat function
load("/Users/llls2012163.com/TCGA/TCGA_BRCA_FPKM/pcExpdata.Rdata")
load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/ICIScore validata/metabric/metabricRNAseq.Rdata")
exp <- pcExpdata$Exp[,113:ncol(pcExpdata$Exp)]
load("~/project2022/ImmuenGeneList.Rdata")
id <- rownames(exp)
exp <- apply(exp,2,as.numeric)
exp <- log2(exp+1)
exp <- as.data.frame(exp)
rownames(exp) <- id
id1 <- rownames(MetaBric_Corhort)
id2 <- intersect(id1,id)
Metacohort <- cbind(exp[id2,],MetaBric_Corhort[id2,])
library(sva)
batch <- paste0('batch',rep(c(1,2),c(1110,1904)))
s <- ComBat(dat=Metacohort, batch=batch, mod=NULL, par.prior=TRUE,  prior.plots=FALSE)
save(s,file = 'noneBatch.Rdata')
###pca plot to see batch effects
Metacohort <- na.omit(Metacohort)
pcab <- prcomp(t(Metacohort[1:100,]),center = T, scale. = T)
pcab.data <- as.data.frame(pcab$x)
platform <- c(rep('TCGA',1110),rep('METABRIC',1904))
pcab.data <- cbind(platform,pcab.data)
library(ggplot2)
library(ggprism)
library(patchwork)
p1 <- ggplot(pcab.data,aes(PC1,PC2,color=platform))+geom_point(size=0.3)+theme_prism()
s <- na.omit(s)
pcaa <- prcomp(t(s[1:100,]),center = T,scale. = T)
pcaa.data <- as.data.frame(pcaa$x)
pcaa.data <- cbind(platform,pcaa.data)
p2 <- ggplot(pcaa.data,aes(PC1,PC2,color=platform))+geom_point(size=0.3)+theme_prism()
p1+p2+plot_layout(ncol = 2,nrow = 2)
### here we perform consensus cluster
d <- s[gene,]
d <- na.omit(d)
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:1234],]
#sweep函数减去中位数进行标准化
d = sweep(d,1, apply(d,1,median,na.rm=T))
library(ConsensusClusterPlus)
d <- as.matrix(d)
p <- ConsensusClusterPlus(d=d,maxK = 8,reps =1000,
  pFeature = 1,clusterAlg = 'hc',title = 'Immunecluster',
  innerLinkage = 'average',finalLinkage = 'average',distance = 'pearson',
  ml=NULL,tmyPal = NULL,writeTable = TRUE,weightsItem = NULL,
  weightsFeature = NULL,verbose = F,corUse = 'everything',
  plot = 'png')
calcICL(res=p,title = 'Immunecluster',writeTable = FALSE)
rm(list = ls())
#4 Survival analysis
`Immunecluster.k=2.consensusClass` <- read.csv("~/project2022/Immunecluster/Immunecluster.k=2.consensusClass.csv", header=FALSE)
clu <- as.data.frame(`Immunecluster.k=2.consensusClass`)
load("~/project2022/clinicTCGA.Rdata")
load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/ICIScore validata/metabric/(MetaBric_CorhortClin.Rdata")
rownames(clinical) <- clinical$ID
colnames(clu) <- c('ID','IS')
library(limma)
id <- clu$ID
id <- strsplit2(id,'-')[,1:3]
id1 <- paste(id[,1],'-',id[,2],sep = '')
id <- paste(id1,'-',id[,3],sep = '')
clu$ID[1:1110] <- id[1:1110]
id3 <- intersect(clu$ID,clinical$ID)
id4 <- intersect(clu$ID,MetaBric_CorhortClin$Clinical$PATIENT_ID)
clu <- clu[!duplicated(clu$ID),]
rownames(clu) <- clu$ID
clu <- clu[c(id3,id4),]
clinical <- clinical[id3,]
clinical$Day <- as.numeric(clinical$Day)
clinical$statu <- ifelse(clinical$statu=='Alive',0,1)
s1 <- cbind(ID=clinical$ID,statu=clinical$statu,Day=clinical$Day)
Metc <- MetaBric_CorhortClin$Clinical[!duplicated(MetaBric_CorhortClin$Clinical),]
rownames(Metc) <- Metc$PATIENT_ID
Metc <- Metc[id4,]
Metc1 <- Metc
Metc1$VITAL_STATUS <- ifelse(Metc1$VITAL_STATUS=='Living',0,1)
Metc1$OS_MONTHS <- 30*as.numeric(Metc1$OS_MONTHS)
s2 <- cbind(ID=Metc1$PATIENT_ID,statu=Metc1$OS_STATUS,Day=Metc1$OS_MONTHS)
Metacohort <- rbind(s1,s2)
preLog <- cbind(Metacohort,clu)
preLog$IS <- ifelse(preLog$IS==1,'IS1',
  ifelse(preLog$IS==2,'IS2',
    ifelse(preLog$IS==3,'IS3','IS4')))
preLog$statu <- ifelse(preLog$statu=='Alive',0,1)
preLog$Day <- as.numeric(preLog$Day)
library(survminer)
library(survival)
library(ggprism)
fit <- survfit(Surv(Day, statu) ~ IS, data = preLog)
ggsurvplot(fit,
  pval = TRUE, conf.int = FALSE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_prism(), # Change ggplot2 theme
  palette = c("red", "blue","green")
)
fit1 <- survfit(Surv(Day, statu) ~ IS, data = preLog[1:1069,])
ggsurvplot(fit1,
  pval = TRUE, conf.int = FALSE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_prism(), # Change ggplot2 theme
  palette = c("red", "blue","green")
)
fit2 <- survfit(Surv(Day, statu) ~ IS, data = preLog[1070:nrow(preLog),])
ggsurvplot(fit2,
  pval = TRUE, conf.int = FALSE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_prism(), # Change ggplot2 theme
  palette = c("red", "blue","green")
)
save(preLog,file = 'preLogIS.Rdata')
#5 Relationships between IS and TNM stage.
load("~/project2022/clinicTCGA.Rdata")
load("~/project2022/preLogIS.Rdata")
brca_metabric_clinical_data <- read.delim("~/project2022/brca_metabric_clinical_data.tsv")
TractClinicaldata <- function(x){
  m <- c("case_submitter_id","age_at_index","gender","vital_status","days_to_last_follow_up",
    "days_to_death","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_t","ajcc_pathologic_stage")
  x <- x[,m]
  x <- x[!duplicated(x$case_submitter_id),]
  Day <- ifelse(x$days_to_death >0,x$days_to_death,x$days_to_last_follow_up)
  x <- cbind(x[,1:6],Day,x[,7:ncol(x)])
  x <- x[x$Day >0,]
  colnames(x) <-c('ID','Age','Gender','statu','Follow','days','Day','M','N',"T",'Stage')
  return(x)
}
clinical <- TractClinicaldata(clinical1)
save(clinical,file = 'clinicTCGA.Rdata')
load("~/project2022/clinicTCGA.Rdata")
TNMtcga <- cbind(ID=clinical$ID,'T'=clinical$'T',N=clinical$N,M=clinical$M,S=clinical$Stage)
TNMtcga <- as.data.frame(TNMtcga)
TNMtcga <- TNMtcga[!TNMtcga$'T' <0,]
TNMtcga <- TNMtcga[!TNMtcga$S <0,]
TNMtcga <- TNMtcga[!TNMtcga$N <0,]
TNMtcga <- TNMtcga[!TNMtcga$M <0,]
TNMtcga$'T' <- ifelse(
  TNMtcga$'T'=='T1' | TNMtcga$'T'=='T1a' | TNMtcga$'T'=='T1b' | TNMtcga$'T'=='T1c','T1',ifelse(
    TNMtcga$'T'=='T2' | TNMtcga$'T'=='T2a' | TNMtcga$'T'== 'T2b', 'T2',
    ifelse(TNMtcga$'T'== 'T3'| TNMtcga$'T'=='T3a', 'T3',
      ifelse(TNMtcga$'T'=='T4'|TNMtcga$'T'=='T4b' | TNMtcga$'T'=='T4d','T4','TX'))
  )
)
TNMtcga$N <- ifelse(TNMtcga$N=='N0' | TNMtcga$N=='N0 (i-)' | TNMtcga$N=='N0 (i+)' | TNMtcga$N=='N0 (mol+)','N0',
  ifelse(TNMtcga$N=='N1' | TNMtcga$N=='N1a' | TNMtcga$N=='N1b' |TNMtcga$N=='N1c' | TNMtcga$N=='N1mi','N1',
    ifelse(TNMtcga$N=='N2'| TNMtcga$N=='N2a','N2',ifelse(
      TNMtcga$N=='N3' | TNMtcga$N=='N3a' | TNMtcga$N=='N2b' | TNMtcga$N=='N3c','N3','NX'
    ))))
TNMtcga$M <- ifelse(TNMtcga$M=='M0' | TNMtcga$M=='cM0 (i+)','MO',ifelse(
  TNMtcga$M=='M1','M1','MX'
))
TNMtcga$S <- ifelse(TNMtcga$S=='Stage I'| TNMtcga$S== 'Stage IA' | TNMtcga$S== 'Stage IB','Stage I',
  ifelse(TNMtcga$S=='Stage II' | TNMtcga$S== 'Stage IIA' | TNMtcga$S== 'Stage IIB','Stage II',
    ifelse(TNMtcga$S== 'Stage III' | TNMtcga$S== 'Stage IIIA' | TNMtcga$S== 'Stage IIIB' | TNMtcga$S== 'Stage IIIC','Stage III',
      ifelse(TNMtcga$S== 'Stage IV','Stage IV','Stage X'))))
TNMtcga <- TNMtcga[!TNMtcga$'T'=='TX',]
TNMtcga <- TNMtcga[!TNMtcga$N=='NX',]
TNMtcga <- TNMtcga[!TNMtcga$M=='MX',]
TNMtcga <- TNMtcga[!TNMtcga$S=='Stage X',]
a1 <- preLog[,4:5]
rownames(a1) <- a1$ID
a1 <- a1[TNMtcga$ID,]
s <- cbind(a1,TNMtcga)
s[c(2,4:7)] <- apply(s[c(2,4:7)],2,as.factor)
library(tidyverse)
library(reshape2)
library(ggplot2)
s1 <- s[,-3]
s1$IS <- ifelse(s1$IS=='IS1',1,2)
#1$T <- ifelse(s1$'T'=='T1',1,ifelse(
  s1$T=='T2',2,ifelse(s1$T=='T3',3,4)
))
#s1$N <- ifelse(s1$N=='N0',1,ifelse(s1$N=='N1',2,ifelse(
  s1$N=='N2',3,4
)))
#s1$M <-  ifelse(s1$M=='MO',0,1)
#s1$S <- ifelse(s1$S=='Stage I',1,ifelse(
  s1$S=='Stage II',2,ifelse(
    s1$S=='Stage III',3,4
  )
))
save(s1,file = 'fig4corplot.Rdata')
load("~/project2022/fig4corplot.Rdata")
s2 <- s1[,c(2,4)]
a <- s1[s1$M=='MO',]$IS
b <- s1[s1$M=='M1',]$IS
t.test(a,b)
wilcox.test(a,b)
wilcox.test(s1[s1$T=='T3',]$IS,s1[s1$T=='T2',]$IS)
## finally the relationships between TNMS stage and IS were detected by wilicox.test after shapiro.test
## to detected the normal distribution one by one.
M <- cbind(M0=c(1,0.7516),M1=c(0.7516,1))
rownames(M) <- c('M0','M1')
N <- cbind(N0=c(1,0.006991,0.4857,0.545),N1=c(0.006991,1,0.2599,0.4465),N2=c(0.4857,0.2599,1,0.9384),
  N3=c(0.545,0.4465,0.9384,1))
rownames(N) <- c('N0','N1','N2','N3')
Ts <- cbind(T1=c(1,0.0132,0.8601,0.5535),T2=c(0.01232,1,0.05104,0.1091),
  T3=c(0.8601,0.05104,1,0.6511),T4=c(0.5535,0.1091,0.6511,1))
rownames(Ts) <- c('T1','T2','T3','T4')
ss <- cbind('Stage I'=c(1,0.05284,0.6889,0.4435),'Stage II'=c(0.05284,1,0.1099,0.9461),
  'Stage III'=c(0.6889,0.1099,1,0.549),'Stage IV'=c(0.4435,0.9461,0.549,1))
rownames(ss) <- c('Stage I','Stage II','Stage III','Stage IV')
library(pheatmap)
pheatmap(ss, 
  show_colnames = TRUE,   # 是否显示列名
  show_rownames=TRUE,     # 是否显示行名
  fontsize=10,             # 字体大小
  color = colorRampPalette(c('#87CEFA','#ffffff','#4682B4'))(100), # 指定热图的颜色
  annotation_legend=TRUE, # 是否显示图例
  border_color=NA,        # 边框颜色 NA表示没有
  scale="none",           # 指定归一化的方式。"row"按行归一化，"column"按列归一化，"none"不处理
  cluster_rows = FALSE,    # 是否对行聚类
  cluster_cols = FALSE     # 是否对列聚类
)
Fig4cor <- list(s1,Ts,N,M,ss)
save(Fig4cor,file = 'Fig4Cor.Rdata')
# stack plot for percentage of IS to TNMS stage.
q <- Fig4cor[[1]]
q$IS <- ifelse(q$IS==1,'IS1','IS2')
Ms <- cbind(M=c('M0','M0','M1','M1'),IS=c('IS1','IS2','IS1','IS2'),Freq=c(580,295,10,6),label=c("66%", "34%", "63%", "37%"),
  precent = c(0.66,0.34,0.63,0.37))
Ms <- as.data.frame(Ms)
library(ggplot2)
p <- ggplot(Ms,aes(M,precent,fill=IS))+geom_bar(stat = 'identity',position = position_stack())+
  scale_fill_manual(values = c('#87CEFA','#4682B4'))+theme_prism(base_size = 12)+theme(legend.position = 'top')+
  labs(y='Precentage',x='')+geom_text(aes(label=label),vjust=2,size=6,color='black')
p

Ts <- cbind('T'=c('T1','T1','T2','T2','T3','T3','T4','T4'),IS =c(rep(c('IS1','IS2'),4)),Frq=c(168,47,328,200,71,27,23,7),
  precent=c(0.71,0.29,0.62,0.38,0.72,0.28,0.77,0.23),label=c('71%','29%','62%','38%','72%','28%','77%','23%'))
Ts <- as.data.frame(Ts)
p1 <- ggplot(Ts,aes(T,precent,fill=IS))+geom_bar(stat = 'identity',position = position_stack())+
  scale_fill_manual(values = c('#87CEFA','#4682B4'))+theme_prism(base_size = 12)+theme(legend.position = 'top')+
  labs(y='Precentage',x='')+geom_text(aes(label=label),vjust=2,size=6,color='black')
p1
Ns <- cbind(N=c('N0','N0','N1','N1','N2','N2','N3','N3'),IS =c(rep(c('IS1','IS2'),4)),Frq=c(273,165,215,84,68,35,34,17),
  precent=c(0.62,0.38,0.72,0.28,0.67,0.33,0.72,0.28),label=c('62%','38%','72%','28%','67%','33%','72%','72%'))
Ns <- as.data.frame(Ns)
p2 <- ggplot(Ns,aes(N,precent,fill=IS))+geom_bar(stat = 'identity',position = position_stack())+
  scale_fill_manual(values = c('#87CEFA','#4682B4'))+theme_prism(base_size = 12)+theme(legend.position = 'top')+
  labs(y='Precentage',x='')+geom_text(aes(label=label),vjust=2,size=6,color='black')
p2

Ss <- cbind(S=c('Stage I','Stage I','Stage II','Stage II','Stage III','Stage III','Stage IV ','Stage IV '),IS =c(rep(c('IS1','IS2'),4)),
  Frq=c(114,45,330,191,136,59,10,6),precent=c(0.72,0.28,0.63,0.37,0.70,0.30,0.63,0.37),
  label=c('72%','28%','63%','37%','70%','30%','63%','37%'))
Ss <- as.data.frame(Ss)
p3 <- ggplot(Ss,aes(S,precent,fill=IS))+geom_bar(stat = 'identity',position = position_stack())+
  scale_fill_manual(values = c('#87CEFA','#4682B4'))+theme_prism(base_size = 12)+theme(legend.position = 'top')+
  labs(y='Precentage',x='')+geom_text(aes(label=label),vjust=2,size=6,color='black')
p3
Fig4cor$Ms <- Ms
Fig4cor$Ns <- Ns
Fig4cor$Ts <- Ts
Fig4cor$Ss <- Ss
save(Fig4cor,file = 'Fig4cor.Rdata')
rm(list = ls())
# TMB, Mutation number, and mutation landscape between IS1 and IS2.
load("/Users/llls2012163.com/project2022/Figure4/preLogIS.Rdata")
library(maftools)
id <- '/Users/llls2012163.com/TCGA/TCGA_BRCA/gdc_download_20210318_073529.382183/6c93f518-1956-4435-9806-37185266d248/TCGA_BRCA.maf.gz'
tmb1 <- read.maf(maf = id)
TCGAtmb <- tmb(tmb1)[,1:3]
Metabrictmb <- brca_metabric_clinical_data[,c(3,23,36)]
colnames(Metabrictmb) <- c('ID','count','TMB')
colnames(TCGAtmb) <- c('ID','count','TMB')
library(limma)
s <- TCGAtmb$ID
s <- strsplit2(s,'-')
s <- s[,1:3]
s1 <- paste(s[,1],'-',s[,2],sep = '')
s <- paste(s1,'-',s[,3],sep = '')
TCGAtmb$ID <- s
MetaTMB <- rbind(TCGAtmb,Metabrictmb)
MetaTMB <- MetaTMB[!duplicated(MetaTMB$ID),]
MetaTMB <- as.data.frame(MetaTMB)
rownames(MetaTMB) <- MetaTMB$ID
id1 <- rownames(preLog)
MetaTMB <- as.data.frame(MetaTMB)
MetaTMB <- MetaTMB[id1,]
preTMB <- cbind(preLog,MetaTMB)
preTMB <- na.omit(preTMB)
save(preTMB,file = 'preTMB.Rdata')
write.csv(preTMB,file = 'preTMB.csv')
preTMB <- preTMB[,-c(4,6)]
preTMB <- as.data.frame(preTMB)
shapiro.test(preTMB[preTMB$IS=='IS1',]$count)
shapiro.test(preTMB[preTMB$count,]$count)
shapiro.test(preTMB[preTMB$IS=='IS1',]$TMB)
shapiro.test(preTMB[preTMB$count,]$TMB)
wilcox.test(preTMB[preTMB$IS=='IS1',]$count,preTMB[preTMB$count,]$count)
wilcox.test(preTMB[preTMB$IS=='IS1',]$TMB,preTMB[preTMB$count,]$TMB)
rm(list = ls())
# The result visulization were performed by Prism 8 software 
# and the data analyszes were performed by R software.
# It seems that there were no mutation data in the Metabric cohort, so we showed the 
# top ten mutated genes landscape between IS1 and IS2 group of TCGA cohort.
library(maftools)
id <- '/Users/llls2012163.com/TCGA/TCGA_BRCA/gdc_download_20210318_073529.382183/6c93f518-1956-4435-9806-37185266d248/TCGA_BRCA.maf.gz'
load("/Users/llls2012163.com/project2022/Figure4/preLogIS.Rdata")
tmb <- read.maf(maf = id)
s <- preLog[,c(1,5)]
library(limma)
s1 <- as.data.frame(tmb@clinical.data)
s2 <- strsplit2(s1[,1],'-')[,1:3]
s3 <- paste(s2[,1],s2[,2],sep = '-')
s4 <- paste(s3,s2[,3],sep = '-')
rownames(s1) <- s4
s <- s[s4,]
clinical <- cbind(s1,s)
clinical <- as.data.frame(clinical)
rownames(clinical) <- 1:961
tmb@clinical.data$ID <- s4 
tmb@clinical.data$IS <- clinical$IS
tmb@clinical.data <- na.omit(tmb@clinical.data)
annocor <- list(IS=c(IS1='#87CEFA',IS2='#4682B4'))
#此处使用RColorBrewer的颜色，当然也可以使用任意颜色
vc_cols = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
oncoplot(maf = tmb,top = 10,clinicalFeatures = c('IS'),
  sortByAnnotation = TRUE,annotationColor = annocor,showTitle = TRUE,
  legend_height = 3,anno_height = 2,legendFontSize = 0.8,fill = TRUE,
  drawRowBar = TRUE,drawBox = FALSE,colbar_pathway = TRUE,colors = vc_cols)
# Relationships between IS and CYT, the geometrical mean of PRF1 and GZMA mRNA expression level
# Shapiro-Wilk normality test were used to identify the normality and the wilcoxon test was used to 
# compare the data.
load("~/project2022/noneBatch.Rdata")
load("/Users/llls2012163.com/project2022/Figure4/preLogIS.Rdata")
s <- as.data.frame(s)
Target <- s[c('PRF1','GZMA'),]
Target <- t(Target)
Target <- as.data.frame(Target)
Target$CYT <- Target$PRF1*Target$GZMA 
Target$CYT <- sqrt(Target$CYT)  
id <- rownames(Target)[1:1110]
TreatTCGAid <- function(x){
  library(limma)
  a <- strsplit2(x,'-')[,1:3]
  b <- paste(a[,1],'-',a[,2],sep ='')
  c <- paste(b,'-',a[,3],sep = '')
  return(c)
}
id1 <- TreatTCGAid(id)
id2 <- c(id1,rownames(Target)[1111:nrow(Target)])
Target$ID <- id2
Target <- Target[!duplicated(Target$ID),]
rownames(Target) <- Target$ID
Target <- Target[rownames(preLog),]
CYT <- cbind(ID=Target$ID,PRF1=Target$PRF1,GZMA=Target$GZMA,CYT=Target$CYT,IS=preLog$IS)
CYT <- as.data.frame(CYT)
write.csv(CYT,file = 'CYTis.csv')
CYT[,2:4] <- apply(CYT[,2:4],2,as.numeric)
shapiro.test(CYT[CYT$IS=='IS1',]$CYT)
shapiro.test(CYT[CYT$IS=='IS2',]$CYT)
wilcox.test(CYT[CYT$IS=='IS1',]$CYT,CYT[CYT$IS=='IS2',]$CYT)
rm(list = ls())
####Here we performed consensus cluster in TCGA and Metabric cohort respectively.
load("~/project2022/noneBatch.Rdata")
load("~/project2022/ImmuenGeneList.Rdata")
exp <- as.data.frame(s)[gene,]
exp <- exp[!duplicated(exp),]
exp <- na.omit(exp)
expTCGA <- exp[,1:1110]
expMetabric <- exp[,1111:ncol(exp)]
d=expTCGA
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:1222],]
#sweep函数减去中位数进行标准化
d = sweep(d,1, apply(d,1,median,na.rm=T))
library(ConsensusClusterPlus)
d <- as.matrix(d)
p <- ConsensusClusterPlus(d=d,maxK = 8,reps =1000,
  pFeature = 1,clusterAlg = 'hc',title = 'ImmuneclusterTCGA',
  innerLinkage = 'average',finalLinkage = 'average',distance = 'pearson',
  ml=NULL,tmyPal = NULL,writeTable = TRUE,weightsItem = NULL,
  weightsFeature = NULL,verbose = F,corUse = 'everything',
  plot = 'png')
calcICL(res=p,title = 'Immunecluster',writeTable = FALSE)
###LogRank test
load("~/project2022/Figure4/preLogIS.Rdata")
`ImmuneclusterTCGA.k=2.consensusClass` <- read.csv("~/project2022/ImmuneclusterTCGA/ImmuneclusterTCGA.k=2.consensusClass.csv", header=FALSE)
TreatTCGAid <- function(x){
  library(limma)
  a <- strsplit2(x,'-')[,1:3]
  b <- paste(a[,1],'-',a[,2],sep ='')
  c <- paste(b,'-',a[,3],sep = '')
  return(c)
}
a <- as.data.frame(`ImmuneclusterTCGA.k=2.consensusClass`)
id <- TreatTCGAid(a$V1)
a$V1 <- id
a <- a[!duplicated(a$V1),]
rownames(a) <- a$V1
id2 <- intersect(a$V1,preLog$ID)
a <- a[id2,]
preLog <- preLog[id2,]
colnames(a) <- c('id','IS1')
preLog1 <- cbind(preLog,a)
preLog1$IS1 <- ifelse(preLog1$IS1==1,'IS1',ifelse(preLog1$IS1==2,'IS2','IS3'))
library(survminer)
library(survival)
library(ggprism)
fit <- survfit(Surv(Day, statu) ~ IS1, data = preLog1)
ggsurvplot(fit,
  pval = TRUE, conf.int = FALSE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_prism(), # Change ggplot2 theme
  palette = c("red", "blue","green")
) # The logRank test showed no significance.
#### Metabric cohort 
d=expMetabric
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:1222],]
#sweep函数减去中位数进行标准化
d = sweep(d,1, apply(d,1,median,na.rm=T))
library(ConsensusClusterPlus)
d <- as.matrix(d)
p <- ConsensusClusterPlus(d=d,maxK = 8,reps =1000,
  pFeature = 1,clusterAlg = 'hc',title = 'ImmuneclusterMetabric',
  innerLinkage = 'average',finalLinkage = 'average',distance = 'pearson',
  ml=NULL,tmyPal = NULL,writeTable = TRUE,weightsItem = NULL,
  weightsFeature = NULL,verbose = F,corUse = 'everything',
  plot = 'png')
calcICL(res=p,title = 'Immunecluster',writeTable = FALSE)
####LogRank test
`ImmuneclusterMetabric.k=2.consensusClass` <- read.csv("~/project2022/ImmuneclusterMetabric/ImmuneclusterMetabric.k=2.consensusClass.csv", header=FALSE)
a <- as.data.frame(`ImmuneclusterMetabric.k=2.consensusClass`)
rownames(a) <- a$V1
colnames(a) <- c('ID','IS2')
preLog1 <- preLog[a$ID,]
a$IS2 <- ifelse(a$IS2==1,'IS1','IS2')
preLog1 <- cbind(preLog1,a)
library(survminer)
library(survival)
library(ggprism)
fit <- survfit(Surv(Day, statu) ~ IS2, data = preLog1)
ggsurvplot(fit,
  pval = TRUE, conf.int = FALSE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_prism(), # Change ggplot2 theme
  palette = c("red", "blue"))
rm(list = ls())
#########p=0.000015  
## Here, we identify the relationships between RFS and IS in the metacohort.
load("~/project2022/Figure4/preLogIS.Rdata")
brca_metabric_clinical_data <- read.delim("~/project2022/brca_metabric_clinical_data.tsv")
s <- brca_metabric_clinical_data[,c(2,30:31)]
colnames(s) <- c('ID','RFS','Statu')
s$Statu <- strsplit2(s$Statu,':')[,1]
s$Statu <- as.numeric(s$Statu)
rownames(s) <- s$ID
id <- intersect(s$ID,preLog$ID)
s<- s[id,]
preLog <- preLog[id,]
preLog1 <- cbind(preLog,s)
library(survminer)
library(survival)
library(ggprism)
fit <- survfit(Surv(RFS, Statu) ~ IS, data = preLog1)
ggsurvplot(fit,
  pval = TRUE, conf.int = FALSE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_prism(), # Change ggplot2 theme
  palette = c("red", "blue"))
rm(list = ls())
### Now we compare the ICD and ICP
load("~/project2022/noneBatch.Rdata")
`Immunecluster.k=2.consensusClass` <- read.csv("~/project2022/Immunecluster/Immunecluster.k=2.consensusClass.csv", header=FALSE)
ICD <- c('ANXA1','CALR','CXCL10','EIF2A','EIF2AK1','EIF2AK2','EIF2AK3','EIF2AK4',
  'FPR1','HGF','HMGB1','IFNAR1','IFNAR2','IFNE','LRP1','MET','P2RX7','P2RY2',
  'PANX1','TLR3','TLR4')
ICP <- c('ADORA2A','BTLA','CD160','CD200','CD200R1','CD244','CD27','CD274','CD276',
  'CD28','CD40','CD40LG','CD44','CD48','CD70','CD80','CD86','CTLA4','HAVCR2','HHLA2','ICOS',
  'ICOSLG','IDO1','LAG3','LAIR1','LGALS9','NRP1','PDCD1','PDCD1LG2','TIGIT','TMIGD2',
  'TNFRSF14','TNFRSF18','TNFRSF25','TNFRSF4','TNFRSF8','TNFRSF9','TNFSF14','TNFSF15','TNFSF18',
  'TNFSF4','TNFSF9','VTCN1')
s <- as.data.frame(s)
ICDexp <- s[ICD,]
ICPexp <- s[ICP,]
ICDexp <- na.omit(ICDexp)
ICPexp <- na.omit(ICPexp)
ICDexp <- t(ICDexp)
ICPexp <- t(ICPexp)
Cluster <- as.data.frame(`Immunecluster.k=2.consensusClass`)
colnames(Cluster) <- c('ID','IS')
Cluster$IS <- ifelse(Cluster$IS==1,'IS1','IS2')
ICDexp <- cbind(ICDexp,IS=Cluster$IS)
ICPexp <- cbind(ICPexp,IS=Cluster$IS)
ICDexptcga <- as.data.frame(ICDexp[1:1110,])
ICDexptmetabric <- as.data.frame(ICDexp[1111:3014,])
ICDexptcga1<-tidyr::gather(
          data=ICDexptcga,
          key="Gene",
          value="value",
          ANXA1:TLR3
  )##################重点代码
ICDexptcga1$value <- as.numeric(ICDexptcga1$value)
library(ggplot2)
library(ggprism)
p <- ggplot(data = ICDexptcga1,mapping = aes(x = Gene,y= value,fill=IS))+geom_boxplot()+
     theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 30),
       legend.position = 'top')+scale_fill_manual(values = c('#87CEFA','#4682B4'))+
     labs(x='',y='Expression value')
p
ICDexpmeta<-tidyr::gather(
  data=ICDexptmetabric,
  key="Gene",
  value="value",
  ANXA1:TLR3
)
ICDexpmeta$value <- as.numeric(ICDexpmeta$value)
p1 <- ggplot(data = ICDexpmeta,mapping = aes(x = Gene,y= value,fill=IS))+geom_boxplot()+
  theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 30),
    legend.position = 'top')+scale_fill_manual(values = c('#87CEFA','#4682B4'))+
  labs(x='',y='Expression value')
p1
ICPTCGA <- ICPexp[1:1110,]
ICPmeta <- ICPexp[1111:3014,]
ICpexpmeta <- as.data.frame(ICPmeta)
ICpexpmeta<-tidyr::gather(
  data=ICpexpmeta,
  key="Gene",
  value="value",
  ADORA2A:VTCN1)
ICpexpmeta$value <- as.numeric(ICpexpmeta$value)
p3 <- ggplot(data = ICpexpmeta,mapping = aes(x = Gene,y= value,fill=IS))+geom_boxplot()+
  theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 30),
    legend.position = 'top')+scale_fill_manual(values = c('#87CEFA','#4682B4'))+
  labs(x='',y='Expression value')
p3
ICPTCGA <- as.data.frame(ICPTCGA)
ICpexpTCGA<-tidyr::gather(
  data=ICPTCGA,
  key="Gene",
  value="value",
  ADORA2A:VTCN1)
ICpexpTCGA$value <- as.numeric(ICpexpTCGA$value)
p4 <- ggplot(data = ICpexpTCGA,mapping = aes(x = Gene,y= value,fill=IS))+geom_boxplot()+
  theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 30),
    legend.position = 'top')+scale_fill_manual(values = c('#87CEFA','#4682B4'))+
  labs(x='',y='Expression value')
p4
ICDICP <- list(ICDtcga=ICDexptcga,ICDmetabric=ICDexptmetabric,ICPtcga=ICPTCGA,ICPmetabric=ICPmeta)
save(ICDICP,file = 'ICDICP.Rdata')
rm(list = ls())
####compare the data of ICD and ICP in IS1 and IS2 groups.
load("~/project2022/ICDICP.Rdata")
ICDtcga <- as.data.frame(ICDICP$ICDtcga)
ICDtcga1 <-cbind(as.data.frame(apply(ICDtcga[,1:19],2,as.numeric)),IS=ICDtcga$IS)
ICDtcga1 <- as.data.frame(ICDtcga1)
rownames(ICDtcga1) <- rownames(ICDtcga)
output <- vector()
gene <- colnames(ICDtcga1)[1:19]
gene <- as.vector(gene)
for(i in gene){
  a1=ICDtcga1[ICDtcga1$IS=='IS1',i]
  a2=ICDtcga1[ICDtcga1$IS=='IS2',i]
  a3=wilcox.test(a1,a2)$p.value
  a4 <- as.vector(a3)
  if(length(output)==0){
    output=a3
  }else{
    output=c(output,a3)
  }
}
ICDtcgaP = cbind(ICD=gene,p=output)
ICDmeta <- ICDICP$ICDmetabric
ICDmeta1 <-cbind(as.data.frame(apply(ICDmeta[,1:19],2,as.numeric)),IS=ICDmeta$IS)
ICDtmeta1 <- as.data.frame(ICDmeta1)
output <- vector()
for(i in gene){
  a1=ICDmeta1[ICDmeta1$IS=='IS1',i]
  a2=ICDmeta1[ICDmeta1$IS=='IS2',i]
  a3=wilcox.test(a1,a2)$p.value
  a4 <- as.vector(a3)
  if(length(output)==0){
    output=a3
  }else{
    output=c(output,a3)
  }
}
ICDmetaP = cbind(ICD=gene,p=output)
ICPTCGA <- ICDICP$ICPtcga
gene <- colnames(ICPTCGA)[1:40]
ICPTCGA1 <-cbind(as.data.frame(apply(ICPTCGA[,1:40],2,as.numeric)),IS=ICPTCGA$IS)
ICPTCGA1<- as.data.frame(ICPTCGA1)
output <- vector()
for(i in gene){
  a1=ICPTCGA1[ICPTCGA1$IS=='IS1',i]
  a2=ICPTCGA1[ICPTCGA1$IS=='IS2',i]
  a3=wilcox.test(a1,a2)$p.value
  a4 <- as.vector(a3)
  if(length(output)==0){
    output=a3
  }else{
    output=c(output,a3)
  }
}
ICPTCGApvalue <- cbind(ICP=gene,p=output)
ICPmeta <- ICDICP$ICPmetabric
ICPmeta <- as.data.frame(ICPmeta)
ICPmeta1 <-cbind(as.data.frame(apply(ICPmeta[,1:40],2,as.numeric)),IS=ICPmeta$IS)
output <- vector()
for(i in gene){
  a1=ICPmeta1[ICPmeta1$IS=='IS1',i]
  a2=ICPmeta1[ICPmeta1$IS=='IS2',i]
  a3=wilcox.test(a1,a2)$p.value
  a4 <- as.vector(a3)
  if(length(output)==0){
    output=a3
  }else{
    output=c(output,a3)
  }
}
ICPmetapvalue <- cbind(ICP=gene,p=output)
ICDICP$ICDtcgaPvalue <- ICDtcgaP
ICDICP$ICDmetabricpvalue <- ICDmetaP
ICDICP$ICPtcgaPvalue <- ICPTCGApvalue
ICDICP$ICPmetabricPvalue <- ICPmetapvalue
ICDICP$anno <- c('all of the data analysis were identified by wilcox test')
save(ICDICP,file = 'Figure6.Rdata')
rm(list = ls())
## Here we visulized the immune cell abundance of Metabric and TCGA cohort among IS1 and IS2 groups by heatmap.
load("/Users/llls2012163.com/project2022/MetaImmunecellabun.Rdata")
`Immunecluster.k=2.consensusClass` <- read.csv("~/project2022/Immunecluster/Immunecluster.k=2.consensusClass.csv", header=FALSE)
s <- as.data.frame(`Immunecluster.k=2.consensusClass`)
colnames(s) <- c('ID','IS')
s$IS <- ifelse(s$IS==1,'IS1','IS2')
library(pheatmap)
s$cohort <- c(rep('TCGA',1110),rep('Metabric',1904))
rownames(s) <- s$ID
annocor <- list(IS=c(IS1='#87CEFA',IS2='#4682B4'),cohort=c(TCGA='red',Metabric='green'))
s <- s[,-1]
s <- s[order(s$cohort,s$IS),]
id <- rownames(s)
ssgsva <- ssgsva[,id]
bk <- c(seq(-3,0,by=0.01),seq(0.01,3,by=0.01))
pheatmap(ssgsva,scale='row',show_colnames = FALSE,show_rownames = TRUE,cluster_rows = TRUE,cluster_cols = FALSE,
  annotation_col = s,annotation_colors = annocor,
  color=c(colorRampPalette(colors=c("blue","white"))(length(bk)/2),colorRampPalette(color=c("white","red"))(length(bk)/2)),
  legend_breaks=seq(-3,3,1),breaks = bk,fontsize_row=10,gaps_col = c(1904,1219,2642))
## We also visulized the immune cell abundance of Metabric and TCGA cohort among IS1 and IS2 groups by boxplot.
ssgsva <- as.data.frame(ssgsva)
ssgsva <- t(ssgsva)
q <- s[rownames(ssgsva),]
ssgsva1 <- cbind(ssgsva,q)
ssgsvaTCGA <- ssgsva1[ssgsva1$cohort=='TCGA',][,-30]
ssgsvaMeta <- ssgsva1[ssgsva1$cohort=='Metabric',][,-30]
ssgsvaTCGA<-tidyr::gather(
  data=ssgsvaTCGA,
  key="cell",
  value="abundance",
  "Activated CD8 T cell":"Neutrophil"
)##################重点代码
library(ggplot2)
library(ggprism)
p <- ggplot(data =ssgsvaTCGA ,mapping = aes(x = cell,y= abundance,fill=IS))+geom_boxplot()+
  theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 45),
    legend.position = 'top')+scale_fill_manual(values = c('#87CEFA','#4682B4'))+
  labs(x='',y='density')
p
ssgsvaMeta<-tidyr::gather(
  data=ssgsvaMeta,
  key="cell",
  value="abundance",
  "Activated CD8 T cell":"Neutrophil"
)
p2 <- ggplot(data =ssgsvaMeta ,mapping = aes(x = cell,y= abundance,fill=IS))+geom_boxplot()+
  theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 45),
    legend.position = 'top')+scale_fill_manual(values = c('#87CEFA','#4682B4'))+
  labs(x='',y='density')
p2
ssgsvaTCGA1 <- ssgsva1[ssgsva1$cohort=='TCGA',][,-30]
ssgsvaMeta1 <- ssgsva1[ssgsva1$cohort=='Metabric',][,-30]
output <- vector()
cell <- colnames(ssgsvaTCGA1)[-29]
for(i in cell){
  a = ssgsvaTCGA1[ssgsvaTCGA1$IS=='IS1',i]
  b = ssgsvaTCGA1[ssgsvaTCGA1$IS=='IS2',i]
  c = wilcox.test(a,b)$p.value
  if(length(output)==0){
    output=c
  }else{
    output=c(output,c)
  }
}
TCGAssgsvaPvalue <- cbind(cell=cell,pvalue=output)
for(i in cell){
  a = ssgsvaMeta1[ssgsvaMeta1$IS=='IS1',i]
  b = ssgsvaMeta1[ssgsvaMeta1$IS=='IS2',i]
  c = wilcox.test(a,b)$p.value
  if(length(output)==0){
    output=c
  }else{
    output=c(output,c)
  }
}
MetassgsvaPvalue <- cbind(cell=cell,pvalue=output)
Fig7 <- list(ssgsvaTCGA=ssgsvaTCGA1,ssgsvaMeta=ssgsvaMeta1,Tcgapvalue=TCGAssgsvaPvalue,MetassgsvaPvalue=MetassgsvaPvalue)
save(Fig7,file = 'Fig7.Rdata')
######DDR tree
library(monocle)
load("~/project2022/noneBatch.Rdata")
load("~/project2022/ImmuenGeneList.Rdata")
gene <- intersect(gene,rownames(s))
`Immunecluster.k=2.consensusClass` <- read.csv("~/project2022/Immunecluster/Immunecluster.k=2.consensusClass.csv", header=FALSE)
pData <- as.data.frame(`Immunecluster.k=2.consensusClass`)
colnames(pData) <- c('cell','IS')
rownames(pData) <- pData$cell
pData <- as.data.frame(pData)
pData$IS <- ifelse(pData$IS==1,'IS1','IS2')
fDat <- cbind(ensembl_gene_id=rownames(s),gene_short_name=rownames(s))
rownames(fDat) <- fDat[,1]
fDat <- as.data.frame(fDat)
fDat <- fDat[gene,]
s <- s[gene,]
cds <- newCellDataSet(s,phenoData = Biobase::AnnotatedDataFrame(pData),featureData = Biobase::AnnotatedDataFrame(fDat))
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- reduceDimension(cds,max_components = 2,method='DDRTree')
cds <- orderCells(cds)
plot_cell_trajectory(cds = cds,color_by = "State")
plot_cell_trajectory(cds = cds,color_by = "State")+# state 1 6
  scale_color_manual(values = c('#FF6600','#778899','#778899','#778899','#778899','#33CCCC','#778899','#778899','#778899','#778899','#778899'))
a <- cds@phenoData@data
a1 <- a[a$State==1 | a$State==6| a$State==8,]
TreatTCGAid <- function(x){
  library(limma)
  a <- strsplit2(x,'-')[,1:3]
  b <- paste(a[,1],'-',a[,2],sep ='')
  c <- paste(b,'-',a[,3],sep = '')
  return(c)
}
a1$q <- 1:1372
id <- a1$cell[1:524]
id <- TreatTCGAid(id)
id <- c(id,a1$cell[525:1372])
a1$ID <- id
a1 <- a1[!duplicated(a1$ID),]
rownames(a1) <- a1$ID
load("~/project2022/Figure4/preLogIS.Rdata")
preLog1 <- preLog[a1$ID,]
preLog1 <- na.omit(preLog1)
a1 <- a1[rownames(preLog1),]
preLog1 <- cbind(ID=rownames(preLog1),statu=preLog1$statu,Day=preLog1$Day,IS=preLog1$IS,State=a1$State)
preLog1 <- as.data.frame(preLog1)
preLog1$Day <- as.numeric(preLog1$Day)
preLog1$statu <- as.numeric(preLog1$statu)
library(survminer)
library(survival)
library(ggprism)
fit <- survfit(Surv(Day, statu) ~ State, data = preLog1)
ggsurvplot(fit,
  pval = TRUE, conf.int = FALSE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_prism(), # Change ggplot2 theme
  palette = c("red", "blue",'green'))
Fig8 <- list(cds=cds,log168=preLog1)
save(Fig8,file = 'Fig8.Rdata')
rm(list = ls())
### The IS subgroups indentification and wilicox test were used to compare the ICI density between different 
# IS subgroups.
load("~/project2022/Figure8/Fig8.Rdata")
library(monocle)
cds <- Fig8$cds
cds@phenoData@data$Cluster <- ifelse(cds@phenoData@data$State==6,'IS2A',
  ifelse(cds@phenoData@data$State==5 | cds@phenoData@data$State==11,'IS2B',
    ifelse(cds@phenoData@data$State==3 | cds@phenoData@data$State==1 | cds@phenoData@data$State==4,'IS1A','IS1B')))
plot_cell_trajectory(cds = cds,color_by = "State")+
  scale_color_manual(values = c('#778899','#778899','#778899','#778899','#FF6600','#FF6600','#FF6600','#FF6600','#FF6600','#FF6600','#778899'))
plot_cell_trajectory(cds = cds,color_by = 'IS')+scale_color_manual(values = c('#87CEFA','#4682B4'))
plot_cell_trajectory(cds = cds,color_by = "Cluster")+scale_color_manual(values = c('#FF6633','#99CC33','#00CCCC','#9966CC'))
#### We also compared and visulized the immune cell abundance between different IS subgroups.
load("~/project2022/MetaImmunecellabun.Rdata")
cellID <- rownames(ssgsva)
cellID[c(3,6)] <- c("Effector memory CD8 T cell",'Effector memory CD4 T cell')
rownames(ssgsva) <- cellID
d <- cds@phenoData@data
IS1 <- d[d$Cluster=='IS1A' | d$Cluster=='IS1B',][,c(1,7)]
idIS1 <- rownames(d[d$Cluster=='IS1A' | d$Cluster=='IS1B',])
IS1ssgsva <- t(ssgsva[,idIS1])
IS1ssgsva <- as.data.frame(IS1ssgsva)
IS1ssgsva$Cluster <- IS1$Cluster
IS1ssgsva<-tidyr::gather(
  data=IS1ssgsva,
  key="cell",
  value="abundance",
  "Activated CD8 T cell":"Neutrophil"
)##################重点代码
library(ggplot2)
library(ggprism)
library(patchwork)
p1 <- ggplot(data =IS1ssgsva ,mapping = aes(x = cell,y= abundance,fill=Cluster))+geom_boxplot()+
  theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 45),
    legend.position = 'top')+scale_fill_manual(values = c('#FF6633','#99CC33'))+
  labs(x='',y='density')
p1
IS2 <- d[d$Cluster=='IS2A' | d$Cluster=='IS2B',][,c(1,7)]
idIS2 <- rownames(d[d$Cluster=='IS2A' | d$Cluster=='IS2B',])
IS2ssgsva <- t(ssgsva[,idIS2])
IS2ssgsva <- as.data.frame(IS2ssgsva)
IS2ssgsva$Cluster <- IS2$Cluster
IS2ssgsva<-tidyr::gather(
  data=IS2ssgsva,
  key="cell",
  value="abundance",
  "Activated CD8 T cell":"Neutrophil"
)##################重点代码
p2 <- ggplot(data =IS2ssgsva ,mapping = aes(x = cell,y= abundance,fill=Cluster))+geom_boxplot()+
  theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 45),
    legend.position = 'top')+scale_fill_manual(values = c('#00CCCC','#9966CC'))+
  labs(x='',y='density')
p2
p1+p2
#####wilicox test
output <- vector()
cell <- rownames(ssgsva)
for(i in cell){
  a = IS1ssgsva[IS1ssgsva$Cluster=='IS1A',i]
  b = IS1ssgsva[IS1ssgsva$Cluster=='IS1B',i]
  c = wilcox.test(a,b)$p.value
  if(length(output)==0){
    output = c
  }else{
    output = c(output,c)
  }
}
IS1subP <- cbind(cell=cell,p=output)
for(i in cell){
  a = IS2ssgsva[IS2ssgsva$Cluster=='IS2A',i]
  b = IS2ssgsva[IS2ssgsva$Cluster=='IS2B',i]
  c = wilcox.test(a,b)$p.value
  if(length(output)==0){
    output = c
  }else{
    output = c(output,c)
  }
}
IS2subP <- cbind(cell=cell,p=output)
subgroups <- list(IS1subP=IS1subP,IS2subP)
save(subgroups,file = 'ISsubgroups.Rdata')
###### correlationship between PC1/PC2 and immune cell abundance.
load("~/project2022/Figure8/Fig8.Rdata")
a2 <- Fig8$cds@reducedDimS#########pca value
load("~/project2022/MetaImmunecellabun.Rdata")
pc1 <- data.frame()
cell <-rownames(ssgsva)
for(i in cell){
  s1 =a2[1,]
  s2 =ssgsva[i,]
  s3 =cor.test(s1,s2)
  cor = s3$estimate
  p = s3$p.value
  s4 = cbind(cor=cor,pvalue=p)
  if(dim(pc1)[1]==0){
    pc1 = s4
  }else{
    pc1 = rbind(pc1,s4)
  }
}
rownames(pc1) <- cell
pc2 <- data.frame()
for(i in cell){
  s1 =a2[2,]
  s2 =ssgsva[i,]
  s3 =cor.test(s1,s2)
  cor = s3$estimate
  p = s3$p.value
  s4 = cbind(cor=cor,pvalue=p)
  if(dim(pc2)[1]==0){
    pc2 = s4
  }else{
    pc2 = rbind(pc2,s4)
  }
}
rownames(pc2) <- cell
pc <- rbind(pc1,pc2)
pc <- cbind(pc,cell=rep(cell,2),pc=c(rep('pc1',28),rep('pc2',28)))
pc <- as.data.frame(pc)
pc[,1:2] <- apply(pc[,1:2],2,as.numeric)
library(ggnewscale)
library(corrplot)
library(ggcorrplot)
geom_rectriangle <- function(mapping = NULL, data = NULL,
  stat = "identity", position = "identity",
  ...,
  linejoin = "mitre",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRectriangle,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}
GeomRectriangle <- ggproto(
  "GeomRectriangle", Geom,
  default_aes = aes(r = 1, colour = "grey35", fill = NA, size = 0.25, linetype = 1,
    alpha = NA,type = "upper"),
  required_aes = c("x", "y"),
  draw_panel = function(self, data, panel_params, coord, linejoin = "mitre",type = "upper") {
    aesthetics <- setdiff(names(data), c("x", "y"))
    
    polys <- lapply(split(data, seq_len(nrow(data))), function(row) {
      rectriangle <- point_to_rectriangle(row$x, row$y, row$r, row$type)
      aes <- new_data_frame(row[aesthetics])[rep(1, 4), ]
      GeomPolygon$draw_panel(cbind(rectriangle, aes), panel_params, coord)
    })
    
    ggplot2:::ggname("geom_rectriangle", do.call("grobTree", polys))
  },
  draw_key = draw_key_polygon
)
library(vctrs)
library(grid)
geom_rectriangle <- function(mapping = NULL, data = NULL,
  stat = "identity", position = "identity",
  ...,
  linejoin = "mitre",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRectriangle,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      linejoin = linejoin,
      na.rm = na.rm,
      ...
    )
  )
}

GeomRectriangle <- ggproto(
  "GeomRectriangle", Geom,
  default_aes = aes(r = 1, colour = "grey35", fill = NA, size = 0.25, linetype = 1,
    alpha = NA,type = "upper"),
  required_aes = c("x", "y"),
  draw_panel = function(self, data, panel_params, coord, linejoin = "mitre",type = "upper") {
    aesthetics <- setdiff(names(data), c("x", "y"))
    
    polys <- lapply(split(data, seq_len(nrow(data))), function(row) {
      rectriangle <- point_to_rectriangle(row$x, row$y, row$r, row$type)
      aes <- new_data_frame(row[aesthetics])[rep(1, 4), ]
      GeomPolygon$draw_panel(cbind(rectriangle, aes), panel_params, coord)
    })
    
    ggplot2:::ggname("geom_rectriangle", do.call("grobTree", polys))
  },
  draw_key = draw_key_polygon
)

point_to_rectriangle <- function(x, y, r, type = type) {
  r <- 0.5 * sign(r) * sqrt(abs(r))
  #r0 = 0.5
  xmin <- - r + x
  xmax <- r + x
  ymin <- - r + y
  ymax <- r + y 
  if(type == "upper"){
    df = new_data_frame(list(
      y = c(ymax, ymax, ymin, ymax),
      x = c(xmin, xmax, xmin, xmin)
    ))
  }else if(type == "lower"){
    df = new_data_frame(list(
      y = c(ymax, ymin, ymin, ymax),
      x = c(xmax, xmax, xmin, xmax) 
    ))
  }
  df
}
pc$label = ifelse(pc$pvalue < 0.0001,'****',ifelse(
  pc$pvalue < 0.001,'***',ifelse(pc$pvalue < 0.01,'**',ifelse(
    pc$pvalue < 0.05,'*',''
  ))
))
ggplot()+
  geom_rectriangle(data = pc, aes(pc, cell, fill = -log10(pvalue)),type = "upper", r = 1)+
  scale_fill_gradient(low = "white", high = "#66a3b1")+
  ggnewscale::new_scale_fill()+
  geom_rectriangle(data = pc, aes(pc, cell, fill = cor),type = "lower", r = 1)+
  scale_fill_gradient2(high = "red", mid = "white",low = "blue")+
  labs(x = "", y = "")+
  geom_text(data = pc,aes(pc,cell,label=label),
    nudge_y = 0.1)+ 
  theme(axis.text.x = element_text(angle = 30, 
    vjust = 0.1))
###Now we perform the WGCNA to identify the immune gene co-expression modules of BRCA ---
rm(list = ls())  
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table) #多核读取文件
### 启用WGCNA多核计算
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 

################################# 0.输入数据准备 ################################
### 读取表达矩阵
load("~/project2022/noneBatch.Rdata")
load("~/project2022/ImmuenGeneList.Rdata")
s <- as.data.frame(s)
geneid <- intersect(gene,rownames(s))
s <- s[geneid,]
data <- s
### 筛选MAD前1203的基因
keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:1203],]

### 创建datTraits，包含分组、表型等信息
datTraits <- data.frame(row.names = colnames(data),group=colnames(data))
fix(datTraits)

### 给分组加上编号
grouptype <- data.frame(group=sort(unique(datTraits$group)),
  groupNo=1:length(unique(datTraits$group)))
# fix(grouptype)
datTraits$groupNo = "NA"
for(i in 1:nrow(grouptype)){
  datTraits[which(datTraits$group == grouptype$group[i]),'groupNo'] <- grouptype$groupNo[i]}
datTraits

### 转置
datExpr0 <- as.data.frame(t(keep_data))

############################## 1.判断数据质量 ################################

### 判断数据质量--缺失值
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
      collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
      paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK

### 绘制样品的系统聚类树
if(T){
  #针对样本做聚类树
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
    cex.axis = 1, cex.main = 1,cex.lab=1)
  ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
    colors = rainbow(length(table(datTraits$group))), 
    signed = FALSE)
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1),cex=0.8)
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6)
  plotDendroAndColors(sampleTree, sample_colors,
    groupLabels = "trait",
    cex.dendroLabels = 0.8,
    marAll = c(1, 4, 3, 1),
    cex.rowText = 0.01,
    main = "Sample dendrogram and trait" )
  ## Plot a line to show the cut
  # abline(h = 23500, col = "red") #根据实际情况而定
}

##若存在显著离群点；剔除掉
if(F){
  clust <- cutreeStatic(sampleTree, cutHeight = 23500, minSize = 10) # cutHeight根据实际情况而定
  table(clust)
  keepSamples <- (clust==1)
  datExpr0 <- datExpr0[keepSamples, ]
  datTraits <- datTraits[keepSamples,]
  dim(datExpr0) 
}

### 判断数据质量 : PCA进行分组查看
rm(list = ls())  
load("step1_input.Rdata")
group_list <- datTraits$group
dat.pca <- PCA(datExpr, graph = F) 
pca <- fviz_pca_ind(dat.pca,
  title = "Principal Component Analysis",
  legend.title = "Groups",
  geom.ind = c("point","text"), #"point","text"
  pointsize = 2,
  labelsize = 4,
  repel = TRUE, #标签不重叠
  col.ind = group_list, # 分组上色
  axes.linetype=NA,  # remove axeslines
  mean.point=F#去除分组中心点
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) #坐标轴的纵横比
pca
ggsave(pca,filename= "step1_Sample PCA analysis.pdf", width = 8, height = 8)

##保存数据
datExpr <- datExpr0
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits,file="step1_input.Rdata")

############################### 2.挑选最佳阈值power ###################################
rm(list = ls())  
load("step1_input.Rdata")
R.sq_cutoff = 0.8  #设置R^2 cut-off，默认为0.85
if(T){
  # Call the network topology analysis function
  #设置power参数选择范围
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
    networkType = "unsigned",
    powerVector = powers, 
    RsquaredCut = R.sq_cutoff,  
    verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  pdf("step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: 寻找拐点，确认最终power取值
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate  #查看估计的最佳power
# power = sft$powerEstimate
power = 2

# 若无向网络在power小于15或有向网络power小于30内，没有一个power值使
# 无标度网络图谱结构R^2达到0.8且平均连接度在100以下，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if(is.na(power)){
  # 官方推荐 "signed" 或 "signed hybrid"
  # 为与原文档一致，故未修改
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
    ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
      ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
        ifelse(type == "unsigned", 6, 12))      
    )
  )
}

save(sft, power, file = "step2_power_value.Rdata")




##################### 3.一步法构建加权共表达网络，识别基因模块 ####################
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 30,    ##越大模块越少
    mergeCutHeight = 0.25, ##越大模块越少
    numericLabels = TRUE, 
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors) 
  # power: 上一步计算的软阈值
  # maxBlockSize:计算机能处理的最大模块的基因数量(默认5000),16G内存可以处理2万个，
  # 计算资源允许的情况下最好放在一个block里面。
  # corType：计算相关性的方法；可选pearson(默认)，bicor。后者更能考虑离群点的影响。
  # networkType:计算邻接矩阵时，是否考虑正负相关性；默认为"unsigned",可选"signed", "signed hybrid"
  # TOMType：计算TOM矩阵时，是否考虑正负相关性；默认为"signed",可选"unsigned"。但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义。
  # numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  # saveTOMs：最耗费时间的计算，可存储起来供后续使用，
  # mergeCutHeight: 合并模块的阈值，越大模块越少,一般为0.25
  # minModuleSize: 每个模块里最少放多少个基因，设定越大模块越少
  # 输出结果根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
  # **0 (grey)**表示**未**分入任何模块的基因。
}

## 模块可视化，层级聚类树展示各个模块
if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "step3_genes_modules.Rdata")


#####################  分布法完成网络构建，一般不用 #################################
if(F){
  ## 构建加权共表达网络分为两步：
  ## 1. 计算邻近值，也是就是两个基因在不同样品中表达量的表达相关系数(pearson correlation rho)，
  ## 2. 计算topology overlap similarity (TOM)。 用TOM表示两个基因在网络结构上的相似性，即两个基因如果具有相似的邻近基因，这两个基因更倾向于有相互作用。
  
  ###(1)网络构建 Co-expression similarity and adjacency 
  adjacency = adjacency(datExpr, power = power) 
  
  ###(2) 邻近矩阵到拓扑矩阵的转换，Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM
  
  ###(3) 聚类拓扑矩阵 Clustering using TOM
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)
  ## 这个时候的geneTree与一步法的 net$dendrograms[[1]] 性质类似，但是还需要进行进一步处理
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
    labels = FALSE, hang = 0.04)
  
  ###(4) 聚类分支的修整 dynamicTreeCut 
  ################# set the minimum module size ##############################
  minModuleSize = 30
  ####
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
    deepSplit = 2, pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize)
  table(dynamicMods)
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram and module colors")
  
  ###(5) 聚类结果相似模块的融合 Merging of modules whose expression profiles are very similar
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  #一般选择 height cut 为0.25,对应于有75%相关性，进行融合
  ###################### set  Merging height cut  ################################
  MEDissThres = 0.5
  ####
  # Plot the result
  plot(METree, main = "Clustering of module eigengenes",
    xlab = "", sub = "")
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  # 统计mergedmodule
  table(mergedColors)
  
  ### (6) plot the gene dendrogram 
  pdf(file = "step3_stepbystep_DynamicTreeCut_genes-modules.pdf", width = 16,height = 12)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
    c("Dynamic Tree Cut", "Merged dynamic"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  ### 保存数据
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  # Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, geneTree, 
    file = "step3_stepByStep_genes_modules.Rdata")
  
}

####################### 4.关联基因模块与表型 #####################################
rm(list = ls())  
load(file = "step1_input.Rdata")
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")

## 模块与表型的相关性热图
if(F){
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0+datTraits$group)
  colnames(design) <- levels(datTraits$group) #get the group
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
    signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("step4_Module-trait-relationship_heatmap.pdf",
    width = 2*length(colnames(design)), 
    height = 0.6*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) #留白：下、左、上、右
  labeledHeatmap(Matrix = moduleTraitCor,
    xLabels = colnames(design),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = F,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = F,
    cex.text = 0.5,
    zlim = c(-1,1), 
    main = "Module-trait relationships")
  dev.off()
  save(design, file = "step4_design.Rdata")
}
### 模块与IS的相关性boxplot图
IS <- read.csv("~/project2022/Immunecluster/Immunecluster.k=2.consensusClass.csv", header=FALSE)
IS <- as.data.frame(IS)
IS$IS <- ifelse(IS$V2==1,'IS1','IS2')
datTraits <- IS[,-2]
MEs$Id <- rownames(MEs)
mes_group <- merge(MEs,datTraits,by.x ="Id",by.y = 'V1') 
library(gplots)
library(ggpubr)
library(grid)
library(gridExtra)
library(ggprism)
mes_group1 <- mes_group[,-4]
mes_group1<-tidyr::gather(
  data=mes_group1,
  key="module",
  value="value",
  MEblue:MEturquoise
)##################重点代码
save(mes_group1,MEs,file = 'MEs.Rdata')
load("~/project2022/MEs.Rdata")
pdf(file = 'eig.pdf',width = 3.27,height = 5.14)
p <- ggplot(data =mes_group1 ,mapping = aes(x = module,y= value,fill=IS))+geom_boxplot()+
  theme_prism()+theme(axis.text.x = element_text(vjust =1,hjust = 1,angle = 45),
    legend.position = 'top')+scale_fill_manual(values = c('#00CCCC','#9966CC'))+
  labs(x='',y='module eigengenes')
p
wilcox.test(mes_group1[mes_group1$IS=='IS1' & mes_group1$module=='MEblue',]$value,mes_group1[mes_group1$IS=='IS2' & mes_group1$module=='MEblue',]$value)
wilcox.test(mes_group1[mes_group1$IS=='IS1' & mes_group1$module=='MEturquoise',]$value,mes_group1[mes_group1$IS=='IS2' & mes_group1$module=='MEturquoise',]$value)
rm(list = ls())

#################### 6. 选择感兴趣基因模块进行GO分析 ####################
setwd("/Users/llls2012163.com/project2022/Figure9")
rm(list = ls())  
load(file = 'step1_input.Rdata')
load(file = "step2_power_value.Rdata")
load(file = "step3_genes_modules.Rdata")
load(file = "step4_design.Rdata")

### 条件设置
OrgDb = "org.Hs.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"
table(moduleColors)
choose_module <- c("blue","turquoise")

if(T){
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gene_module <- data.frame(gene=colnames(datExpr),
    module=moduleColors)
  write.csv(gene_module,file = "step6_gene_moduleColors.csv",row.names = F, quote = F) 
  tmp <- bitr(gene_module$gene,fromType = genetype,  # "SYMBOL"   "ENSEMBL"
    toType = "ENTREZID",
    OrgDb = OrgDb )
  gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  ###run go analysis
  Goblue <- enrichGO(choose_gene_module_entrz[choose_gene_module_entrz$module=='blue',]$ENTREZID,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',pvalueCutoff = 1,
    pAdjustMethod = "BH",minGSSize = 4,maxGSSize = 100,ont = 'ALL')
  Goturquoise <- enrichGO(choose_gene_module_entrz[choose_gene_module_entrz$module=='turquoise',]$ENTREZID,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',pvalueCutoff = 1,
    pAdjustMethod = "BH",minGSSize = 4,maxGSSize = 100,ont = 'ALL')
  picktop5 <- function(x){
    a <- x[x$ONTOLOGY=='BP',]
    a1 <- a[order(a$pvalue),][1:5,]
    b <- x[x$ONTOLOGY=='CC',]
    b1 <- b[order(b$pvalue),][1:5,]
    c <- x[x$ONTOLOGY=='MF',]
    c1 <- c[order(c$pvalue),][1:5,]
    d <- rbind(a1,b1,c1)
    return(d)
  }
  s1 <- picktop5(Goblue)
  s2 <- picktop5(Goturquoise)
  library(ggplot2)
  p1 <- ggplot(s1,aes(x=GeneRatio,y=Description,color=-log10(pvalue)))+geom_point(aes(size=Count))+
    facet_grid(ONTOLOGY~.,scales = 'free_y')+theme_light()+
    scale_color_gradient(low = 'green',high = 'red')
  p2 <- ggplot(s2,aes(x=GeneRatio,y=Description,color=-log10(pvalue)))+geom_point(aes(size=Count))+
    facet_grid(ONTOLOGY~.,scales = 'free_y')+theme_light()+
    scale_color_gradient(low = 'green',high = 'red')
  p1
  p2
}
########模块特征向量与IS主成分之间的相关性
rm(list = ls())
load("/Users/llls2012163.com/project2022/Figure9/step3_stepByStep_genes_modules.Rdata")
load("/Users/llls2012163.com/project2022/Figure8/Fig8.Rdata")
pca <- t(Fig8$cds@reducedDimS)
colnames(pca) <- c('PCA1','PCA2')
s <- cbind(MEs,pca)
library(ggplot2)
library(ggprism)
p1 <- ggplot(s,aes(MEblue,PCA1))+geom_point(size=0.8,color='black')+ theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),panel.border = element_rect(fill = NA,colour = 'black',size = 2))+
  geom_smooth(method = "lm")
p1
p2 <- ggplot(s,aes(MEturquoise,PCA1))+geom_point(size=0.8,color='black')+ theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),panel.border = element_rect(fill = NA,colour = 'black',size = 2))+
  geom_smooth(method = "lm")
p2
cor.test(s$MEblue,s$PCA1)
cor.test(s$MEturquoise,s$PCA1)
#Here we identify the prognostic value of feature vector with blue and turquoise module.
rm(list = ls())
load("/Users/llls2012163.com/project2022/Figure9/step3_stepByStep_genes_modules.Rdata")
load("/Users/llls2012163.com/project2022/Figure4/preLogIS.Rdata")
TreatTCGAid <- function(x){
  library(limma)
  a <- strsplit2(x,'-')[,1:3]
  b <- paste(a[,1],'-',a[,2],sep ='')
  c <- paste(b,'-',a[,3],sep = '')
  return(c)
}
id <- rownames(MEs)[1:1110]
id <- TreatTCGAid(id)
id <- c(id,rownames(MEs)[1111:3014])
MEs <- cbind(id,MEs)
MEs <- MEs[!duplicated(MEs$id),]
rownames(MEs) <- MEs$id
id1 <- intersect(id,rownames(preLog))
preLog <- preLog[id1,]
MEs <- MEs[id1,]
preLog1 <- cbind(MEs,preLog)[,-c(5,8)]
###初始化用于输出结果的表格
library(survival)
outTab <- data.frame()
s <- colnames(preLog1)[2:4]
for(i in s){
  expr=preLog1[,i]
  cox=coxph(Surv(Day,statu)~ expr,preLog1)
  coxsummary=summary(cox)
  outTab=rbind(outTab,cbind(gene=i,HR=round(coxsummary$coefficients[,'exp(coef)'],2),
    ###HR表示风险比
    z=round(coxsummary$coefficients[,'z'],2),##z值
    '95%CI'=paste(round(coxsummary$conf.int[,3],2),
      round(coxsummary$conf.int[,4],2),sep = '-'),##95%CI
    pvalue=round(coxsummary$coefficients[,'Pr(>|z|)'],2)))
  ##p值
}### No significant difference with anly module.
## log rank to test the prognostic value
preLog1$blue <- ifelse(preLog1$MEblue > mean(preLog1$MEblue),'High','Low')
preLog1$turq <- ifelse(preLog1$MEturquoise > mean(preLog1$MEturquoise),'High','Low')
library(survminer)
library(survival)
library(ggprism)
fit <- survfit(Surv(Day, statu) ~ turq###can be change to blue, data = preLog1)
ggsurvplot(fit,
  pval = TRUE, conf.int = FALSE,
  risk.table = TRUE, # Add risk table
  risk.table.col = "strata", # Change risk table color by groups
  linetype = "strata", # Change line type by groups
  surv.median.line = "hv", # Specify median survival
  ggtheme = theme_prism(), # Change ggplot2 theme
  palette = c("red", "blue"))
save(preLog1,file = 'moduleprognosis.Rdata')
rm(list = ls())
#Here we detected the prognosis of T N M S PCA IS age------------------------------------------------------
load("/Users/llls2012163.com/project2022/Figure4/Fig4cor.Rdata")
load("/Users/llls2012163.com/project2022/Figure9/step3_stepByStep_genes_modules.Rdata")
load("/Users/llls2012163.com/project2022/Figure9/moduleprognosis.Rdata")
load("/Users/llls2012163.com/project2022/clinicTCGA.Rdata")
load("/Users/llls2012163.com/R-scrpit/ImmuneFilither/ICIScore validata/metabric/(MetaBric_CorhortClin.Rdata")
preLog1 <- preLog1[,-10]
tnm <- Fig4cor[[1]]
tnm1 <- preLog1[rownames(tnm),]
tnm <- cbind(tnm,tnm1)
tnm <- tnm[,-c(2,7,8,9,10)]
tnm$M <- ifelse(tnm$M=='MO','M0','M1')
rownames(clinical) <- clinical$ID
clinical <- clinical[rownames(tnm),]
tnm <- cbind(tnm,Age=clinical$Age)
tnm$'T' <- ifelse(tnm$'T'=='T1' | tnm$'T'=='T2',0,1)
tnm$N <- ifelse(tnm$N=='N0',0,1)
tnm$M <- ifelse(tnm$M=='M0',0,1)
tnm$S <- ifelse(tnm$S=='Stage I' | tnm$S=='Stage II',0,1)
metaclin <- MetaBric_CorhortClin$Clinical
rownames(metaclin) <- metaclin$PATIENT_ID
id <- intersect(metaclin$PATIENT_ID,preLog1$id)
metaclin <- metaclin[id,]
metaclin <- metaclin[,c(1,13)]
preLog2 <- preLog1[id,]
preLog2 <- cbind(preLog2,Age=metaclin$AGE_AT_DIAGNOSIS)
load("/Users/llls2012163.com/project2022/clinicTCGA.Rdata")
rownames(clinical) <- clinical$ID
id1 <- intersect(clinical$ID,preLog1$id)
clinical <- clinical[id1,]
preLog3 <- preLog1[id1,]
preLog3 <- cbind(preLog3,Age=clinical$Age)
preLog4 <- rbind(preLog2,preLog3)
##prognostic value of TNMS
library(survival)
outTab <- data.frame()
s <- c('T','N','M','S')
for(i in s){
  expr=tnm[,i]
  cox=coxph(Surv(Day,statu)~ expr,tnm)
  coxsummary=summary(cox)
  outTab=rbind(outTab,cbind(gene=i,HR=round(coxsummary$coefficients[,'exp(coef)'],2),
    ###HR表示风险比
    z=round(coxsummary$coefficients[,'z'],2),##z值
    '95%CI'=paste(round(coxsummary$conf.int[,3],2),
      round(coxsummary$conf.int[,4],2),sep = '-'),##95%CI
    pvalue=round(coxsummary$coefficients[,'Pr(>|z|)'],2)))}
##prognostic value of PCA IS age
  preLog4$Age <- ifelse(preLog4$Age < 65,0,1)
  outTab1 <- data.frame() 
s1 <- c("IS","blue", "turquoise","Age")  
for(i in s1){
  expr=preLog4[,i]
  cox=coxph(Surv(Day,statu)~ expr,preLog4)
  coxsummary=summary(cox)
  outTab1=rbind(outTab1,cbind(gene=i,HR=round(coxsummary$coefficients[,'exp(coef)'],2),
    ###HR表示风险比
    z=round(coxsummary$coefficients[,'z'],2),##z值
    '95%CI'=paste(round(coxsummary$conf.int[,3],2),
      round(coxsummary$conf.int[,4],2),sep = '-'),##95%CI
    pvalue=round(coxsummary$coefficients[,'Pr(>|z|)'],2)))}
##multcox regression test
mult <- coxph(Surv(Day,statu)~IS+turquoise+Age,preLog4)
summMuticox <- summary(mult)
outTab3 <- cbind(HR=round(summMuticox$coefficients[,'exp(coef)'],4),
  z=round(summMuticox$coefficients[,'z'],4),
  '95%CI'=paste(round(summMuticox$conf.int[,3],4),
    round(summMuticox$conf.int[,4],4),sep = '-'),
  pvalue=round(summMuticox$coefficients[,'Pr(>|z|)'],4))
# All of those results including unvariate and multivariate Cox regression seems unsuit to be demostrated
# in the paper. T N M S were not prognostic factors in the univariate cox regression. The IS turquoise module
# were prognostic factors in the univariate cox regression. Only the age were the prognostic factor in the 
# multivariate cox regression. So it seems that those results were unsuit to be demostrated in the paper.
# We finally finished the study here. The last works were checked those results and prepare the manuscript.