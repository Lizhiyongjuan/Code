
setwd('D://newanalysis/EMTA_TME/')

#install.packages('data.table')
library(data.table)

expr= fread(file='TCGA-LIHC.htseq_fpkm.tsv',header = T,sep = '\t')
library(tibble)                 
expr=column_to_rownames(expr,var = 'Ensembl_ID')
max(expr)
##
expr = 2^expr -1
max(expr)
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
exprSet <- as.data.frame(apply(expr,2,fpkmToTpm))
library(dplyr)
library(tibble)
exprSet = rownames_to_column(exprSet, var ='gene_id')
library(tidyr)
exprSet <- exprSet %>% 
  tidyr::separate(gene_id,into = c("gene_id"),sep="\\.")


load('gtf_df.Rdata')



mRNA_exprSet <- gtf_df %>% 
  dplyr::filter(type=="gene",gene_biotype=="protein_coding") %>% 
  dplyr::select(c(gene_name,gene_id,gene_biotype)) %>% 
  dplyr::inner_join(exprSet,by ="gene_id") %>% 
  tidyr::unite(gene_id,gene_name,gene_id,gene_biotype,sep = " | ")

exprSet_mRNA <-mRNA_exprSet %>%
  tidyr::separate(gene_id,c('gene_name','gene_id','gene_biotype'),
                  sep = " \\| ")
exprSet_mRNA<- exprSet_mRNA[,-c(2,3)]
exprSet_mRNA<-exprSet_mRNA[exprSet_mRNA$gene_name != 'NA',]
rownames(exprSet_mRNA)=NULL

exprSet_tcga_mRNA<-aggregate(x=exprSet_mRNA[,2:(ncol(exprSet_mRNA))],by=list(exprSet_mRNA$gene_name),FUN=mean)
exprSet_tcga_mRNA <- column_to_rownames(exprSet_tcga_mRNA, var='Group.1')
exprSet_tcga_mRNA=log2(exprSet_tcga_mRNA+1)
max(exprSet_tcga_mRNA)
library(limma)
boxplot(exprSet_tcga_mRNA)
exprSet_tcga_mRNA <- normalizeBetweenArrays(exprSet_tcga_mRNA)
boxplot(exprSet_tcga_mRNA)

save(exprSet_tcga_mRNA,file = 'LIHC_tpm.Rdata')



############CIBERSORT
#install.packages('e1071')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library("limma")         
load('LIHC_tpm.Rdata')

data=exprSet_tcga_mRNA
data=data[rowMeans(data)>0,]

write.table(data,'LIHC_tpm_for_cibersort.txt',quote = F,sep = '\t',col.names = NA)

source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "LIHC_tpm_for_cibersort.txt", perm=100, QN=TRUE)
write.table(results,"CIBERSORT-Results.txt",sep = "\t")





setwd("D://newanalysis/EMTA_TME/tcga/")

library(boot)

library(survival)      
pFilter= 0.05         

library(limma)               


load('TCGA.Rdata')

gpr=read.table('MFEM.txt',sep = "\t",header = T)

gpr=gpr$X



data=TCGA[rownames(TCGA) %in% gpr,]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group_list=ifelse(group=="0",'tumor','normal')
group_list=factor(group_list,levels = c('normal','tumor'))
library(limma)
design=model.matrix(~ group_list)

fit=lmFit(data,design)
fit=eBayes(fit) 
allDiff=topTable(fit,adjust='fdr',coef=2,number=Inf,p.value=0.05) 


write.csv(allDiff,file ='CRGDiff.csv',quote = F)



data=data[rownames(allDiff),]
gpr <- as.data.frame(gpr)
rownames(gpr) <- gpr$gpr
data=data[rownames(gpr),]

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]


keep <- rowSums(data>0) >= floor(0.75*ncol(data))
table(keep)

data<- data[keep,]

data=as.data.frame(t(data))


suv=read.table('TCGA-LIHC.survival.tsv',row.names = 1,header = T,check.names = F)
cli=dplyr::select(suv,'OS.time','OS')
colnames(cli)=c("futime", "fustat")


sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]

out=cbind(cli,data)
out=cbind(id=row.names(out),out)

write.table(out,file="expTime.txt",sep="\t",row.names=F,quote=F)

rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件



outTab=data.frame()
sigGenes=c("futime","fustat")
for(gene in colnames(rt[,3:ncol(rt)])){
  set.seed(123456)
  cox=coxph(Surv(futime, fustat) ~ rt[,gene], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP<pFilter){
    sigGenes=c(sigGenes,gene)
    outTab=rbind(outTab,
                 cbind(gene=gene,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxP))
    print(coxP)
  }
}



write.table(outTab,file="uniCox_CRG.txt",sep="\t",row.names=F,quote=F)
surSigExp=rt[,sigGenes]
surSigExp=cbind(id=row.names(surSigExp),surSigExp)
write.table(surSigExp,file="uniSigExp_CRG.txt",sep="\t",row.names=F,quote=F)



load('TCGA.Rdata')
gene=read.table('uniCox_CRG.txt',header = T)
gene=gene$gene
data=TCGA[gene,]
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)

conNum=length(group[group==1])       
treatNum=length(group[group==0])     

sampleType=ifelse(group=='1',1,2)

identical(colnames(data),colnames(TCGA))


sigVec=c()
allDiff=read.csv('CRGDiff.csv',header = T,row.names = 1)
alldiff_cox=allDiff[gene,]
pvalue=alldiff_cox$adj.P.Val
Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
sigVec=paste0(gene, Sig)



compare=data

row.names(compare)=sigVec


normal=compare[,sampleType==1]
tumor=compare[,sampleType==2]
compare=cbind(normal,tumor)


Type=c(rep("Normal",conNum), rep("Tumor",treatNum))
names(Type)=colnames(compare)
Type=as.data.frame(Type)
library(pheatmap)
pheatmap::pheatmap(compare,
                   annotation=Type,
                   breaks = c(seq(-3,3,length.out = 100)),
                   cluster_cols =F,
                   cluster_rows =T,
                   scale="row",
                   show_colnames=F,
                   show_rownames=T,
                   fontsize=12,
                   fontsize_row=10,
                   fontsize_col=6)


gene=read.table('uniCox_CRG.txt',header = T)
gene=gene$gene


################lasso,bootstrap_multicox#####################
gene=read.table('uniCox_CRG.txt',header = T)
gene=gene$gene
gene <- c("FCN3","AP3M2","CDH19")

#install.packages('survival')
library(survival)
library(survminer)
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt=rt[,c('futime','fustat',gene)]





set.seed(2)   
x=as.matrix(rt[,c(3:ncol(rt))]) 
y=data.matrix(Surv(rt$futime,rt$fustat)) 


library(glmnet)
fit=glmnet(x, y, family = "cox", alpha = 1) 
plot(fit, xvar = "lambda", label = F)

cvfit = cv.glmnet(x, y, family="cox",nfolds = 10,alpha=1) 
plot(cvfit) 

abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")

coef =coef(fit,s = cvfit$lambda.min)
index = which(coef !=0)
actCoef = coef[index] 
lassoGene = row.names(coef)[index] 
geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
geneCoef   


gene=read.table('uniCox_CRG.txt',header = T)
gene=gene[gene$gene %in% geneCoef[,1],]
write.table(gene,file = 'uniCox_lasso_CRG.txt',quote = F,sep = '\t',row.names = F)



########bootstrap_multicox####
#install.packages('boot')
library(boot)
gene=read.table('uniCox_lasso_CRG.txt',header = T)
gene <- read.table('NEWGENE.txt',sep = "\t")
# boot_coef=coef/Boot_sd
gene=gene$gene
library(survival)
library(survminer)
gene=gene$V1
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)     #读取输入文件
rt=rt[,c('futime','fustat',gene)]

# HR
cox=coxph(Surv(futime, fustat) ~.,data = rt)
ggforest(cox)

#install.packages('boot')
library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

## bootstrap，
set.seed(123456)
boot_results <- boot(data=rt, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ .)


print(boot_results)


coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)


ratio=coef/sd

MFEM=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(MFEM,file= 'CRG_coef.csv',quote = F)


load('TCGA.Rdata')
gene=read.table('uniCox_lasso_CRG.txt',header = T)
gene=gene$gene
gene <- read.table('NEWGENE.txt',sep = "\t")

gene=gene$V1
data=TCGA[gene,]


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]


coef=read.csv('CRG_coef.csv',header = T,check.names = F)
coef=coef$Coef..boot_SD
CRG_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  CRG_score=c(CRG_score,score)
}

data=as.data.frame(t(data))
data$CRG_score=CRG_score


suv=read.table('TCGA-LIHC.survival.tsv',row.names = 1,header = T,check.names = F)
cli=dplyr::select(suv,'OS.time','OS')
colnames(cli)=c("futime", "fustat")


sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]



rt=cbind(cli,data)
rt$futime=rt$futime/30

library(survival)
library(survminer)


Type=ifelse(data[,'CRG_score']<= median(rt$CRG_score), "Low", "High")
data=rt
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='PCG_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = "lancet",
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    
ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           palette = "lancet", 
           legend.labs = c("Low", "High"), 
           size = 1,
           xlim = c(0,120), 
           break.time.by = 20, 
           legend.title = "",
           surv.median.line = "hv", 
           ylab = "Survival probability (%)", 
           xlab = "Time (Months)", 
           #ncensor.plot = TRUE, 
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
save(data,file ='CRGscore_and_group.Rdata')





## CIBERSORT
immune=read.table('CIBERSORT-Results.txt',sep = '\t',header = T,check.names = F,row.names = 1)
immune=immune[,-c(23:25)]

immune=as.data.frame(t(immune))


group=sapply(strsplit(colnames(immune),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
immune=immune[,group==0]

data=as.data.frame(t(immune))


suv=read.table('TCGA-LIHC.survival.tsv',row.names = 1,header = T,check.names = F)
cli=dplyr::select(suv,'OS.time','OS')
colnames(cli)=c("futime", "fustat")


sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]


rt=cbind(cli,data)
rt$futime=rt$futime/30

immune1=immune[-c(2,5,18,20,21),]
rt1=rt[,-c(4,7,20,22,23)]

immune_p=c()
immune_figure=list()
library(survival)
library(survminer)

#install.packages('cowplot')
library(cowplot)
dir.create('k_m')
for (i in rownames(immune1)) {
  res.cut=surv_cutpoint(rt1, time="futime", event="fustat", variables=i)
  cutoff=as.numeric(res.cut$cutpoint[1])
  print(cutoff)
  Type=ifelse(data[,i]<= cutoff, "Low", "High")
  data=rt1
  data$group=Type
  data$group=factor(data$group, levels=c("Low", "High"))
  diff=survdiff(Surv(futime, fustat) ~ group, data = data)
  length=length(levels(factor(data[,"group"])))
  pValue=1-pchisq(diff$chisq, df=length-1)
  immune_p=c(immune_p,pValue)
  if(pValue<0.05){
    fit <- survfit(Surv(futime, fustat) ~ group, data = data)
    bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
    bioCol=bioCol[1:length]
    
    p=ggsurvplot(fit, 
                 data=data,
                 conf.int=F,
                 pval=pValue,
                 pval.size=6,
                 legend.title=i,
                 legend.labs=levels(factor(data[,"group"])),
                 legend = c(0.88, 0.9),
                 font.legend=12,
                 xlab="Time(Months)",
                 palette = bioCol,
                 surv.median.line = "hv",
                 risk.table=F,
                 cumevents=F,
                 risk.table.height=.25)
    ggsave2(filename = paste0('./k_m/',i,'.pdf'),width = 4,height = 4)
  }
}





library(survival)
library(survminer)
 cox = coxph(Surv(futime, fustat) ~`Mast cells resting` +
                `T cells CD8` + `Plasma cells`+
               `Macrophages M1` + `Monocytes`, data = rt)

ggforest(cox)

############
#install.packages('boot')
library(boot)
rsq <- function(formula, data, indices) { 
  d <- data[indices,] 
  fit <- coxph(formula, data=d) 
  return(fit$coefficients) 
} 

set.seed(123456)
## 1000
# boot_results <- boot(data=rt, statistic=rsq, 
#                      R=1000, formula=Surv(futime, fustat) ~ `Macrophages M1` +`Mast cells resting` +
#                        `Plasma cells` + `T cells CD4 memory resting` + `T cells CD8` + 
#                        `T cells gamma delta`)
boot_results <- boot(data=rt, statistic=rsq, 
                     R=1000, formula=Surv(futime, fustat) ~ `Mast cells resting` +
                        `T cells CD8` + `Plasma cells`+
                       `Macrophages M1` + `Monocytes`)


cox$coefficients

print(boot_results)


coef=boot_results$t0
sd=as.data.frame(boot_results$t)
sd=apply(sd, 2, sd)


ratio=coef/sd

TME=data.frame('Coef'=coef,'boot_SD'=sd,'Coef\\/boot_SD'=ratio)
write.csv(TME,file= 'TME_coef.csv',quote = F)








###########TME_score

immune=read.table('CIBERSORT-Results.txt',sep = '\t',header = T,check.names = F,row.names = 1)
immune=immune[,-c(23:25)]

immune=as.data.frame(t(immune))


group=sapply(strsplit(colnames(immune),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
immune=immune[,group==0]

data=as.data.frame(t(immune))

 data=data[,c("Mast cells resting",
              "T cells CD8",
              "Plasma cells","Macrophages M1",
              "Monocytes")]

 
# data=data[,c("Mast cells resting","T cells CD4 memory resting","T cells CD8")]

data=as.data.frame(t(data))

coef=read.csv('TME_coef.csv',header = T,check.names = F)
coef=coef$Coef..boot_SD

TME_score=c()
for (i in 1:ncol(data)) {
  score=sum(as.numeric(data[,i])*coef)
  TME_score=c(TME_score,score)
}

data=as.data.frame(t(data))
data$TME_score=TME_score


suv=read.table('TCGA-LIHC.survival.tsv',row.names = 1,header = T,check.names = F)
cli=dplyr::select(suv,'OS.time','OS')
colnames(cli)=c("futime", "fustat")

sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]



rt=cbind(cli,data)
rt$futime=rt$futime/30

library(survival)
library(survminer)


Type=ifelse(rt[,'TME_score'] <= median(rt$TME_score), "High", "Low")
data=rt
data$group=Type
data$group=factor(data$group, levels=c("Low", "High"))
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)

fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='TME_score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.88, 0.9),
             font.legend=12,
             xlab="Time(Months)",
             palette = "lancet",
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)

p    


data_immune=data
save(data_immune,file ='TMEscore_and_group.Rdata')



load('TCGA.Rdata')
gene=read.table('uniCox_lasso_CRG.txt',header = T)
gene <- read.table("NEWGENE.txt",sep = "\t")
gene=gene$V1
data=TCGA[gene,]


group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

data=as.data.frame(t(data))


immune=read.table('CIBERSORT-Results.txt',sep = '\t',header = T,check.names = F,row.names = 1)
immune=immune[,-c(23:25)]

immune=as.data.frame(t(immune))

group=sapply(strsplit(colnames(immune),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
immune=immune[,group==0]
immune=as.data.frame(t(immune))

identical(rownames(immune),rownames(data))

ss=intersect(rownames(immune),rownames(data))
immune=immune[ss,]
colnames(immune)
 immune=immune[,c("Mast cells resting",
                  "Macrophages M1","T cells CD8",
                  "Plasma cells","Monocytes")]

data=data[ss,]

rt=cbind(data,immune)


#install.packages('corrplot')
#install.packages('ggcorrplot')

library(corrplot)
M <- cor(rt)
res1 <- cor.mtest(rt, conf.level = .95) 
library(ggcorrplot)
ggcorrplot(
  M,
  hc.order = F,
  #type = "lower",
  outline.color = NA,
  ggtheme = ggplot2::theme_gray,
  colors = c("blue", "white", "red")
)

setwd("D://newanalysis/EMTA_TME/tcga/")

library(survival)
library(survminer)
load('CRGscore_and_group.Rdata')
load('TMEscore_and_group.Rdata')


data$PCG_group=ifelse(data$group=='High','PCG_high','PCG_low')
data$TME_group=ifelse(data_immune$group=='High','TME_high','TME_low')
data$group=paste0(data$PCG_group,'+',data$TME_group)
table(data$group)
data$group[data$group=='PCG_high+TME_high']='Mix'
data$group[data$group=='PCG_low+TME_low']='Mix'
diff=survdiff(Surv(futime, fustat) ~ group, data = data)
length=length(levels(factor(data[,"group"])))
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
  pValue="p<0.001"
}else{
  pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ group, data = data)
bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]

p=ggsurvplot(fit, 
             data=data,
             conf.int=F,
             pval=pValue,
             pval.size=6,
             legend.title='PCG-TME score',
             legend.labs=levels(factor(data[,"group"])),
             legend = c(0.8, 0.8),
             font.legend=12,
             xlab="Time(Months)",
             palette = "lancet",
             surv.median.line = "hv",
             risk.table=T,
             cumevents=F,
             risk.table.height=.25)
p    
table(data$group)
rt_roc=data[data$group %in% c('PCG_high+TME_low','PCG_low+TME_high'),]

rt_roc$class=ifelse(rt_roc$group=='PCG_high+TME_low',1,0)
pROC::plot.roc(rt_roc$fustat,rt_roc$class,print.auc=T)

#install.packages('timeROC')
library(timeROC)
library(survival)
library(survivalROC)

data=rt_roc

time_roc_res <- timeROC(
  T = data$futime,
  delta = data$fustat,
  marker = data$class,
  cause = 1,
  weighting="marginal",
  times = c(3 * 12, 5 * 12,7*12),
  ROC = TRUE,
  iid = TRUE
)


time_ROC_df <- data.frame(
  TP_3year = time_roc_res$TP[, 1],
  FP_3year = time_roc_res$FP[, 1],
  TP_5year = time_roc_res$TP[, 2],
  FP_5year = time_roc_res$FP[, 2],
  TP_7year = time_roc_res$TP[, 3],
  FP_7year = time_roc_res$FP[, 3]
)


bioCol=c("#0073C2","#EFC000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")

library(ggplot2)
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1.3, color = "#223D6C") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1.3, color = "#D20A13") +
  geom_line(aes(x = FP_7year, y = TP_7year), size = 1.3, color = "#11AA4D") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.7, y = 0.25, size = 3.2,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#223D6C"
  ) +
  annotate("text",
           x = 0.7, y = 0.15, size = 3.2,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#D20A13"
  ) +
  annotate("text",
           x = 0.7, y = 0.05, size = 3.2,
           label = paste0("AUC at 7 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#11AA4D"
  ) +
  labs(x = "False positive rate", y = "True positive rate") 







