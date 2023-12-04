setwd("XXX")
library(GEOquery)
gset=getGEO('XXX',destdir = '.',getGPL = F)
exprSet=exprs(gset[[1]])
exprSet=as.data.frame(exprSet)
pdata=pData(gset[[1]])
pdata$group <- c(rep("XXX",XXX),rep("XXX",XXX))
group_list=pdata$group
group_list2=factor(group_list,levels = c('XXX','XXX'))
library(limma)
gpl="XXX.txt"
a=1#平台文件的ID列号
b=11#平台文件的SYMBOL列号
gene_name="DDR1"#目标基因SYMBOL
out_name="data_out.txt"#输出文件名（记得加拓展名）
library(pacman)
library(GEOquery)#加载数据包
library(limma)
library(affy)
library(tidyverse)
library(stringr)
GPL_file=read.table(gpl,header=T,
                    quote="",sep="\t",dec=".",
                    comment.char="#",na.strings =c("NA"),fill=T )#读取平台文件（GPL）
gpl_file=GPL_file[,c(a,b)]#将平台文件的ID列和SYMBOL列取出
#??str_split
e <- apply(gpl_file,1,
           function(x){
             paste(x[1],
                   str_split(x[2],'///',simplify=T),
                   sep = "...")
           })
x = tibble(unlist(e))
colnames(x) <- "f" 
write.table(exprSet,"exprSet.txt",sep = "\t")
exprSet <- read.table("exprSet.txt",sep="\t",header = T)
gpl_file <- separate(x,f,c("ID","symbol"),sep = "\\...")
exp<-as.data.frame(exprSet)
colnames(exp)[1]="ID"
exp_symbol<-merge(exp,gpl_file,by="ID")
exp_symbol[exp_symbol==""]<-NA
exp_symbol<-na.omit(exp_symbol)
exp_symbol[,grep("symbol", colnames(exp_symbol))]=
  trimws(exp_symbol[,grep("symbol", colnames(exp_symbol))])#去除数据头尾空格
if("DDR1"%in% exp_symbol[,grep("symbol", colnames(exp_symbol))]== F)
{ errorCondition("XXX",class = NULL,call = NULL)}#检测是否含有目标基因
table(duplicated(exp_symbol[,ncol(exp_symbol)]))#检测重复基因数
d=data.frame(duplicated(exp_symbol[,ncol(exp_symbol)]))#将布尔值写入数据框
exp_symbol=avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$symbol)#对重复值取平均后混合
table(duplicated(rownames(exp_symbol)))#再次检测重复基因数
write.table(exp_symbol,out_name,sep="\t",#写出文件
            quote=F,
            col.names=T)
save(exp,exp_symbol,gpl_file,GPL_file,group_list1,file="ID.Rdata")
save(group_list2,file = "group.Rdata")



##DEG
setwd("XXX")
load("group.Rdata")
exp <- read.table("data_out.txt",header = T,sep="\t")
max(exp)
#exprSet1=log2(exp+1)
#max(exprSet1)
boxplot(exp)
exprSet=normalizeBetweenArrays(exp)
boxplot(exprSet)
designtable<-model.matrix(~0+factor(group_list2))
colnames(designtable)<-levels(factor(group_list2))
rownames(designtable)<-colnames(exprSet)
contrast.matrix=makeContrasts(HCC-control,levels=designtable)
#limma DEG
fit<-lmFit(exprSet,designtable)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)

DEG<-topTable(fit2,coef="HCC - control",n=Inf)
write.table(DEG,'DEG.txt',sep = "\t")


allDiff_HCC2_up=DEG[DEG$logFC>0&DEG$P.Value<0.05,]
allDiff_HCC2_down=DEG[DEG$logFC< 0&DEG$P.Value<0.05,]
exprSet2 <- exprSet
save(exprSet2,group_list2,allDiff_HCC2_up,allDiff_HCC2_down,file='XXX.Rdata')
load('XXX.Rdata')
write.table(exprSet,file='XXX.txt',quote = F,sep = '\t',col.names = NA)
write.table(allDiff_HCC2_up,file='XXX.txt',quote = F,sep = '\t',col.names = NA)
write.table(allDiff_HCC2_down,file='XXX.txt',quote = F,sep = '\t',col.names = NA)
