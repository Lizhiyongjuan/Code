setwd("D://analysis/GSVA/")

load('GS.Rdata')

library(GSVA)
library(Biobase)

data=read.table('xxx.txt',header = T,sep="\t",check.names = F,row.names = 1)

uni_matrix=data

list= list


gsva_matrix<- gsva(as.matrix(uni_matrix), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

##前5个字符
rownames(gsva_matrix)=sub('^......','',rownames(gsva_matrix))

save(gsva_matrix,file = 'gsva_arrayHCC.Rdata')
load('gsva_arrayHCC.Rdata')
metabolism=gsva_matrix[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                         "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                         "HALLMARK_TGF_BETA_SIGNALING",
                         "HALLMARK_COMPLEMENT",
                         "HALLMARK_XENOBIOTIC_METABOLISM",
                         "HALLMARK_UV_RESPONSE_UP",
                         "HALLMARK_INFLAMMATORY_RESPONSE",
                         "HALLMARK_P53_PATHWAY",
                         "HALLMARK_KRAS_SIGNALING_UP",
                         "HALLMARK_IL2_STAT5_SIGNALING"),]





gene=read.csv('xxx.csv')
gene=gene$x
rt=data[gene,]


load("xxx.Rdata")
#### 先看HCC中
a=grep('xxxx',colnames(rt))
b=grep('xxx',colnames(rt))
c=grep('xxx',colnames(rt))
rt_pr=rt[,c(a,b,c)]
group1=ifelse(group_list2=='control','Normal','HCC')
group2=ifelse(group_list2=='control','Normal','HCC' )
group3=ifelse(group_list3=='control','Normal','HCC' )
group=c(group1,group2,group3)

anno=data.frame(row.names = colnames(rt_pr),group=group)

# 要取交集
metabolism_pr=metabolism[,rownames(anno)]
metabolism_pr=t(metabolism_pr)

library(Hmisc)

nc =cbind(metabolism_pr,t(rt_pr))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(metabolism_pr),(ncol(nc)-length(gene)+1):ncol(nc)]
rownames(m)=stringr::str_remove(rownames(m),pattern = 'HALLMARK_')
p = rcorr(nc)$P[1:ncol(metabolism_pr),(ncol(nc)-length(gene)+1):ncol(nc)]
library(dplyr)
library(gtools)
p <- stars.pval(p)
tmp <- as.data.frame(p)

library(pheatmap)

pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#0072b5", "white", "#bc3c29"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)



load('gsva_array.Rdata')

setwd("xxx")

load('GS.Rdata')

library(GSVA)
library(Biobase)

data=read.table('merge.normalzieHF.txt',header = T,sep="\t",check.names = F,row.names = 1)

uni_matrix=data

list= list


gsva_matrix<- gsva(as.matrix(uni_matrix), list,
                   method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

##
rownames(gsva_matrix)=sub('^......','',rownames(gsva_matrix))

save(gsva_matrix,file = 'gsva_arrayHF.Rdata')
load('gsva_arrayHF.Rdata')
metabolism=gsva_matrix[c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                         "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                         "HALLMARK_TGF_BETA_SIGNALING",
                         "HALLMARK_COMPLEMENT",
                         "HALLMARK_XENOBIOTIC_METABOLISM",
                         "HALLMARK_UV_RESPONSE_UP",
                         "HALLMARK_INFLAMMATORY_RESPONSE",
                         "HALLMARK_P53_PATHWAY",
                         "HALLMARK_KRAS_SIGNALING_UP",
                         "HALLMARK_IL2_STAT5_SIGNALING"),]



#write.table(gsva_matrix,file="HCCGSVA通路.txt",sep="\t")
## 
gene=read.csv('genes.csv')
gene=gene$x
#data=read.table('merge.normalzie.txt',header = T,sep <- "\t",check.names = F,row.names = 1)
rt=data[gene,]

#
load('HF573.Rdata')
load('HF29819.Rdata')
load('HF5406.Rdata')
#### 
a=grep('GSE57338',colnames(rt))
b=grep('GSE29819',colnames(rt))
c=grep('GSE5406',colnames(rt))
rt_MS=rt[,c(a,b,c)]
group1=ifelse(group_list3=='control','Normal','HF')
group2=ifelse(group_list1=='control','Normal','HF' )
group3=ifelse(group_list3=='control','Normal','HF' )
group=c(group1,group2,group3)

anno=data.frame(row.names = colnames(rt_MS),group=group)


metabolism_MS=metabolism[,rownames(anno)]
metabolism_MS=t(metabolism_MS)

nc =cbind(metabolism_MS,t(rt_MS))
nc=as.matrix(nc)
m = rcorr(nc)$r[1:ncol(metabolism_MS),(ncol(nc)-length(gene)+1):ncol(nc)]
rownames(m)=stringr::str_remove(rownames(m),pattern = 'HALLMARK_')
p = rcorr(nc)$P[1:ncol(metabolism_MS),(ncol(nc)-length(gene)+1):ncol(nc)]
library(dplyr)

library(gtools)
p <- stars.pval(p)
tmp <- as.data.frame(p)

library(pheatmap)



pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#20854E", "white", "#E18727"))(100),
         border_color = "white",
         cellwidth = 20, 
         cellheight = 20,
         width = 7, 
         height=9.1,
         treeheight_col = 0,
         treeheight_row = 0)

