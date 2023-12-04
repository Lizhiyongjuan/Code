setwd("xxx")
load("TCGA.rdata")
TCGA <- t(TCGA)
TCGA <- as.data.frame(TCGA)
TCGA$sample <- rownames(TCGA)
suv <- read.table("TCGAsurvival.tsv",header = T,sep = "\t")
TCGASUV <- merge(TCGA,suv,by='sample')


#####
library(survival)
library(survminer)
setwd("D://analysis/课题计划/TCGA/survival/")
surv <- read.table("5GENE.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
surv <- TCGASUV
surv$OS.time <- as.numeric(surv$OS.time)
surv$OS <- as.numeric(surv$OS)

surv$OS.time <- surv$OS.time/30
# surv$OS.time <- ceiling(surv$OS.time)
# surv$OS.time <- floor(surv$OS.time)
#surv=subset(surv,surv$OS.time<60)
res.cut <- surv_cutpoint(surv, 
                         time = "OS.time",
                         event = "OS", 
                         variables = c("FCN3","CDH19","AP3M2","MAP2K1"))
res.cut
# cutpoint statistic
# cutpoint statistic
# FCN3   0.5560135  3.215261
# HMGN2  5.4514344  3.587267
# COL1A2 2.5627251  1.818004
# ID1    2.3470943  1.652753
# HMGB2  4.0057690  3.980841
# FCN3   0.5560135  3.938122
# MAP2K1 3.8514342  1.602456
# AP3M2  0.8392200  2.460216
# CDH19  0.1009096  2.171856
# CCNA2 2.9518820  4.132955
#dataframe$group <- cutree(cutpoint$points, dataframe$Var1)
#median
#median(surv$ZNF318)
surv$group <- ifelse(surv$FCN3 > 1.3593486 ,"High","Low")
surv$group <- ifelse(surv$HMGB2 > 5.3956620 ,"High","Low")
surv$group <- ifelse(surv$COL15A1 > 2.5737815 ,"High","Low")
surv$group <- ifelse(surv$PDGFRL > 0.5182962 ,"High","Low")
surv$group <- ifelse(surv$DDX17 > 6.5984700 ,"High","Low")

surv$group <- factor(surv$group, levels = c("Low","High")) 
class(surv$group)
table(surv$group)
#install.packages("survival")
library(survival)
fitd <- survdiff(Surv(OS.time, OS) ~ group,
                 data = surv,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)

#2.2 
fit <- survfit(Surv(OS.time, OS)~ group, data = surv)
summary(fit)

#3. 

plot(fit, conf.int = T,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Time(Months)",
     ylab = "Survival probablity(%)"
)

legend("topright",
       title = "Group",
       c("Low", "High"),
       lwd = 2, lty = 1,
       col = c("blue", "red"))

p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
text(25, 0.2, p.lab)
dev.off()



library(survminer)

ggsurvplot(fit,
           data = surv,
           pval = p.lab,
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           palette = "lancet", #
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
dev.off()




####cox####

rm(list = ls())
setwd("xxx")

library(survival)
library(forestplot)
library(tidyverse)
library(survminer)

#xena：https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

surv = read.table(file = 'diaease.tsv', sep = '\t', header = TRUE) 

surv$sample <- gsub("-",".",surv$sample) #
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]

expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]

res_deseq2 <- as.data.frame(res)%>% 
  arrange(pvalue) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, pvalue < 0.05)

deg_expr <- expr[rownames(res_deseq2),] %>% t() %>% as.data.frame()
surv.expr <- cbind(surv,deg_expr)
#write.table(surv.expr,"disease.txt",sep = "\t")



# 导入数据集，假设为dataframe类型，其中含有生存时间、生存状态、以及需要计算的变量

# 计算最优分割点
surv <- read.table("5gene新.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
res.cut <- surv_cutpoint(surv, #数据集
                         time = "OS.time", #生存状态
                         event = "OS", #生存时间
                         variables = c("FCN3", "MAP2K1", "AP3M2","CDH19"))
res.cut
# cutpoint statistic
# FCN3   0.5560135  3.215261
# MAP2K1 3.8514342  1.942675
# AP3M2  1.6043188  2.697498
# CDH19  0.1030765  2.167039
#res.cut

surv$group <- ifelse(surv$CDH19  > 0.1030765 ,"High","Low")
surv$group <- factor(surv$group, levels = c("Low","High"))

# 进行Cox回归分析
cox_model <- coxph(Surv(OS.time,OS) ~ group, data = surv)


summary(cox_model)

topgene <- read.table("5genecox.txt",stringsAsFactors = F,header = T)
topgene <- read.table("HR.txt",stringsAsFactors = F,header = T)

tabletext <- cbind(c("Gene",topgene$Genes),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$Lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$Upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$Pvalue),3),nsmall = 3)))
tabletext <- topgene

forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$Lower)), 
           upper=c(NA,as.numeric(topgene$Upper)),
           graph.pos=5,# 
           graphwidth = unit(.25,"npc"),# 
           fn.ci_norm="fpDrawCircleCI",#
           col=fpColors(box="red", lines="black", zero = "black"),
           
           boxsize=0.1,# 
           lwd.ci=1.5,
           ci.vertices.height = 0.1,ci.vertices=T,
           zero=1,#
           lwd.zero=1,
           xticks = c(0.5,1,1.5,2),
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), 
                           "2" = gpar(lwd=1.5, col="black"), 
                           "6" = gpar(lwd=2, col="black")), 
           lineheight = unit(.75,"cm"),
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)



DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter( padj < 0.05)
write.table(DEG,"HCCDEG.txt",sep = "\t")














