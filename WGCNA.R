rm(list = ls())
setwd("XXX")#
library(WGCNA)
#BiocManager::install(c("openxlsx","WGCNA","tidyverse"))
# BiocManager::install('edgeR')
# pacman::p_load("openxlsx","WGCNA","tidyverse")

options(stringsAsFactors = FALSE)
enableWGCNAThreads()
collectGarbage()
exp <- read.table("XXX.txt",header = T,sep = "\t")
rownames(exp) = exp[,1] 
exp <- exp[,-1]
library(edgeR)
library(limma)
# row.names(exp) <- exp[,1]
# exp <- exp[,-1]
#exp1 <- log2(cpm(exp)+1)
all_exp <- exp


m.mad <- apply(all_exp,1,mad)
dataExprVar <- all_exp
a=sort(m.mad,decreasing=T)
a <- a[1:24428]
dataExprVar <- all_exp[which(a > max(quantile(a, probs=seq(0, 1, 0.25))[3],0)),]

#dataExprVar <- all_exp[which(m.mad <- a),]
#dataExprVar <- all_exp[which(m.mad <- a),]
# Check the value "NA"
datExpr0 = as.data.frame(t(dataExprVar[,]));  #
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
table(gsg$allOK)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
table(gsg$allOK)

# Hierarchical cluster of sample

sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)

par(cex = 0.6)	
par(mar =c(0,5,2,0))	
plot(sampleTree, main ="Sample clustering to detectoutliers",sub="", xlab="", cex.lab = 1.5,
     cex.axis= 1.5, cex.main = 2)
dev.off()
?plot

# Remove outliers
#h=7.2*10^4
h=98
cutHeight=h
abline(h = h, col = "red");
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
clust = cutreeStatic(sampleTree, cutHeight = h, minSize = 1)
table(clust)
keepSamples = (clust<2)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# Hierarchical cluster of sample
sampleTree = hclust(dist(datExpr), method = "average");
sizeGrWindow(12,9)

# pdf(file="Plots/sampleClustering.pdf",width=12,height=9);
par(cex = 0.6)	
par(mar =c(0,5,2,0))	
plot(sampleTree, main ="Sample clustering to detectoutliers",sub="", xlab="", cex.lab = 1.5,
     cex.axis= 1.5, cex.main = 2)



# 构建表型数据
dataExprVar = as.data.frame(t(datExpr[,]))
traitData = data.frame(sample=colnames(dataExprVar[,1:ncol(dataExprVar)]),
                       HCC=c(rep(0,15),rep(1,174)),
                       NHCC=c(rep(1,15),rep(0,174)))

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr);
traitRows = match(Samples, traitData$sample)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage();   #

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()
save(datExpr, datTraits, file = "XXX")

rm(list = ls())
lnames = load(file = "XXX");
lnames
powers = c(seq(1,20,by=1))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
par(cex = 0.6)	
par(mar =c(2,5,2,1))	
par(mfrow = c(1,2));
cex1 = 1.1;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.85,col="red",lty=4,lwd=1)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

power = sft$powerEstimate
power=17 
minModuleSize=30;
net = blockwiseModules(datExpr, power = power,
                       TOMType = "unsigned", 
                       minModuleSize = minModuleSize,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.15,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "XXX", 
                       verbose = 3)

# open a graphics window
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)

plotDendroAndColors(net$dendrograms[[1]], 
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)

dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "XXX.RData")

##################################  Step-by-step method  #######################
library(WGCNA)
rm(list = ls())
options(stringsAsFactors = FALSE);
collectGarbage()
# Load the data saved in the first part
lnames = load(file = "XXX.RData");
#The variable lnames contains the names of loaded variables.
lnames

power=4;
minModuleSize = 50;
softPower = power;
adjacency = adjacency(datExpr, power = softPower);

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, 
     xlab="", 
     sub="", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, 
     hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, 
                            distM = dissTOM,
                            deepSplit = 2, 
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
library(openxlsx)
write.xlsx(as.data.frame(table(dynamicColors)),file = "XXX.xlsx")

# Plot the dendrogram and colors underneath

plotDendroAndColors(geneTree, 
                    dynamicColors, 
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

par(cex = 0.8)	
par(mar =c(1,5,2,1))	
plot(METree,
     xlab = "", 
     sub = "",
     main = "Clustering of module eigengenes")
dev.off()

MEDissThres = 0.25

abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, 
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, 
                    hang = 0.03,
                    addGuide = TRUE, 
                    guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "D://analysis/112790WGCNA/03-networkConstruction-stepByStep.RData")


rm(list = ls())
library(WGCNA)
options(stringsAsFactors = FALSE);

# Load the expression and trait data saved in the first part
lnames = load(file = "XXX.RData");


# Load network data saved in the second part.
lnames = load(file = "XXX.RData");

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p"); 
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples); 


textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)



par(mar = c(3, 9.5, 2.7, 4));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(300),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               xLabelsAngle = 45,
               yLabelsPosition = 'left',
               x.adj.lab.y = 1
)
dev.off()
?labeledHeatmap
# Define variable weight containing the weight column of datTrait
InterestedModule = datTraits$HCC
InterestedModule = as.data.frame(InterestedModule);
names(InterestedModule) = "InterestedModule"
module <- ("brown")
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, InterestedModule, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(InterestedModule), sep="");
names(GSPvalue) = paste("p.GS.", names(InterestedModule), sep="");

names(datExpr)
names(datExpr)[moduleColors==module]
