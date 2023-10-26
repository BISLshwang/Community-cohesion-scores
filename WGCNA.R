getwd();
workingDir = ".";
setwd(workingDir);
install.packages("xfun");
library(WGCNA);
library(rhdf5);
options(stringsAsFactors = FALSE);

# -----------------------------------------------------------------------
#
# WGCNA.R
# Last modified: 2022.07
# This code is for to run WGCNA algorithm.
#
# Reference:
# Zhang, Bin, and Steve Horvath. "A general framework for weighted gene co-expression network analysis." 
# Statistical applications in genetics and molecular biology 4.1 (2005).
# Langfelder, Peter, and Steve Horvath. "Tutorials for the WGCNA Package." UCLA. Los Ageles (2014).
#
# -----------------------------------------------------------------------


# ========================================================================
# load reference samples
# row: samples, columns:gene
exp = read.csv("./example/reference_sample_liver_normalized.csv");
num = dim(exp)[1];

# split reference samples into training set and test 
trainNum = round(num*0.8);
testNum = round(num*0.2);

datTrain = as.data.frame(exp[1:trainNum, -c(1)]);
rownames(datTrain) = exp$sample[1:trainNum]
dim(datTrain);


datTest = as.data.frame(exp[trainNum+1:testNum, -c(1)]);
rownames(datTest) = exp$sample[trainNum+1:testNum]
dim(datTest)

gsg = goodSamplesGenes(datTrain, verbose=3);
gsg$allOK;

# ========================================================================
# determine parameters (Soft threshold)
powers = c(c(1:10), seq(from=12, to=20, by=2))
sft = pickSoftThreshold(datTrain, powerVector = powers, verbose=5)

sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

beta = sft$fitIndices[,];


# ========================================================================
# calculate adjacency matrix
adjacency = adjacency(datTrain, power=10);

# calculate weight matrix
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM;

# save weight matrix
write.csv(TOM, file="./example/df_tom_similarity.csv")

# ========================================================================
# detect co-expressed modules
geneTree = hclust(as.dist(dissTOM), method = "average");
minModuleSize = 50;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)

table(dynamicColors)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

genes = colnames(datTrain)
df_modules = data.frame(genes, dynamicColors)

# save module detection results
write.csv(df_modules, file = "./example/df_modules.csv")

# ========================================================================
# test module preservation using test set
setLabels = c("Train", "Test");
multiExpr = list(Train = list(data = datTrain), Test = list(data = datTest));
multiColor = list(Train = dynamicColors);

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );


ref = 1
test = 2

statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])

statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

temp = cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))

# save module preservation test results
write.csv(temp, file="./example/df_zsummary.csv")

print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

modColors = rownames(mp$preservation$observed[[ref]][[test]])

moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
text = modColors[plotMods];
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
mains = c("Preservation Median rank", "Preservation Zsummary");

sizeGrWindow(10, 5);
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  #labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=15, col = "darkgreen", lty = 2)
  }
}