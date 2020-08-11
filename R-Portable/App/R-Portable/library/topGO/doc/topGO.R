### R code from vignette source 'topGO.Rnw'

###################################################
### code chunk number 1: topGO.Rnw:40-43
###################################################
options(width = 95)
## we better keep the data in data frames as strings
options(stringsAsFactors = FALSE)


###################################################
### code chunk number 2: topGO.Rnw:67-73
###################################################
x <- topGO:::.algoComp[, -8]
x <- x[, colSums(x) > 0]
yesPic <- "\\includegraphics[width=3mm]{green_ckmark.png}"
noPic <- "\\includegraphics[width=3mm]{red_ckmark.png}"
x[x == 1] <-  yesPic
x[x == "0"] <-  noPic


###################################################
### code chunk number 3: topGO.Rnw:78-83
###################################################
if(require(xtable))
  print(xtable(x, align =  c("l", rep("c", ncol(x)))), 
        type = "latex", sanitize.text.function = function(x) return(x),
        floating = FALSE)



###################################################
### code chunk number 4: topGO.Rnw:112-115 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install()


###################################################
### code chunk number 5: topGO.Rnw:120-123 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("topGO")


###################################################
### code chunk number 6: topGO.Rnw:165-169
###################################################
library(topGO)
library(ALL)
data(ALL)
data(geneList)


###################################################
### code chunk number 7: topGO.Rnw:179-181
###################################################
affyLib <- paste(annotation(ALL), "db", sep = ".")
library(package = affyLib, character.only = TRUE)


###################################################
### code chunk number 8: topGO.Rnw:191-192
###################################################
sum(topDiffGenes(geneList))


###################################################
### code chunk number 9: topGO.Rnw:199-204
###################################################
sampleGOdata <- new("topGOdata", 
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)


###################################################
### code chunk number 10: topGO.Rnw:217-218
###################################################
sampleGOdata


###################################################
### code chunk number 11: topGO.Rnw:239-240
###################################################
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")


###################################################
### code chunk number 12: topGO.Rnw:245-246
###################################################
resultFisher


###################################################
### code chunk number 13: topGO.Rnw:252-254
###################################################
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")


###################################################
### code chunk number 14: topGO.Rnw:272-275
###################################################
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, 
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)


###################################################
### code chunk number 15: topGO.Rnw:285-287
###################################################
if(require(xtable))
  print(xtable(apply(allRes, 2, as.character)), floating = FALSE)


###################################################
### code chunk number 16: topGO.Rnw:302-306
###################################################
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}


###################################################
### code chunk number 17: topGO.Rnw:309-318 (eval = FALSE)
###################################################
## pValue.classic <- score(resultKS)
## pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
## 
## gstat <- termStat(sampleGOdata, names(pValue.classic))
## gSize <- gstat$Annotated / max(gstat$Annotated) * 4
## gCol <- colMap(gstat$Significant)
## 
## plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
##      pch = 19, cex = gSize, col = gCol)


###################################################
### code chunk number 18: topGO.Rnw:325-339
###################################################
pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]

gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)

par(mfcol = c(1, 2), cex = 1)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

plot(pValue.classic, pValue.elim, log = "xy", xlab = "log(p-value) classic", ylab = "log(p-value) elim",
     pch = 19, cex = gSize, col = gCol)



###################################################
### code chunk number 19: topGO.Rnw:356-360
###################################################
sel.go <- names(pValue.classic)[pValue.elim < pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
      elim = pValue.elim[sel.go],
      classic = pValue.classic[sel.go])


###################################################
### code chunk number 20: topGO.Rnw:376-377 (eval = FALSE)
###################################################
## showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')


###################################################
### code chunk number 21: topGO.Rnw:380-381
###################################################
printGraph(sampleGOdata, resultKS.elim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)


###################################################
### code chunk number 22: topGO.Rnw:413-416
###################################################
library(topGO)
library(ALL)
data(ALL)


###################################################
### code chunk number 23: topGO.Rnw:424-426
###################################################
BPterms <- ls(GOBPTerm)
head(BPterms)


###################################################
### code chunk number 24: topGO.Rnw:439-442
###################################################
library(genefilter)
selProbes <- genefilter(ALL, filterfun(pOverA(0.20, log2(100)), function(x) (IQR(x) > 0.25)))
eset <- ALL[selProbes, ]


###################################################
### code chunk number 25: topGO.Rnw:541-543
###################################################
geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
str(head(geneID2GO))


###################################################
### code chunk number 26: topGO.Rnw:565-567
###################################################
GO2geneID <- inverseList(geneID2GO)
str(head(GO2geneID))


###################################################
### code chunk number 27: topGO.Rnw:581-583
###################################################
geneNames <- names(geneID2GO)
head(geneNames)


###################################################
### code chunk number 28: topGO.Rnw:589-593
###################################################
myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)


###################################################
### code chunk number 29: topGO.Rnw:607-609
###################################################
GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO)


###################################################
### code chunk number 30: topGO.Rnw:623-624
###################################################
GOdata


###################################################
### code chunk number 31: topGO.Rnw:647-649
###################################################
y <- as.integer(sapply(eset$BT, function(x) return(substr(x, 1, 1) == 'T')))
table(y)


###################################################
### code chunk number 32: topGO.Rnw:658-659
###################################################
geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")


###################################################
### code chunk number 33: topGO.Rnw:672-677
###################################################
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)
}
x <- topDiffGenes(geneList)
sum(x) ## the number of selected genes


###################################################
### code chunk number 34: topGO.Rnw:683-691
###################################################
GOdata <- new("topGOdata", 
              description = "GO analysis of ALL data; B-cell vs T-cell",
              ontology = "BP",
              allGenes = geneList,
              geneSel = topDiffGenes,
              annot = annFUN.db,
              nodeSize = 5,
              affyLib = affyLib)


###################################################
### code chunk number 35: topGO.Rnw:723-731
###################################################
allProb <- featureNames(ALL)
groupProb <- integer(length(allProb)) + 1
groupProb[allProb %in% genes(GOdata)] <- 0
groupProb[!selProbes] <- 2
groupProb <- factor(groupProb, labels = c("Used", "Not annotated", "Filtered"))

tt <- table(groupProb)
tt


###################################################
### code chunk number 36: topGO.Rnw:753-757 (eval = FALSE)
###################################################
## pValue <- getPvalues(exprs(ALL), classlabel = y, alternative = "greater")
## geneVar <- apply(exprs(ALL), 1, var)
## dd <- data.frame(x = geneVar[allProb], y = log10(pValue[allProb]), groups = groupProb)
## xyplot(y ~ x | groups, data = dd, groups = groups)


###################################################
### code chunk number 37: topGO.Rnw:761-784
###################################################
pValue <- getPvalues(exprs(ALL), classlabel = y, alternative = "greater")
geneVar <- apply(exprs(ALL), 1, var)
dd <- data.frame(x = geneVar[allProb], y = log10(pValue[allProb]), groups = groupProb)

library(lattice)
trellis.device(device = pdf, theme = col.whitebg(), file = "whichProbe.pdf", width = 9, height = 7)
legendLab <- paste(names(table(groupProb)), " (#", table(groupProb), ")", sep = "")
pP <- xyplot(y ~ x | groups, data = dd, groups = groups,
             xlab = "Variance", ylab = "Log of p-values",
             layout = c(2, 2),
             key = list(text = list(lab = legendLab),
               points = list(pch = 20, cex = 2, 
                 col = Rows(trellis.par.get("superpose.symbol"), 1:3)$col),
               size = 7, padding.text = 3,
               x = .65, y = .7, corner = c(0, 0), border = TRUE, cex = 1),
             panel = function(x, y, ...) {
               selY <- y <= -2
               panel.xyplot(x[selY], y[selY], pch = 2, ...)
               panel.xyplot(x[!selY], y[!selY], pch = 20, ...)
               panel.abline(h = -2, lty = 2, col = "black")
             })
print(pP)
dev.off()


###################################################
### code chunk number 38: topGO.Rnw:807-810
###################################################
description(GOdata)
description(GOdata) <- paste(description(GOdata), "Object modified on:", format(Sys.time(), "%d %b %Y"), sep = " ")
description(GOdata)


###################################################
### code chunk number 39: topGO.Rnw:816-819
###################################################
a <- genes(GOdata) ## obtain the list of genes
head(a)
numGenes(GOdata)


###################################################
### code chunk number 40: topGO.Rnw:826-829
###################################################
selGenes <- sample(a, 10)
gs <- geneScore(GOdata, whichGenes = selGenes) 
print(gs)


###################################################
### code chunk number 41: topGO.Rnw:833-838
###################################################
gs <- geneScore(GOdata, whichGenes = selGenes, use.names = FALSE)
print(gs)

gs <- geneScore(GOdata, use.names = FALSE)
str(gs)


###################################################
### code chunk number 42: topGO.Rnw:842-845
###################################################
sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)


###################################################
### code chunk number 43: topGO.Rnw:852-856
###################################################
.geneList <- geneScore(GOdata, use.names = TRUE)
GOdata ## more available genes
GOdata <- updateGenes(GOdata, .geneList, topDiffGenes)
GOdata ## the available genes are now the feasible genes


###################################################
### code chunk number 44: topGO.Rnw:863-867
###################################################
graph(GOdata) ## returns the GO graph

ug <- usedGO(GOdata)
head(ug) 


###################################################
### code chunk number 45: topGO.Rnw:872-879
###################################################
sel.terms <- sample(usedGO(GOdata), 10)

num.ann.genes <- countGenesInTerm(GOdata, sel.terms) ## the number of annotated genes
num.ann.genes

ann.genes <- genesInTerm(GOdata, sel.terms) ## get the annotations
head(ann.genes)


###################################################
### code chunk number 46: topGO.Rnw:884-888
###################################################
ann.score <- scoresInTerm(GOdata, sel.terms)
head(ann.score)
ann.score <- scoresInTerm(GOdata, sel.terms, use.names = TRUE)
head(ann.score)


###################################################
### code chunk number 47: topGO.Rnw:894-895
###################################################
termStat(GOdata, sel.terms)


###################################################
### code chunk number 48: topGO.Rnw:975-979
###################################################
goID <- "GO:0044255"
gene.universe <- genes(GOdata)
go.genes <- genesInTerm(GOdata, goID)[[1]]
sig.genes <- sigGenes(GOdata)


###################################################
### code chunk number 49: topGO.Rnw:985-991
###################################################
my.group <- new("classicCount", testStatistic = GOFisherTest, name = "fisher",
                allMembers = gene.universe, groupMembers = go.genes,
                sigMembers = sig.genes)

contTable(my.group)
runTest(my.group)


###################################################
### code chunk number 50: topGO.Rnw:1012-1020
###################################################
set.seed(123)
elim.genes <- sample(go.genes, length(go.genes) / 4)
elim.group <- new("elimCount", testStatistic = GOFisherTest, name = "fisher",
                allMembers = gene.universe, groupMembers = go.genes,
                sigMembers = sig.genes, elim = elim.genes)

contTable(elim.group)
runTest(elim.group)


###################################################
### code chunk number 51: topGO.Rnw:1052-1054
###################################################
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
resultFisher <- getSigGroups(GOdata, test.stat)


###################################################
### code chunk number 52: topGO.Rnw:1060-1061
###################################################
resultFisher


###################################################
### code chunk number 53: topGO.Rnw:1069-1071
###################################################
test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
resultKS <- getSigGroups(GOdata, test.stat)


###################################################
### code chunk number 54: topGO.Rnw:1084-1086 (eval = FALSE)
###################################################
## test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
## resultElim <- getSigGroups(GOdata, test.stat)


###################################################
### code chunk number 55: topGO.Rnw:1091-1093
###################################################
test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
resultWeight <- getSigGroups(GOdata, test.stat)


###################################################
### code chunk number 56: topGO.Rnw:1151-1152
###################################################
resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")


###################################################
### code chunk number 57: topGO.Rnw:1160-1165 (eval = FALSE)
###################################################
## weight01.fisher <- runTest(GOdata, statistic = "fisher")
## weight01.t <- runTest(GOdata, algorithm = "weight01", statistic = "t")
## elim.ks <- runTest(GOdata, algorithm = "elim", statistic = "ks")
## 
## weight.ks <- runTest(GOdata, algorithm = "weight", statistic = "ks") #will not work!!!


###################################################
### code chunk number 58: topGO.Rnw:1172-1174
###################################################
whichTests()
whichAlgorithms()


###################################################
### code chunk number 59: topGO.Rnw:1208-1211
###################################################
pvalFis <- score(resultFis)
head(pvalFis)
hist(pvalFis, 50, xlab = "p-values")


###################################################
### code chunk number 60: topGO.Rnw:1218-1219
###################################################
head(score(resultWeight))


###################################################
### code chunk number 61: topGO.Rnw:1228-1231
###################################################
pvalWeight <- score(resultWeight, whichGO = names(pvalFis))
head(pvalWeight)
cor(pvalFis, pvalWeight)


###################################################
### code chunk number 62: topGO.Rnw:1238-1239
###################################################
geneData(resultWeight)


###################################################
### code chunk number 63: topGO.Rnw:1251-1253
###################################################
allRes <- GenTable(GOdata, classic = resultFis, KS = resultKS, weight = resultWeight,
                   orderBy = "weight", ranksOf = "classic", topNodes = 20)


###################################################
### code chunk number 64: topGO.Rnw:1263-1265
###################################################
if(require(xtable))
  print(xtable(apply(allRes, 2, as.character)), floating = FALSE)


###################################################
### code chunk number 65: topGO.Rnw:1286-1288
###################################################
goID <- allRes[1, "GO.ID"]
print(showGroupDensity(GOdata, goID, ranks = TRUE))


###################################################
### code chunk number 66: topGO.Rnw:1311-1313
###################################################
goID <- allRes[10, "GO.ID"]
gt <- printGenes(GOdata, whichTerms = goID, chip = affyLib, numChar = 40)


###################################################
### code chunk number 67: topGO.Rnw:1318-1320
###################################################
if(require(xtable))
  print(xtable(gt), floating = FALSE)


###################################################
### code chunk number 68: topGO.Rnw:1344-1346 (eval = FALSE)
###################################################
## showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
## showSigOfNodes(GOdata, score(resultWeight), firstSigNodes = 5, useInfo = 'def')


###################################################
### code chunk number 69: topGO.Rnw:1349-1351
###################################################
printGraph(GOdata, resultFis, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)
printGraph(GOdata, resultWeight, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "def", pdfSW = TRUE)


###################################################
### code chunk number 70: topGO.Rnw:1381-1383 (eval = FALSE)
###################################################
## printGraph(GOdata, resultWeight, firstSigNodes = 10, resultFis, fn.prefix = "tGO", useInfo = "def")
## printGraph(GOdata, resultElim, firstSigNodes = 15, resultFis, fn.prefix = "tGO", useInfo = "all")


###################################################
### code chunk number 71: topGO.Rnw:1392-1393
###################################################
toLatex(sessionInfo())


