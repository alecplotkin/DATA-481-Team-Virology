# -----------------------------------------------------------------------------------------------------
# Dirk Dittmer on 04082014
# -----------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------
# STEP 0.a:    Initial install of the libraries from bioconductor
# ----------------------------------------------------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")
# qvalue also need TCL/TK installed from the tool section of R CRAN

# -----------------------------------------------------------------------------------------------------
# STEP 0.a:    Read in the libraries.
# -----------------------------------------------------------------------------------------------------
library(gdata)
# library(lattice)
# library(reshape)
library(reshape2)
library(MASS)
# library(DAAG)
library(RColorBrewer)
library(psych)
library(FactoMineR)
library(heatmap.plus)
# library(ggplot2) #function "qplot" in ggplot2 is the same as "qplot" in qvalue
# library(beeswarm)
# library(vioplot)
library(DMwR)
library(Hmisc)
library(vegan)
# library(qvalue)
library(xlsx)


# -----------------------------------------------------------------------------------------------------
# STEP 0.b:    Define digits
# -----------------------------------------------------------------------------------------------------
options(digits=3)


# -----------------------------------------------------------------------------------------------------
# STEP 0.c:    Define working directory
# -----------------------------------------------------------------------------------------------------
 setwd("/Users/dirk/DittmerSync/paper 2014 array Malawi Array/R Try1 heatmap n16")


# -----------------------------------------------------------------------------------------------------
# Define functions
#
# -----------------------------------------------------------------------------------------------------
# FUNCTION read.rain
# read in multiple excel spreadsheets each with primers and one or more experiments
# the experiments are in  columns
#
read.rain <- function(x)
{
	filnamelistlength <- length(x)
	 
	mir <- read.delim(x[1], header = TRUE, sep = "\t")
	# mir <- na.omit(mir)
	# mir <- drop.levels(mir)
	mir$filename <- factor(x[1])
	rain <- mir
	print(nrow(rain))
	
	if (filnamelistlength< 2) return (rain)
	
	for (i in 2:filnamelistlength)
	{
	mir <- read.delim(x[i], header = TRUE, sep = "\t")
	# mir <- na.omit(mir)
	# mir <- drop.levels(mir)
	mir$filename <- factor(x[i])
	
	rain <- rbind(rain,mir)
	print(nrow(rain))
	}	
	return (rain)
}


# -----------------------------------------------------------------------------------------------------
# STEP 1.a: 	READ IN THE FILES 
#
# -----------------------------------------------------------------------------------------------------
filenames <- c("dCT.txt")
rain <- read.rain(filenames)

summary(rain)
head(rain)
table(rain$filename)


# -----------------------------------------------------------------------------------------------------
# STEP 1.b.:	SAFE COPY OF DATASET
#
# -----------------------------------------------------------------------------------------------------
rainall <- rain
write.table(rainall,file = "b_CompleteDataSet_dCT.txt", sep = "\t", row.names = FALSE)

#  This will always re-create the initial data set
rain <- rainall


# -----------------------------------------------------------------------------------------------------
# STEP 1:    DATA CLEANING
#
#																				OUT:		rain
# -----------------------------------------------------------------------------------------------------
head(rain)
colnames(rain)

#drop out filename column
rain <- rain[,-ncol(rain)]  


# -------------------------------------------------------
# STEP 1:    format variables
# -------------------------------------------------------
summary(rain[,1:6])
rain$samples <- factor(rain$samples)
rain$name <- factor(rain$name)
rownames(rain) <- rain$name


# -----------------------------------------------------------------------------------------------------
# STEP:    DATA TRANSFORMATION 
#
# -----------------------------------------------------------------------------------------------------
# STEP: drop out known bad experiments and rows
# -------------------------------------------------------
colnames(rain)[1:2]
rain <- rain[,-1]  
rain <- rain[,-1]  

# find NA
nrow(rain)
nrow(is.na(rain))
nrow(rain) == nrow(is.na(rain))

# find complete cases
nrow(rain)
length(complete.cases(rain))
nrow(rain) == length(complete.cases(rain))

# find duplicate cases 
nrow(rain)
rain[duplicated(rain[,"name"]),]


# convert to data frame with samples as columns
summary(rain)
rain.wide <- na.omit(rain)
rain.wide <- as.matrix(rain.wide[,2:ncol(rain.wide)])
rain.wide <- t(rain.wide)
rain.wide <- as.data.frame(rain.wide)

head(rain.wide)
summary(rain.wide)


# -----------------------------------------------------------------------------------------------------
# SCREEN densityplot of real primers and NTC
#
colnames(rain.wide[1:5])
par(mfrow = c(1,1))
par(bty = "n")
par(cex=1.8)
# samples
mean <- rowMeans(rain.wide, na.rm = TRUE)
plot(density(mean), main = "", col = "red", lwd = 2)


# -----------------------------------------------------------------------------------------------------
# cut out primers with rowMean > 0
# actin and b2m
#
rain.wide[rowMeans(rain.wide, na.rm = TRUE)>0,]
rain.wide <- rain.wide[rowMeans(rain.wide, na.rm = TRUE)<0,]
rain.wide <- drop.levels(rain.wide)


# -----------------------------------------------------------------------------------------------------
# cut out experiment 4428 because of low RNA (see prior figure)
#
colnames(rain.wide)
rain.wide <- rain.wide[,-12]
summary(rain.wide)


# -------------------------------------------------------
# PANEL A: ECDF
#
png(file = "b_PrimerPerformanceNohouse.png", units = "in", width = 10, height = 6, res = 300)
{
par(mfrow = c(1,2))
par(cex=1.7)
par(bty="l")
plot(ecdf(apply( rain.wide ,1, FUN = mean, na.rm = TRUE)), col = "darkblue", verticals = TRUE, do.points = FALSE, lwd = 2, main = "a.", xlab = "mean dCT' / primer", ylab = "Culm. density", las = 1)    
histSpike(apply(rain.wide,1, FUN = mean, na.rm = TRUE), add=TRUE, col = "darkred", frac = 0.2, lwd = 4)
abline(h = .5, lty =3, col = "darkgray")

# -------------------------------------------------------
# PANEL B:    explore primer means  and  SD
#
par(cex=1.7)
par(bty = "l")
plot( 
apply(rain.wide,1, mean, na.rm = TRUE), 
apply(rain.wide,1, sd, na.rm = TRUE),
pch = 16,
cex=1,
# ylim = c(0,8),
# xlim = c(10, 40),
col = "darkred",  
lwd = 2, 
las = 1,
main = "b.", 
ylab = "sd dCT' / primer", 
xlab = "mean dCT' / primer")

abline(h = 1, lty =2, col = "lightgray")
abline(h = 2, lty =2, col = "lightgray")
abline(h = 3, lty =2, col = "lightgray")
abline(h = 4, lty =2, col = "lightgray")
abline(h = 5, lty =2, col = "lightgray")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
#  STEP:		convert data frame into an all numeric data frame
#
# -----------------------------------------------------------------------------------------------------
# delete primer column and label rownames with primer names
# rain <- rain[,8:(ncol(rain))] 
head(rain.wide)
colnames(rain.wide)
GLOBAL.allnumericStart = 1

rain.wide <- as.matrix(rain.wide[,GLOBAL.allnumericStart:ncol(rain.wide)])
rain.wide <- t(rain.wide)
head(rain.wide)


# NOTE THIS REASSIGN RAIN GLOBALLY
rain <- as.data.frame(rain.wide)


# -------------------------------------------------------
# write out dCT to file
#
write.table(rain,file = "b_dCT'.txt", sep = "\t",row.names = FALSE)



# -----------------------------------------------------------------------------------------------------
# FIGURE:		Comparison of sample sets: Positive, NTC, samples
#
rownames(rain)
png(file = "b_AlteredGenesBasedonSDofdCT'.png", units = "in", width = 8, height = 6, res = 300)
{
par(mfrow = c(1,1))
par(cex = 1.5)
par(bty = "l")
# positive
plot(ecdf(apply(rain[,GLOBAL.allnumericStart:ncol(rain)],2, sd)), col = "darkred", verticals = TRUE, do.points = FALSE, lwd = 2, main = "n = 15 samples", xlab = "sd of dCT'", lty =1, cex = 2, las = 1)    
# apply: 2 is for column, mad is for median absolute deviation
histSpike(apply(rain[,GLOBAL.allnumericStart:ncol(rain)],2, sd), add=TRUE, col = "darkgray", frac = 0.3, lwd = 5)
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# STEP: Evaluate the entire dCT distribution of the samples
# -----------------------------------------------------------------------------------------------------
colnames(rain)
rownames(rain)
head(rain)
rainCluster	<- melt(rain)
rainCluster <- na.omit(rainCluster)
rainCluster <- drop.levels(rainCluster)
colnames(rainCluster) <- c("name","value")
summary(rainCluster)
nrow(rainCluster)
head(rainCluster)


# -------------------------------------------------------
# Test for Normality
IQR(rainCluster$value, na.rm = TRUE)
fivenum(rainCluster$value, na.rm = TRUE)
shapiro.test(rainCluster$value)


# -----------------------------------------------------------------------------------------------------
# FIGURE: qqPLOT of entire dCT distribution
#
png(file = "b_qqPlot Of All dCT'.png", units = "in", width = 8, height = 8, res = 300)
{
par(mfrow = c(1,1))
par(cex = 1.5)
par(bty = "o")
qqnorm(rainCluster[,"value"], 
pch = 19,
cex = 0.5, 
las = 1,
col = "darkred",
main = "Samples")
qqline(rainCluster[,"value"], col = "darkblue", lwd = 2, lty = 3)
histSpike(rainCluster[,"value"], add=TRUE, col = "darkgray", frac = 0.5, lwd = 5, side =2)
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# FIGURE:  	HEATMAP
#					IN: rain.wide
# -----------------------------------------------------------------------------------------------------
summary(rain[,1:4])
GLOBAL.allnumericStart
rownames(rain)

# cutting out houskeeping mRNAs
colnames(rain)
rain[,76:78]
rain <- rain[,-c(76:78)]


# cutting out the NTC and pos component
try2m <- as.matrix(rain) 
nrow(try2m)
nrow(try2m[!complete.cases(try2m),])
rownames(try2m)
head(try2m)
nrow(try2m)


# -----------------------------------------------------------------------------------------------------
# FIGURE:   variant #1 "orange"
png(file = "b_heatmap dCT' orange.png", units = "in", width = 21, height = 7, res = 300)
{
rgb.palette <- colorRampPalette(c("white", "orange", "red"), space = "rgb", bias = 0.8)
heatmap.plus(try2m, na.rm = T, scale = "none", col = rgb.palette(25), hclustfun=function(m) hclust(m, method="ward"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,15), Colv = T, Rowv = T)
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# FIGURE:	variant #2 "color scheme"
#
png(file = "b_heatmap dCT' alternativeColor.png", units = "in", width = 21, height = 7, res = 300)
{
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")),space = "rgb", bias = 1.0, interpolate = "linear")
heatmap.plus(try2m, na.rm = T, scale = "none", col = brewer.palette (25), hclustfun=function(m) hclust(m, method="ward"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,15), Colv = T, Rowv = T, main = "dCT", xlab = "CDS", ylab = "samples")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# FIGURE:	variant #3 "color scheme" scaled by row=sample
#
png(file = "heatmap dCT'alternativeColor scaled by sample.png", units = "in", width = 21, height = 7, res = 300)
{
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")),space = "rgb", bias = 1, interpolate = "linear")
heatmap.plus(try2m, na.rm = T, scale = "row", col = brewer.palette (25), 
hclustfun=function(m) hclust(m, method="ward"), 
distfun = function(x) dist(x, method = "manhattan"), 
margins = c(12,24), 
Colv = T, Rowv = T, 
cexRow = 1, cexCol = 0.5,
main = "dCT' scaled by sample", xlab = "CDS", ylab = "samples")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# Now setting aside the dCT for other calculations
#
GLOBAL.tryDCT <- try2m
summary(GLOBAL.tryDCT)
head(GLOBAL.tryDCT)


# -----------------------------------------------------------------------------------------------------
# STEP:  	CENTER by primer
#
# -----------------------------------------------------------------------------------------------------
#experimentByprimer <- t(as.matrix(rain))
experimentByprimer <- try2m
head(experimentByprimer)

# center by primer
experimentByprimer  <- scale(experimentByprimer , center = TRUE, scale = FALSE)
head(experimentByprimer)

primer.cor <- cor(experimentByprimer, use="complete.obs",method="pearson")
head(primer.cor)
cor.plot(primer.cor, main = "all genes", zlim = c(-1,1), cex.axis = 0.8, n = 12)

# sorted
primer.cor <- mat.sort(primer.cor)
head(primer.cor)
cor.plot(primer.cor, main = "dCT'", zlim = c(-1,1), cex.axis = 0.8, n = 8)


# -----------------------------------------------------------------------------------------------------
# heatmap of centerd but not scaled primers
# 
head(experimentByprimer[,1:4])
png(file = "b_heatmap dCT' centered By Primer.png", units = "in", width = 21, height = 7, res = 300)
{
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")),space = "rgb", bias = .7, interpolate = "linear")
heatmap.plus(experimentByprimer, na.rm = T, scale = "none", col = brewer.palette (25),
hclustfun=function(m) hclust(m, method="ward"), 
distfun = function(x) dist(x, method = "manhattan"), 
margins = c(12,24), 
Colv = T, Rowv = T, 
cexRow = 1, cexCol = 0.5,
main = "dCT' centered by primer", xlab = "primer", ylab = "samples")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# heatmap of centerd primers PuOr
# 
head(experimentByprimer)
png(file = "b_heatmap dCT' centered By Primer purple.png", units = "in", width = 21, height = 7, res = 300)
{
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "PuOr")),space = "rgb", bias = 0.8, interpolate = "spline")
heatmap.plus(experimentByprimer, na.rm = T, scale = "none", col = brewer.palette (25),
hclustfun=function(m) hclust(m, method="ward"), 
distfun = function(x) dist(x, method = "manhattan"), 
margins = c(12,24), 
Colv = T, Rowv = T, 
cexRow = 1, cexCol = 0.5,
main = "dCT centered by primer", xlab = "primer", ylab = "samples")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# STEP:  PCA
#  
# -----------------------------------------------------------------------------------------------------
#  
head(try2m)
for.pca <- t(try2m)
for.pca <- unique(for.pca)
colnames(for.pca) <- rownames(try2m)
rownames(for.pca) <- colnames(try2m)
head(for.pca)

mir.pca <- PCA(for.pca, scale.unit = T, ncp = 5, graph = F)
par(mfrow = c(2,2))
plot.PCA(mir.pca, cex=1.0, choix ="ind", axes = c(1,2), title = "Samples", col.ind = "darkred")
plot.PCA(mir.pca, cex=1.0, choix ="ind", axes = c(2,3), title = "Samples", col.ind = "darkred")
plot.PCA(mir.pca, cex=1.0, choix ="ind", axes = c(1,3), title = "Samples", col.ind = "darkred")

dimdesc(mir.pca)
barplot(mir.pca$eig[,1], main = "Eigenvalues", names.arg = paste ("Dim", 1:nrow(mir.pca$eig), sep = ""))


# -----------------------------------------------------------------------------------------------------
# STEP:  PCA  & cluster tree
#
mir.hcpc <- HCPC(mir.pca)
mir.hcpc$desc.var 
mir.hcpc$desc.axe
mir.hcpc$desc.ind 


# -----------------------------------------------------------------------------------------------------
# STEP:  PCA on samples
#
head(try2m)
for.pca <- try2m
for.pca <- unique(for.pca)
head(for.pca)

mir.pca <- PCA(for.pca, scale.unit = F, ncp = Inf, graph = F)

png(file = "b_PCA of samples.png", units = "in", width = 8, height = 8, res = 300)
{
par(mfrow = c(2,2))
par(cex = 1.2)
plot.PCA(mir.pca, cex=0.5, choix ="ind", axes = c(1,2), title = "Samples", col.ind = "darkblue")
plot.PCA(mir.pca, cex=0.5, choix ="ind", axes = c(2,3), title = "Samples", col.ind = "darkblue")
plot.PCA(mir.pca, cex=0.5, choix ="ind", axes = c(1,3), title = "Samples", col.ind = "darkblue")
dimdesc(mir.pca)
barplot(mir.pca$eig[,1], main = "Eigenvalues", names.arg = paste ("Dim", 1:nrow(mir.pca$eig), sep = ""), col = "darkgray", las = 2, cex.names = 0.5)
}
dev.off()

mir.hcpc <- HCPC(mir.pca)
mir.hcpc$desc.var 
mir.hcpc$desc.axe
print(mir.hcpc$desc.ind)


# -----------------------------------------------------------------------------------------------------
# STEP:  ddCT waterfallplot of ratio using PLOT
#
head(GLOBAL.tryDCT)

waterfall <- t(GLOBAL.tryDCT)
waterfall <- data.frame(waterfall)
head(waterfall)

# normalizing to control by group
control <- waterfall[,"X8239"] 
waterfall <- sweep(waterfall,1,control,FUN = "-")

# order by del, which is in colum 1
waterfall <- waterfall[order(waterfall[,"X8239"], decreasing = F),]

colnames(waterfall)
par(mfrow = c(1,1))
par(cex=0.8)
par(bty="n")
barplot(waterfall[,1], 
main = "Expression", 
pch = 19, 
col = "black", 
# type = "s", 
names.arg = rownames(waterfall), 
ylab = "ddCT' ", 
xaxs = "r",yaxs = "r", 
xaxt ="s",
las = 3
# cex = 0.2, lwd = 2,  
#ylim = c(0,5)
)



abline(h = median(waterfall[,1]), lwd = 1, col = "gray", lty = 2)
abline(h = (median(waterfall[,1])-sd(waterfall[,1])), lwd = 1, col = "gray", lty = 2)
abline(h = (median(waterfall[,1])+sd(waterfall[,1])), lwd = 1, col = "gray", lty = 2)



head(waterfall)
tail(waterfall)
write.table(waterfall,file = "dCTordered.txt", sep = "\t")