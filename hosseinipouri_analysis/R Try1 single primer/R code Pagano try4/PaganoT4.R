# -----------------------------------------------------------------------------------------------------
# Dirk Dittmer on 12052012
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
library(lattice)
library(reshape)
library(reshape2)
library(MASS)
library(DAAG)
library(RColorBrewer)
library(psych)
library(FactoMineR)
library(heatmap.plus)
# library(ggplot2) function "qplot" in ggplot2 is the same as "qplot" in qvalue
library(beeswarm)
library(vioplot)
library(DMwR)
library(Hmisc)
library(vegan)
library(qvalue)


# -----------------------------------------------------------------------------------------------------
# STEP 0.b:    Define working directory
# -----------------------------------------------------------------------------------------------------
setwd("/Users/dirk/Desktop/try4")


# -----------------------------------------------------------------------------------------------------
# STEP 0.c:    Define functions
#
# -----------------------------------------------------------------------------------------------------
# read in multiple spreadsheets each with primers and one or more experiments
# the experiments are in  columns
#
read.rain <- function(x)
{
	filnamelistlength <- length(x)
	 
	mir <- read.delim(x[1], header = TRUE, sep = "\t")
	mir <- na.omit(mir)
	mir <- drop.levels(mir)
	mir$filename <- factor(x[1])
	rain <- mir
	print(nrow(rain))
	
	if (filnamelistlength< 2) return (rain)
	
	for (i in 2:filnamelistlength)
	{
	mir <- read.delim(x[i], header = TRUE, sep = "\t")
	mir <- na.omit(mir)
	mir <- drop.levels(mir)
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
filenames <- c("maribavir2012.txt")
rain <- read.rain(filenames)

summary(rain)
head(rain)
table(rain$filename)


# -----------------------------------------------------------------------------------------------------
# STEP 1.b.:	SAFE COPY OF DATASET
#
# -----------------------------------------------------------------------------------------------------
rainall <- rain
write.table(rainall,file = "CompleteDataSet_CT.txt", sep = "\t")

#  This will always re-create the initial data set
rain <- rainall


# -----------------------------------------------------------------------------------------------------
# STEP 1.c:    DATA CLEANING
#
#																				OUT:		rain
# -----------------------------------------------------------------------------------------------------
head(rain)
colnames(rain)

#drop out filename column
rain <- rain[,-ncol(rain)]  


# -----------------------------------------------------------------------------------------------------
# STEP 2:    DATA TRANSFORMATION 
#
# -----------------------------------------------------------------------------------------------------
# STEP 2.a: drop out known bad experiments
# -------------------------------------------------------
colnames(rain)[15]
rain <- rain[,-15]  


# -------------------------------------------------------
# STEP 2.b:    explore primer means  and  remove mean>35
rain$mean <- rowMeans(rain[,2:ncol(rain)])

# densityplot using the lattice function
trellis.par.set("fontsize", list(text=18))
trellis.par.set("plot.symbol", list(pch=16))
densityplot(~mean, data = rain, main = "primer mean CT distribution")

# -------------------------------------------------------
# ecdf distribution
# apply: 1 is for row
par(cex=1.8)
par(bty="l")
plot(ecdf(apply(rain[,2:(ncol(rain)-1)],1, mean)), col = "darkblue", verticals = TRUE, do.points = FALSE, lwd = 2, main = "mean CT", xlab = "mean CT/ primer", ylab = "Culm. density")    
histSpike(apply(rain[,2:(ncol(rain)-1)],1, mean), add=TRUE, col = "darkred", frac = 0.2, lwd = 4)
abline(v = 35, lty =1, col = "darkgreen")
text(35,0,".", col = "darkgreen")

# -------------------------------------------------------
# cut out genes that are not expressed mean > 45
subset(rain,rain$mean>35)
rain <- subset(rain,rain$mean<35)

# -------------------------------------------------------
# STEP 2.b:    explore primer means  and  SD
par(cex=1.5)
par(bty = "o")
plot( 
apply(rain[,2:(ncol(rain))],1, mean), 
apply(rain[,2:(ncol(rain))],1, sd),
pch = 1,
cex=2,
col = "darkblue",  
lwd = 2, 
main = "M/A plot", 
ylab = "SD CT / primer", 
xlab = "mean CT / primer",
xlim = c(20,35),
ylim = c(1,9)
)
abline(h = 2, lty =2, col = "lightgray")
abline(h = 4, lty =2, col = "lightgray")
abline(h = 6, lty =2, col = "lightgray")
abline(h = 8, lty =2, col = "lightgray")

rain <- rain[,-ncol(rain)]  #drop out mean and restore original data frame
summary(rain)


# -----------------------------------------------------------------------------------------------------
#  STEP:		convert data frame into an all numeric data frame
# -----------------------------------------------------------------------------------------------------
# delete primer column and label rownames with primer names
table(rain$primer)
rownames(rain) <- rain$primer
rain <- rain[,-1] 
head(rain)


# -----------------------------------------------------------------------------------------------------
# STEP 3:  dCT CALCULATION
#  
# -----------------------------------------------------------------------------------------------------
# evaluate housekeeping genes
house <- as.matrix(rain[c("ACTIN","HPRT","house","GAPDH","B-ACT"),])
house <- t(house)
head(house)

# VARIANT 1:		pairwise comparison using SPLOM
splom(house, col = "darkred", pch =19, pscales=2, type=c("g","p"), main = "housekeeping genes")

# -------------------------------------------------------
# VARIANT 2:		pairwise comparisons using PAIRS.PANELS
png(file = "Housekeeping.png", units = "in", width = 8, height = 8, res = 600)
{
pairs.panels (house, 
lm = TRUE, 
ellipses = FALSE, 
hist.col = "lightgray",
lty = 1, 
pch = 1,
lwd = 1, 
method = "pearson", 
cex = 2, 
col = "darkred", 
cex.labels = 2,font.labels = 2,
gap =1.2, 
scale = TRUE,
digits = 4,
main = "housekeeping genes")
}
dev.off()

# -------------------------------------------------------
# VARIANT 3:		heatmap of pariwaise correlation coefficients using COR.PLOT
#unsorted
house.cor <- cor(house,use="complete.obs",method="pearson")
cor.plot(house.cor, main = "Housekeeping genes", zlim = c(0.8,1), cex.axis = 1.0, n = 32)

# sorted
house.cor <- mat.sort(house.cor)
head(house.cor)
cor.plot(house.cor, main = "Housekeeping genes", zlim = c(0.9,1), cex.axis = 1.0, n = 32)


# --------------------------------------------------------------------------------------------------------
# dCT calculation  based on house (house was calcuated as geomean in excel)
# --------------------------------------------------------------------------------------------------------
for (i in 1:nrow(rain))
	{
	rain[i,] <- 		rain["house",] - rain[i,]
	}
summary(rain)

# -------------------------------------------------------
# write out dCT to file
write.table(rain,file = "dCT.txt", sep = "\t")


# -----------------------------------------------------------------------------------------------------
# STEP 4:  DELETE NON-CHANGED GENES based on m.a.d distribution
#  
# -----------------------------------------------------------------------------------------------------

# VARIANT 1:		subsetting based on mad, i.e. taking out non-expressed genes
png(file = "AlteredGenesMAD.png", units = "in", width = 8, height = 8, res = 600)
{
par(mfrow = c(1,2))
par(cex = 1.5)
plot(ecdf(apply(rain,1, mad)), col = "darkblue", verticals = TRUE, do.points = FALSE, lwd = 2, main = "initial", xlab = "CT")    
# apply: 1 is for row, mad is for median absolute deviation
histSpike(apply(rain,1, mad), add=TRUE, col = "darkred", frac = 0.2, lwd = 3)

# set mad > 0.1
rain <- subset(rain, apply(rain,1, mad) > 0.1)

plot(ecdf(apply(rain,1, mad)), col = "darkblue", verticals = TRUE, do.points = FALSE, lwd = 2, main = "changed", xlab = "CT")    
# apply: 1 is for row, mad is for median absolute deviation
histSpike(apply(rain,1, mad), add=TRUE, col = "darkred", frac = 0.2, lwd = 3)
nrow(rain)
}
dev.off()

# -----------------------------------------------------------------------------------------------------
# Plot the entire dCT distribution
# -----------------------------------------------------------------------------------------------------
rainCluster	<- melt(rain)
rainCluster <- na.omit(rainCluster)
rainCluster <- drop.levels(rainCluster)
rainCluster <- rainCluster[,c("variable","value")]
summary(rainCluster)
nrow(rainCluster)
head(rainCluster)

# -------------------------------------------------------
# Test for Normality
IQR(rainCluster$value, na.rm = TRUE)
fivenum(rainCluster$value, na.rm = TRUE)
shapiro.test(rainCluster$value)

# -------------------------------------------------------
# qqPLOT of entire dCT distribution
png(file = "qqPlotOfAlldCThouse.png", units = "in", width = 8, height = 8, res = 600)
{
par(mfrow = c(1,1))
par(cex = 1.5)
par(bty = "o")
qqnorm(rainCluster$value, 
pch = 19,
cex = 0.5, 
col = "darkblue",
main = "dCT distribution")
qqline(rainCluster$value, col = "darkred", lwd = 2, lty = 3)
histSpike(rainCluster$value, add=TRUE, col = "darkgray", frac = 0.3, lwd = 3, side =2)
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# STEP 5:  	HEATMAP
#					IN: rain
# -----------------------------------------------------------------------------------------------------
summary(rain)

# -------------------------------------------------------
# dataframe to matrix conversion
try2m <- as.matrix(rain)
nrow(try2m)
nrow(try2m[!complete.cases(try2m),])
head(try2m)

png(file = "heatmapAllScalebyRow.png", units = "in", width = 8, height = 16, res = 600)
{
rgb.palette <- colorRampPalette(c("white", "orange", "red"), space = "rgb", bias = 0.5)
heatmap.plus(try2m, na.rm = T, scale = "row", col = rgb.palette(25), hclustfun=function(m) hclust(m, method="ward"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,10), Colv = T, Rowv = T)
}
dev.off()

# -------------------------------------------------------
# variant color scheme
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "PuOr")))
heatmap.plus(try2m, na.rm = T, scale = "row", col = brewer.palette (15), hclustfun=function(m) hclust(m, method="ward"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,10), Colv = T, Rowv = T)


# -----------------------------------------------------------------------------------------------------
# STEP 6:  	CENTER by primer
#
# -----------------------------------------------------------------------------------------------------
experimentByprimer <- t(as.matrix(rain))
head(experimentByprimer)

#scale by primer
experimentByprimer  <- scale(experimentByprimer , center = TRUE, scale = TRUE)
head(experimentByprimer)

primer.cor <- cor(experimentByprimer, use="complete.obs",method="pearson")
head(primer.cor)
cor.plot(primer.cor, main = "all genes", zlim = c(-1,1), cex.axis = 0.5, n = 32)

# sorted
primer.cor <- mat.sort(primer.cor)
head(primer.cor)
cor.plot(primer.cor, main = "dCT", zlim = c(0,1), cex.axis = 0.8, n = 32)
primer.kmeans <- kmeans(primer.cor, nrow(primer.cor)/3, iter.max=1000, nstart=10000)
primer.kmeans.assignment <- as.data.frame(primer.kmeans$cluster)
primer.kmeans.assignment$primer <- rownames(primer.kmeans.assignment)
# primer.kmeans.assignment[,order(primer.kmeans.assignment$cluster)]

# -------------------------------------------------------
# kmeans plot
primer_dist <- dist(primer.cor)
cmd <- cmdscale(primer_dist)
groups <- levels(factor(primer.kmeans$cluster))
ordiplot(cmd)
cols <- c("steelblue", "darkred", "darkgreen", "pink","green","red","blue","yellow")
for(i in seq_along(groups)){
  points(cmd[factor(primer.kmeans$cluster) == groups[i], ], col = cols[i], pch = 16)
}
# add spider and hull
ordispider(cmd, factor(primer.kmeans$cluster), label = TRUE)
ordihull(cmd, factor(primer.kmeans$cluster), lty = "dotted")

# -------------------------------------------------------
# heatmap of scaled primers
heatmap.plus(experimentByprimer, na.rm = T, scale = "none", col = rgb.palette (25), hclustfun=function(m) hclust(m, method="ward"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,10), Colv = T, Rowv = T)


# -----------------------------------------------------------------------------------------------------
# STEP 7:  subset of only induced cells
#  					(take out t0   then take out ACV treatment)
# -----------------------------------------------------------------------------------------------------
summary(rain)
rownames(rain)
colnames(rain)
# take out t0
rain <- rain[,-c(1:6)]
# take out ACV
rain <- rain[,-c(3,6,8,11,14,17)]
summary(rain)
head(rain)

# -------------------------------------------------------
# dataframe to matrix conversion
try2m <- as.matrix(rain)
nrow(try2m)
nrow(try2m[!complete.cases(try2m),])
head(try2m)

# -----------------------------------------------------------------------------------------------------
# STEP 8:  PCA
#  
# -----------------------------------------------------------------------------------------------------
# STEP 8.b:  PCA  
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
# STEP 8.b:  PCA  & cluster tree

mir.hcpc <- HCPC(mir.pca)
mir.hcpc$desc.var 
mir.hcpc$desc.axe
mir.hcpc$desc.ind 

# -----------------------------------------------------------------------------------------------------
# STEP 8.c:  PCA on primer

for.pca <- t(for.pca)
for.pca <- unique(for.pca)
head(for.pca)

mir.pca <- PCA(for.pca, scale.unit = T, ncp = Inf, graph = F)

par(mfrow = c(2,2))
plot.PCA(mir.pca, cex=1.0, choix ="ind", axes = c(1,2), title = "Primers", col.ind = "darkblue")
plot.PCA(mir.pca, cex=1.0, choix ="ind", axes = c(2,3), title = "Primers", col.ind = "darkblue")
plot.PCA(mir.pca, cex=1.0, choix ="ind", axes = c(1,3), title = "Primers", col.ind = "darkblue")
dimdesc(mir.pca)
barplot(mir.pca$eig[,1], main = "Eigenvalues", names.arg = paste ("Dim", 1:nrow(mir.pca$eig), sep = ""), col = "darkblue")

mir.hcpc <- HCPC(mir.pca)
mir.hcpc$desc.var 
mir.hcpc$desc.axe
print(mir.hcpc$desc.ind)


# -----------------------------------------------------------------------------------------------------
# STEP:		SAFE COPY OF REDUCED DATASET
#
# -----------------------------------------------------------------------------------------------------
raindCT <- rain
write.table(raindCT,file = "dCT_reducedSetOfExperiments.txt", sep = "\t")
rain <- raindCT

# -----------------------------------------------------------------------------------------------------
# STEP:		heatmap of REDUCED DATASET
#
# -----------------------------------------------------------------------------------------------------
head(rain)
rownames(rain)
# take out the housekeeping genes
rain <- rain[-c(1:2,72,73),]

# -------------------------------------------------------
# dataframe to matrix conversion
try2m <- as.matrix(rain)
nrow(try2m)
nrow(try2m[!complete.cases(try2m),])
head(try2m)

png(file = "heatmapReducedExperiments.png", units = "in", width = 8, height = 16, res = 600)
{
rgb.palette <- colorRampPalette(c( "white","orange", "red" ), space = "rgb", bias = 0.5)
heatmap.plus(try2m, na.rm = T, scale = "none", col = rgb.palette (25), hclustfun=function(m) hclust(m, method="ward"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,10), Colv = T, Rowv = T)
}
dev.off()

# -------------------------------------------------------
# variant color scheme
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "PuOr")))
heatmap.plus(try2m, na.rm = T, scale = "none", col = brewer.palette (15), hclustfun=function(m) hclust(m, method="ward"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,10), Colv = T, Rowv = T)


# -----------------------------------------------------------------------------------------------------
# STEP:		Single primer statistics
#
# -----------------------------------------------------------------------------------------------------
# test
filenames <- c("factor.txt")
rain <- read.rain(filenames)
summary(rain)

#drop out filename column
rain <- rain[,-ncol(rain)]  

# -----------------------------------------------------------------------------------------------------
#  STEP:		convert data frame into an all numeric data frame
# delete primer column and label rownames with primer names
colnames(rain)
rownames(rain) <- rain[,1]
rain <- rain[,-1] 
head(rain)

# -------------------------------------------------------
# simple linear model
	a <- lm(PK ~ time + MBV + rain[,4], data = rain, na.action = na.omit)
	name <- (colnames(rain[4]))
	p <- anova(a)$Pr[3]
b <- as.data.frame(cbind(name, p))
b

# -------------------------------------------------------
# single variate analysis for each primer
for (i in 5:ncol(rain))
{
	a <- lm(PK ~ time + MBV +  rain[,i], data = rain, na.action = na.omit)
	p <- anova(a)$Pr[3]
	name <- (colnames(rain[i]))
	p <- as.numeric(p)
	b <- rbind(b, cbind(name, p))
}
summary(b)

b$p <- as.character(b$p)
b$p <- as.numeric(b$p)

# -------------------------------------------------------
# q-value adjustment
summary(b)
head(b)
rainq <- qvalue(b$p, lambda =0 , robust = FALSE, fdr.level = TRUE, pi0.method="bootstrap")
qsummary(rainq)
qwrite(rainq, filename = "adjustedp values.txt")
b$q <- rainq$qvalues
summary(b)
write.table(b,file = "qAndpValuesOfreducedSetOfExperiments.txt", sep = "\t",row.names = FALSE)

png(file = "qplotOfANOVAfit.png", units = "in", width = 8, height = 8, res = 600)
{
qplot(rainq, rng = c(0.0, 0.8) )
}
dev.off()

# -------------------------------------------------------
# Densityplot of s
png(file = "SingleANOVA.png", units = "in", width = 8, height = 8, res = 600)
{
par(mfrow = c(2,1))
par(cex = 1.2)
plot(density(log10(b$p)), main = "p-value", xlab = "lg(p-value)", col = "darkred", lwd = 2, cex = 2, pch =1)
abline(v = log10(0.001), lty = 2, col = "darkblue")

plot(density(log10(b$q)), main = "q-value", xlab = "lg(q-value)", col = "darkred", lwd = 2, cex = 2, pch =1)
abline(v = log10(0.01), lty = 2, col = "darkblue")
}
dev.off()


# -------------------------------------------------------
# looking at the best ones by general linear regression
b[log10(b$p) < -3,]

rain$PK <- as.factor(rain$PK)
fit <- glm(PK ~ time + MBV + rain[,"LMP.2A_59"], data = rain, family = binomial())
summary(fit)
exp(confint(fit))
predict(fit, type = "response")
cdplot(PK ~ rain[,"LMP.2A_59"],data = rain )

