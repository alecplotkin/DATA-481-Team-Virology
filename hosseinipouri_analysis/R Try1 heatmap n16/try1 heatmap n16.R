# -----------------------------------------------------------------------------------------------------
# Dirk Dittmer on 05032014
# -----------------------------------------------------------------------------------------------------
# R version 3.1.0 (2014-04-10) -- "Spring Dance"
# [R.app GUI 1.63 (6734) x86_64-apple-darwin13.1.0]

# -----------------------------------------------------------------------------------------------------
# STEP:    Initial install of the libraries from bioconductor 
# 
# source("http://bioconductor.org/biocLite.R")
# biocLite("qvalue")
# qvalue under less than R 3.1.0 and less than OsX Maverick
# also need TCL/TK installed from the tool section of R CRAN

# -----------------------------------------------------------------------------------------------------
# STEP:    Read in the R libraries.
#
library(gdata)
library(reshape2)
library(MASS)
library(RColorBrewer)
library(psych)
library(FactoMineR)
library(heatmap.plus)
# library(beeswarm)
# library(vioplot)
library(DMwR)
library(Hmisc)
library(vegan)
# library(qvalue)
library(xlsx)


# -----------------------------------------------------------------------------------------------------
# STEP:    Define digits
#
options(digits=3)


# -----------------------------------------------------------------------------------------------------
# STEP: Define working directory
# 
setwd("/Users/dirk/DittmerSync/paper 2014 array Malawi Array/R Try1 heatmap n16")


# -----------------------------------------------------------------------------------------------------
# STEP: Define functions
#

# -----------------------------------------------------------------------------------------------------
# FUNCTION read.exceldata
# read in multiple excel spreadsheets each with primers and one or more experiments
# the experiments are in  columns
#
read.exceldata <- function(x)
{
	filnamelistlength <- length(x)
	 
	mir <- read.xlsx(x[1], sheetIndex=1)
	# mir <- na.omit(mir)
	# mir <- drop.levels(mir)
	mir$filename <- factor(x[1])
	rain <- mir
	print(nrow(rain))
	
	if (filnamelistlength< 2) return (rain)
	
	for (i in 2:filnamelistlength)
	{
	mir <- read.xlsx(x[i], sheetIndex=1)
	# mir <- na.omit(mir)
	# mir <- drop.levels(mir)
	mir$filename <- factor(x[i])
	
	rain <- rbind(rain,mir)
	print(nrow(rain))
	}	
	return (rain)
}


# -----------------------------------------------------------------------------------------------------
# STEP: 	Read in the files
#
filenames <- c("a.xlsx","b.xlsx")
rain <- read.exceldata(filenames)

summary(rain)
head(rain)
table(rain$filename)


# -----------------------------------------------------------------------------------------------------
# STEP:	 Safe a copy of the data set
#
rainall <- rain
write.table(rainall,file = "CompleteDataSet_CT.txt", sep = "\t", row.names = FALSE)

#  This will always re-create the initial data set
rain <- rainall


# -----------------------------------------------------------------------------------------------------
# STEP:    DATA CLEANING
#
head(rain)
colnames(rain)

colnames(rain) <- c("Primer","629","731","1314","1700","2572","2741","2858","2910","2969","3544","3937","4428","4678","5485","7121","8239","filename" )

#drop out filename column
rain <- rain[,-ncol(rain)]  


# -----------------------------------------------------------------------------------------------------
# STEP:    format variables
#
summary(rain)
#rain$Col <- factor(rain$Col)


# -----------------------------------------------------------------------------------------------------
# STEP:    DATA TRANSFORMATION 
# -----------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------
# STEP: drop out known bad experiments and rows
#
# colnames(rain)[15]
# rain <- rain[,-15]  

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
rain[duplicated(rain[,"Primer"]),]
rain <- rain[-144,]

summary(rain)
rain.wide <- na.omit(rain)
rain.wide <- drop.levels(rain)

# -----------------------------------------------------------------------------------------------------
# STEP: Generate unique Row names
#
#crain.wide$uniqueRowName <- interaction(rain.wide$Primer, drop.levels = TRUE)
#vtable(rain.wide$uniqueRowName)
rownames(rain.wide) <- rain.wide$Primer
rownames(rain.wide)
# rain.wide <- rain.wide[,-ncol(rain.wide)]

head(rain.wide)
summary(rain.wide)


# -----------------------------------------------------------------------------------------------------
# SCREEN: Densityplot of real primers and NTC
#
colnames(rain.wide)
par(mfrow = c(1,1))
par(bty = "n")
par(cex=1.8)
# samples
mean <- rowMeans(rain.wide[,2:ncol(rain.wide)], na.rm = TRUE)
plot(density(mean), main = "", col = "red", lwd = 2)
par(new = FALSE)


# -----------------------------------------------------------------------------------------------------
# STEP Cut out genes that are not expressed at all, i.e. mean > 45
# apply: 2 is for column
# this should not be used if a single sample may be interesting
#
colnames(rain.wide)
rain.wide$mean <- rowMeans(rain.wide[,2:ncol(rain.wide)])
rain.wide[rain.wide$mean>=45,]
rain.wide <- subset(rain.wide,rain.wide$mean<=45)
rain.wide <- rain.wide[,-ncol(rain.wide)]
colnames(rain.wide)


# -----------------------------------------------------------------------------------------------------
# FIGURE: ECDF
#
png(file = "PrimerPerformance.png", units = "in", width = 10, height = 6, res = 300)
{
par(mfrow = c(1,2))
par(cex=1.7)
par(bty="l")
plot(ecdf(apply( rain.wide[,c(2:2:ncol(rain.wide))] ,1, FUN = mean, na.rm = TRUE)), col = "darkblue", verticals = TRUE, do.points = FALSE, lwd = 2, main = "a.", xlab = "mean CT' / primer", ylab = "Culm. density", las = 1, xlim = c(20, 50))    
histSpike(apply(rain.wide[,c(2:2:ncol(rain.wide))] ,1, FUN = mean, na.rm = TRUE), add=TRUE, col = "darkred", frac = 0.2, lwd = 4)

abline(v = 21, lty =3, col = "darkgreen")
text(25,0,".", col = "darkgreen")
text(25,1.0,".", col = "darkgreen")

# -------------------------------------------------------
# PANEL B:    explore primer means  and  SD
#
par(cex=1.7)
par(bty = "l")
plot( 
apply(rain.wide[,c(2:ncol(rain.wide))],1, mean, na.rm = TRUE), 
apply(rain.wide[,c(2:ncol(rain.wide))],1, sd, na.rm = TRUE),
pch = 16,
cex=1,
# ylim = c(0,8),
xlim = c(20, 50),
col = "darkred",  
lwd = 2, 
las = 1,
main = "b.", 
ylab = "sd CT' / primer", 
xlab = "mean CT' / primer")

abline(h = 2, lty =2, col = "lightgray")
abline(h = 4, lty =2, col = "lightgray")
abline(h = 6, lty =2, col = "lightgray")
abline(h = 8, lty =2, col = "lightgray")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
#  STEP:		convert data frame into an all numeric data frame
#
# delete primer column and label rownames with primer names
# rain <- rain[,8:(ncol(rain))] 
head(rain.wide)
colnames(rain.wide)
GLOBAL.allnumericStart = 2

rain.wide <- as.matrix(rain.wide[,GLOBAL.allnumericStart:ncol(rain.wide)])
rain.wide <- t(rain.wide)
head(rain.wide)


# -----------------------------------------------------------------------------------------------------
# NOTE THIS REASSIGN RAIN GLOBALLY
rain <- as.data.frame(rain.wide)
head(rain[,1:4])


# -----------------------------------------------------------------------------------------------------
# STEP:  dCT CALCULATION
#  
# -----------------------------------------------------------------------------------------------------
# remove NTC samples   NEEDS TESTING
rain2 <- rain
# rain2[rain2$NTC == "yes",-ncol(rain2)] <- 40

# evaluate housekeeping genes
head(rain)
colnames(rain)
house <- rain[,c("BACT_2","B2M_2","ACTIN_2")]


# -----------------------------------------------------------------------------------------------------
# remove indiviudal experiments with too low RNA:   In this case sample 4438
#
splom(house, col = "darkred", pch =19, pscales=2, type=c("g","p"), main = "housekeeping genes")

house$sample <-rownames(house)
house
house <- house[house$sample != "4428",]

rain2$sample <- rownames(rain2)

rain2 <- rain2[rain2$sample != "4428",]
rain2 <- rain2[,-ncol(rain2)]


house <- as.matrix(house[,-ncol(house)])
house

nrow(rain)
nrow(rain2)
rain <- rain2

# -----------------------------------------------------------------------------------------------------
# FIGURE:		pairwise comparisons using PAIRS.PANELS
#
png(file = "HousekeepingGenes.png", units = "in", width = 8, height = 8, res = 300)
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
cex.labels = 1,font.labels = 2,
gap =1.2, 
las = 1,
scale = TRUE,
digits = 3,
main = "Housekeeping mRNAs"
# xlim = c(5,35), ylim = c(5,35)
)
}
dev.off()


# --------------------------------------------------------------------------------------------------------
# SCREEN:		heatmap of paIwise correlation coefficients using COR.PLOT
#
# unsorted
house.cor <- cor(house,use="complete.obs",method="pearson")
cor.plot(house.cor, main = "Housekeeping genes", zlim = c(0.8,1), cex.axis = 1.0, n = 32)

# sorted
house.cor <- mat.sort(house.cor)
head(house.cor)
cor.plot(house.cor, main = "Housekeeping Genes", zlim = c(0.5,1), cex.axis = 1.0, n = 32)
house


# --------------------------------------------------------------------------------------------------------
# dCT calculation  based on geometric mean of house 
# 

# Base on the above plots all or only a subset of samples are used to calculate dCT
house <- as.matrix(rain[,c("BACT_2","B2M_2","ACTIN_2")])
house
geo <- as.data.frame(apply(house, 1, geometric.mean))
colnames(geo) <- c("house")
geo$name <- rownames(house)
geo

head(rain)
rain$name <- as.factor(rownames(rain))
rain <- cbind(geo,rain)
rain <- rain[,-ncol(rain)]
rain$name <- as.factor(rain$name)
colnames(rain)

# modify global
GLOBAL.allnumericStart <- GLOBAL.allnumericStart + 1
GLOBAL.allnumericStart
head(rain[,1:4])


# -------------------------------------------------------
# calculate dCT 
# still need to convert to apply or sweep
#
for (i in GLOBAL.allnumericStart:ncol(rain))
	{
	rain[,i] <- 		rain[,"house"] - rain[,i]
	}

summary(rain)


# -------------------------------------------------------
# write out dCT to file
#
write.table(rain,file = "dCT.txt", sep = "\t",row.names = FALSE)


# -----------------------------------------------------------------------------------------------------
# FIGURE:		Comparison of sample sets: Positive, NTC, samples
#
rownames(rain)
png(file = "AlteredGenesBasedonSDofdCT'.png", units = "in", width = 6, height = 6, res = 300)
{
par(mfrow = c(1,1))
par(cex = 1.6)
par(bty = "l")
# positive
plot(ecdf(apply(rain[,GLOBAL.allnumericStart:ncol(rain)],2, sd)), col = "darkred", verticals = TRUE, do.points = FALSE, lwd = 2, main = "", xlab = "s.d. of dCT'", ylab = "Culm. Density", xlim = c(0,10), lty =1, cex = 2, las = 1)    
# apply: 2 is for column, mad is for median absolute deviation
histSpike(apply(rain[,GLOBAL.allnumericStart:ncol(rain)],2, sd), add=TRUE, col = "darkgray", frac = 0.3, lwd = 5)
}
dev.off()

head(rain[,1:4])


# -----------------------------------------------------------------------------------------------------
# STEP 4:  DELETE NON-CHANGED GENES based on sd distribution
#  
head(rain)

# set sd < 1
subset(rain[,GLOBAL.allnumericStart:ncol(rain)], apply(rain[,GLOBAL.allnumericStart:ncol(rain)],1, sd) < 1)
rain.2 <- subset(rain[,GLOBAL.allnumericStart:ncol(rain)], apply(rain[,GLOBAL.allnumericStart:ncol(rain)],1, sd) > 1)
ncol(rain)
ncol(rain.2)
head(rain.2[,1:4])

png(file = "SD filter.png", units = "in", width = 8, height = 8, res = 600)
{
par(mfrow = c(1,2))
par(cex = 1.5)
plot(ecdf(apply(rain[,GLOBAL.allnumericStart:ncol(rain)],2, sd)), col = "darkblue", verticals = TRUE, do.points = FALSE, lwd = 2, main = ncol(rain), xlab = "sd CT' ")    
histSpike(apply(rain[,GLOBAL.allnumericStart:ncol(rain)],2, sd), add=TRUE, col = "darkred", frac = 0.2, lwd = 3)

plot(ecdf(apply(rain.2[,GLOBAL.allnumericStart:ncol(rain.2)],2, sd)), col = "darkblue", verticals = TRUE, do.points = FALSE, lwd = 2, main = ncol(rain.2), xlab = "CT' ")    
histSpike(apply(rain.2[,GLOBAL.allnumericStart:ncol(rain.2)],2, sd), add=TRUE, col = "darkred", frac = 0.2, lwd = 3)
}
dev.off()

# GLOBALLY assign rain.2 to rain
rain <- rain.2
GLOBAL.allnumericStart <- 1


# -----------------------------------------------------------------------------------------------------
# STEP: Evaluate the entire dCT distribution of the samples
# -----------------------------------------------------------------------------------------------------
colnames(rain)
rownames(rain)
GLOBAL.allnumericStart
head(rain[,1:5])

rainCluster	<- melt(rain)
rainCluster <- na.omit(rainCluster)
rainCluster <- drop.levels(rainCluster)
colnames(rainCluster) <- c("Primer","dCT")
summary(rainCluster)
nrow(rainCluster)
head(rainCluster)


# -------------------------------------------------------
# Formal test for Normality
#
IQR(rainCluster$dCT, na.rm = TRUE)
fivenum(rainCluster$dCT, na.rm = TRUE)
shapiro.test(rainCluster$dCT)


# -----------------------------------------------------------------------------------------------------
# FIGURE: qqPLOT of entire dCT distribution
#
png(file = "qqPlot Of All dCT.png", units = "in", width = 8, height = 8, res = 300)
{
par(mfrow = c(1,1))
par(cex = 1.8)
par(bty = "o")
qqnorm(rainCluster[,"dCT"], 
pch = 19,
cex = 0.5, 
las = 1,
col = "darkred",
main = "Samples")
qqline(rainCluster[,"dCT"], col = "darkblue", lwd = 2, lty = 3)
histSpike(rainCluster[,"dCT"], add=TRUE, col = "darkgray", frac = 0.5, lwd = 5, side =2)

abline(h = 0, lty =2, col = "lightgray")
abline(v = 0, lty =2, col = "lightgray")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# FIGURE:  	HEATMAP
#					IN: rain.wide
#
summary(rain[,1:4])
GLOBAL.allnumericStart
rownames(rain)


# -----------------------------------------------------------------------------------------------------
# remove housekeeping genes
#
colnames(rain)
rain[,53:55]
rain <- rain[,-c(53:55)]


# -----------------------------------------------------------------------------------------------------
# remove indiviudal experiments
#
rownames(rain)


try2m <- as.matrix(rain[,GLOBAL.allnumericStart:ncol(rain)]) 
nrow(try2m)
nrow(try2m[!complete.cases(try2m),])
rownames(try2m)
head(try2m[,1:5])
nrow(try2m)


# -----------------------------------------------------------------------------------------------------
# FIGURE:   variant #1 "orange"
png(file = "heatmap dCT' orange.png", units = "in", width = 21, height = 7, res = 300)
{
rgb.palette <- colorRampPalette(c("white", "orange", "red"), space = "rgb", bias = 1)
heatmap.plus(try2m, na.rm = T, scale = "none", col = rgb.palette(25), hclustfun=function(m) hclust(m, method="ward.D"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,15), Colv = T, Rowv = T)
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# FIGURE:	variant #2 "color scheme"
#
png(file = "heatmap dCT' alternativeColor.png", units = "in", width = 21, height = 7, res = 300)
{
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")),space = "rgb", bias = 1.2, interpolate = "linear")
heatmap.plus(try2m, na.rm = T, scale = "none", col = brewer.palette (55), hclustfun=function(m) hclust(m, method="ward.D"), distfun = function(x) dist(x, method = "manhattan"), margins = c(10,15), Colv = T, Rowv = T, main = "dCT", xlab = "CDS", ylab = "samples")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# FIGURE:	variant #3 "color scheme" scaled by row
#
png(file = "heatmap dCT'alternativeColor scaled by sample.png", units = "in", width = 21, height = 7, res = 300)
{
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")),space = "rgb", bias = 1, interpolate = "linear")
heatmap.plus(try2m, na.rm = T, scale = "row", col = brewer.palette (55), 
hclustfun=function(m) hclust(m, method="ward.D"), 
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
cor.plot(primer.cor, main = "all genes", zlim = c(0.5,1), cex.axis = 0.8, n = 12)

# sorted
primer.cor <- mat.sort(primer.cor)
head(primer.cor)
cor.plot(primer.cor, main = "dCT", zlim = c(-1,1), cex.axis = 0.8, n = 12)


# -----------------------------------------------------------------------------------------------------
# heatmap of centerd but not scaled primers
# 
head(experimentByprimer[,1:4])
png(file = "heatmap dCT' centered By Primer.png", units = "in", width = 21, height = 7, res = 300)
{
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")),space = "rgb", bias = 0.8, interpolate = "linear")
heatmap.plus(experimentByprimer, na.rm = T, scale = "none", col = brewer.palette (55),
hclustfun=function(m) hclust(m, method="ward.D"), 
distfun = function(x) dist(x, method = "manhattan"), 
margins = c(12,24), 
Colv = T, Rowv = T, 
cexRow = 1, cexCol = 0.5,
main = "dCT' centered by primer", xlab = "primer", ylab = "samples")
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# heatmap of centerd primers scaled by experiment
# 
head(experimentByprimer)
png(file = "heatmap dCT' centered By Primer purple.png", units = "in", width = 21, height = 7, res = 300)
{
brewer.palette <- colorRampPalette(rev(brewer.pal(11, "PuOr")),space = "rgb", bias = 0.8, interpolate = "spline")
heatmap.plus(experimentByprimer, na.rm = T, scale = "none", col = brewer.palette (55),
hclustfun=function(m) hclust(m, method="ward.D"), 
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

png(file = "PCA of samples.png", units = "in", width = 8, height = 8, res = 300)
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
waterfall <- waterfall[order(waterfall[,1], decreasing = F),]

colnames(waterfall)
par(mfrow = c(1,1))
par(cex=1.0)
par(bty="n")
barplot(1.8^waterfall[,1], 
main = "K1 deletion", 
pch = 19, 
col = "gray", 
# type = "s", 
names.arg = rownames(waterfall), 
ylab = "fold induction", 
xaxs = "r",yaxs = "r", 
xaxt ="s",
las = 3,
# cex = 0.2, lwd = 2,  
ylim = c(0,5)
)


abline(h = median(waterfall[,1]), lwd = 1, col = "gray", lty = 2)
abline(h = (median(waterfall[,1])-sd(waterfall[,1])), lwd = 1, col = "gray", lty = 2)
abline(h = (median(waterfall[,1])+sd(waterfall[,1])), lwd = 1, col = "gray", lty = 2)


head(waterfall)
tail(waterfall)
write.table(waterfall,file = "dCTordered.txt", sep = "\t")