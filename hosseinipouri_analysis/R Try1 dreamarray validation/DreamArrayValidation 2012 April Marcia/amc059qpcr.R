# -----------------------------------------------------------------------------------------------------
# Dirk Dittmer on 04242014
# -----------------------------------------------------------------------------------------------------
#
# -----------------------------------------------------------------------------------------------------
# STEP 0:    Read in the libraries.
# -----------------------------------------------------------------------------------------------------
#
library(gdata)
library(MASS)
library(DAAG)
library(RColorBrewer)
library(FactoMineR)
library(heatmap.plus)
library(ggplot2)
library(DMwR)
library(reshape2)


# -----------------------------------------------------------------------------------------------------
# STEP 1:    Read in the data set
#
# -----------------------------------------------------------------------------------------------------
#
# Define read files function (needs > 1 filename)
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
# Read in the data     		must be in the folder as defined by workspace
#                       				must be at least two file names
#                  
filenames = c("PP1_a.txt","PP1_b.txt","PP2_a.txt","PP2_b.txt")
rain <- read.rain(filenames)
summary(rain)
table(rain$filename)


# -----------------------------------------------------------------------------------------------------
# Safe a copy of the data set on disk
# Safe a copy of the data set as variable rainall
# RECOVER with  rain <- rainall
#
rainall <- rain
write.table(rainall,file = "CompleteDataSet.txt", sep = "\t")


# -----------------------------------------------------------------------------------------------------
# STEP 2:    Data cleaning
#
# -----------------------------------------------------------------------------------------------------
#
head(rain)
colnames(rain)
colnames(rain) <- c("Row","Col","Well","Primer","Sample" , "CT","Tm1","Tm2","Notes","x","x1","Filename")

# -----------------------------------------------------------------------------------------------------
# Use a subset of Variables (columns)
#
rain <-rain[,c("Row","Col","Well","Primer","Sample" , "CT","Tm1","Tm2","Filename")] 
summary(rain)


# -----------------------------------------------------------------------------------------------------
# STEP 3:    Data transformation 
#
# -----------------------------------------------------------------------------------------------------
#
rain$Col	<- factor(rain$Col)

levels(rain$Filename)
levels(rain$Filename) <- list(PP1 = c("PP1_a.txt","PP1_b.txt"), PP2 = c("PP2_a.txt","PP2_b.txt"))

rain$Primer2 <- interaction(rain$Filename, rain$Primer)
table(rain$Primer)

rain[is.na(rain$CT),] 
rain[is.na(rain$CT),"CT"] <- 40


# -----------------------------------------------------------------------------------------------------
# adjustment of Tms

rain[is.na(rain$Tm1),"Tm1"] <- 99.999
rain[is.na(rain$Tm2),"Tm2"] <- 99.999

rain[rain$CT>39.9,"Tm1"] <- 99.999
rain[rain$CT>39.9,"Tm2"] <- 99.999

rain[rain$Tm1<1,"Tm1"] <- 99.999
rain[rain$Tm2<1,"Tm2"] <- 99.999

temp.Tm1 <- rain[rain$Tm2 < rain$Tm1, "Tm1"]
rain[rain$Tm2 < rain$Tm1, "Tm1"] <- rain$Tm2
rain[rain$Tm2 < rain$Tm1, "Tm2"] <- temp.Tm1 

# -----------------------------------------------------------------------------------------------------
# drop out NAs
rain[is.na(rain$Primer),]
rain[is.na(rain$Sample),]
rain <- na.omit(rain)
rain <- drop.levels(rain)


summary(rain)
head(rain)
table(rain$Primer2)


# -----------------------------------------------------------------------------------------------------
# STEP 4:  Unsupervised clustering
#  
# -----------------------------------------------------------------------------------------------------
#
rainCluster	<- rain[,c("Primer2","Sample","CT")]


# -----------------------------------------------------------------------------------------------------
# STEP 5:  figure  "input"   All the data without clustering
#  
# -----------------------------------------------------------------------------------------------------
#
png(file = "Input.png", units = "in", width = 5, height = 20, res = 600)
{
a <- qplot(Sample, Primer2, col = cut(CT,4), data = rain, size = I(2), alpha = I(1/2), geom= c("point"))
a <- a + opts(strip.background = theme_rect(colour = "white")) + opts(legend.position = "right") 
a <- a + opts(axis.text.x = theme_text(size = 8, angle = 90)) + opts(axis.text.y = theme_text(size = 8, angle = 0)) + opts(axis.title.y = theme_text(size = 14, angle = 90)) + opts(plot.title = theme_text(size = 18)) + opts(axis.title.x = theme_text(size = 14)) + opts(strip.text.x = theme_text(size = 14))
a <- a + opts(panel.grid.minor = theme_blank())+ opts(panel.grid.major = theme_blank())
a
}
dev.off()

# convert to long format with CAST
#
rainCluster <- dcast(rainCluster, Primer2 ~ Sample, mean)
summary(rainCluster)
nrow(rainCluster)
head(rainCluster)
colnames(rainCluster)


# -----------------------------------------------------------------------------------------------------
# STEP 6:		subsetting to identify those entries with any variation 
#					as measured by m.a.d > 1
# 					evaluate using ecdf
# -----------------------------------------------------------------------------------------------------
#
# define function to compute m.a.d. with na's
#
mad.na <- function(x)
{
mad(x, na.rm = TRUE)
}


# -----------------------------------------------------------------------------------------------------
# Figure: Comparison of  m.a.d > 1 based filtering
#
png(file = "mad.png", units = "in", width = 5, height = 5, res = 600)
{
par(mfrow = c(1,2))
par(cex=1.4)
plot(ecdf(apply(rainCluster[,2:ncol(rainCluster)],1, mad.na)), main = nrow(rainCluster), col = "darkred", cex = 0.5, xlab = "m.a.d in CT units")

rainCluster.mad <- subset(rainCluster, apply(rainCluster[,2:ncol(rainCluster)],1, mad.na) > 0)

plot(ecdf(apply(rainCluster.mad[,2:ncol(rainCluster)],1, mad.na)), main = nrow(rainCluster.mad), col = "darkblue", cex = 0.5, xlab = "m.a.d in CT units")
nrow(rainCluster)
}
dev.off()

# convert data frame to matrix of number
#
try2m <- as.matrix(rainCluster[,2:ncol(rainCluster)])
nrow(try2m)
nrow(try2m[!complete.cases(try2m),])
rownames(try2m) <- rainCluster[,1]
colnames(try2m) <- colnames(rainCluster[,2:ncol(rainCluster)])
head(try2m)


# -----------------------------------------------------------------------------------------------------
# Figure1 :  heatmap save to file
# -----------------------------------------------------------------------------------------------------
# Set up file and colors
#

# heatmap and clustering
#
#	method = average
# distance = manhattan
# 

png(file = "heatmap.png", units = "in", width = 4, height = 20, res = 300)
{
rgb.palette <- colorRampPalette(c("red", "orange", "white"), space = "rgb", bias = 2)
heatmap.plus(try2m, 
na.rm = T, scale = "none", 
col = rgb.palette(20), 
hclustfun=function(m) hclust(m, method="average"), 
distfun = function(x) dist(x, method = "manhattan"), 
margins = c(10,10), Colv = NA, Rowv = T)
}
dev.off()

# -----------------------------------------------------------------------------------------------------
# Figure 2
# 
head(try2m)



png(file = "Figure_ALL_heatmap centeredPrimer.png", units = "in", width = 4, height = 20, res = 300)
{
brewer.palette <- colorRampPalette(brewer.pal(11, "RdYlBu"),space = "rgb", bias = 1, interpolate = "linear")
heatmap.plus(try2m, na.rm = T, scale = "row", col = brewer.palette (55),
hclustfun=function(m) hclust(m, method="ward"), 
distfun = function(x) dist(x, method = "manhattan"), 
margins = c(10,10), 
Colv = NA, Rowv = T, 
cexRow = 0.7, cexCol = 1.4,
main = "", xlab = "sample", ylab = "Primer")
}
dev.off()

head(try2m)



# -----------------------------------------------------------------------------------------------------
# Optimization
# 


