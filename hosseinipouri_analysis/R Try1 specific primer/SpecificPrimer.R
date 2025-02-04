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
library(xlsx)
library(Hmisc)
library(psych)
library(ggplot2)

# library(FactoMineR)
# library(heatmap.plus)
# library(beeswarm)
# library(vioplot)
# library(DMwR)

# library(vegan)
# library(qvalue)



# -----------------------------------------------------------------------------------------------------
# STEP:    Define digits
#
options(digits=3)


# -----------------------------------------------------------------------------------------------------
# STEP: Define working directory
# 
setwd("/Users/dirk/DittmerSync/paper 2014 array Malawi Array/R Try1 specific primer")


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
filenames <- c("dCT.xlsx")
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

colnames(rain) <- c(
"sample","orf","position","B","B.1","B.2","A" ,"B.3" ,"A.1"  , "A.2" ,"B.4" , "A.3" ,"A.4" , "A.5" ,"A.6", "B.5","B.6"     ,
"A.7" ,"B.7" ,"B.8" ,"B.9" ,"A.8" ,"B.10","A.9" , "B.11" , "B.12" , "B.13" , "B.14" ,"A.10" , "A.11" ,"A.12" ,"A.13", "B.15", "A.14",    
 "A.15" ,"B.16" , "B.17","B.18" ,  "filename"
)

#drop out filename column
rain <- rain[,-ncol(rain)]  
summary(rain)

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
colnames(rain)
rain[duplicated(rain),]
summary(rain)


# -----------------------------------------------------------------------------------------------------
# Define data frame rain.wide
#
rain.wide <- na.omit(rain)
rain.wide <- drop.levels(rain)

head(rain.wide)
summary(rain.wide[,1:11])
GLOBAL.allnumericStart = 4

# Density plot on screen
par(mfrow = c(1,1))
par(bty = "n")
par(cex=1.8)
# samples
mean <- rowMeans(rain.wide[,GLOBAL.allnumericStart:ncol(rain.wide)], na.rm = TRUE)
plot(density(mean), main = "", col = "red", lwd = 2)


# -----------------------------------------------------------------------------------------------------
# FIGURE: ECDF
#
png(file = "PrimerPerformance.png", units = "in", width = 10, height = 6, res = 300)
{
par(mfrow = c(1,2))
par(cex=1.7)
par(bty="l")
plot(ecdf(apply( rain.wide[,GLOBAL.allnumericStart:ncol(rain.wide)] ,1, FUN = mean, na.rm = TRUE)), col = "darkblue", verticals = TRUE, do.points = FALSE, lwd = 2, main = "a.", xlab = "mean dCT / primer", ylab = "Cumulative density", las = 1)    
histSpike(apply(rain.wide[,GLOBAL.allnumericStart:ncol(rain.wide)] ,1, FUN = mean, na.rm = TRUE), add=TRUE, col = "darkred", frac = 0.2, lwd = 4)

abline(v = 21, lty =3, col = "darkgreen")
text(25,0,".", col = "darkgreen")
text(25,1.0,".", col = "darkgreen")

# -------------------------------------------------------
# PANEL B:    explore primer means  and  SD
#
par(cex=1.7)
par(bty = "l")
plot( 
apply(rain.wide[,c(GLOBAL.allnumericStart:ncol(rain.wide))],1, mean, na.rm = TRUE), 
apply(rain.wide[,c(GLOBAL.allnumericStart:ncol(rain.wide))],1, sd, na.rm = TRUE),
pch = 16,
cex=1,
# ylim = c(0,8),
# xlim = c(20, 60),
col = "darkblue",  
lwd = 2, 
las = 1,
main = "b.", 
ylab = "sd CT / primer", 
xlab = "mean CT / primer")

abline(h = 2, lty =2, col = "lightgray")
abline(h = 4, lty =2, col = "lightgray")
abline(h = 6, lty =2, col = "lightgray")
abline(h = 8, lty =2, col = "lightgray")
}
dev.off()


summary(rain.wide)


# -----------------------------------------------------------------------------------------------------
# Change to long format
#
colnames(rain.wide)
rain <- melt(data = rain.wide, id.var = c("sample","orf","position"), variable.name = "cluster", value.name = "dCT")
summary(rain)

tail(rain)
rain$clustername <- as.character(rain$cluster)
rain$clustername <- substring(rain$clustername,1,1)
rain$clustername <- as.factor(rain$clustername)

fivenum(rain$dCT)[2]
rain[rain$dCT< -29,"dCT"] <- jitter(-30, amount = 1)

rain$foldhouse <- 1.8^(rain$dCT)

attach (rain)
rain <- rain[order(rain$position),]
detach(rain)

# -----------------------------------------------------------------------------------------------------
# FIGURE ggplot
#
head(rain)


png(file = "FigureByPosition.png", units = "in", width = 20, height = 6, res = 300)
{
a <-ggplot(aes(position, dCT, colour = clustername), data = rain) # + geom_boxplot(fill = "lightgray", col="black") 
a <- a + geom_point (size = I(2), alpha = I(0.5)) 
a <- a + theme_bw() + xlab ("KSHV") + ylab ("dCT") + labs(title = "mRNA levels")
a <- a # + ylim(10^5,10^10) 
a <- a # + scale_y_log10()
a <- a + scale_colour_manual(values = c("red","blue"))
a <- a # + facet_grid( . ~ isotype)  
a <- a + theme(legend.position = "none") + theme(strip.background = element_rect(colour = "white", fill = "white"))   
a <- a + theme(axis.text.x = element_text(size = 16, angle = 0)) + theme(axis.text.y = element_text(size = 16, angle = 0)) + theme(axis.title.y = element_text(size = 16, angle = 90)) + theme(plot.title = element_text(size = 16)) + theme(axis.title.x = element_text(size = 16)) + theme(strip.text.x = element_text(size = 16))
a <- a + theme(panel.grid.minor = element_blank()) #+ theme(panel.grid.major = element_blank())
a 
}
dev.off()


head(rain)
png(file = "FigureByOrf.png", units = "in", width = 20, height = 6, res = 300)
{
a <-ggplot(aes(x=reorder(orf, dCT, FUN=sum), y = dCT, fill = clustername), data = rain)  + geom_boxplot(notch = FALSE, outlier.size = 0.1, outlier.colour = "white", alpha = I(0.5), color = "white") 
a <- a # + geom_point (size = I(2), alpha = I(0.5)) 
a <- a + theme_bw() + xlab ("KSHV") + ylab ("dCT") + labs(title = "mRNA levels")
a <- a # + ylim(10^5,10^10) 
a <- a # + scale_y_log10()
a <- a + scale_fill_manual(values = c("red","blue"))
a <- a # + facet_grid( . ~ isotype)  
a <- a + theme(legend.position = "none") + theme(strip.background = element_rect(colour = "white", fill = "white"))   
a <- a + theme(axis.text.x = element_text(size = 12, angle = 90)) + theme(axis.text.y = element_text(size = 16, angle = 0)) + theme(axis.title.y = element_text(size = 16, angle = 90)) + theme(plot.title = element_text(size = 16)) + theme(axis.title.x = element_text(size = 16)) + theme(strip.text.x = element_text(size = 16))
a <- a + theme(panel.grid.minor = element_blank()) #+ theme(panel.grid.major = element_blank())
a 
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# aggregate
#
summary(rain)


# calculate median
#
by(rain[,"dCT"], INDICES = list(rain$clustername, rain$orf), median)

# calculate median
#
tapply(rain$dCT, list(rain$orf, rain$clustername), FUN=median)

# calculate median
#
dcast(rain, clustername ~ orf, value.var="dCT", fun.aggregate=median)


# -----------------------------------------------------------------------------------------------------
# calculate statistics by group
#
head(rain)
names <- levels(rain$orf)
all.p <- as.data.frame(names)
head(all.p)
nrow(all.p)

# LOOP
#
i <- 1
for (i in 1:nrow(all.p))
{
a <- as.data.frame (split(rain, rain$orf)[i])
head(a)
nrow(a)
a.p <- as.numeric((wilcox.test(a[,5] ~ a[,6])$p.value))
print(a.p)
all.p[i,"ap"] <- a.p
all.p
}
all.p


# example of one case
#
a <- as.data.frame (split(rain, rain$orf)[6])
wilcox.test(a[,5] ~ a[,6])


# output
#
write.table(all.p,file = "pvalues.txt", sep = "\t", row.names = FALSE)
