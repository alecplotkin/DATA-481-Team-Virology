# -------------------------------------------------------------------------------------------------
# Preprocessing 
# by DPD 04182014
#

# -------------------------------------------------------------------------------------------------
# read in libraries
#
library(gdata)
library(xlsx)
library(reshape2)
library(ggplot2)
library(psych)
library(beeswarm)
library(vioplot)

# -----------------------------------------------------------------------------------------------------
# Define digits
# 
options(digits=3)


# -----------------------------------------------------------------------------------------------------
# Define working directory
# 
 setwd("/Users/dirk/DittmerSync/paper 2014 array Malawi Array/R Try1  calculate CTprime")


# -----------------------------------------------------------------------------------------------------
# Evaluate CT training file 1 :  "CT_training_PP1.xlsx"
# 
data<- read.csv('toR_calcCTPrime_45max_04182014a.csv')
head(data)

#  Check unique primer names
summary(data[,1:10])


# Corrective action No. 1
# this primer was inadvertantly duplicated. It is deleted 
data[data$Primer =="K4 (21778)",]
nrow(data[data$Primer =="K4 (21778)",])
data <- data[-23,]
summary(data[,1:10])


# Corrective action No. 2
# this primername was duplicated. However, the data are different. Hence a new primername is made
data[data$Primer =="Orf73 (LANA) (124002,2)",]
data[147,]
duplicate.name <- data[147,]
data <- data[-147,]
duplicate.name[,"Primer"] <- as.factor("Orf73 (LANA) (124002,2)a")
duplicate.name<- drop.levels(duplicate.name)
summary(duplicate.name)

nrow(data)
data <- rbind(data,duplicate.name)
nrow(data) 

# Corrective action No. 3
# this primer was inadvertantly duplicated. It is merged
summary(data[,1:4])
data[data$Primer =="ORF70 (20979)",]
duplicate.name <- data[141,]
colnames(data)
ncol(data)
duplicate.name [,10:ncol(data)] <- colMeans(data[data$Primer =="ORF70 (20979)",10:ncol(data)], na.rm = TRUE)
duplicate.name

data[data$Primer =="ORF70 (20979)",]
data <- subset(data, data$Primer !="ORF70 (20979)")
data <- drop.levels(data)
summary(data[,1:10])

duplicate.name<- drop.levels(duplicate.name)
summary(duplicate.name)

nrow(data)
data <- rbind(data,duplicate.name)
nrow(data) 

summary(data[,1:10])



# -----------------------------------------------------------------------------------------------------
# Replace >40 with NA
#
head(data)

GLOBAL.numeric <- 10	#start of the numeric data
i <- GLOBAL.numeric
for (i in GLOBAL.numeric:ncol(data))
{
data[is.na(data[,i]),i] <- 99	
data[data[,i]>40,i] <- NA
}
head(data)


# -----------------------------------------------------------------------------------------------------
# Melt to Count NA
summary(data[,1:10])
data$Date <- as.factor(data$Date)
data <- melt(data)
summary(data)

# -----------------------------------------------------------------------------------------------------
# Count NA using the psych package
#

# Number of primers to human mRNAs
human <- data[data$Type == "human",]
human <- drop.levels(human)
no.humanP <- nlevels(human$Primer)
no.humanP

# Number of primers to KSHV mRNAs
kshv <- data[data$Type != "human",]
kshv <- drop.levels(kshv)
no.kshvP <- nlevels(kshv$Primer)
no.kshvP

# ----------------------------------------------------------------------------------------------------- 
# generate a summary table
#
data.s <- describeBy(	data$value, 	list(data$Type,data$variable),	skew=FALSE,ranges=FALSE, mat = TRUE)
data.s <- as.data.frame(data.s[,-1])  # describeBy adds a column name as new variable
head(data.s)
summary(data.s)

data.s <- data.s[,-3]
colnames(data.s)
colnames(data.s) <- c("Type","Sample","n","mean", "sd", "se")


# -----------------------------------------------------------------------------------------------------
# calculate percent present
#
data.s$percent <- 0
data.s[data.s$Type == "human","percent"] <- (no.humanP - data.s[data.s$Type == "human","n"]) / no.humanP
data.s[data.s$Type != "human","percent"] <- (no.kshvP -data.s[data.s$Type != "human","n"]) / no.kshvP

# -----------------------------------------------------------------------------------------------------
# Plot
#
plot(data.s$Type,data.s$percent, )

png(file = "b_DetectableBee.png", units = "in", width = 8, height = 8, res = 300)
{
summary(data.s)
par(cex = 2)
par(bty ="o")
beeswarm(percent ~ Type, data = data.s, 
col = "darkred",
cex = 1.3,
pch = 16,
method = "center",
las =1,
ylab = "% absent")
abline(h = 0.25, lty = 2)
abline(h = 0.5, lty = 2)
abline(h = 0.75, lty = 2)
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# Plot
#
png(file = "b_Detectable.png", units = "in", width = 8, height = 8, res = 300)
{
summary(data.s)
zp1 <- ggplot(data.s, aes(y = percent, color = Type, factor(Type)))
zp1 <- zp1 + geom_violin(scale = "width") #+ stat_smooth(se = TRUE, method = "rlm", alpha = 0.15, level = 0.90, lwd = 1)
zp1 <- zp1 + geom_point()
zp1 <- zp1 + scale_color_manual(values = c("blue","red"))
zp1 <- zp1 # + scale_x_log10(breaks = c(0.0001,0.001,0.01,0.1,1)) 
zp1 <- zp1 # + facet_wrap( ~ variable) #, scales = "free_y"
zp1 <- zp1 + theme_bw()
zp1 <- zp1 + theme(legend.position = "top")  + theme(strip.background = element_rect(colour = "white", fill = "white"))
zp1 <- zp1 + theme(axis.text.x = element_text(size = 20, angle = 0))  + theme(axis.text.y = element_text(size = 14, angle = 0))  + theme(axis.title.y = element_text(size = 20, angle = 90))  + theme(plot.title = element_text(size = 20))  + theme(axis.title.x = element_text(size = 14))  + theme(strip.text.x = element_text(size = 20))
zp1 <- zp1 + ggtitle (" ") + xlab("Group") + ylab("% absent")
print(zp1) 
}
dev.off()


# -----------------------------------------------------------------------------------------------------
# Table
#
data.s <- data.s[order(data.s[,"Type"]),]
write.csv(data.s, file = "b_Table_data.csv",row.names=FALSE)

