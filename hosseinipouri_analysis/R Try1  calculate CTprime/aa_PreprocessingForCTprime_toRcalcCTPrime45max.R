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
data<- read.csv('toR_calcCTPrime_no60.csv')
head(data)

#  Check unique primer names
summary(data[,1:10])


# Corrective action No. 1
# this primer was inadvertantly duplicated. It is deleted 
data[data$Primer =="K4 (21778)",]
nrow(data[data$Primer =="K4 (21778)",])
nrow(data)
data <- data[-23,]
nrow(data)
summary(data[,1:10])


# Corrective action No. 2
# this primername was duplicated. However, the data are different. Hence a new primername is made
data[data$Primer =="Orf73 (LANA) (124002,2)",]
data[147,]
duplicate.name <- data[147,]
nrow(data)
data <- data[-147,]
nrow(data)
duplicate.name[,"Primer"] <- as.factor("Orf73 (LANA) (124002,2)a")
duplicate.name<- drop.levels(duplicate.name)
summary(duplicate.name)

nrow(data)
data <- rbind(data,duplicate.name)
nrow(data) 

# Corrective action No. 3
# this primer  is present on both plates
summary(data[,1:4])
nrow(data)
data[data$Primer =="ORF70 (20979)",]
# duplicate.name <- data[141,]
# data <- data[-141,]
# colnames(data)
# ncol(data)
# duplicate.name [,10:ncol(data)] <- colMeans(data[data$Primer =="ORF70 (20979)",10:ncol(data)], na.rm = TRUE)
# duplicate.name<- drop.levels(duplicate.name)
# duplicate.name

# data[data$Primer =="ORF70 (20979)",]
# data <- subset(data, data$Primer !="ORF70 (20979)")
# data <- drop.levels(data)
# summary(data[,1:10])

# duplicate.name<- drop.levels(duplicate.name)
# summary(duplicate.name)

# nrow(data)
# data <- rbind(data,duplicate.name)
# nrow(data) 

summary(data[,1:10])

write.csv(data, file = "aa_calcCTPrime_no60_04182014a.csv",row.names=FALSE)


# -----------------------------------------------------------------------------------------------------
# Evaluate CT training file 1 :  "CT_training_PP1v04182014.xlsx"
#
data=read.csv('aa_calcCTPrime_no60_04182014a.csv')
head(data)
summary(data[,1:10])
data$Date <- as.factor(data$Date)

# -----------------------------------------------------------------------------------------------------
# Replace NA with jittered 45
#
head(data)
GLOBAL.numeric <- 10	#start of the numeric data
i <- GLOBAL.numeric
for (i in GLOBAL.numeric:ncol(data))
{
data[is.na(data[,i]),i] <- 55
data[data[,i]>=45,i] <- jitter(data[data[,i]>=45,i], amount = 1)
}
head(data)

write.csv(data, file = "aa_calcCTPrime_45max_04182014a.csv",row.names=FALSE)


# FIGURE:  Density distribution of data
data <- melt(data)
summary(data)

 
par(cex = 2)
plot(density(data$value, na.rm = TRUE), col = "blue", lwd = 2, main = "CT values")

png(file = "aa_CTdistribution.png", units = "in", width = 8, height = 8, res = 300)
{
summary(data)
zp1 <- ggplot(data, aes(x = value, color = Type))
zp1 <- zp1 + geom_density() #+ stat_smooth(se = TRUE, method = "rlm", alpha = 0.15, level = 0.90, lwd = 1)
zp1 <- zp1 + scale_colour_manual(values = c("blue","red"))
zp1 <- zp1 # + scale_x_log10(breaks = c(0.0001,0.001,0.01,0.1,1)) 
zp1 <- zp1 + facet_wrap( ~ variable) #, scales = "free_y"
zp1 <- zp1 + theme_bw()
zp1 <- zp1 + theme(legend.position = "top")  + theme(strip.background = element_rect(colour = "white", fill = "white"))
zp1 <- zp1 + theme(axis.text.x = element_text(size = 20, angle = 0))  + theme(axis.text.y = element_text(size = 14, angle = 0))  + theme(axis.title.y = element_text(size = 20, angle = 90))  + theme(plot.title = element_text(size = 20))  + theme(axis.title.x = element_text(size = 14))  + theme(strip.text.x = element_text(size = 20))
zp1 <- zp1 + ggtitle ("CT distribution") + xlab("CT")
print(zp1) 
}
dev.off()



