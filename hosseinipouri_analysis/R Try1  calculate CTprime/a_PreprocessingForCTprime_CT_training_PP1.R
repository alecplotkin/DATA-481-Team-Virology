# -------------------------------------------------------------------------------------------------
# Preprocessing CT'
# by DPD 04182014
#

# -------------------------------------------------------------------------------------------------
# read in libraries
#
library(xlsx)
library(reshape2)
library(ggplot2)
library(gdata)

# -----------------------------------------------------------------------------------------------------
# Define digits
# 
options(digits=3)


# -----------------------------------------------------------------------------------------------------
# Define working directory
# 
# setwd("/Users/dirk/DittmerSync/paper 2014 array Malawi Array/R Try1  calculate CTprime")


# -----------------------------------------------------------------------------------------------------
# Evaluate CT training file 1 :  "CT_training_PP1.xlsx"
# 
data=read.xlsx('CT_training_PP1.xlsx',sheetIndex=1)
head(data)

#  Check unique primer names
summary(data)
data[data$Primer =="K4 (21778)",]
nrow(data[data$Primer =="K4 (21778)",])

# this primer was inadvertantly duplicated. It is deleted new file written
data <- data[-7,]
write.xlsx(data, file = "a_CT_training_PP1v04182014.xlsx",row.names=FALSE)


# -----------------------------------------------------------------------------------------------------
# Evaluate CT training file 1 :  "CT_training_PP1v04182014.xlsx"
data=read.xlsx('a_CT_training_PP1v04182014.xlsx',sheetIndex=1)
head(data)
summary(data)

# FIGURE:  Density distribution of data
data <- melt(data)
summary(data)

 
par(cex = 2)
plot(density(data$value, na.rm = TRUE), col = "blue", lwd = 2, main = "CT values")

summary(data)
zp1 <- ggplot(data, aes(x = value, fill = variable))
zp1 <- zp1 + geom_density() #+ stat_smooth(se = TRUE, method = "rlm", alpha = 0.15, level = 0.90, lwd = 1)
zp1 <- zp1 # + scale_colour_manual(values = c("blue","red"))
zp1 <- zp1 # + scale_x_log10(breaks = c(0.0001,0.001,0.01,0.1,1)) #, scales = "free_y"
zp1 <- zp1 # + facet_wrap( ~ cancerdx)
zp1 <- zp1 + theme_bw()
zp1 <- zp1 + theme(legend.position = "right")  + theme(strip.background = element_rect(colour = "white"))
zp1 <- zp1 + theme(axis.text.x = element_text(size = 24, angle = 0))  + theme(axis.text.y = element_text(size = 24, angle = 0))  + theme(axis.title.y = element_text(size = 24, angle = 90))  + theme(plot.title = element_text(size = 24))  + theme(axis.title.x = element_text(size = 24))  + theme(strip.text.x = element_text(size = 24))
print(zp1) 

