# -----------------------------------------------------------------------------------------------------
# read in libraries
# -----------------------------------------------------------------------------------------------------
library(foreign)
library(gdata)
library(MASS)
library(DAAG)

# -----------------------------------------------------------------------------------------------------
# enrollment
# -----------------------------------------------------------------------------------------------------
mir <- read.delim("days.txt", header = TRUE, sep = "\t")
mir <- na.omit(mir)
mir <- drop.levels(mir)
nrow(mir)
summary(mir)

png(file = "Enrollment.png", units = "in", width = 10, height = 6, res = 300)
{
par(mfrow = c(1,1))
par(bty="l")
par(cex=1.5)
plot(mir$days, mir$step, main = "Enrollment into AMC S001", pch = 19, col = "red", type = "S", xlab = "days since opening", ylab = "No. of subjects", xaxs = "i",yaxs = "i", cex = 1.2, las = 1, lwd = 4, ylim = c(0,70),yaxp = c(0,70,2), xlim = c(0,700), xaxp = c(0,700,7))
rug(mir$days, col = "black", lwd = 3, cex = 1.0, side = 1)
abline(h = 35, col = "darkgray", lty = 4, lwd = 2)
}
dev.off()


