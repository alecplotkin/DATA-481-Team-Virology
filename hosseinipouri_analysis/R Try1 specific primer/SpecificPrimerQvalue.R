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
library(qvalue)


# -----------------------------------------------------------------------------------------------------
# STEP:    Read in p value file
#
b <- read.delim("ddCT_pvalues.txt", header = TRUE, sep = "\t",col.names =c("gene","p"))


# -----------------------------------------------------------------------------------------------------
# q-value adjustment
summary(b)
head(b)
rainq <- qvalue(b$p, robust = TRUE, fdr.level = .5, pi0.method="bootstrap")
## qsummary(rainq)
## qwrite(rainq, filename = "adjustedp values.txt")
b$q <- rainq$qvalues
summary(b)
## write.table(b,file = "qAndpValuesOfreducedSetOfExperiments.txt", sep = "\t",row.names = FALSE)


# -----------------------------------------------------------------------------------------------------
#
png(file = "qplotOfWilcox.png", units = "in", width = 8, height = 8, res = 600)
{

    plot(rainq, rng = c(0.0, 0.8) )

}
dev.off()`


# -----------------------------------------------------------------------------------------------------
# Densityplot of s
png(file = "SinglePvalue.png", units = "in", width = 8, height = 8, res = 600)
{
par(mfrow = c(2,1))
par(cex = 1.2)
plot(density(log10(b$p)), main = "p-value", xlab = "lg(p-value)", col = "darkred", lwd = 2, cex = 2, pch =1)
abline(v = log10(0.001), lty = 2, col = "darkblue")

plot(density(log10(b$q)), main = "q-value", xlab = "lg(q-value)", col = "darkred", lwd = 2, cex = 2, pch =1)
abline(v = log10(0.001), lty = 2, col = "darkblue")
}
dev.off()


