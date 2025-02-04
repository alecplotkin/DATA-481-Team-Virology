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
