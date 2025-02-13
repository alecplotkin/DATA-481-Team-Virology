# -----------------------------------------------------------------------------------------------------
# read in R libraries
# -----------------------------------------------------------------------------------------------------
library(foreign)
library(gdata)
library(MASS)
library(DAAG)
library(RColorBrewer)
library(reshape2)
library(psych)
library(FactoMineR)
library(heatmap.plus)
library(ggplot2)
library(SiZer)


# -----------------------------------------------------------------------------------------------------
# read in the STATA data
# -----------------------------------------------------------------------------------------------------
here::here()
s001 <- read.dta("AMC 001 for AMC stats stata verion 9.dta")
summary(s001)
head(s001)
nrow(s001)


# -----------------------------------------------------------------------------------------------------
# calculate lytic index
# -----------------------------------------------------------------------------------------------------
s001$lyticindex <- s001$orf36 + s001$orf21
s001$lyticindex <- factor(s001$lyticindex, ordered = TRUE)


# -----------------------------------------------------------------------------------------------------
# encode the factors
# -----------------------------------------------------------------------------------------------------
s001$age_cat <- factor(s001$age_cat, ordered = TRUE)
s001$weightloss <- factor(s001$weightloss, labels = c("no","yes"))
s001$bpmmhg <- factor(s001$bpmmhg)
s001$lymphnode <- factor(s001$lymphnode, labels = c("no","yes"))
s001$karnovisky <- factor(s001$karnovisky, ordered = TRUE)
s001$t_stage<- factor(s001$t_stage, ordered = TRUE)
s001$tumor1edema<- factor(s001$tumor1edema, labels = c("no","yes"))
s001$numberoflessions<- factor(s001$numberoflessions, ordered = TRUE)

s001$tumor1oral <- factor(s001$tumor1oral)
s001$tumor1gastro <- factor(s001$tumor1gastro)
s001$tumor1viscera <- factor(s001$tumor1viscera)
s001$gender <- factor(s001$gender, labels = c("female","male"))

s001$orf36 <- factor(s001$orf36, labels = c("no","yes"))
s001$Kover70 <- factor(s001$Kover70, labels = c("no","yes"))
s001$cd4under200 <- factor(s001$cd4under200, labels = c("no","yes"))

s001$bmi_low <- factor(s001$bmi_low, labels = c("no","yes"))
s001$KS_nd <- factor(s001$KS_nd, labels = c("no","yes"))
s001$orf21 <- factor(s001$orf21, labels = c("no","yes"))


# -----------------------------------------------------------------------------------------------------
# composite factor: lytic
# -----------------------------------------------------------------------------------------------------
s001$lytic <- interaction(s001$orf21, s001$orf36)

# -----------------------------------------------------------------------------------------------------
# composite factor: lytic
# -----------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------
# log transformation
# -----------------------------------------------------------------------------------------------------
s001$logcpskshv_<- log10(s001$cpskshv_+1)
s001$loghivrna <- log10(s001$hivrna+1)

# -----------------------------------------------------------------------------------------------------
# eliminate incomplete records
# -----------------------------------------------------------------------------------------------------
s001 <- na.omit(s001)
s001 <- drop.levels(s001)


# -----------------------------------------------------------------------------------------------------
# eliminate dulicate records
# -----------------------------------------------------------------------------------------------------
nrow(s001)
s001 <- unique(s001)
nrow(s001)


# -----------------------------------------------------------------------------------------------------
# fix temperture recording outlier
# -----------------------------------------------------------------------------------------------------
summary(s001)
s001[(s001$temperaturec > 100),"temperaturec"] <- s001[(s001$temperaturec > 136),"temperaturec"] - 100


# -----------------------------------------------------------------------------------------------------
# subset : drop out rare factors
# -----------------------------------------------------------------------------------------------------

# drop out: "no visceral involvement"
s001small <- subset(s001, s001$tumor1viscera != 1)
s001small <- drop.levels(s001small)
summary(s001small)

# drop out: "oral involvement"
s001small <- subset(s001small, s001small$tumor1oral != 1)
s001small <- drop.levels(s001small)
summary(s001small)

# drop out: "KSHV not detectable"
plot(density(s001small$logcpskshv_))
s001small <- subset(s001small, s001small$KS_nd == "no")
s001small <- drop.levels(s001small)
summary(s001small)


# -----------------------------------------------------------------------------------------------------
# subset: PCA on factors only
# -----------------------------------------------------------------------------------------------------

colnames(s001small)
s001factor <- s001small[,c("t_stage","tumor1edema","gender","weightloss","lymphnode","orf36","Kover70","cd4under200","bmi_low","numberoflessions","age_cat","orf21","lytic","logcpskshv_")]
summary(s001factor)
rownames(s001factor) <- s001small$pid
colnames(s001factor)
ncol(s001factor)

res.mca <- MCA(s001factor, quanti.sup = c(14))
plot.MCA(res.mca, invisible = c("var"))
plot.MCA(res.mca, invisible = c("ind"), habillage = "quali", title = "AMC S001")
dimdesc(res.mca)


# -----------------------------------------------------------------------------------------------------
# tables
# -----------------------------------------------------------------------------------------------------
attach(s001small)


table(t_stage, lymphnode)
fisher.test(t_stage, lymphnode,
            or = 1,
            alternative = "two.sided",
            conf.int = TRUE,
            conf.level = 0.95)

tapply(cd4count, cd4under200, median)
wilcox.test(cd4count ~ cd4under200)

tapply(cd4count, lymphnode, median)
wilcox.test(cd4count ~ lymphnode)

wilcox.test(cd4count ~ orf36)
wilcox.test(cd4count ~ orf21)
summary(aov(cd4count ~ lytic))
summary(aov(cd4count ~ lyticindex))

table(lytic, t_stage)
fisher.test(t_stage,lytic,
            or = 1,
            alternative = "two.sided",
            conf.int = TRUE,
            conf.level = 0.95)

summary(aov(cd4count ~ lytic,projections = TRUE))
summary(aov(cd4count ~ lyticindex,projections = TRUE))

summary(aov(loghivrna ~ lytic,projections = TRUE))
summary(aov(loghivrna ~ lyticindex))

summary(aov(logcpskshv_ ~ lytic,projections = TRUE))
summary(aov(logcpskshv_ ~ lyticindex))

summary(aov(age ~ lytic,projections = TRUE))
summary(aov(age ~ lyticindex,projections = TRUE))

detach(s001small)


# -----------------------------------------------------------------------------------------------------
# exploratory graphics
# -----------------------------------------------------------------------------------------------------
summary(s001small)
colnames(s001small)

# no optical evidence of correlation
lattice::splom(s001small[,c("t_stage","tumor1edema","age","gender","weightloss","temperaturec","bpmmhg","pulsebpm","resp","karnovisky","lymphnode","cd4count","orf36","bmi","numberoflessions","age_cat","orf21","lytic","logcpskshv_","loghivrna")])

lattice::splom(s001small[,c("t_stage","tumor1edema","gender","karnovisky","lymphnode","cd4count","orf36","bmi","numberoflessions","orf21","lytic","lyticindex","logcpskshv_","loghivrna")])


# -----------------------------------------------------------------------------------------------------
# broken stick regression
# -----------------------------------------------------------------------------------------------------
# standard linear regression
cfar <- lm(loghivrna ~ cd4count, data = s001small)
anova(cfar)
summary(cfar)

par(mfrow = c(1,2))
par(cex=1.2)
plot(s001small$cd4count, s001small$loghivrna,
	pch = 16,
	xlim = c(0, 800),
	#ylim = c(0,0.01),
	xlab = "CD4 count",
	ylab = "log(HIV Vl)",
	main = "AMC S001")
abline(cfar,
	lty = "dotted",
	lwd = 1,
	col = "darkred")
# mtext(expression(r^2:0.3139), side = 4, col = "darkred")

plot(cfar,2, xlim = c(-2,2))


# variants of broken stick regression of CD4 < 500
ss001small <- subset(s001small, s001small$cd4count < 500)
ss001small <- drop.levels(ss001small)
ss001small <- subset(ss001small, ss001small$cd4count >= 10)
ss001small <- drop.levels(ss001small)

cfar <- rlm(loghivrna ~ cd4count, data = ss001small)
model <- piecewise.linear(ss001small$cd4count, ss001small$loghivrna, CI = FALSE)
print(model)
model.bent <- bent.cable(ss001small$cd4count, ss001small$loghivrn, grid.size=20)
print(model.bent)

x.grid <- seq(min(ss001small$cd4count), max(ss001small$cd4count), length=20)

# -----------------------------------------------------------------------------------------------------
# Figure
png(file = "Figure_HIVvsCD4.png",
    units = "in",
    width = 10,
    height = 5,
    res = 600)

par(mfrow = c(1,3))
par(cex = 1.2)
par(bty = "l")

# PLOT 1
plot(s001small$cd4count, s001small$loghivrna,
	pch = 16,
	xlim = c(0, 500),
	#ylim = c(0,0.01),
	xlab = "CD4 count",
	ylab = "log(HIV vl)",
	main = "robust linear")
abline(cfar,
	lty = "dotted",
	lwd = 2,
	col = "darkred")
mtext(expression(r^2:0.06333), side = 4, col = "darkred")

# PLOT 2
plot(s001small$cd4count, s001small$loghivrna,
	pch = 16,
	xlim = c(0, 500),
	# ylim = c(0,0.01),
	xlab = "CD4 count",
	ylab = "log(HIV vl)",
	main = "broken stick")
lines(x.grid, predict(model, x.grid),
		lty = "dotted",
		lwd = 2,
		col = "darkred")
mtext(round(model$change.point), side = 4, col = "darkred")


# PLOT 3
plot(s001small$cd4count, s001small$loghivrna,
	pch = 16,
	xlim = c(0, 500),
	# ylim = c(0,0.01),
	xlab = "CD4 count",
	ylab = "log(HIV vl)",
	main = "bent cable")
lines(x.grid, predict(model.bent, x.grid),
		lty = "dotted",
		lwd = 2,
		col = "darkred")
mtext(round(model.bent$alpha), side = 4, col = "darkred")

dev.off()
