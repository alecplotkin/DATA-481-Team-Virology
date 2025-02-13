# -----------------------------------------------------------------------------------------------------
# read in libraries
# -----------------------------------------------------------------------------------------------------
library(foreign)
library(gdata)
library(MASS)
library(DAAG)
library(RColorBrewer)
library(psych)
library(FactoMineR)
library(heatmap.plus)
library(ggplot2)


# -----------------------------------------------------------------------------------------------------
# read in the STATA data
# -----------------------------------------------------------------------------------------------------
s001 <- read.dta("AMC 001 for AMC stats stata verion 9.dta")
summary(s001)
head(s001)
nrow(s001)


# -----------------------------------------------------------------------------------------------------
# lytic index
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

s001$tumor1oral<- factor(s001$tumor1oral)
s001$tumor1gastro<- factor(s001$tumor1gastro)
s001$tumor1viscera<- factor(s001$tumor1viscera)
s001$gender <- factor(s001$gender, labels = c("female","male"))

s001$orf36<- factor(s001$orf36, labels = c("no","yes"))
s001$Kover70<- factor(s001$Kover70, labels = c("no","yes"))
s001$cd4under200<- factor(s001$cd4under200, labels = c("no","yes"))

s001$bmi_low<- factor(s001$bmi_low, labels = c("no","yes"))
s001$KS_nd<- factor(s001$KS_nd, labels = c("no","yes"))
s001$orf21<- factor(s001$orf21, labels = c("no","yes"))

summary(s001)

# -----------------------------------------------------------------------------------------------------
# composite factor
# -----------------------------------------------------------------------------------------------------
s001$lytic<- interaction(s001$orf21, s001$orf36)

# -----------------------------------------------------------------------------------------------------
# eliminate incomplete records
# -----------------------------------------------------------------------------------------------------
s001 <- na.omit(s001)
s001 <- drop.levels(s001)
nrow(s001)
head(s001)

summary(s001)



# -----------------------------------------------------------------------------------------------------
# eliminate dulicate records
# -----------------------------------------------------------------------------------------------------
nrow(s001)
s001 <- unique(s001)
nrow(s001)
summary(s001)


# -----------------------------------------------------------------------------------------------------
# subset : drop out rare factors
# -----------------------------------------------------------------------------------------------------

# no visceral involvement
s001small <- subset(s001, s001$tumor1viscera != 1)
s001small <- drop.levels(s001small)
summary(s001small)

# no oral involvement
s001small <- subset(s001small, s001small$tumor1oral != 1)
s001small <- drop.levels(s001small)
summary(s001small)

head(s001small)

# -----------------------------------------------------------------------------------------------------
# subset: factors only

s001factor <- s001small[,c("t_stage", "age_cat","tumor1edema","gender","weightloss","Kover70","lymphnode","cd4under200","bmi_low","KS_nd","lyticindex", "cpskshv_")]
summary(s001factor)
rownames(s001factor) <- s001small$pid 
colnames(s001factor)
ncol(s001factor)


# -----------------------------------------------------------------------------------------------------
# MCA analysis
res.mca <- MCA(s001factor, quanti.sup = c(12))
plot.MCA(res.mca, invisible=c("var"))
plot.MCA(res.mca, invisible=c("ind"),habillage = "quali", title = "AMC S001")
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

tapply(cd4count, lymphnode, median)
wilcox.test(cd4count ~ lymphnode)
wilcox.test(cd4count ~ orf26)
wilcox.test(cd4count ~ orf21)

table(lyticindex, t_stage)
fisher.test(t_stage,lyticindex, 
            or = 1, 
            alternative = "two.sided",
            conf.int = TRUE, 
            conf.level = 0.95)

summary(aov(cd4count ~ lyticindex,projections = TRUE))
summary(aov(hivrna ~ lyticindex,projections = TRUE))
summary(aov(cpskshv_ ~ lyticindex,projections = TRUE))
summary(aov(age ~ lyticindex,projections = TRUE))
