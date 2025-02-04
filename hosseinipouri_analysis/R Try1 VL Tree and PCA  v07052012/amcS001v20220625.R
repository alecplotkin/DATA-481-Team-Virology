# STEP 0:    Read in the libraries ----
library(foreign)
library(MASS)
library(DAAG)
library(gdata)
library(Hmisc)
library(RColorBrewer)
library(FactoMineR)
#library(heatmap.plus)
library(reshape2)
library(ggplot2)
#library(DMwR)
library(klaR)
library(party)
library(SiZer)
library(rpart)
library(magrittr)
library(dplyr)
library(lattice)

# OPTIONS ----
options(digits = 3)

# FUNCTIONS ----

#' plotPNG
#' @param GGplotObject    - An object of type GGplot object
#'                        - An object of type Plot object
#' @param NameOfPlot      The filename of the plot file to be generates
#'                        MUST contain ".png"
#' @return TRUE
#' @export SIDE_EFFECTS a png file SIDE_EFFECTS as the file is printed to disk
#' @examples plotPNG()
#'
plotPNG = function(PlotObject = plot(x = seq(1:10),
                                     y = seq(1:10)),
                   NameOfPlot = "UnitTestFigureFromR.png"){
  png(file = NameOfPlot,
      units = "in",
      width = 5,
      height = 4,
      res = 300)
  print(PlotObject)
  dev.off()
  return(TRUE)
}

theme_dirk <- function(base_size = 22,
                       base_family = "",
                       base_line_size = base_size/22,
                       base_rect_size = base_size/22, time_stamp = FALSE){

  obj <- theme_classic(base_size = base_size,
                       base_family = base_family,
                       base_line_size = base_line_size,
                       base_rect_size = base_rect_size) %+replace%
    theme(
      axis.text = element_text(size = base_size),
      axis.text.x = element_text(vjust = 0.5),
      axis.title.x = element_text(margin = margin(t = .8 * base_size), vjust = 0),
      axis.title.x.top = element_text(margin = margin(b = .8 * base_size), vjust = 1),
      axis.title.y = element_text(margin = margin(r = .8 * base_size), vjust = 1, angle = 90),
      axis.title.y.right = element_text(margin = margin(l = .8 * base_size), vjust = 0, angle = 90),
      strip.text = element_text(size = base_size),
      strip.background = element_rect(colour = "white", fill = "white"),
      panel.grid = element_blank(),
      panel.background = element_rect(colour = "black", fill = "white"),
      plot.caption = element_text(size = rel(.5), hjust = 1, vjust = 0),
      plot.caption.position = "plot",
      complete = T
    )

  if (time_stamp) {
    obj <- list(obj, labs(caption = Sys.time()))
  }

  return(obj)

}

figurefromU54 = function(df_tb_vl = s001small){
vars <- c("CD4",
          "CD4_sqrt",
          "log_viral_load",
          "CopyNumber")
df_sub <- df_tb_vl[, vars]
df_sub$log_KSH <- log10(df_sub$CopyNumber)

df_sub = as_tibble(df_sub)

df_sub = df_sub %>%
  mutate(log_viral_load = if_else(log_viral_load < 1.61, 1, log_viral_load))

ggplot(data = df_sub, aes(x = log_viral_load)) +
  geom_density(col = "darkred", lty = 1) +
  geom_rug() +
  xlab("log10 HIV VL") +
  theme_dirk()

ggplot(data = df_sub, aes(x = log_KSH)) +
  geom_density(col = "darkred", lty = 1) +
  geom_rug() +
  xlab("log10 KSHV VL") +
  theme_dirk()

ggplot(data = df_sub, aes(x = sqrt(CD4))) +
  geom_density(col = "darkred", lty = 1) +
  geom_rug() +
  xlab(expression(sqrt("CD4 cells/µl"))) +
  theme_dirk()

lm <- lm(log_viral_load ~ sqrt(CD4), data = df_sub)
summary(lm)

ggplot(data = df_sub, aes(x = sqrt(CD4), y = log_viral_load)) +
  geom_point(color = "darkgray", size = 3.5) +
  geom_vline(xintercept = sqrt(200), col = "darkgray", lty = 2) +
  geom_vline(xintercept = sqrt(400), col = "darkgray", lty = 2) +
  geom_vline(xintercept = sqrt(600), col = "darkgray", lty = 2) +
  geom_smooth(method = 'lm', col = "darkred", lty = 1) +
  ylab("log10 HIV VL") +
  theme_dirk()

shapiro.test(df_sub$log_viral_load)
shapiro.test(df_sub$CD4)

df_sub$CD4group <- as.factor(cut(df_sub$CD4,
                                 ordered_result = TRUE,
                                 breaks = c(0,100,200,400,800),
                                 labels = c("less100","less200","less400","over400")))
table(df_sub$CD4group)

ggplot(data = df_sub[complete.cases(df_sub),],
       aes(x = log_KSH,
           y = log_viral_load,
           color = CD4group,
           shape = CD4group)) +
  geom_point(size = 4.5) +
  scale_color_manual(values = c("darkred", "darkgray","darkgray", "darkblue"), labels = c("<=100", "100-200", "200-400", ">=400")) +
  scale_shape_manual(values = c(16,16,17,17), labels = c("<=100", "100-200", "200-400", ">=400")) +
  xlim(0,6) +
  ylim(0,6) +
  geom_vline(xintercept = 2.0, col = "black", lty = 2) +
  geom_vline(xintercept = 4.0, col = "black", lty = 2) +
  geom_hline(yintercept = 1.5, col = "black", lty = 2) +
  xlab("log10 KSHV VL") +
  ylab("log10 HIV VL") +
  theme_dirk()

}

# -----------------------------------------------------------------------------------------------------
# STEP 1: Read in the data     		must be in the folder as defined by workspace
#                       							must be at least two file names
# -----------------------------------------------------------------------------------------------------
#
here::here()
rain <- read.delim(here::here("R Try1 VL Tree and PCA  v07052012",
                              "AMC001revised05092012_KT.txt"),
                   header = TRUE,
                   sep = "\t")
summary(rain)


# -----------------------------------------------------------------------------------------------------
# Safe a copy of the data set on disk
# Safe a copy of the data set as variable rainall
# RECOVER with  rain <- rainall
#
rainall <- rain

## write.table(rainall,file = "CompleteDataSet.txt", sep = "\t")


# STEP 2:    Data cleaning

head(rain)
colnames(rain)

# eliminate dulicate records
nrow(rain)
rain <- unique(rain)
nrow(rain)


# -----------------------------------------------------------------------------------------------------
# Use a subset of Variables (columns)
#
rain <- rainall[,c(
  "t_stage",
  "tumor1edema",
  "cpskshv_",
  "KS_nd",
  "orf21",
  "orf36",
  "numberoflessions",
  "Kover70",
  "karnovisky",
  "hivrna",
  "cd4count",
  "cd4under200",
  "gender",
  "age",
  "weightkg",
  "heightcm",
  "bpmmhg",
  "pulsebpm",
  "resp",
  "bmi",
  "bmi_low",
  "weightloss",
  "tumor1oral",
  "tumor1gastro",
  "tumor1viscera",
  "lymphnode")]
summary(rain)


# -----------------------------------------------------------------------------------------------------
# STEP 3:    Data transformation
#
# -----------------------------------------------------------------------------------------------------
#
rain$cpskshv_ <- as.numeric(rain$cpskshv_)
rain$logKS	 <- log10(rain$cpskshv_ + 1)			#allways add +1 for log transformation
rain$KS_nd	<- factor(rain$KS_nd)
levels(rain$KS_nd) <- c("no", "yes")
table(rain$KS_nd)

rain$logHIV	<- log10(rain$hivrna + 1)				#standard
rain$log2HIV	<- log2(rain$hivrna + 1)				#as in the original Mellors et al. 1996

rain$sqrtcd4 <- sqrt(rain$cd4count)
rain$cd4under200	<- factor(rain$cd4under200)
levels(rain$cd4under200) <- c("no", "yes")
table(rain$cd4under200)

rain$bmi_low	<- factor(rain$bmi_low)
levels(rain$bmi_low) <- c("no", "yes")
table(rain$bmi_low)

rain$weightloss	<- factor(rain$weightloss)
levels(rain$weightloss) <- c("no", "yes")
table(rain$weightloss)

rain$lymphnode	<- factor(rain$lymphnode	)
levels(rain$lymphnode) <- c("no", "yes")
table(rain$lymphnode)

rain$numberoflessions	<- factor(rain$numberoflessions)
rain$Kover70	<- factor(rain$Kover70)
levels(rain$Kover70) <- c("no", "yes")
table(rain$Kover70)

# 3 levels of lyticness
rain$lyticscore <- factor(rain$orf36 + rain$orf21)
levels(rain$lyticscore) <- c("0","+","*")
table(rain$lyticscore)

# 2 levels of lyticness
rain$verylytic <- factor(rain$orf36 * rain$orf21)
rain$orf21	<- factor(rain$orf21, ordered = TRUE)
rain$orf36	<- factor(rain$orf36, ordered = TRUE)
table(rain$verylytic)

# 4 levels of lyticness
rain$lyticscorefactor <- interaction(rain$orf36,rain$orf21)
table(rain$lyticscorefactor)

rain$gender	<- factor(rain$gender)
levels(rain$gender) <- c("female", "male")
table(rain$gender)

rain$tumor1oral	<- factor(rain$tumor1oral)
rain$tumor1gastro	<- factor(rain$tumor1gastro)
rain$tumor1viscera	<- factor(rain$tumor1viscera)
rain$tumor1edema 	<- factor(rain$tumor1edema)
rain$t_stage 	<- factor(rain$t_stage)

# this is possible after imputation
# rain <- na.omit(rain)
# rain <- drop.levels(rain)

rain$CD4group <- as.factor(cut(rain$cd4count,
                                 ordered_result = TRUE,
                                 breaks = c(0,100,200,400,800),
                                 labels = c("less100","less200","less400","over400")))
table(rain$CD4group)

## write.csv(rain, "2008_2010_hivvl_kshvvl.csv", row.names = FALSE)

ggplot(data = rain,
       aes(x = logKS,
           y = logHIV,
           color = CD4group,
           shape = CD4group)) +
  geom_point(size = 4.5) +
  scale_color_manual(values = c("darkred", "darkgray","darkgray", "darkblue"), labels = c("<=100", "100-200", "200-400", ">=400")) +
  scale_shape_manual(values = c(16,16,17,17), labels = c("<=100", "100-200", "200-400", ">=400")) +
  xlim(0,6) +
  ylim(0,6) +
  geom_vline(xintercept = 2.0, col = "black", lty = 2) +
  geom_vline(xintercept = 4.0, col = "black", lty = 2) +
  geom_hline(yintercept = 1.5, col = "black", lty = 2) +
  xlab("log10 KSHV VL") +
  ylab("log10 HIV VL") +
  theme_dirk()
# converting to an optimized data set
summary(rain)

# the working dataset
malawi <- rain
## readr::write_csv(malawi, paste0("AMCS001",
##                          stringr::str_sub(lubridate::today()),
##                          ".csv"))

# -----------------------------------------------------------------------------------------------------
# STEP5 :  test for Normality:  logHIV, logKS, sqrtCD4
#
colnames(malawi)
malawi <- malawi[,c("logHIV",
                    "log2HIV",
                    "cpskshv_",
                    "logKS",
                    "cd4count",
                    "sqrtcd4",
                    "cd4under200",
                    "KS_nd")]
malawi <- na.omit(malawi)
malawi <- drop.levels(malawi)

# adding the level of detection ceiling at 200 copies per ml
malawi[(malawi$cpskshv_ < 200),"logKS"] <- 2.30
summary(malawi)

attach(malawi)

# Test for Normality
IQR(logHIV, na.rm = TRUE)
fivenum(logHIV, na.rm = TRUE)
shapiro.test(logHIV)

IQR(logKS, na.rm = TRUE)
fivenum(logKS, na.rm = TRUE)
shapiro.test(logKS)
shapiro.test(malawi[malawi$logK > 2.3,"logKS"])

IQR(cd4count, na.rm = TRUE)
fivenum(cd4count, na.rm = TRUE)
shapiro.test(cd4count)
shapiro.test(sqrtcd4)


# -----------------------------------------------------------------------------------------------------
# Figure
png(file = here::here("R Try1 VL Tree and PCA  v07052012",
                      "Norm.png"),
    units = "in",
    width = 10,
    height = 6,
    res = 600)

par(mfrow = c(1,3))
par(cex = 1.2)

qqnorm(logHIV,
pch = 19,
cex = 1.0,
col = "darkgray",
main = "log HIV /ml"
)
qqline(logHIV,
       col = "darkred",
       lwd = 2,
       lty = 3)
histSpike(logHIV,
          add = TRUE,
          col = "darkblue",
          frac = 0.2,
          lwd = 2,
          side = 2)

qqnorm(logKS,
pch = 19,
cex = 1.0,
col = "darkgray",
main = "log KSHV /ml"
)
qqline(logKS,
       col = "darkred",
       lwd = 2,
       lty = 3)
histSpike(logKS,
          add = TRUE,
          col = "darkblue",
          frac = 0.2,
          lwd = 2,
          side = 2)

qqnorm(sqrtcd4,
       pch = 19,
       cex = 1.0,
       col = "darkgray",
       main = expression(sqrt("CD4 cells /µl"))
       )
qqline(sqrtcd4,
       col = "darkred",
       lwd = 2,
       lty = 3)
histSpike(sqrtcd4,
          add = TRUE,
          col = "darkblue",
          frac = 0.2,
          lwd = 2,
          side = 2)

dev.off()


# -----------------------------------------------------------------------------------------------------
# alternative plot ecdf
par(mfrow = c(1,3))
par(cex = 1.3)

plot(
  ecdf(logHIV),
  main = "",
  col = "darkred",
  cex = 0.5,
  xlab = "log10 HIV RNA copies /ml")
histSpike(logHIV,
          add=TRUE,
          col = "darkblue",
          frac = 0.2,
          lwd = 2)

plot(
  ecdf(logKS),
  main = "",
  col = "darkred",
  cex = 0.5,
  xlab = "log10 KSHV DNA copies /ml")
histSpike(logKS,
          add = TRUE,
          col = "darkblue",
          frac = 0.2,
          lwd = 2)

plot(
  ecdf(sqrtcd4),
  main = "",
  col = "darkred",
  cex = 0.5,
  xlab = expression(sqrt("CD4 cells/µl"))
  )
histSpike(sqrtcd4,
          add = TRUE,
          col = "darkblue",
          frac = 0.2,
          lwd = 2)

detach(malawi)


# -----------------------------------------------------------------------------------------------------
# STEP 4:    Introducing CD4 segments and HIV segments
#
png(file = here::here("R Try1 VL Tree and PCA  v07052012",
                      "Discretization.png"),
    units = "in",
    width = 8,
    height = 8,
    res = 600)


par(mfrow = c(2,2))
par(cex = 1.1)

par(bty = "o")
hist(malawi$sqrtcd4,
     col = "gray",
     main = "Distribution of  CD4 count",
     xlab = expression(sqrt("CD4 cells/µl")),
     breaks = 7,
     labels = TRUE,
     ylim = c(0,30),
     include.lowest = TRUE)

plot(
  ecdf(malawi$sqrtcd4),
  main = "",
  col = "darkred",
  cex = 0.5,
  xlab = expression(sqrt("CD4 cells/µl"))
                  )
histSpike(malawi$sqrtcd4,
          add=TRUE,
          col = "darkblue",
          frac = 0.2,
          lwd = 2)

hist(malawi$logHIV,
     col = "gray",
     main = "Distribution of  HIV Vl",
     xlab = "log10 HIV RNA copies /ml",
     breaks = 7,
     labels = TRUE,
     ylim = c(0,30),
     include.lowest = TRUE
)

plot(
  ecdf(malawi$logHIV),
  main = "",
  col = "darkred",
  cex = 0.5,
  xlab = "log10 HIV RNA copies /ml")
histSpike(malawi$logHIV,
          add=TRUE,
          col = "darkblue",
          frac = 0.2,
          lwd = 2)


dev.off()


# NOTE	this modifies rain ----

# CD4
# see also CD4under200
#labels = c("less150","less300","less600","over600")
rain$CD4group <- cut(rain$cd4count,
                     ordered_result = TRUE,
                     breaks = c(0,150,300,600,800),
                     labels = c("less150","less300","less600","over600"))
table(rain$CD4group )

# HIV
#labels = c("less4","less4.7","less5.2","over5.2")
rain$HIVgroup4 <- cut(rain$logHIV,
                      ordered_result = TRUE,
                      breaks = 4,
                      labels = c("less4","less4.7","less5.2","over5.2"))
table(rain$HIVgroup4)

#labels = c("low","high")
rain$HIVunder4.4 <- cut(rain$logHIV,
                        ordered_result = TRUE,
                        breaks = c(0,4.4,6),
                        labels = c("yes","no"))
table(rain$HIVunder4.4)

# -----------------------------------------------------------------------------------------------------
# standard linear regression
s001small <- rain
summary(s001small)

cfar <- lm(logHIV ~ sqrtcd4, data = s001small)
anova(cfar)
summary(cfar)

par(mfrow = c(1,2))
par(cex=1.2)
plot(s001small$sqrtcd4,
     s001small$logHIV,
     pch = 16,
     xlim = c(0, sqrt(800)),
     #ylim = c(0,0.01),
     xlab = expression(sqrt("CD4 cells/µl")),
     ylab = "log(HIV VL)",
     main = "AMC S001")
abline(cfar,
	lty = "dotted",
	lwd = 1,
	col = "darkred")
mtext(expression(r^2:0.00918),
      side = 4,
      col = "darkred")
abline(v = 24,
       col = "darkblue",
       lwd = 1, lty = 22.4)

plot(cfar,2,
     xlim = c(-2,2))
s001small[13,]
s001small[45,]

# Figure ggplot version  ----

s001small = tibble::tibble(s001small)
s001small = s001small %>%
  dplyr::select(cd4count, sqrtcd4, cpskshv_, logHIV) %>%
  dplyr::mutate(CD4 = cd4count) %>%
  dplyr::mutate(CD4_sqrt = sqrtcd4) %>%
  dplyr::mutate(log_viral_load = logHIV) %>%
  dplyr::mutate(CopyNumber = cpskshv_)

##descriptive stats for copy number
s001small %>%
    mutate(log10cp = log10(CopyNumber)) %>%
    mutate(mean = mean(log10cp, na.rm = TRUE), sd = sd(log10cp, na.rm = TRUE)) %>%
    select(mean, sd) %>% unique()


png(file = here::here("R Try1 VL Tree and PCA  v07052012","S001HIVvsKSHV.png"),
    units = "in",
    width =8,
    height = 10,
    res = 300)

print(figurefromU54(df_tb_vl = s001small))

dev.off()



# broken stick regression ----

# variants of broken stick regression of CD4 < 500
ss001small <- subset(s001small, s001small$sqrtcd4< 22.4)
ss001small <- drop.levels(ss001small)
summary(ss001small)
ss001small <- drop.levels(ss001small)
ss001small <- subset(ss001small, ss001small$cd4count >= 10)
ss001small <- drop.levels(ss001small)

cfar <- lm(logHIV~ sqrtcd4,
           data = ss001small)
summary(cfar)

model <- piecewise.linear(ss001small$sqrtcd4, ss001small$logHIV, CI = FALSE)
print(model)

model.bent <- bent.cable(ss001small$sqrtcd4, ss001small$logHIV, grid.size=20)
print(model.bent)

x.grid <- seq(min(ss001small$sqrtcd4),
              max(ss001small$sqrtcd4),
              length=20)

# -----------------------------------------------------------------------------------------------------
# Figure
png(file = here::here("R Try1 VL Tree and PCA  v07052012",
                      "Figure_HIVvsCD4.png"),
    units = "in",
    width = 10,
    height = 5,
    res = 600)

par(mfrow = c(1,3))
par(cex = 1.2)
par(bty = "l")

# PLOT 1
plot(s001small$sqrtcd4, s001small$logHIV,
	pch = 16,
	xlim = c(0, sqrt(500)),
	#ylim = c(0,0.01),
	xlab = expression(sqrt("CD4 cells/µl")),
	ylab = "log(HIV vl)",
	main = "robust linear")
abline(cfar,
	lty = "dotted",
	lwd = 2,
	col = "darkred")
mtext(expression(r^2:0.0636), side = 4, col = "darkred")

# PLOT 2
plot(s001small$sqrtcd4, s001small$logHIV,
	pch = 16,
	xlim = c(0, sqrt(500)),
	# ylim = c(0,0.01),
	xlab = expression(sqrt("CD4 cells/µl")),
	ylab = "log(HIV vl)",
	main = "broken stick")
lines(x.grid, predict(model, x.grid),
		lty = "dotted",
		lwd = 2,
		col = "darkred")
mtext(round(model$change.point), side = 4, col = "darkred")


# PLOT 3
plot(s001small$sqrtcd4, s001small$logHIV,
	pch = 16,
	xlim = c(0, sqrt(500)),
	# ylim = c(0,0.01),
	xlab = expression(sqrt("CD4 cells/µl")),
	ylab = "log(HIV vl)",
	main = "bent cable")
lines(x.grid, predict(model.bent, x.grid),
		lty = "dotted",
		lwd = 2,
		col = "darkred")
mtext(round(model.bent$alpha), side = 4, col = "darkred")

dev.off()



# -----------------------------------------------------------------------------------------------------
# STEP5 :  Categorical 2x2 comparisons
#

malawi <- rain
summary(malawi)
colnames(malawi)

fisher.test(malawi$KS_nd, malawi$cd4under200)
plot(malawi$KS_nd,
     malawi$cd4under200)
table (malawi$KS_nd, malawi$cd4under200)

fisher.test(malawi$t_stage, malawi$KS_nd)

table(malawi$KS_nd, malawi$CD4group)
fisher.test(malawi$KS_nd, malawi$CD4group)

#labels = c("low","high")
malawi$CD400low400 <- cut(malawi$cd4count,
                          ordered_result = F,
                          breaks = c(0,400,600),
                          labels = c("yes","no"))
table(malawi$CD400low400)

# -----------------------------------------------------------------------------------------------------
# STEP 6:    bivariate CD4, HIV, KSHV correlations
#


malawi <- rain
summary(malawi)


# SIMPLE
malawi.lm <- lm(logHIV ~ cd4count,
                data = subset(malawi, malawi$cd4count < 600))
summary(malawi.lm)
malawi.lm <- rlm(logHIV ~ cd4count,
                 data = subset(malawi, malawi$cd4count < 600))
summary(malawi.lm)


# VARIANT 1:		Using the parameters of Mellors 1996: log2
malawi.Mellors <- lm(log2HIV ~ cd4count,
                     data = subset(malawi, malawi$cd4count < 600))
summary(malawi.Mellors)
malawi.Mellors <- rlm(logHIV ~ cd4count,
                      data = subset(malawi, malawi$cd4count < 600))
summary(malawi.Mellors)

# VARIANT2: 		VARIANT 1 & Restrict to CD4 < 600
malawiLowCD4 <- subset(malawi,
                       malawi$cd4count < 600)
malawiLowCD4 <- na.omit(malawiLowCD4)
malawiLowCD4 <- drop.levels(malawiLowCD4)
summary(malawiLowCD4)
cor(malawiLowCD4$log2HIV,
    malawiLowCD4$cd4count,
    method = "spearman")
cor.test(malawiLowCD4$log2HIV,
         malawiLowCD4$cd4count,
         method = "spearman")

# VARIANT 3:		now using sqrt(CD4 count)
cor(malawiLowCD4$log2HIV,
    malawiLowCD4$sqrtcd4, method = "spearman")
cor.test(malawiLowCD4$log2HIV,
         malawiLowCD4$sqrtcd4, method = "spearman")

malawi.Mellorssqrt <- lm(log2HIV ~ sqrtcd4, data = malawiLowCD4)
summary(malawi.Mellorssqrt)

# Figure
par(mfrow = c(2,2))
plot(malawi.Mellorssqrt)

plot(malawiLowCD4$sqrtcd4,malawiLowCD4$log2HIV,
     col = "darkgray",
     pch = 16, main = "",
     xlab = "sqrt(CD4 cells / µl)",
     ylab = "lg HIV cps/ml",
     cex =1.4)
abline(v = sqrt(150) , col = "darkblue", lty = 2)
abline(v = sqrt(300), col = "darkblue", lty = 2)
abline(v = sqrt(600), col = "darkblue", lty = 2)
abline(malawi.Mellorssqrt, col = "darkred", lty = 1)


# VARIANT 4
colnames(malawi)
malawi.lm4 <- rlm(sqrtcd4 ~ log2HIV, data = malawi)
summary(malawi.lm4)

colnames(malawi)
malawi.lm5 <- lm(sqrtcd4 ~ HIVunder4.4 + log2HIV + logKS, data = malawi)
summary(malawi.lm5)

malawi.lm5 <- lm(logKS ~ sqrtcd4 + log2HIV, data = malawi)
summary(malawi.lm5)


# STEP 6:    Figure VL
malawi <- rain
colnames(malawi)

malawi.lm <- lm(logHIV ~ cd4count,
                data = subset(malawi, malawi$cd4count < 600)
                )

png(file = here::here("R Try1 VL Tree and PCA  v07052012",
                      "VL.png"),
    units = "in",
    width = 8,
    height = 10)

attach(malawi)

par(mfrow = c(2,1))
par(cex =1.5)
par(bty="o")

plot(cd4count,logHIV, col = "darkgray", pch = 16, main = "", xlab = "CD4 cells / µl", ylab = "lg HIV cps/ml", ylim = c(3.0,6), cex =1.4)
abline(v = 150 , col = "darkblue", lty = 2)
abline(v = 300, col = "darkblue", lty = 2)
abline(v = 600, col = "darkblue", lty = 2)
abline(malawi.lm, col = "darkred", lty = 1)


low <- subset(malawi, malawi$CD4group == "less150")
par(bty="n")
plot(low$logKS, low$logHIV,
col = "darkred",
pch = 1,
main = "",
xlab = "",
ylab = "",
ylim = c(3.0,6),
xlim = c(2,5),
xaxt = "n",
yaxt = "n",
cex = 1.0)


par(new = TRUE)
medium <- subset(malawi, malawi$CD4group == "less300")
par(bty="n")
plot(medium$logKS, medium$logHIV,
col = "darkblue",
pch = 1,
main = "",
xlab = "",
ylab = "",
ylim = c(3.0,6),
xlim = c(2,5),
xaxt = "n",
yaxt = "n",
cex = 1.0)

par(new = TRUE)
high <- subset(malawi, malawi$CD4group == "less600")
par(bty="n")
plot(high$logKS, high$logHIV,
col = "darkgreen",
pch = 3,
main = "",
xlab = "",
ylab = "",
ylim = c(3.0,6),
xlim = c(2,5),
xaxt = "n",
yaxt = "n",
cex = 1.0)

par(new = TRUE)
highhigh <- subset(malawi, malawi$CD4group == "over600")
par(bty="o")
plot(highhigh$logKS, highhigh$logHIV,
col = "black",
pch = 2,
main = "",
xlab = "lg KSHV cps/ml",
ylab = "lg HIV cps/ml",
ylim = c(3.0,6),
xlim = c(2,5),
xaxt = "s",
yaxt = "s",
cex = 1.0)


abline(v = 2.3, col = "black", lty = 2, lwd =1)
abline(h = 4.4, col = "black", lty = 2, lwd =1)

detach(malawi)

dev.off()


# --- STEP 7:    Unsupervised models
#

malawi <- rain
summary(malawi)

# -----------------------------------------------------------------------------------------------------
# remove all NAs
#
malawi <- na.omit(malawi)
malawi<- drop.levels(malawi)

# -----------------------------------------------------------------------------------------------------
# adding the level of detection ceiling at 200 copies per ml
#
malawi[(malawi$logKS < 2.3),"logKS"] <- 2.30
summary(malawi)

# -----------------------------------------------------------------------------------------------------
# select the number of factors/variables
#
malawi <- malawi[,c(

# quantiative measures
	"logKS",
	"logHIV",
	"sqrtcd4",
	"age",
	"weightkg",
	"heightcm",
	"pulsebpm",
	"resp",
	"bmi",

# kinase stage
	"lyticscore","verylytic","lyticscorefactor",

# qualitative factores
	"t_stage",
	"tumor1edema",
	"KS_nd",
	"orf21","orf36",
	"numberoflessions",
	"Kover70",
	"cd4under200","HIVunder4.4",
	"gender","bpmmhg","bmi_low","weightloss",
	"tumor1oral", "tumor1viscera","lymphnode")]

colnames(malawi)
ncol(malawi)


# -----------------------------------------------------------------------------------------------------
# PCA
#
summary(malawi)
colnames(malawi[,1:10])
ncol(malawi)

malawi.pca <- PCA(malawi, scale.unit = T, ncp = 5, graph = F, quanti.sup = 4:9, quali.sup = 10:ncol(malawi))

png(file = here::here("R Try1 VL Tree and PCA  v07052012",
      "PCAvar.png"),
    units = "in",
    width = 10,
    height = 8,
    res = 600)

par(mfrow = c(2,2))
colnames(malawi)
malawi.colorcode <- 11

plot(malawi.pca, cex=0.8, choix ="ind", axes = c(1,2), habillage = malawi.colorcode, invisible = "quali")
plot.PCA(malawi.pca, cex=0.8, choix ="ind", axes = c(1,3), habillage = malawi.colorcode, invisible = "quali")
plot.PCA(malawi.pca, cex=0.8, choix ="ind", axes = c(2,3), habillage = malawi.colorcode, invisible = "quali")

plot.PCA(malawi.pca, cex=0.8, choix ="var", title = "", lim.cos2.var = 0.1)

dev.off()

malawi[12,]
malawi[54,]
malawi[64,]

dimdesc(malawi.pca)
barplot(malawi.pca$eig[,1],
        main = "Eigenvalues",
        names.arg = paste ("Dim", 1:nrow(malawi.pca$eig), sep = ""))


# Building a tree cluster on top of it
malawi.hcpc = HCPC(malawi.pca)


# Using the standard R function to derive loadings and numbers
summary(malawi)
colnames(malawi[,1:10])
malawi.pca <- princomp(scale(malawi[,c("logKS","logHIV","sqrtcd4")]))
summary(malawi.pca)
plot(malawi.pca)
malawi.pca$loadings


# -----------------------------------------------------------------------------------------------------
# follow-up bivariate analyses
#
summary(malawi)

malawi.bi <- lm(logKS~ lyticscorefactor +sqrtcd4 + logHIV + numberoflessions + Kover70 + tumor1viscera,
data = malawi)
summary(malawi.bi)
anova(malawi.bi)
cor(malawi$logKS, malawi$bmi, method = "spearman")
cor.test(malawiLowCD4$log2HIV, malawiLowCD4$cd4count, method = "spearman")


# -----------------------------------------------------------------------------------------------------
#  LDA
#
colnames(malawi)

# explicit lab parameters only
malawi.lda <- lda(lyticscore ~ logKS + logHIV + sqrtcd4, data = malawi, method = "t")
malawi.lda
par(bty = "o")
par(cex = 1.5)
plot(malawi.lda, col = "darkred", main ="all")


# -----------------------------------------------------------------------------------------------------
#  QDA
#
colnames(malawi[,1:10])
partimat(lyticscore ~ logKS + logHIV + sqrtcd4, data = malawi, method = "qda", plot.matrix = TRUE)

# -----------------------------------------------------------------------------------------------------
#
# TREEs
#
# -----------------------------------------------------------------------------------------------------
#  regression tree 1: lyticscorefactor - logKS + logHIV + sqrtcd4
#
colnames(malawi[,1:10])

malawi.rpart <- rpart(lyticscorefactor ~ logKS + logHIV + sqrtcd4, data = malawi, cp = 0.001)
summary(malawi.rpart)


par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(malawi.rpart) # visualize cross-validation results
printcp(malawi.rpart)

# now pruning
malawi.rpart <- prune.rpart(malawi.rpart, cp = 0.033)
summary(malawi.rpart)


par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(malawi.rpart) # visualize cross-validation results
printcp(malawi.rpart)

post(malawi.rpart)


png(file =  here::here("R Try1 VL Tree and PCA  v07052012",
      "TREE_lyticscorefactor.png"),
    units = "in", width = 6, height = 10, res = 600)

par(mfrow=c(1,1))
plot(malawi.rpart)
text(malawi.rpart, use.n=TRUE, cex = 1.3, col = "darkblue")

dev.off()


# Now using the party package
colnames(malawi)
malawi.ctree <- ctree(lyticscorefactor ~ logKS + logHIV + sqrtcd4, data = na.omit(malawi))
plot(malawi.ctree)
summary(malawi.ctree)


#
colnames(malawi)
malawi.rpart2 <- rpart(lyticscorefactor ~ logKS + logHIV + sqrtcd4, data = subset(malawi, malawi$cd4under200 == "yes"), cp = 0.001)
summary(malawi.rpart2)
printcp(malawi.rpart2)
plot(malawi.rpart2)
text(malawi.rpart2, use.n=TRUE, cex = 1, col = "darkblue")

# -----------------------------------------------------------------------------------------------------
#  regression tree 3: lyticscore - logKS + logHIV + sqrtcd4
#
colnames(malawi)
malawi.rpart3 <- rpart(lyticscore ~ logKS + logHIV + sqrtcd4, data = malawi, cp = 0.001)
summary(malawi.rpart3)

png(file = here::here("R Try1 VL Tree and PCA  v07052012","TREE.png"),
    units = "in", width = 6, height = 10, res = 600)

plot(malawi.rpart3)
text(malawi.rpart3,
     use.n=TRUE,
     cex = 1.5,
     col = "darkblue")

dev.off()

# -----------------------------------------------------------------------------------------------------
#  regression tree 3: lyticscore - logKS + logHIV + sqrtcd4
#
summary(malawi)

malawi.rpart <- rpart(lyticscore ~ logKS + logHIV + sqrtcd4, data = subset(malawi, malawi$t_stage == 1))
summary(malawi.rpart)
plot(malawi.rpart)
text(malawi.rpart, use.n=TRUE, cex = 0.5, col = "darkgreen")


# -----------------------------------------------------------------------------------------------------
#  regression tree 4: lyticscore - logKS + logHIV + sqrtcd4
#
summary(malawi)

malawi.rpart <- rpart(lyticscore ~ logKS + logHIV + sqrtcd4, data = subset(malawi, malawi$ tumor1oral  == "0"))
summary(malawi.rpart)
plot(malawi.rpart)
plotcp(malawi.rpart)
printcp(malawi.rpart)
text(malawi.rpart, use.n=TRUE, cex = 0.8, col = "black")


# -----------------------------------------------------------------------------------------------------
# further subsetting
#
# -----------------------------------------------------------------------------------------------------
malawi <- rain
summary(malawi)

# Where are the NA's
apply(malawi,2, function(x) sum (is.na(x)))
nacases <- apply(malawi,1, function(x) sum (is.na(x)))
malawi[nacases > 0,]
malawi <- na.omit(malawi)
malawi <- drop.levels(malawi)
summary(malawi)
nrow(malawi)

# No oral lesions
# malawi <- subset(malawi, (malawi$tumor1oral != "1"))

# No visceral involvement
# malawi <- subset(malawi, (malawi$tumor1viscera != "1"))

# No assay failures for KSHV
# malawi <- subset(malawi, (malawi$KS_nd == "no")

#  template
malawi <-malawi[,c(

# quantiative measures
	"logKS",
	"logHIV",
	"sqrtcd4",
	"age", "weightkg","heightcm",
	"pulsebpm","resp","bmi",

# kinase stage
	"lyticscore","verylytic","lyticscorefactor",

# qualitative factores
	"t_stage",
	"KS_nd",
	"tumor1edema",
	"orf21","orf36",
	"numberoflessions",
	"Kover70",
	"cd4under200","HIVunder4.4","CD4group","HIVgroup4",
	"gender","bpmmhg","bmi_low","weightloss",
	"tumor1oral","tumor1viscera",
	"lymphnode")]


#  slection
malawi <-malawi[,c(

# quantiative measures
	"logKS",
	"logHIV",
	"sqrtcd4",

# kinase stage
	"lyticscorefactor",

# qualitative factores
	"tumor1edema",
	"tumor1oral","tumor1viscera"
	)]




colnames(malawi)
ncol(malawi)
par(mfrow=c(1,1))
splom(malawi)

malawi.rpart <- rpart(lyticscorefactor ~ ., data = malawi)
summary(malawi.rpart)
plot(malawi.rpart)
text(malawi.rpart, use.n=TRUE, cex = 1.0, col = "darkblue")

par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(malawi.rpart) # visualize cross-validation results
printcp(malawi.rpart)

malawi.rpart <- prune.rpart(malawi.rpart, cp = 0.04)
par(mfrow=c(1,2)) # two plots on one page
rsq.rpart(malawi.rpart) # visualize cross-validation results
printcp(malawi.rpart)


malawi$CD4group <- as.factor(cut(malawi$cd4count,
                                 ordered_result = TRUE,
                                 breaks = c(0,100,200,400,800),
                                 labels = c("less100","less200","less400","over400")))
table(malawi$CD4group)

ggplot(data = malawi[complete.cases(malawi),],
       aes(x = logKS,
           y = logHIV,
           color = CD4group,
           shape = CD4group)) +
  geom_point(size = 4.5) +
  scale_color_manual(values = c("darkred", "darkgray","darkgray", "darkblue"), labels = c("<=100", "100-200", "200-400", ">=400")) +
  scale_shape_manual(values = c(16,16,17,17), labels = c("<=100", "100-200", "200-400", ">=400")) +
  xlim(0,6) +
  ylim(0,6) +
  geom_vline(xintercept = 2.0, col = "black", lty = 2) +
  geom_vline(xintercept = 4.0, col = "black", lty = 2) +
  geom_hline(yintercept = 1.5, col = "black", lty = 2) +
  xlab("log10 KSHV VL") +
  ylab("log10 HIV VL") +
    theme_dirk()

