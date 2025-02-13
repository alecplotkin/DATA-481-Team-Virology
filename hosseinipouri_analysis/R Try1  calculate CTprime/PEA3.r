IdentifyUnreliable = function(Data, Dilution, Display = FALSE){

X <- log2(Dilution);
Width <- length(Data)
Length <- length(Data[,1])

Slope <- c()
LowerBnd <- c()
UpperBnd <- c()
Yintercept <- c()
Rsquare <- c()
for(i in 1:Length){
	Slope[i] <- lm(t(Data[i,3:Width]) ~ X)$coefficients[2]
	LowerBnd[i]<-confint(lm(t(Data[i,3:Width])~X),"X", 0.90)[1]
	UpperBnd[i]<-confint(lm(t(Data[i,3:Width])~X),"X", 0.90)[2]
	Rsquare[i] <- summary(lm(t(Data[i,3:Width]) ~ X))$r.squared  
}

Insignificant <- Data$primer[UpperBnd>0]

PurgedData<-subset(Data, UpperBnd<0);
Slope <- subset(Slope, UpperBnd < 0);
LowerBnd <- subset(LowerBnd, UpperBnd < 0);
Rsquare <- subset(Rsquare, UpperBnd < 0);
UpperBnd <- subset(UpperBnd, UpperBnd < 0);

Keff <- 2^(-1/Slope);
LowerKeff <- 2^(-1/LowerBnd);
UpperKeff <- 2^(-1/UpperBnd);

KeffLength <- UpperKeff-LowerKeff
Unreliable <- PurgedData$primer[KeffLength > quantile(KeffLength, .75)+1.5*IQR(KeffLength)]
PurgedData <- subset(PurgedData, KeffLength < quantile(KeffLength, .75)+1.5*IQR(KeffLength))

if(Display == TRUE){
	plot(Keff[KeffLength < 3], KeffLength[KeffLength < 3], xlab = "Efficiency", ylab = "Length of Transformed CI", main = "Efficiency Vs Transformed CI", pch = 19, cex = 1.7)
	points(Keff[KeffLength > quantile(KeffLength, .75)+1.5*IQR(KeffLength)], KeffLength[KeffLength > quantile(KeffLength, .75)+1.5*IQR(KeffLength)], pch = 19, col = 'red', cex = 1.7) 
}
return(list(Purged = PurgedData, Insignificant = Insignificant, Unreliable = Unreliable))
}


##############################################
ClusterEffs <- function(PurgedData, Dilution, Display = FALSE){
#PurgedData = NewData$Purged
#Dilution = c(1,.2,.04,.008,.0016,.00032)
#Display = TRUE

Length <- length(PurgedData[,1])
Width <- length(PurgedData[1,])
X <- log2(Dilution)

Keff <- c()
for(i in 1:Length){
	Slope <- lm(PurgedData[i,3:Width] ~ X)$coefficients[2]
	Keff[i] <- 2^(-1/Slope)
}
PurgedData.Sorted <- PurgedData[order(Keff),]
Y <- c()
for(i in 3:Width){
    Y <- c(Y, PurgedData.Sorted[,i])}
newX <- rep(X, each = Length)
#Primer <- rep(PurgedData.Sorted$primer, length(Dilution))
Primer <- rep(PurgedData.Sorted[,2], length(Dilution))

Groups <- matrix(nrow=Length, ncol=Length)
BestRsquare <- 0
Groups[,1] <- sample(1,Length, TRUE)
for(k in 2:Length){
 for(i in 1:(k-1)){  
   for(j in 1:(sum(Groups[,k-1]==i)-1)){
	Dummy <- Groups[,k-1]
	Dummy[Groups[,k-1]==i] <- c(rep(i,j), rep(k,sum(Groups[,k-1]==i)-j))
	Rsquare <- summary(lm(Y ~ newX:factor(rep(Dummy,length(Dilution))) + factor(Primer)))$r.squared
	if(Rsquare>BestRsquare){
		BestRsquare <- Rsquare
		Groups[,k]<-Dummy
      }
  }}
 if(k==2)
    p <- anova(lm(Y ~ newX+ factor(Primer)), lm(Y ~ newX:factor(rep(Groups[,k],length(Dilution))) + factor(Primer)))$Pr[2]
 else
    p <- anova(lm(Y ~ newX:factor(rep(Groups[,k-1],length(Dilution))) + factor(Primer)), lm(Y ~ newX:factor(rep(Groups[,k],length(Dilution))) + factor(Primer)))$Pr[2]
 if((p > 0.05)){
    Group <- Groups[,k-1] 
    NumClusts <- k-1
    break; 
 }}

cat("Found ")
cat(NumClusts)
cat(" clusters. \n")

if(Display == TRUE){
  ClusColors <- rep(c("red", "green", "blue", "yellow", "black", "purple", "brown", colors()[79], colors()[450], colors()[81], colors()[209]), 1+NumClusts/12)
  stripchart(sort(Keff), col = "white", xlab = "Efficiency", main = "Primer Efficiency Clusters")
  title(main = "Primer Efficiency Clusters")
  for(i in 1:NumClusts)
  	stripchart(sort(Keff)[Group==i], method = "jitter", jitter=.1, add=TRUE, col = ClusColors[i], pch=".", cex = 9)
}
Group_Unsorted <- Group[rank(Keff)]
Efficiencies <- c()
for(i in 1:NumClusts){
	Efficiencies[Group_Unsorted==i] <- mean(sort(Keff)[Group==i])
}

 return(Efficiencies)
}



#############################################################
AdjustCt <- function(Data, Efficiencies){
	Width <- length(Data)
	Adjusted <- Data
	Adjusted[, 3:Width] <- log2(Efficiencies^Data[,3:Width]) 
	return(Adjusted)
}



