# -------------------------------------------------------------------------------------------------
# this code is to calculate CT' based on CT values
# Efficiencies are calculated from training set with the same set of primers
# Notice that for different primers in training set, 
# the number of CT value measured are different

# -------------------------------------------------------------------------------------------------
# read in libraries
library(xlsx)

# -------------------------------------------------------------------------------------------------
# read data from excel file
# Note the data have to be in the following format
# "Plate"  "Primer"    "dilution1"    ... "dilution 7"
# the highest CT is on the left

data=read.xlsx('a_CT_training_PP1v04182014.xlsx',sheetIndex=1)
head(data)
colnames(data)

# -------------------------------------------------------------------------------------------------
# store primer number = rownumber
length=dim(data)[1]
length
#width=dim(data)[2] #==7

# -------------------------------------------------------------------------------------------------
# missing values in some of the higher folds (righthand side of matrix)
CT5=matrix(nrow=length,ncol=7) # subset with 5 CTs measured
CT6=matrix(nrow=length,ncol=8) # subset with 6 CTs measured
CT7=matrix(nrow=length,ncol=9) # subset with 7 CTs measured
i5=1
i6=1
i7=1
ind5=rep(0,length)
ind6=rep(0,length)
ind7=rep(0,length)

# -------------------------------------------------------------------------------------------------
# Missing values in some of the higher folds (MUST BE righthand side of matrix)
# At most 3 NA's are handles by this code
# This sets up 3 submatrices
for(i in 1:length){
	n=sum(is.na(data[i,]))
	if(n==2){
		CT5[i5,]=as.matrix(data[i,1:7])
		ind5[i5]=i
		i5=i5+1
		}
	if(n==1){
		CT6[i6,]=as.matrix(data[i,1:8])
		ind6[i6]=i
		i6=i6+1
		}
	if(n==0){
		CT7[i7,]=as.matrix(data[i,1:9])
		ind7[i7]=i
		i7=i7+1
		}
	}
	
# deletes empty rows in each submatrix
head(CT5)
CT5=na.omit(CT5)
CT6=na.omit(CT6)
CT7=na.omit(CT7)

# counts empty columns in each submatrix
head(ind5)
ind5=ind5[ind5!=0]
ind6=ind6[ind6!=0]
ind7=ind7[ind7!=0]

# now we have complete matrices with no NAs
head(CT5)


# -------------------------------------------------------------------------------------------------
# use Lock's code to calculate beta1/efficienct related to each primer
# setwd('/Users/Tam/Desktop')
source('PEA3.r')

# Note how the list c() hardcodes the dilutions
par(mfrow = c(3,1))
eff5=ClusterEffs(CT5, c(16,64,256,1024,4096), Display = T)
eff6=ClusterEffs(CT6, c(4,16,64,256,1024,4096), Display = T)
eff7=ClusterEffs(CT7, c(1,4,16,64,256,1024,4096), Display = T)
head(eff5)


# -------------------------------------------------------------------------------------------------
# combine estimated efficiencies together in original order
ind=c(ind5,ind6,ind7)
primer=c(CT5[,2],CT6[,2],CT7[,2])
eff=c(eff5,eff6,eff7)
eff.matrix=cbind(primer,eff)
eff.sorted=eff.matrix[order(ind),]
eff.sorted
head(eff.sorted)


# -------------------------------------------------------------------------------------------------
# summary plot
eff.sorted.frame <- as.data.frame(eff.sorted)
eff.sorted.frame[,"eff"] <- as.character(eff.sorted.frame[,"eff"])
eff.sorted.frame[,"eff"] <- as.numeric(eff.sorted.frame[,"eff"])
summary(eff.sorted.frame)

par(mfrow = c(1,1))
plot(density(eff.sorted.frame[,"eff"] ))


# -------------------------------------------------------------------------------------------------
# calculate new CT'
new.data=read.csv("aa_calcCTPrime_45max_04182014a.csv")
index=(new.data[,8]=='good') & ((new.data[,9]=='PP1'))
new.data=new.data[index,]
dim(new.data)
CT=new.data[,10:25]
scale=c(log2(as.numeric(eff.sorted[,2])))
CT.prime=CT*scale
new.sheet=cbind(eff.sorted[,1],CT.prime)
write.xlsx(new.sheet,file='c_121412_OUTPUT-CT prime_PP1_no60.xlsx', row.names = FALSE)



#--------------------------------------------------------------------------------------
data=read.xlsx('a_CT_training_PP2v04182014.xlsx',sheetIndex=1)
data
length=dim(data)[1]
#width=dim(data)[2] #==7
CT5=matrix(nrow=length,ncol=7) # subset with 5 CTs measured
CT6=matrix(nrow=length,ncol=8) # subset with 6 CTs measured
CT7=matrix(nrow=length,ncol=9) # subset with 7 CTs measured
i5=1
i6=1
i7=1
ind5=rep(0,length)
ind6=rep(0,length)
ind7=rep(0,length)
for(i in 1:length){
	n=sum(is.na(data[i,]))
	if(n==2){
		CT5[i5,]=as.matrix(data[i,1:7])
		ind5[i5]=i
		i5=i5+1
		}
	if(n==1){
		CT6[i6,]=as.matrix(data[i,1:8])
		ind6[i6]=i
		i6=i6+1
		}
	if(n==0){
		CT7[i7,]=as.matrix(data[i,1:9])
		ind7[i7]=i
		i7=i7+1
		}
	}
CT5=na.omit(CT5)
CT6=na.omit(CT6)
CT7=na.omit(CT7)
ind5=ind5[ind5!=0]
ind6=ind6[ind6!=0]
ind7=ind7[ind7!=0]

# use Lock's code to calculate beta1/efficienct related to each primer
# setwd('/Users/Tam/Desktop')
source('PEA3.r')
eff5=ClusterEffs(CT5, c(4096, 1024, 256, 64, 16), Display = T)
eff6=ClusterEffs(CT6, c(4096, 1024, 256, 64, 16, 4), Display = T)
eff7=ClusterEffs(CT7, c(4096, 1024, 256, 64, 16, 4, 1), Display = T)

# combine estimated efficiencies together in original order
ind=c(ind5,ind6,ind7)
primer=c(CT5[,2],CT6[,2],CT7[,2])
eff=c(eff5,eff6,eff7)
eff.matrix=cbind(primer,eff)
eff.sorted=eff.matrix[order(ind),]
eff.sorted


# calculate new CT'
new.data=read.csv("aa_calcCTPrime_45max_04182014a.csv")
index=(new.data[,8]=='good') & ((new.data[,9]=='PP2'))
new.data=new.data[index,]
dim(new.data)
CT=new.data[,10:25]
scale=c(log2(as.numeric(eff.sorted[,2])))
CT.prime=CT*scale
new.sheet=cbind(eff.sorted[,1],CT.prime)
write.xlsx(new.sheet,file='c_121412_OUTPUT-CT prime_PP2_no60.xlsx', row.names = FALSE)

