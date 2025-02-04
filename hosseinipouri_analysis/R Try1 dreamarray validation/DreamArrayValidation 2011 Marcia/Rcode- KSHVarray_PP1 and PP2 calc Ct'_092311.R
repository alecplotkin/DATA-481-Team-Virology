#PP1
source(("PEA3.r"))
DataPP1 <- read.csv("PP1forCtPrime2.csv", header = TRUE)
names(DataPP1)[2] = 'primer'
DataPP1


##perform unreliable-primer check##
Dilution <- c(80, 400, 2000, 10000)
NewData = IdentifyUnreliable(DataPP1, Dilution=c(80, 400, 2000, 10000), Display=TRUE)
NewData

##################################################################

Efficiencies <- ClusterEffs(NewData$Purged, Dilution= c(80, 400, 2000, 10000), Display=TRUE)

##################################################################
#calc new adjusted Cts
Data_AdjustedPP1 = AdjustCt(NewData$Purged,Efficiencies)
Data_AdjustedPP1

write.table(Data_AdjustedPP1, file= "102411Malawi_CTprime_PP1.csv", sep= "\t")

##################################################################
##################################################################
#################################################################

#PP2
source(("PEA3.r"))
DataPP2 <- read.csv("PP2forCtPrime_3.csv", header = TRUE)
names(DataPP2)[2] = 'primer'
DataPP2


##perform unreliable-primer check##
Dilution <- c(312.5, 625, 1250, 2500, 5000, 10000)
NewData = IdentifyUnreliable(DataPP2, Dilution=c(312.5, 625, 1250, 2500, 5000, 10000), Display=TRUE)
NewData

##################################################################

Efficiencies <- ClusterEffs(NewData$Purged, Dilution= c(312.5, 625, 2500, 2500, 5000, 10000), Display=TRUE)

##################################################################
#calc new adjusted Cts
Data_AdjustedPP2 = AdjustCt(NewData$Purged,Efficiencies)
Data_AdjustedPP2

write.table(Data_AdjustedPP2, file= "101411Malawi_CTprime_PP2.csv", sep= "\t")

##################################################################
####################################################################################################################################
##################################################################
#
#calc dCT'
#open Malawi_CTprime_allData.csv with excel, calc geometric mean for reference genes, make new table with primers in columns and samples going down rows.  save as Malawi-fordCT2.csv



Data = read.csv("Malawi-fordCT2.csv", header = TRUE)

#calc dCT for actin minus gene

dCt = Data
for(i in 2:92)
	dCt[,i] = Data[,2]-Data[,i]
	
dCt	

write.table(dCt, file= "050411_Malawi_dCT.csv", sep= "\t")
names(dCt)
