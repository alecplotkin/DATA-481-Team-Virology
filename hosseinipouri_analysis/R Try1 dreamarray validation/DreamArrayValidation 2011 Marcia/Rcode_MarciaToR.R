
source(("PEA3.r"))
DataPP1 <- read.csv("MarciaToR.csv", header = TRUE)
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

write.table(Data_AdjustedPP1, file= "MarciaToR-AdjustedCT'.csv", sep= "\t")

##################################################################
##################################################################
#