#-- this script reads in the same data as the variance decomposition --#
#-- coxme models --#

filepath <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/"

relocdata <- read.csv(paste(filepath,
											 "RelocsWithNetworkMeasures/FullEweDataAllSummerRelocs_MinEdge.1_18Sept2013.csv", sep = ""), header = T)
relocdata$component.ind <- paste(relocdata$Pop, "_", relocdata$Year, "_",
																 relocdata$res.component, sep = "")
relocdata$popear.ind <- paste(relocdata$Pop, "_", relocdata$Year, sep = "")

lambdata <- read.csv(paste(filepath,
													 "CleanLambSurvData/FullEweRelocsLambSurvDat_MinEdge.1_allewerelocs_18Sept2013.csv",
													 sep = ""), header = T)

lambdata$LastYearKnownCompMort <- lambdata$LastYearLambStatus <-
	lambdata$LastYearCompLambDiedOrNoLamb <- lambdata$ThisYearCompKnownMort <-
		lambdata$ThisYearCompLambDiedOrNoLamb <-
			lambdata$ThisPopyrLambDiedOrNoLamb <- lambdata$ThisPopyrKnownMort <- rep(NA, dim(lambdata)[1])

for(i in 1:dim(lambdata)[1]){
	ewedat <- subset(relocdata, as.character(EWEID) ==
																	as.character(lambdata$EWEID)[i] &
																	as.numeric(as.character(Year)) ==
																	as.numeric(as.character(lambdata$Year[i])) -
																	1)[1, ] 
	EwesLastYearComponent <- as.character(ewedat$component.ind)
	if(length(EwesLastYearComponent) == 0){
		lambdata$LastYearKnownCompMort[i] <- lambdata$LastYearLambStatus[i] <-
			lambdata$LastYearCompLambDiedOrNoLamb[i] <- NA
	} else {
	LastYearCompEwes <- levels(factor(subset(relocdata,
											 as.character(component.ind) ==
											 as.character(EwesLastYearComponent))$EWEID))
	LastYearCompLambs <- subset(lambdata, as.numeric(as.character(Year)) == as.numeric(as.character(lambdata$Year[i])) - 1 & as.character(EWEID) %in% LastYearCompEwes)
	lambdata$LastYearKnownCompMort[i] <- sum(LastYearCompLambs$CENSOR2) /
	dim(LastYearCompLambs)[1]
	lambdata$LastYearCompLambDiedOrNoLamb[i] <- sum(LastYearCompLambs$CENSOR2) /
	length(LastYearCompEwes)
	lambdata$LastYearLambStatus[i] <- ifelse(as.character(ewedat$HasLamb) ==
																					 "NoLamb", "NoLamb",
																					 ifelse(ewedat$CENSOR2 == 0,
																									"LambDied",
																									ifelse(ewedat$CENSOR2 == 1,
																												 "LambSurvived", NA))) 
	}

	ewedatnow <- subset(relocdata, as.character(EWEID) ==
																	as.character(lambdata$EWEID)[i] &
																	as.numeric(as.character(Year)) ==
																	as.numeric(as.character(lambdata$Year[i])))[1, ] 
	EwesThisYearComponent <- as.character(ewedatnow$component.ind)
	ThisYearCompEwes <- levels(factor(subset(relocdata,
											 as.character(component.ind) ==
											 as.character(EwesThisYearComponent))$EWEID))
	ThisYearCompLambs <- subset(lambdata, as.numeric(as.character(Year)) ==
															as.numeric(as.character(lambdata$Year[i])) &
															as.character(EWEID) %in% ThisYearCompEwes)
	lambdata$ThisYearKnownCompMort[i] <- sum(ThisYearCompLambs$CENSOR2) /
	dim(ThisYearCompLambs)[1]
	lambdata$ThisYearCompLambDiedOrNoLamb[i] <- sum(ThisYearCompLambs$CENSOR2) /
	length(ThisYearCompEwes)

	#-- extract mort levels for this popyear --#
	EwesThisPopyr <- as.character(ewedatnow$popear.ind)
	ThisPopyrEwes <- levels(factor(subset(relocdata,
											 as.character(popear.ind) ==
											 as.character(EwesThisPopyr))$EWEID))
	ThisPopyrLambs <- subset(lambdata, as.numeric(as.character(Year)) ==
															as.numeric(as.character(lambdata$Year[i])) &
															as.character(EWEID) %in% ThisPopyrEwes)
	lambdata$ThisPopyrKnownMort[i] <- sum(ThisPopyrLambs$CENSOR2) /
	dim(ThisPopyrLambs)[1]
	lambdata$ThisPopyrLambDiedOrNoLamb[i] <- sum(ThisPopyrLambs$CENSOR2) /
	length(ThisPopyrEwes)
}

write.path <- "~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/LambSurvDatWithLastYearCompCovs/"
write.csv(lambdata, paste(write.path,
													"LambDataWithLastYearCompCovs_19Sept2013.csv", sep =
													""))
