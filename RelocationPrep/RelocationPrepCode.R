filepath <- "~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/"
pop.in <- "BlackButte"

new.dat <- read.csv(paste(filepath, pop.in, "Relocs_Through2012_Clean.csv", sep
													=""), header = T, sep = "\t")

#-- match individual IDs to sexes and extract females --#
studysheep <- read.csv(paste(filepath,
														 "Study_sheep_toothage_original_012612.csv", sep =
														 ""), header = T)
eweids <- levels(factor(subset(studysheep, SEX == "F")$ID))
imnewerelocs <- subset(new.dat, as.character(ID) %in% as.character(eweids))

#-- cut down to only summer (May 1 to Sept 30) observations --#
sumimnewes <- subset(imnewerelocs, Month %in% c(5:9))

#-- pull in lamb data corresponding to each ewe in each year --#
lambs <- read.csv(paste(filepath, "MergedLambData_080512_krmmod.csv", sep =
												""), header = T)

#-- build vectors to stick lamb data in --#
sumimnewes$EWEID <- sumimnewes$PNYear <- sumimnewes$LAMBID <- sumimnewes$CENSOR2 <- sumimnewes$SurvDays <- sumimnewes$HasLamb <- rep(NA, dim(sumimnewes)[1])
for(i in 1:dim(sumimnewes)[1]){
	k <- subset(lambs, as.character(EWEID) == as.character(sumimnewes$ID[i]) &
							as.numeric(as.character(YEAR))
							== as.numeric(as.character(sumimnewes$Year[i])))
	sumimnewes$HasLamb[i] <- ifelse(dim(k)[1] ==0, "NoLamb", "HadLamb")
	sumimnewes$CENSOR2[i] <- ifelse(dim(k)[1] == 0, NA,
																 as.numeric(as.character(k$CENSOR2)))
	sumimnewes$SurvDays[i] <- ifelse(dim(k)[1] == 0, NA,
																	as.numeric(as.character(k$SURV_DAYS)))
	sumimnewes$PNYear[i] <- ifelse(dim(k)[1] == 0, NA,
																 (as.character(k$PNYear)))
	sumimnewes$LAMBID[i] <- ifelse(dim(k)[1] == 0, NA,
																 (as.character(k$LAMBID)))
#	sumimnewes$EWEID[i] <- ifelse(dim(k)[1] == 0, NA,
#																as.numeric(as.character(k$EWEID)))
}

sumimnewes$Pop <- as.character(toupper(sumimnewes$Herd))
sumimnewes$EWEID <- as.character(sumimnewes$ID)
write.csv(sumimnewes, paste(filepath, "BBRelocsPrepped.csv", sep = ""))
