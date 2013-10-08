#-- script to build lamb survival data from ewe network-level relocation data --#
#-- KRM 14 Feb 2013; revised 17 Sept 2013 --#

#-- read in data relocation data (constructed with Main datacleaning script)
filepath <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/RelocsWithNetworkMeasures/"

#data <- read.csv(paste(filepath,
#											 "FullEweDataAllSummerRelocs_16Sept2013.csv",sep =
#											 ""), header = T) 

data <- read.csv(paste(filepath,
											 "FullEweDataAllSummerRelocs_MinEdge.1_18Sept2013.csv",sep =
											 ""), header = T) 

#-- read in merged.lamb.data, and subset down to a manageable size --#
#--(cut out pops not included in networks) --#
lambpath <- "~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/"
merged.lamb.data <- read.csv(paste(lambpath, "MergedLambData_080512_krmmod.csv", sep =
																	 ""), header = T, sep = ",")

#-- build factor covariates indicating component --#
data$component.ind <- paste(data$Pop,"_",data$Year,"_",data$res.component,sep="")
data$popyear.ind <- paste(data$Pop, "_", data$Year, sep = "")

#-- check to be sure degree never exceeds resident component size --#
node.component.size.mismatch <- subset(data,node.degree>=res.component.size)

#-- subset data down to only include ewes that currently have lambs --#
#-- (to extract lambs for lamb survival analysis) --#
data.withlambs <- subset(data,HasLamb == "HadLamb")
withlambs <- subset(data, HasLamb == "HadLamb")

for(i in 1:dim(withlambs)[1]){
  k <- subset(merged.lamb.data, EWEID == as.character(withlambs$EWEID[i]) &
							YEAR == withlambs$Year[i])
	withlambs$LAMBID[i] <- ifelse(dim(k)[1] == 0, NA, as.character(k$LAMBID))
}

#-- now, subset withlambs to just one row per lamb --#
#-- (no info lost here, since lamb survival data is replicated in all ewe rows --#
lamb.data <- vector("list",length(levels(factor(withlambs$LAMBID))))
for(i in 1:length(levels(factor(withlambs$LAMBID)))){
  k <- subset(withlambs,LAMBID == as.character(levels(factor(withlambs$LAMBID)))[i])
  out.vec <- k[1, ]
  lamb.data[[i]] <- out.vec
}

#-- rbind all elements in the "lamb.data" list together --#
lamb.data.out <- do.call(rbind, lamb.data)

#-- do some post-processing to clean up the lamb dataframe --# 
names(lamb.data.out) <- names(withlambs)
lamb.data <- lamb.data.out
lamb.data$SURV_DAYS <- as.numeric(as.character(lamb.data$SurvDays))
lamb.data$CENSOR2 <- as.numeric(as.character(lamb.data$CENSOR2))
lamb.data$Pop <- as.character(levels(data$Pop)[lamb.data$Pop])

#-- write out the newly processed lamb.data as a csv. --#
#-- path for all edges --#
#path <-
#"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/CleanLambSurvData/FullEweRelocsLambSurvDat_allewerelocs_17Sept2013.csv"

#-- path for edges with weights >= .1 --#
path <- "~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/CleanLambSurvData/FullEweRelocsLambSurvDat_MinEdge.1_allewerelocs_18Sept2013.csv"

write.csv(lamb.data, paste(path))
