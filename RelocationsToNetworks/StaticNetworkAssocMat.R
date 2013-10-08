AssocTimeVect<-function(AssocData){
  ind <- unique(factor(AssocData$EWEID))
  n <-  length(ind) # number of individuals in the dataset
  
  dyads <-  n*(n-1)/2
  TimesTogether<-DyadID<-Ind1<-Ind2<-TimesInd1Separate<-TimesInd2Separate<-TotalInd1<-TotalInd2<-rep(NA,length(dyads))
  Avec <- rep(NA, dyads*5)
  dim(Avec) <- c(dyads,5)
  cooccur <- vector("list",dyads)
  
  i=1
  j=2  
  for(k in 1:dyads){
    Avec[k,1] <- ind[i]  # individual 1
    Avec[k,2] <- ind[j]  # individual 2
    ind1<-ind[i]
    ind2<-ind[j]
    temp.a<-subset(AssocData,as.character(AssocData$EWEID) ==
																					as.character(ind[i]))
    temp.2a<-subset(AssocData,as.character(AssocData$EWEID) == as.character(ind[j]))
    
    cooccur[[k]]<-merge(temp.a,temp.2a,by = c("DATE","UTME","UTMN")) #-- HERE'S THE PROBLEM. 
    TimesTogether[k]<-ifelse(dim(cooccur[[k]])[1]==0,0,dim(cooccur[[k]])[1])
    TimesInd1Separate[k]<-ifelse(dim(temp.a)[1]==0,0,dim(temp.a)[1])-TimesTogether[k]
    TimesInd2Separate[k]<-ifelse(dim(temp.2a)[1]==0,0,dim(temp.2a)[1])-TimesTogether[k]
    TotalInd1[k]<-ifelse(dim(temp.a)[1]==0,0,dim(temp.a)[1])
    TotalInd2[k]<-ifelse(dim(temp.2a)[1]==0,0,dim(temp.2a)[1])
    DyadID[k]<-paste(ind[i],"_",ind[j],sep="")
    Ind1[k]<-as.character(ind1)
    Ind2[k]<-as.character(ind2)
#    cooccur[[k]]$Ind1<-rep(ind[i],dim(cooccur[[k]])[1])
#    cooccur[[k]]$Ind2<-rep(ind[j],dim(cooccur[[k]])[1])
#    cooccur[[k]]$Sex1<-rep(ind1$SEX[1],dim(cooccur[[k]])[1])
#    cooccur[[k]]$Pop1<-rep(ind1$Population[1],dim(cooccur[[k]])[1])
#    cooccur[[k]]$Sex2<-rep(ind2$SEX[1],dim(cooccur[[k]])[1])
#    cooccur[[k]]$Pop2<-rep(ind2$Population[1],dim(cooccur[[k]])[1])
#    cooccur[[k]]$split<-cooccur[[k]]$diff<-cooccur[[k]]$lag<-rep(NA,dim(cooccur[[k]])[1])
#    cooccur[[k]]$LambEndDate1<-cooccur[[k]]$LambEndDate2<-cooccur[[k]]$LambBirthDate1<-cooccur[[k]]$LambBirthDate2<-cooccur[[k]]$next.together<-cooccur[[k]]$next.overall<-cooccur[[k]]$next.overall.date<-cooccur[[k]]$next.together.date<-rep(NA,dim(cooccur[[k]])[1])
#    
#    AllPairRelocs<-data.frame(rbind(temp.a,temp.2a))
#    
#    if(dim(cooccur[[k]])[1]>=1){
#      for(m in 1:dim(cooccur[[k]])[1]){
#        #-- new 30 Oct 2012 --#
#        #-- merge in lamb data --#
#        lamb1<-subset(LambData,as.character(EWEID)==as.character(cooccur[[k]]$Ind1[m]) & YEAR==as.POSIXlt(cooccur[[k]]$DATE3[m])$year+1900)
#        cooccur[[k]]$LambBirthDate1[m]<-ifelse(dim(lamb1)[1]==0,NA,as.character(lamb1$BIRTH[1]))
#        cooccur[[k]]$LambEndDate1[m]<-ifelse(dim(lamb1)[1]==0,NA,as.character(lamb1$END[1]))
#        lamb2<-subset(LambData,as.character(EWEID)==as.character(cooccur[[k]]$Ind2[m]) & YEAR==as.POSIXlt(cooccur[[k]]$DATE3[m])$year+1900)
#        cooccur[[k]]$LambBirthDate2[m]<-ifelse(dim(lamb2)[1]==0,NA,as.character(lamb2$BIRTH[1]))
#        cooccur[[k]]$LambEndDate2[m]<-ifelse(dim(lamb2)[1]==0,NA,as.character(lamb2$END[1]))
#        
#        #-- end 30 Oct 2012 --#
#        
#        FutureMerge<-subset(cooccur[[k]],DATE3>cooccur[[k]]$DATE3[m])
#        future.overall<-subset(AllPairRelocs,DATE3>cooccur[[k]]$DATE3[m])
#        cooccur[[k]]$next.together[m]<-ifelse(dim(FutureMerge)[1]==0,NA,min(as.numeric(FutureMerge$DATE3)))
#        cooccur[[k]]$next.together.date[m]<-as.character(min(FutureMerge$DATE3))
#        cooccur[[k]]$next.overall[m]<-ifelse(dim(future.overall)[1]==0,NA,min(as.numeric(future.overall$DATE3)))
#        cooccur[[k]]$next.overall.date[m]<-as.character(min(future.overall$DATE3))
#        cooccur[[k]]$diff[m]<-ifelse(is.na(cooccur[[k]]$next.together[m])==T,1,cooccur[[k]]$next.together[m]-cooccur[[k]]$next.overall[m])
#        cooccur[[k]]$nextobs[m]<-min(cooccur[[k]]$next.together[m],cooccur[[k]]$next.overall[m])
#        cooccur[[k]]$split[m]<-ifelse(cooccur[[k]]$diff[m]==0,0,1)
#        print(paste(i,"_",j,"_",m,sep=""))
#      }
#    }
    
    i_new<-ifelse(j!=n,i,i+1)
    i<-i_new
    j_new<-ifelse(j==n,i_new+1,j+1)
    j<-j_new
    print(i)
    print(j)
  }
  out.mat<-as.data.frame(cbind(Ind1,Ind2,DyadID,TimesTogether,TimesInd1Separate,TimesInd2Separate,TotalInd1,TotalInd2))
  out.list<-list(ind,out.mat)
  return(out.list)
}
