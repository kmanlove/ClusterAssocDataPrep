#-------------------------------------------------#
#-- Code to build networks based on ewe relocations, and exctract --#
#-- network metrics from the newly-constructed ewe networks --#-
#------------------------------------------#
#------------------------------------------#
#-- this script writes data contained in --#
#-- work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/LambSurviva_Dec2012/CleanData/EweRelocationsWithNetworkMetrics/ --#
#------------------------------------------#

#-- START 12 SEPT 2013 CODE BLOCK--#
filepath <-
	"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/RevisedData_11Sept2013/"
data <- read.csv(paste(filepath, "FourPopRelocs_Through2012.csv", sep = ""), header = T)
#-- END 12 SEPT 2013 CODE BLOCK --#
#-- build pop-year-specific relocation datasets --#
table(data$Pop,data$Year)
pops<-levels(factor(data$Pop))
years<-levels(factor(data$Year))
pop.id<-year.id<-rep(NA,length(pops)*length(years))
data.subsets<-list(NA,length(pops)*length(years))
for(i in 1:length(pops)){
  for(j in 1:length(years)){
    data.subsets[[(i-1)*length(years)+j]]<-subset(data,Pop==pops[i] & Year==years[j])
    pop.id[(i-1)*length(years)+j]<-pops[i]
    year.id[(i-1)*length(years)+j]<-years[j]
  }
}

#-- source in function to build association matrix --#

source("~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Code/DataCleaning/RelocationsToNetworks/StaticNetworkAssocMat.R")
#test<-AssocTimeVect(AssocData=data.subsets[[2]])

edgelists.nozeros<-edgelists<-inds<-output.info<-list(NA,length(pops)*length(years))

#-- clique percolation method (Palla et al., 2005) stolen from this site: --#
#-- igraph.wikidot.com/community-detection-in-r --#

clique.community <- function(graph, k){
	clq <- cliques(graph, min = k, max = k)
	edges <- c()
	for(i in seq_along(clq)){
		for(j in seq_along(clq)){
			if((length(unique(c(clq[[i]], clq[[j]]))) == k + 1)){
#				edges <- c(edges, c(i, j) - 1)
				edges <- c(edges, c(i, j) )
#				edges <- c(edges, c(V(static.graph)$name[i], V(static.graph)$name[j]))
			}
		}
	}
	clq.graph <- simplify(graph(edges))
	V(clq.graph)$name <- seq_len(vcount(clq.graph))
	comps <- decompose.graph(clq.graph)
	lapply(comps, function(x){
	unique(unlist(clq[V(x)$name]))
	})
}
clique.test.fun <- function(x){
	y <- i %in% x 
	return(y)
}

require(igraph)
edgeweight.min <- .1

#for(i in 1:length(data.subsets)){
for(i in 105:length(data.subsets)){
  if(dim(data.subsets[[i]])[1] == 0 | length(levels(factor(data.subsets[[i]]$EWEID))) == 1){
    edgelists[[i]] <- NA
    inds[[i]] <- NA
    output.info[[i]] <- NA
  } else{
    out <- AssocTimeVect(AssocData=data.subsets[[i]])
    inds[[i]] <- out[[1]]
    edgelists[[i]] <- out[[2]]
    edgelists[[i]]$edgeweights <- as.numeric(as.character(edgelists[[i]]$TimesTogether)) / (as.numeric(as.character(edgelists[[i]]$TotalInd1)) + as.numeric(as.character(edgelists[[i]]$TotalInd2)) - as.numeric(as.character(edgelists[[i]]$TimesTogether)))
      #-- subset to remove all edges where edgeweight == 0 --#
#   edgelists.nozeros[[i]] <- subset(edgelists[[i]], edgeweights != 0)
    edgelists.nozeros[[i]] <- subset(edgelists[[i]], edgeweights >=
																		 edgeweight.min)
		#-- included edges with weight = 0.  No good. 
    el <- cbind(as.character(edgelists.nozeros[[i]]$Ind1), as.character(edgelists.nozeros[[i]]$Ind2))
      #-- Were any associations observed? Generate list of isolated nodes --#
			#-- to add to network later --#	
    connected.nodes <- levels(factor(c(el[, 1], el[, 2])))
    isolated.nodes <- inds[[i]][! (inds[[i]] %in% connected.nodes) == T] 
      
    if(dim(el)[1] == 0){  #-- if NO associations observed:   
      res.component <- no.components <- max.component.size <- res.component.size <- node.degree <- cv.edgeweights <- mean.edgeweights <- sd.edgeweights <- local.efficiency <- rep(NA, length(out[[1]]))
      for(m in 1:length(isolated.nodes)){
        res.component[m] <- m
        no.components[m] <- length(isolated.nodes)
        max.component.size[m] <- 1
        node.degree[m] <- cv.edgeweights[m] <- mean.edgeweights[m] <- sd.edgeweights[m] <- local.efficiency[m] <- 0
      }
      individ.ids <- as.character(isolated.nodes)
      #-- leave everything as NAs for now....
    } else {
      static.graph.orig <- graph.edgelist(el, directed=F)
      static.graph.weighted <- set.edge.attribute(static.graph.orig, "weight", value = edgelists.nozeros[[i]]$edgeweights)
      E(static.graph.weighted)$weight<-edgelists.nozeros[[i]]$edgeweights
    if (length(isolated.nodes) >= 1){
      static.graph <- add.vertices(static.graph.weighted, nv = length(isolated.nodes), name = as.character(isolated.nodes))
    } else {
      static.graph <- static.graph.weighted
    }
		graph.clusters <- clusters(static.graph)
		obs.clique.communities <- clique.community(static.graph, k = 3)
		components <- length(levels(factor(graph.clusters$membership)))
		numK3communities <- length(obs.clique.communities)
		k3communities <- res.component <- no.components <- max.component.size <- res.component.size <- node.degree <- cv.edgeweights <- mean.edgeweights <- sd.edgeweights <- rep(NA, length(out[[1]]))
		no.k3com <- max.k3com.size <- no.communities.for.this.node <- rep(NA, length(out[[1]]))
		res.k3com <- vector("list", length(out[[1]]))
		for(m in 1:length(V(static.graph)$name)){
      if(dim(edgelists.nozeros[[i]])[1] == 0){ #-- if the graph has no edges
        node.degree[m] <- cv.edgeweights[m] <- mean.edgeweights[m] <- sd.edgeweights[m] <- local.efficiency[m] <- 0
        no.components[m] <- length(levels(factor(V(static.graph)$name)))
        max.component.size[m] <- 1
        res.component.size[m] <- 1
        res.component[m] <- m
        mean.edgeweights[m] <- 0
        sd.edgeweights[m] <- NA
        cv.edgeweights[m] <- NA
				res.k3com[[m]] <- NA
				no.k3com[m] <- NA
				no.communities.for.this.node[m] <- NA
				max.k3com.size[m] <- 1
				k3communities[m] <- numK3communities
      } else {
			#-- pull all edges linked to node m --#
      ind.edges <- subset(edgelists.nozeros[[i]], Ind1 ==
													as.character(V(static.graph)$name)[m] | Ind2 ==
													as.character(V(static.graph)$name)[m])
      if(dim(ind.edges)[1] == 0){
        node.degree[m] <- cv.edgeweights[m] <- mean.edgeweights[m] <- sd.edgeweights[m] <- 0
        no.components[m] <- length(levels(factor(graph.clusters$membership)))
        max.component.size[m] <- max(graph.clusters$csize)
        res.component.size[m] <- graph.clusters$csize[graph.clusters$membership[m]]
        res.component[m] <- graph.clusters$membership[m]
        mean.edgeweights[m] <- 0
        sd.edgeweights[m] <- NA
        cv.edgeweights[m] <- NA
				res.k3com[[m]] <- NA
				no.communities.for.this.node[m] <- NA
				max.k3com.size[m] <- max(unlist(lapply(obs.clique.communities, length)))
				k3communities[m] <- numK3communities
			} 
      else {
      cv.edgeweights[m] <- (sd(ind.edges$edgeweights)) / (mean(ind.edges$edgeweights))
      mean.edgeweights[m] <- mean(ind.edges$edgeweights)
      sd.edgeweights[m] <- sd(ind.edges$edgeweights)[1]
      node.degree[m] <- degree(static.graph.weighted)[as.character(V(static.graph)$name)[m]]
      no.components[m] <- length(levels(factor(graph.clusters$membership)))
      max.component.size[m] <- max(graph.clusters$csize)
      res.component.size[m] <- graph.clusters$csize[graph.clusters$membership[m]]
      res.component[m] <- graph.clusters$membership[m]
			res.k3com[[m]] <-  which(unlist(lapply(obs.clique.communities, clique.test.fun)) == T)
			no.communities.for.this.node[m] <- length(res.k3com[[m]])
			max.k3com.size[m] <- max(unlist(lapply(obs.clique.communities, length)))
			k3communities[m] <- numK3communities
      }
    }
  }
      individ.ids <- as.character(V(static.graph)$name)
    }
  pop <- rep(paste(pop.id),length(V(static.graph)$name))
  year <- rep(paste(year.id),length(V(static.graph)$name))
  output.info[[i]] <- as.data.frame(cbind(individ.ids,node.degree,
																					cv.edgeweights,mean.edgeweights,
																					sd.edgeweights, no.components,
																					max.component.size,
																					res.component.size, res.component,
																					no.communities.for.this.node,
																					max.k3com.size, k3communities, rep(as.character(pop.id)[i], length(individ.ids)), rep(year.id[i], length(individ.ids))))
  }
}

#-- in FourPopRelocsPrepped, breaks at --#
#-- i = 25, 50, 79, 80, 83, 103, 104
#-- for min edgeweight = .1, breaks at 
#-- i = 23, 25, 50, 89, 103, 104
#-- now, stick output.info information onto "data" --#
#-- unlist output.info and get it into a dataframe --#

full.output.info <- do.call(rbind,output.info)
names(full.output.info)[1] <- "EWEID"
names(full.output.info)[13] <- "HERD"
names(full.output.info)[14] <- "YEAR"

data$node.degree <- data$cv.edgeweights <- data$mean.edgeweights <-
	data$sd.edgeweights <- data$no.components <- data$max.component.size <-
		data$res.component.size <- data$res.component <-
			data$mean.incomponent.edgeweights <- data$sd.incomponent.edgeweights <-
				data$no.communities.for.this.node <- data$numk3communities <- data$max.k3com.size <- rep(NA,dim(data)[1])

for(i in 1:dim(data)[1]){
  k <- subset(full.output.info, EWEID == as.character(data$EWEID)[i] & as.character(YEAR) == as.character(data$Year[i]))
  data$node.degree[i] <- ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$node.degree)))
  data$cv.edgeweights[i] <- ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$cv.edgeweights)))
  data$mean.edgeweights[i] <- ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$mean.edgeweights)))
  data$sd.edgeweights[i] <- ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$sd.edgeweights)))
  data$no.components[i]<-ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$no.components)))
  data$max.component.size[i]<-ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$max.component.size)))
  data$res.component.size[i]<-ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$res.component.size)))
  data$res.component[i]<-ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$res.component)))
  data$mean.incomponent.edgeweights[i] <- ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$mean.incomponent.edgeweights)))
  data$sd.incomponent.edgeweights[i] <- ifelse(dim(k)[1] == 0, NA, as.numeric(as.character(k$sd.incomponent.edgeweights)))
	data$no.communities.for.this.node[i] <- ifelse(dim(k)[1] == 0, NA,
																								 as.numeric(as.character(k$no.communities.for.this.node)))
	data$max.k3com.size[i] <- ifelse(dim(k)[1] == 0, NA,
																	 as.numeric(as.character(k$max.k3com.size)))
	data$numk3communities[i] <- ifelse(dim(k)[1] == 0, NA,
																		 as.numeric(as.character(k$k3communities)))
	data$LAMBID[i] <- ifelse(dim(k)[1] == 0, NA,
																		 (as.character(k$LAMBID)))
}

#write.csv(data,"~/work/Kezia/Research/EcologyPapers/ClustersAssociations/Data/LambSurvival_Dec2012/CleanData/EweRelocationsWithNetworkMetrics/FullEweDataAllSummerRelocs_MinEdge.1_18Sept2013.csv")
