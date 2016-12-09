########################################################
#
#Functions for handling connectivity Map (cMap) data. 
# Guillermo de Anda - Jauregui
# guillermo.deandajaur@med.und.edu
#
########################################################

###########
#Libraries#
###########
library("data.table")
library("stringr")
library("GeneExpressionSignature")
library("igraph")
###########

######
# fread.Ranked.Matrix 
#Read Ranked Matrices
######

fread.Ranked.Matrix <- function(filepath){
  mx<-fread(filepath, data.table = FALSE)
  rownames(mx) = mx[,1]
  colnames(mx) = mx[1,]
  mx = mx[-1,-1]
  return(mx)
}

######
# drug.eset
# make an ExpressionSet of samples 
# perturbed by drugs on a list, 
# based on a file with sample information
######

drug.eset<- function(RankedMatrix, SampleInfoFile, DrugList){
  #Identify samples treated with drugs of interest
  DrugInstances = SampleInfoFile$instance_id[which(SampleInfoFile$normalized_name%in%DrugList)]
  #make a subset of the RankedMatrix, containing only these Ids
  Drug.RankMatrix = RankedMatrix[,as.character(DrugInstances)]
  #dummy eset, to force sample order
  eset = ExpressionSet(assayData = as.matrix(Drug.RankMatrix))
  #make drug annotation object for the ExpressionSet
  DrugData = SampleInfoFile[which(SampleInfoFile$instance_id%in%colnames(exprs(eset))), c(1, 4)] #id, drugname
  rownames(DrugData) = DrugData[,1]
  DrugData = DrugData[,-1, drop = FALSE]
  DrugAFD = new("AnnotatedDataFrame", data = DrugData)
  #make the ExpressionSet
  eset = ExpressionSet(assayData = as.matrix(Drug.RankMatrix), phenoData = DrugAFD)
  return(eset)
}

######
# cMap.eset
# make an ExpressionSet of samples 
# based on a file with cMap instance information
######

cMap.eset<- function(RankedMatrix, SampleInfoFile){
  #dummy eset, to force sample order
  eset = ExpressionSet(assayData = as.matrix(RankedMatrix))
  #make drug annotation object for the ExpressionSet
  cMapData = SampleInfoFile[which(SampleInfoFile$instance_id%in%colnames(exprs(eset))), ] 
  rownames(cMapData) = cMapData[,1]
  cMapData = DrugData[,-1, drop = FALSE]
  cMapAFD = new("AnnotatedDataFrame", data = cMapData)
  #make the ExpressionSet
  eset = ExpressionSet(assayData = as.matrix(RankedMatrix), phenoData = cMapAFD)
  return(eset)
}

######
# read.Drugs
# helper function to read and format a drug list file
# 
######

read.Drugs<- function(DrugList){
  drugList = read.table(file = DrugList)
  drugList = drugList$V1
  drugList = str_to_upper(drugList)
}

#######
# Drug.Gene.Graph
# Takes a Kru-Bor merged set of ranked differentially expressed genes 
# after drug treatment
# and a threshold
# Returns a directed, unweighted network from DRUGS to GENES
#######

Drug.Gene.Graph <- function(MergedDrugEset, Threshold){
  DrugGraph = exprs(MergedDrugEset)
  minz = Threshold
  maxz = max(DrugGraph) - Threshold
  DrugGraph <- ifelse(minz>=DrugGraph | DrugGraph>maxz, 1, 0)
  DrugGraph = graph_from_incidence_matrix(incidence = DrugGraph, 
                                          directed = TRUE, 
                                          mode = "in")
  return(DrugGraph)
}

######################################
# TopBottomComparer
# Compares Top and Bottom Genes
# Which are the neighbors of drugs in
# The bipartite graph
######################################

TopBottomComparer <-function(G1, G2, NodeList, listed = FALSE){
#Make sure NodeList follows a determined order; 
  #not important here,
  #but very important if you are comparing non-Graph objects
    IntersectList<- list()
  for(i in seq_along(NodeList)){
    G1_neighbors = neighbors(G1, v = NodeList[i])$name
    G2_neighbors = neighbors(G2, v = NodeList[i])$name
    q<- intersect(G1_neighbors, G2_neighbors)
    if(listed==TRUE){
      IntersectList<-append(IntersectList, list(q))
      names(IntersectList)[i] <- NodeList[i]
    }
    else{
      IntersectList<-c(IntersectList, length(q))
      names(IntersectList)[i] <- NodeList[i]
    }
  }
    return(IntersectList)
}

######################################
# TopBottomSetDiff
# Compares Top and Bottom Genes
# Which are the neighbors of drugs in
# The bipartite graph
#returns the unique genes of the first graph
######################################

TopBottomSetDiff <-function(G1, G2, NodeList, listed = FALSE){
  #Make sure NodeList follows a determined order; 
  #not important here,
  #but very important if you are comparing non-Graph objects
  IntersectList<- list()
  for(i in seq_along(NodeList)){
    G1_neighbors = neighbors(G1, v = NodeList[i])$name
    G2_neighbors = neighbors(G2, v = NodeList[i])$name
    q<- setdiff(G1_neighbors, G2_neighbors)
    if(listed==TRUE){
      IntersectList<-append(IntersectList, list(q))
      names(IntersectList)[i] <- NodeList[i]
    }
    else{
      IntersectList<-c(IntersectList, length(q))
      names(IntersectList)[i] <- NodeList[i]
    }
  }
  return(IntersectList)
}

#######
# Drug.Gene.Sign.Graph
# Takes a Kru-Bor merged set of ranked differentially expressed genes 
# after drug treatment
# and a threshold
# Returns a directed, weighted network from DRUGS to GENES
#with weight being either "plus" or "minus"
#interpretable as an activation or repression of said GENE expression
#by DRUG
#######

Drug.Gene.Sign.Graph <- function(MergedDrugEset, Threshold){
  DrugGraph = exprs(MergedDrugEset)
  minz = Threshold
  maxz = max(DrugGraph) - Threshold
  
  DrugGraph<-ifelse(minz>=DrugGraph, -1, 
                    ifelse(DrugGraph>maxz, 1,
                           0)
  )
  DrugGraph = graph_from_incidence_matrix(incidence = DrugGraph, 
                                          directed = TRUE, 
                                          mode = "in",
                                          weighted = TRUE)
  return(DrugGraph)
}

#######
# ConnectedComponentMembership
# Takes a graph (and optionally, a components object, default to components(graph))
# and a minimum Component Size (default 1)
#returns a list with the NAMES of NODES in each component.
#######

ConnectedComponentMembership<-function(graph, 
                                       ComponentObject= components(graph), 
                                       ComponentSize = 1){
  q = lapply(seq_along(ComponentObject$csize)[ComponentObject$csize>=ComponentSize],
             function(x) V(graph)$name[q$membership%in%x]
  )
  return(q)
}

#######
# ActivatorsInhibitors
# Takes a graph (and optionally, a set of GENE vertices)
#returns a DataFrame with the number of activators and inhibitors
#For said gene.
#######

ActivatorsInhibitors <- function(DrugGraph, 
                                 nodes=V(DrugGraph)[V(DrugGraph)$type!=TRUE]){
  
  df = data.frame(gene = as.character(),
                  activators = as.numeric(),
                  inhibitors = as.numeric(),
                  stringsAsFactors = FALSE
  )
  
  for(i in seq_along(nodes)){
    gen = nodes[i]$name
    Nhoods = ego(graph = DrugGraph, 
                 order = 1, 
                 nodes = gen, 
                 mode = "in")
    Nhood_sb = induced.subgraph(graph = DrugGraph, 
                                vids = unlist(Nhoods)
    )
    activator = 0
    inhibitor = 0
    for(j in E(Nhood_sb)$weight){
      if(j==1){
        activator=activator+1
      }else
        if(j==-1){
          inhibitor=inhibitor+1
        }
    }
    k=list(gen, activator, inhibitor)
    df[i,]<-k
  }
  return(df)
}
