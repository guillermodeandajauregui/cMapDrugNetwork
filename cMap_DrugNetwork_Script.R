########################################
#
# Networks based on cMap perturbations
#
########################################

########################################
#Libraries
########################################
library("argparser")
library("GeneExpressionSignature")
library("igraph")
source("cMap_functions.R")
########################################
#Argument parsing
########################################
p <- arg_parser('Drug-Gene Perturbations based on cMap')
p <- add_argument(p, "--RankedMatrix",
                  help="filepath to matrix with the rank values of differential expression; samples as columns ", default=7)

p <- add_argument(p, "--SampleInfoFile",
                  help="filepath to Sample Info File")
p <- add_argument(p, "--DrugList",
                  help="filepath to a List of Drugs")
p <- add_argument(p, "--Threshold",
                  help="Number of top and bottom genes to include in network",
                  default = 100)
p <- add_argument(p, "--OutName",
                  help="Graph Output filename")


argv <- parse_args(p)

#########################################
#Define variables
#########################################

MatrixPath = argv$RankedMatrix
SampleInfoPath = argv$SampleInfoFile
DrugFile = argv$DrugList
Threshold = argv$Threshold
odir = argv$Outname
#########################################
#Load files 
#########################################

RankedMatrix = fread.Ranked.Matrix(MatrixPath)
SampleInfoFile = fread(SampleInfoFile, 
                       data.table = FALSE)
DrugList = read.Drugs(DrugFile)

#########################################
#Make ExpressionSet
#########################################

DrugEset = drug.eset(RankedMatrix, SampleInfoFile, DrugList)

#########################################
#Kru-Bor merge samples by drug
#########################################
mergedset = RankMerging(DrugEset)

#########################################
#Graph construction
#########################################
DrugGraph = exprs(mergedset)
uppercut = max(DrugGraph) - Threshold
lowercut = Threshold
DrugGraph <- ifelse(lowercut>=DrugGraph | DrugGraph>uppercut, 
                    1, 
                    0)
G = graph_from_incidence_matrix(incidence = DrugGraph, directed = TRUE, mode = "in")
#########################################
#Graph Writing
#########################################

write.graph(G, 
            file = odir, 
            format = "ncol")
