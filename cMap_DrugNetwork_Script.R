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
p <- add_argument(p, "--OutPath",
                  help="Graph Output Path")


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
# Kruskal Algorithm - Borda Merging 
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
            file = paste0(odir, "graph.txt"), 
            format = "ncol")

#########################################
#Null model construction
#1000 random networks
#Same DRUGS
#With the same number of links to 
#Genes, randomly ranked
#########################################

list_Random = ListShuffleRanks(y, n = 1000)
Random_Networks = mclapply(X = list_Random, mc.cores = 30,
                           FUN = Drug.Gene.Sign.Graph.for.multi, 
                           Threshold = 100, 
                           MaxDrugraph = 12438)
names(Random_Networks)<-lapply(X = 1:1000, FUN = function(x) paste0("G",x))

#Gene Degree
Random_Degree =degree_table_function(Random_Networks)
write.table(Random_Degree, 
            file = paste0(odir, "degree.null"), 
            sep ="\t", quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)
#Activation Gene Degree (links to activator drugs)
ACT_Random_Degree=signed_degree_table_function(Random_Networks, act.or.inhib = "act")
write.table(ACT_Random_Degree, 
            file = paste0(odir, "degree.null"), 
            sep ="\t", quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)
#Inhibition Gene Degree
INHIB_Random_Degree=signed_degree_table_function(Random_Networks, act.or.inhib = "inhib")
write.table(ACT_Random_Degree, 
            file = paste0(odir, "degree.null"), 
            sep ="\t", quote = FALSE, 
            row.names = FALSE, 
            col.names = TRUE)