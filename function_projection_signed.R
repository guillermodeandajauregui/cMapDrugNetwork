#projection function for bipartite, signed graphs
library(igraph)

bipartite.projection.signed <- function(g, type=TRUE){
  #break into plus signed and minus signed networks
  q_pos = igraph::subgraph.edges(graph = g, 
                         eids = E(g)[which(E(g)$weight == 1)], 
                         delete.vertices = FALSE)
  
  q_neg = igraph::subgraph.edges(graph = g, 
                         eids = E(g)[which(E(g)$weight == -1)], 
                         delete.vertices = FALSE)
  
  #project each graph
  q_pos_proj = igraph::bipartite.projection(q_pos, which = type, multiplicity = TRUE)
  q_neg_proj = igraph::bipartite.projection(q_neg, which = type, multiplicity = TRUE)
  
  #make union 
  q_union_proj = igraph::union(q_pos_proj, q_neg_proj)
  
  #change NAs in "weights"
  E(q_union_proj)$weight_1 = ifelse(is.na(E(q_union_proj)$weight_1), 0, E(q_union_proj)$weight_1)
  E(q_union_proj)$weight_2 = ifelse(is.na(E(q_union_proj)$weight_2), 0, E(q_union_proj)$weight_2)
  
  #make final weight 
  E(q_union_proj)$weight = E(q_union_proj)$weight_1 + E(q_union_proj)$weight_2
  
  return(q_union_proj)
}
