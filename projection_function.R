
Projection_Function<-function(graph){

players<-V(graph)[V(graph)$type==TRUE]
Dg = make_empty_graph(directed = FALSE)
Dg = add_vertices(graph = Dg, nv = length(players), name = players)


for(i in seq_along(players)){
  for(j in seq_along(players)){
    if(i != j){
      if(j > i){
        genes = intersect(neighbors(graph = graph, v = players[i], mode = "all")$name,
                          neighbors(graph = graph, v = players[j], mode = "all")$name)
        peso = 0
        print(genes)
        print(peso)
        for(k in genes){
          if(E(g)[V(g)[players[i]]%--%k]$weight == E(g)[V(g)[players[i]]%--%k]$weight){
            peso = peso + 1
            }

          }
        if(peso!=0){
          Dg<-add.edges(graph = Dg,
                        edges = c(V(Dg)[V(Dg)$name==players[i]],
                                  V(Dg)[V(Dg)$name==players[j]]),
                        weight = peso
                        )
          
        }
      }
    }
  }
}
return(Dg)
}

