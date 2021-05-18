# This code here shows how to visually reduced a network as "simplified" version of it. This simplification can be the 
# output of a clustering algorithm or that of a blockmodel

# Required packages
library(networkdata);library(igraph);library(sna)

# Let us pick one network: that of INSNA, for example
INSNA <- as.matrix(as_adj(insna, type = 'both', names = TRUE))

# This is how the original network looks like visually
plot(graph_from_adjacency_matrix(INSNA),
     main='INSNA',vertex.size=7.5,edge.arrow.size=.5,
     layout=layout_with_kk(graph_from_adjacency_matrix(INSNA)))

# Say we wanted a simplified version, for instance, a plot of the relationships between/within blocks in the INSA

# Equivalence cluster detection
INSNA.clust <- equiv.clust(INSNA,equiv.fun="sedist",cluster.method='complete')
plot(INSNA.clust)
(INSNA.block <- blockmodel(INSNA,INSNA.clust,k=5)) # I chose 5 blocks

# The output is transformed into a new igraph object
INSNA.reduced <- graph_from_adjacency_matrix(INSNA.block$block.model, weighted = TRUE)

INSNA.reduced
E(INSNA.reduced)$weight
# Make the with of the ties a function of this weight for the plot 
E(INSNA.reduced)$width <- 10 * E(INSNA.reduced)$weight

plot(INSNA.reduced,
     main='INSNA reduced',vertex.size=7.5,edge.arrow.size=1,
     layout=layout_with_kk(INSNA.reduced))