########################################################################################################################
## SOCIAL NETWORK ANALYSIS (771A23, 2022)
## Lab session on Community detection 
## R script written by Jose Luis Estevez (email: jose.luis.estevez@liu.se)
## Date: March 4th, 2021
########################################################################################################################

# GENERAL INFORMATION

# In this lab, we will learn how  to detect sub-groups or communities within a given a network using different algorithms.  
# For the sake of simplicity, in this session only undirected dichotomous [0,1] networks will be used. Some of these 
# measurements, however, can be extended to directed and even weighted networks. 

########################################################################################################################

# REQUIRED PACKAGES
library(sna);library(igraph)

# Clean the session, and set the directory to wherever you stored the data set
rm(list=ls())
load('lab6_2022.03.04.RData')

########################################################################################################################

# Some overview of the data: attributes
str(attributes) 
head(attributes,10) # show the first ten rows

# Visual overview
par(mfrow=c(1,3))
hist(attributes[attributes$matrix == 1,]$age,main='Age (network 1)',xlim=c(20,60))
hist(attributes[attributes$matrix == 2,]$age,main='Age (network 2)',xlim=c(20,60))
hist(attributes[attributes$matrix == 3,]$age,main='Age (network 3)',xlim=c(20,60))

# Note that I am performing the exact same analyses in all three networks. This means that the script is somewhat 
# redundant in terms of code. However, I preferred to make the code humanly accessible rather than computationally 
# efficient. If you prefer it, you can use loops to simplify the syntax as below: 
for(i in 1:3){ hist(attributes[attributes$matrix == i,]$age,main=paste('Age (network ',i,')',sep=''),xlim=c(20,60)) }

plot(attributes[attributes$matrix == 1,]$gender,main='Gender (network 1)')
plot(attributes[attributes$matrix == 2,]$gender,main='Gender (network 2)')
plot(attributes[attributes$matrix == 3,]$gender,main='Gender (network 3)')

plot(attributes[attributes$matrix == 1,]$position,main='Position (network 1)')
plot(attributes[attributes$matrix == 2,]$position,main='Position (network 2)')
plot(attributes[attributes$matrix == 3,]$position,main='Position (network 3)')

# You can cross variables if you want to
table(attributes$gender,attributes$position) 

########################################################################################################################

# Some overview of networks

# Turn matrices into igraph objects (to be able to call igraph functions)
m1 <- graph_from_adjacency_matrix(matrix_1,mode='undirected',diag=FALSE)
m2 <- graph_from_adjacency_matrix(matrix_2,mode='undirected',diag=FALSE)
m3 <- graph_from_adjacency_matrix(matrix_3,mode='undirected',diag=FALSE)

########################################################################################################################

# SOME DESCRIPTIVE STATS OF THE NETWORK
# Density (number of actual ties in the network over all possible ties)
sna::gden(matrix_1,mode='graph') # using sna's function "gden"
igraph::edge_density(m1) # using igraph's function "edge_density"
gden(matrix_2,mode='graph')
gden(matrix_3,mode='graph')
# These networks are undirected, so reciprocity does not make sense, but you can obtain it (equal to 1)
sna::grecip(matrix_1,measure='edgewise') # using sna's function "grecip"
igraph::reciprocity(m1) # using igraph's function "reciprocity"
# Transitivity
sna::gtrans(matrix_1) # using sna's function "gtrans"
igraph::transitivity(m1) # using igraph's function "transitivity"
gtrans(matrix_2)
gtrans(matrix_3)

# Degree distribution
hist(sna::degree(matrix_1,gmode='graph'),xlim=c(0,13),main='Degree distr. (network 1)')
hist(sna::degree(matrix_2,gmode='graph'),xlim=c(0,13),main='Degree distr. (network 2)')
hist(sna::degree(matrix_3,gmode='graph'),xlim=c(0,13),main='Degree distr. (network 3)')

# Assortativity
# Assortativity (gender)
assortativity(m1,attributes[attributes$matrix==1,]$gender)
assortativity(m2,attributes[attributes$matrix==2,]$gender)
assortativity(m3,attributes[attributes$matrix==3,]$gender)
# Assortativity (degree)
assortativity_degree(m1)
assortativity_degree(m2)
assortativity_degree(m3)

# Connectivity
is_connected(m1) # are all nodes connected
is_connected(m2)
is_connected(m3)

# which nodes are the isolates?
isolates(matrix_1) # it shows the place holder, to see the ID
rownames(matrix_1)[isolates(matrix_1)]
isolates(matrix_2)
rownames(matrix_2)[isolates(matrix_2)]
isolates(matrix_3)
rownames(matrix_3)[isolates(matrix_3)]

# Distance
diameter(m1) # diameter
diameter(m2)
diameter(m3)

mean_distance(m1) # you can also use the function average.path.length()
mean_distance(m2) 
mean_distance(m3)

# Now that we have some descritpive stats of the newtorks, let us visualise them:
# Here I am just saving the Fruchterman-Reingold (fr) layout for future plots
lay1 <- layout_with_fr(m1)
lay2 <- layout_with_fr(m2)
lay3 <- layout_with_fr(m3)

par(mfrow=c(1,3))
plot(m1,
     vertex.color = ifelse(attributes[attributes$matrix==1,]$gender=='female','magenta','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==1,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==1,]$age),
     vertex.label = substr(attributes[attributes$matrix==1,]$ID,4,5),
     main='Network 1',
     layout = lay1)

plot(m2,
     vertex.color = ifelse(attributes[attributes$matrix==2,]$gender=='female','magenta','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==2,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==2,]$age),
     vertex.label = substr(attributes[attributes$matrix==2,]$ID,4,5),
     main='Network 2',
     layout = lay2)

plot(m3,
     vertex.color = ifelse(attributes[attributes$matrix==3,]$gender=='female','magenta','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==3,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==3,]$age),
     vertex.label = substr(attributes[attributes$matrix==3,]$ID,4,5),
     main='Network 3',
     layout = lay3)

# We clearly see that, whereas Network 1 consists of one big component (plus some isolates), Network 3 consists of two
# components (and a couple of isolates again). Network 2 is somewhere in between the other two: there exists one large 
# component, but this is tightly connected in some parts, and loosely connected in others.

########################################################################################################################

# NETWORK COMPONENTS 

# Depending on the nature of your data, sometimes it makes sense to consider different components of a network 
# separately, as "networks" in their own right.
m3_components <- decompose(m3)
table(sapply(m3_components,vcount))
sapply(m3_components,vcount)
m3_components

# This new object has 5 components: three of them are the isolates (actors 51, 62, 74). The other two are the two large
# components we could spot in the plot:
m3_components[[2]] # 12 nodes and 18 ties 
m3_components[[3]] # 10 nodes and 16 ties

# Now, we can analyse the components separately
edge_density(m3_components[[2]])
transitivity(m3_components[[2]])
diameter(m3_components[[2]])
mean_distance(m3_components[[2]])

# I most cases, however, what you have is something akin to Network 2. In such cases, one can try to find "communities"
# (or relatively close-knit communities) rather than components

########################################################################################################################

# CLIQUES

# Cliques are groups of nodes that are maximally connected (i.e., all connected to each other)
cliques(m2,min=5,max=NULL) # find clique of min 5 nodes all connected
cliques(m2,min=4,max=NULL)
cliques(m2,min=3,max=NULL)

largest_cliques(m2) # find the largest clique
clique_num(m2) # tell me the size of the largest clique

########################################################################################################################

# HIERARCHICAL CLUSTERING

# Hierarchical clustering can be applied to the matrix of clique overlaps or to the matrix of geodesic distances. Here,
# I will do the second.

# Geodesic distances (min number of tie a node need to traverse to find another node in the network)
geo_distances <- geodist(matrix_2)
str(geo_distances) # two objects: the important is gdist
rownames(geo_distances$gdist) <- colnames(geo_distances$gdist) <- rownames(matrix_2)
#View(geo_distances$gdist)
# Infinite is obtain whenever a node cannot reach another throughout the network (aka presence of isolates)
# The problem is, we cannot run hierarchical clustering if there are infinite then
h_clustering <- hclust(as.dist(geo_distances$gdist), # don't forget to turn the matrix into a distance object
                       method='complete') # several methods available: maximum or complete, minimum or single, and average

# Solution: assign the number of nodes in the network to those infinite...
geo_distances$gdist[is.infinite(geo_distances$gdist)] <- nrow(geo_distances$gdist)
geo_distances$gdist
# and then divide the geodesic distances by the nember of nodes in the network
geo_distances$gdist <- geo_distances$gdist/nrow(geo_distances$gdist)
geo_distances

h_clustering <- hclust(as.dist(geo_distances$gdist),method='complete')

# Visualisation
par(mfrow=c(1,2))
h_clustering$labels <- colnames(matrix_2) # names need to be added to the plot
plot(h_clustering)

plot(m2,
     vertex.color = ifelse(attributes[attributes$matrix==2,]$gender=='female','magenta','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==2,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==2,]$age),
     vertex.label = substr(attributes[attributes$matrix==2,]$ID,4,5),
     main='Network 2',
     layout = lay2)

# The problem with hierarchical clustering is, where do we cut the tree?
# say we want 5 clusters:
cutree(h_clustering,k=3) # this solution is the large component plus the two isolates

par(mfrow=c(1,1))
plot(h_clustering)
rect.hclust(h_clustering,3,border='green')

# But what if we do not have a pre-defined number of clusters to cut the network? ...
rect.hclust(h_clustering,4,border='red')
rect.hclust(h_clustering,5,border='blue')
rect.hclust(h_clustering,6,border='purple')

########################################################################################################################

# CLUSTERING ALGORITHMS

# GIRVAN-NEWMAN
# The Girvan-Newman method partitions the graph based on edge betweeness. 
# The algorithm works iteratively. First, it calculates edge betweeness, and then it deletes the edge with the highest
# score. After deleting this edge, it recalculates edge betweeeness again and repeats the same process.
ecount(m2) # there are 39 ties
E(m2)
(edge_btw_1 <- edge_betweenness(m2))
max(edge_btw_1)
which(edge_btw_1 == max(edge_btw_1))
delete.edges(m2,19)
# and we keep repeating this process for all 39 ties
gn_clustering <- cluster_edge_betweenness(m2,modularity=TRUE,membership=TRUE)
gn_clustering
gn_clustering$modularity
# Communities are identified as components in the edge-pruned graph. 
# The pruning is set to the level where modularity is the highest.

# WALKTRAP
# The walktrap algorithm finds communities through a series of short random walks. 
# The idea is that these random walks tend to stay within the same community.
walk_clustering <- cluster_walktrap(m2,step=10)
walk_clustering

# MOODY-WHITE
blocks_clustering <- cohesive_blocks(m2)
blocks_clustering$blocks

# Other clustering algorithms you can use are:
cluster_louvain(m2)
cluster_fast_greedy(m2)
cluster_leading_eigen(m2)

# How to recover the groups for node-level analysis
gn_clustering$membership
attributes_m2 <- attributes[attributes$matrix == 2,]
attributes_m2$group <- as.factor(gn_clustering$membership)
attributes_m2

# Also in network form (same-group membership)
outer(attributes_m2$group,attributes_m2$group,'==')*1

# And adding clustering methods to a visualisation is very easy
par(mfrow=c(1,3))
plot(m2,
     vertex.color = ifelse(attributes[attributes$matrix==2,]$gender=='female','magenta','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==2,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==2,]$age),
     vertex.label = substr(attributes[attributes$matrix==2,]$ID,4,5),
     main='Network 2 (Girvan-Newman)',
     mark.groups = cluster_edge_betweenness(m2),
     layout = lay2)

plot(m2,
     vertex.color = ifelse(attributes[attributes$matrix==2,]$gender=='female','magenta','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==2,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==2,]$age),
     vertex.label = substr(attributes[attributes$matrix==2,]$ID,4,5),
     main='Network 2 (Walktrap)',
     mark.groups = cluster_walktrap(m2,step=10),
     layout = lay2)

plot(m2,
     vertex.color = ifelse(attributes[attributes$matrix==2,]$gender=='female','magenta','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==2,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==2,]$age),
     vertex.label = substr(attributes[attributes$matrix==2,]$ID,4,5),
     main='Network 2 (Cohesive blocks)',
     mark.groups = cohesive_blocks(m2)$blocks,
     layout = lay2)

# End of the script