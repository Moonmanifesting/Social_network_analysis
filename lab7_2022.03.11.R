########################################################################################################################

# SOCIAL NETWORK ANALYSIS
# Computational Social Science - Linkoping University
# MATHEMATICAL MODELS FOR NETWORK GRAPHS
# 11 March 2022
# Jose Luis Estevez (jose.luis.estevez@liu.se or jose.luis.a.estevez.navarro@liu)

########################################################################################################################

# Clear the environment
rm(list=ls())
# Load the igraph and networkdata packages
library(igraph);library(networkdata);library(ggplot2)

########################################################################################################################

# 1) "Classical" Random Graph model (Erdos-Renyi)

# In this model, we begin with n isolated nodes. 
# Then, with probability p > 0 each pair of nodes is connected by a link.
# Thus, the network is determined only the number of nodes (n) and the probability of a link (p).
# Alternatively, you can opt for specifying the number of links (m), rather than the probability

prob <- c(0,0.05,0.2,0.5) # Let's pick 4 prob: 0%, 5%, 20% and 50%
ER_graphs <- list() # This is an empty list where results will be stored

for(i in 1:length(prob)){
  # you can use erdos.renyi.game() or random.graph.game()
  ER_graphs[[i]] <- erdos.renyi.game(n=25, # the number of nodes
                                     p=prob[i], # the probability of a tie (here, 0%, 5%, 20% and 50%)
                                     type='gnp', # gnp if you use prob. of a tie. If fixed number of ties, gnm instead
                                     # This notation comes from G(n,p) and G(n,m), respectively
                                     directed=FALSE,loops=FALSE)
}

# Quick look at the networks we just created
ER_graphs 

# Let's visualise them
par(mfrow=c(2,2)) # A 2*2 grid

for(i in seq_along(ER_graphs)){
  plot(ER_graphs[[i]],
       vertex.size = 5,
       vertex.label = NA,
       layout=layout_in_circle(ER_graphs[[i]]), # I will use the layout in circle to facilitate comparison
       main=paste('n=25, p=',prob[i],sep=''))
}

# Random graph models have properties that has been studies extensively

# For example, the expected average degree of the nodes (k) is also determined by n and p
# Specifically, k_{hat} = p*(n-1)
# Let's see this with one example
set.seed(123)
ER_1000 <- erdos.renyi.game(1000, # A network of 1,000 nodes
                            .04, # with prob. of a tie = 4% 
                            type='gnp',directed=FALSE,loops=FALSE)
# In our case, k_{hat} should be approximately (1000-1)*.04, what is 39.96
ER_1000deg <- as.data.frame(table(degree(ER_1000)))
summary(ER_1000deg)
ER_1000deg$Var1 <- as.numeric(as.character(ER_1000deg$Var1)) # define the variable as numeric

# Let's visualise it
ggplot(data=ER_1000deg) +
  geom_point(aes(x=Var1,y=Freq)) 

# let's add the normal curve
ggplot(data=ER_1000deg) +
  geom_point(aes(x=Var1,y=Freq/sum(Freq))) +
  xlab('k')+ylab('p(k)') +
  stat_function(fun = dnorm,
                args = list(mean = mean(degree(ER_1000)), sd = sd(degree(ER_1000))),colour='blue')

# For large values of n, it approaches a Poisson distribution

# Random graph networks display other properties, for example:
# - very small average degree
diameter(ER_1000)
mean_distance(ER_1000)
# - low clustering coefficient, much smaller than that for real world networks with the same density
transitivity(ER_1000)

# Another interesting property (in fact, this is the fundamental result of Erdos and Renyi's original paper in 1960) is 
# that when p increases, most nodes then to be clustered in on giant component, while the rest of the nodes are isolated,
# or clustered in very small components

# According to their results, the structure of G_{ER}(n,p) changes as a function of p = k_{hat}/(n-1), giving raise
# to three stages: 
# 1) sub-critical (k_{hat} < 1; where all components are simple and very small); 
# 2) critical (k_{hat} = 1);
# 3) super-critical (k_{hat} > 1); where there is a giant component and the rest are very small ones

# Let's see this visually
# Say n = 100, and we want three values of k representing each of the 3 stages (say 0.75, 1, and 2.5)
critical_stages <- c(0.75/(100-1), 1/(100-1), 2.5/(100-1))
round(critical_stages,3) # These are the corresponding p 

ER_stages <- list()

set.seed(456)
for(i in 1:length(critical_stages)){
  ER_stages[[i]] <- erdos.renyi.game(n=100, p=critical_stages[i], type='gnp', directed=FALSE, loops=FALSE)
}

# Let's see it visually
par(mfrow=c(1,3)) # A 1*3 grid

for(i in seq_along(ER_stages)){
  plot(ER_stages[[i]],
       vertex.size = 5,
       vertex.label = NA,
       layout=layout_with_fr(ER_stages[[i]]), # I am using the FR algorithm now
       main=paste('n=100, p=',round(critical_stages[i],3),sep=''))
}

########################################################################################################################

# 2) Generalised random graph models

# The most common are random graph models with a fixed degree sequence. 
# Just for fun, let's pick a movie from the networkdata package and use it for obtaining degrees
data(package = 'networkdata')

# I will choose The Lord of the Rings (part 3) because there are sort of parallel stories (and cluster of characters)
LOTR <- movie_440 

# Let's visualise it
par(mfrow=c(1,2))
V(LOTR)$name <- tolower(V(LOTR)$name) # use lowercase letters 
plot(LOTR,
     vertex.size=5,
     layout=layout_with_fr(LOTR), 
     mark.groups=cluster_edge_betweenness(LOTR),
     main='LOTR (part 3)') # Girvan-Newman algorithm for graph partitioning

# Let' see the degree of every character 
degree(LOTR)[order(degree(LOTR),decreasing=TRUE)]

# And use the degrees to create a random graph
set.seed(789)
fake_LOTR <- degree.sequence.game(degree(LOTR),method='vl')
plot(fake_LOTR,
     vertex.size=5,
     vertex.label=NA,
     layout=layout_with_kk(fake_LOTR),
     mark.groups=cluster_edge_betweenness(fake_LOTR),
     main='Random graph with same degrees as LOTR')

# Let's compare some properties
diameter(LOTR)
diameter(fake_LOTR) # the diameter has reduced

transitivity(LOTR)
transitivity(fake_LOTR) # transitivity is almost half

# As with classical random graph models, shorter distance and lower clustering coefficient

########################################################################################################################

# 3) Small world random networks (Watts-Strogatz)

# We have seen that random graph models display small average paths, but they fail in reproducing a large clustering 
# coefficient which is characteristics of most social networks.

# See Milgram's experiments (1967): short average path ('six degree of separation') and large group interconnection.

# In 1998, Watts and Strogatz proposed a model which reproduces networks with short average paths and large clustering 
# coefficients in a very simple way. Let n be the number of nodes, and let k be an even number. 
# 1) Place all nodes in a circle, and connect every node to its first k/2 clockwise nearest neighbours as well as to
#    its k/2 counterclockwise neigbours. This will create a ring, which for k > 2 is full of triangles (therefore,
#     has a large clustering coeff.).

# let's make a network of 25 nodes, and k = 4
sworld_1 <- watts.strogatz.game(1, # number of examples
                                25, # 25 nodes
                                4, # k = 4
                                0, # rewiring prob. = 0
                                loops=FALSE,multiple = FALSE)

# let's visualise it
par(mfrow=c(1,2))

plot(sworld_1,
     vertex.size=5,
     vertex.label=NA,
     layout=layout_in_circle(sworld_1),
     main='n=25, k=4, p=0')

# This network has a high clustering coeff.
transitivity(sworld_1)
# But large distances
diameter(sworld_1) # here 3 steps away, in a network of only 25 (recall the 'six degree of separation')
# In fact, when the network is sparse, as the network size grows the average path lengths tend to n/2
sworld_sparce <- watts.strogatz.game(1,25,1,0, loops=FALSE, multiple = FALSE) # n = 25 and k = 1

plot(sworld_sparce,
     vertex.size=5,
     vertex.label=NA,
     layout=layout_in_circle(sworld_1),
     main='n=25, k=1, p=0')

diameter(sworld_sparce) # the diameter here is exactly n/2

# 2) Watts and Storgatz solved the issue of large distances by considering a probability (p) of rewiring the links in
#    the ring. By doing so, the ave. path length decreases very fast while th clustering coeff. still remains high. 
#    Obviously, as p approaches 1, the network tends to a completetely random graph

# Let's see this with a network of n= 25 and k = 4 using the same probabilities than we used before for random graphs 
# (0, 0.05, 0.2, and 0.5)
prob
SW_graphs <- list()

set.seed(123)
for(i in seq_along(prob)){
  SW_graphs[[i]] <- watts.strogatz.game(1,25,4,prob[i],loops=FALSE, multiple = FALSE)
}

par(mfrow=c(2,2))

for(i in seq_along(SW_graphs)){
  plot(SW_graphs[[i]],
       vertex.size = 5,
       vertex.label = NA,
       layout=layout_in_circle(SW_graphs[[i]]),
       main=paste('n=25, k=4, p=',prob[i],sep=''))
}

# As we see, whereas a random graph could be written as a two-parameter graph G(n,p), a small world network can be 
# written as a three-parameter graph: G(n,k,p)

# What happens to small world graphs is that, while the ave. path length decreases very rapidly as p increases, 
# the clustering coeff. decreases but more slowly. As a consequence of this, there is a region of p values for which
# the networks display relatively large clustering coefficients and small distances. 

# To see this better, let us use a larger network (n = 400) and k = 12

(prob <- seq(0,1,by=.001)) # from 0 to 0.5 by 0.001: 0, 0.001, 0.002, ... 0.999, 1.
SW_400 <- list()
for(i in seq_along(prob)){
  SW_400[[i]] <- watts.strogatz.game(1,400,12,prob[i],loops=FALSE, multiple = FALSE) # same n and k, different p
}

SW_400clus <- SW_400path <- list()

for(i in seq_along(SW_400)){
  SW_400path[[i]] <- mean_distance(SW_400[[i]])
  SW_400clus[[i]] <- transitivity(SW_400[[i]])
}

# Results as a data frame
SW_400results <- data.frame(p = prob,
                            ave_path = unlist(SW_400path),
                            clustering = unlist(SW_400clus))

head(SW_400results,20) # show first 20 lines

# Let's visualise it
SW_400results$p <- log10(as.numeric(SW_400results$p)) # the log10 transformation of p_value
SW_400results$p[1] <- log10(0.0009999) # otherwise it is -Inf
SW_400results$ave_path <- SW_400results$ave_path/max(SW_400results$ave_path) 
SW_400results$clustering <-SW_400results$clustering/max(SW_400results$clustering) 

par(mfrow=c(1,2))
plot(SW_400results$p,SW_400results$ave_path,
     type='l',col='red',xlab='Rewiring probability (log10)',ylab='Mean average path (norm)')
plot(SW_400results$p,SW_400results$clustering,
     type='l',col='blue',xlab='Rewiring probability (log 10)',ylab='Clustering (norm)')

# The same as random graph networks, small world random networks, also display Poissoinian degree distribution
# In this case, the distribution is centered around 2*k

# Let's see this with a network a small world netword of n=1000, k=4, and p=.2
SW_400example <- SW_400[[51]]

SW_400deg <- as.data.frame(table(degree(SW_400example)))
summary(SW_400deg)
SW_400deg$Var1 <- as.numeric(as.character(SW_400deg$Var1)) # define the variable as numeric

SW_400deg

ggplot(data=SW_400deg) +
  geom_point(aes(x=Var1,y=Freq/sum(Freq))) +
  xlab('k')+ylab('p(k)') +
  stat_function(fun = dnorm,
                args = list(mean = mean(degree(SW_400example)), sd = sd(degree(SW_400example))),colour='blue')

########################################################################################################################

# 4) Preferential attachment or 'scale-free' networks (Barabasi-Albert)

# This property of Poissoinian degree distribution, however, does not match many real-world networks. 
# As observed by Barabasi and Albert, many real-world networks are characterised by a few nodes of high degree and a
# large proportion of nodes with relatively few degree.
# The easiest way of conceptualising such topological characteristic is to consider a model in which p(k) ~ k^{-gamma}:
# A model in which the prob. of finding a node with degree k decreases as a power-law of its degree

# Barabasi's and Albert's model follows this procedure:
# Begin with a small number, m_0, of nodes. At each step, add a new node u to the network, and connect it to m <= m_0
# of the existing nodes v with prob.: p_u = {k_v} / {sum_w(k_w)}

# Preferential attachment networks can be written as two-parameter graphs G(n,power), where power is the power of 
# attraction

# Let's see it with three powers: 1, 1.5, and 3
power <- c(1,1.5,3)
PA_graphs <- list()

for(i in seq_along(power)){
  PA_graphs[[i]] <- barabasi.game(400,
                                  power[i],
                                  directed=FALSE)
}

# let's visualise this
par(mfrow=c(1,3))
for(i in seq_along(PA_graphs)){
  plot(PA_graphs[[i]],
       vertex.size=5,
       vertex.label=NA,
       layout=layout_with_fr(PA_graphs[[i]]),
       main=paste('n=400, power=1',power[i],sep=''))
}

# We can see already that the degrees in these graphs do not follow a normal distribution
PAexample <- PA_graphs[[1]]

PAdeg <- as.data.frame(table(degree(PAexample)))
summary(PAdeg)
PAdeg$Var1 <- as.numeric(as.character(PAdeg$Var1)) # define the variable as numeric
PAdeg

ggplot(data=PAdeg) +
  geom_point(aes(x=Var1,y=Freq/sum(Freq))) +
  xlab('k')+ylab('p(k)')
# The distribution is heavy-tailed instead

# We can compare with the small world random network we used before
par(mfrow=c(2,2))
plot(PAexample,
     vertex.size=5,
     vertex.label=NA,
     layout=layout_with_fr(PAexample),
     main='n=400, power=1')
hist(degree(PAexample),xlab='Degree',ylab='Frequency',breaks=25)

plot(SW_400example,
     vertex.size=5,
     vertex.label=NA,
     layout=layout_with_fr(SW_400example),
     main='n=400, k=12, p=0.05')
hist(degree(SW_400example),xlab='Degree',ylab='Frequency',breaks=25)

# Other properties of preferentia attachment networks are
is_connected(PAexample) # one component
diameter(PAexample) # large distances
transitivity(PAexample) # zero transitivity

# End of the script