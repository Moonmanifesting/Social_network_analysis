########################################################################################################################
## SOCIAL NETWORK ANALYSIS (771A23, 2020-2021)
## Lab session 4 (Mathematical models for network graphs and the Quadratic Assignment Procedure)
## R script written by Jose Luis Estevez (Linköping University) (email: jose.luis.estevez@liu.se)
## Date: February 17th, 2021
########################################################################################################################

# GENERAL INFORMATION

# This lab session consists of three parts. First, we will revisit the exercises of the previous lab session devoted to 
# community detection methods. Second, we will learn how to use the quadratic assignment procedure (QAP for short) for 
# network overlap. Finally, we will learn how to use algorithms (i.e., random graph, small world, preferential 
# attachment) to generate well-known network graphs.

########################################################################################################################

# REQUIRED PACKAGES
# If you did not install these yet, first run:
#install.packages(c('sna','igraph'))
library(sna);library(igraph)
#packageDescription('sna')
#packageDescription('igraph')

# Clean the session, and set the directory to wherever you stored the data set (i.e., lab3.RData)
rm(list=ls())
getwd()
setwd("C:/Users/joses32/OneDrive - Linköpings universitet/Desktop/SNA/Lab4")
load('lab3_2021.02.10.RData')

########################################################################################################################

# 1) COMMUNITY DETECTION (EXERCISES)

# EXERCISE: Do a 6-cluster partition of Network 1, and a 5-cluster partition of Network 3 (line 253).
m1 <- graph_from_adjacency_matrix(matrix_1,mode='undirected',diag=FALSE) # matrix to igraph object for visualisation
par(mfrow=c(1,3))
lay1 <- layout_with_fr(m1)
plot(m1,
     vertex.color = ifelse(attributes[attributes$matrix==1,]$gender=='female','yellow','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==1,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==1,]$age),
     vertex.label = substr(attributes[attributes$matrix==1,]$ID,4,5),
     main='Network 1',
     layout = lay1)

geo_distances <- geodist(matrix_1) # geodesic distance between actors 
geo_distances$gdist[is.infinite(geo_distances$gdist)] <- nrow(geo_distances$gdist) # Inf to number of nodes
geo_distances$gdist <- geo_distances$gdist/nrow(geo_distances$gdist) # Divide by the number of nodes

h_clustering <- hclust(as.dist(geo_distances$gdist),method='complete') # Hierarchical clustering
h_clustering$labels <- colnames(matrix_1) # names for the plot
plot(h_clustering)

hcl_solution <- cutree(h_clustering,k=6) # k = 6 (six clusters or communities)
rect.hclust(h_clustering,6,border='red')

# If you want to see the groups visually...
hcl_part <- vector('list',6) # a list with the number of clusters

for(i in seq_along(hcl_part)){
  hcl_part[[i]] <- names(hcl_solution[hcl_solution==i]) 
}

hcl_part

plot(m1,
     vertex.color = ifelse(attributes[attributes$matrix==1,]$gender=='female','yellow','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==1,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==1,]$age),
     vertex.label = substr(attributes[attributes$matrix==1,]$ID,4,5),
     main='Six-cluster solution',
     layout = lay1,
     mark.groups=hcl_part)

# EXERCISE: Apply the Girman-Newman and the Fast Greedy methods to Network 1. Apply one clustering method (up to you 
# which one) to Network 3 (lines 326-327).
par(mfrow=c(1,2))

cluster_edge_betweenness(m1)
plot(m1,
     vertex.color = ifelse(attributes[attributes$matrix==1,]$gender=='female','yellow','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==1,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==1,]$age),
     vertex.label = substr(attributes[attributes$matrix==1,]$ID,4,5),
     main="Girman-Newman's solution",
     layout = lay1,
     mark.groups=cluster_edge_betweenness(m1))

cluster_fast_greedy(m1)
plot(m1,
     vertex.color = ifelse(attributes[attributes$matrix==1,]$gender=='female','yellow','lightblue'),
     vertex.shape = ifelse(attributes[attributes$matrix==1,]$position=='manager','square','circle'),
     vertex.size = 2*sqrt(attributes[attributes$matrix==1,]$age),
     vertex.label = substr(attributes[attributes$matrix==1,]$ID,4,5),
     main='Fast greedy solution',
     layout = lay1,
     mark.groups=cluster_fast_greedy(m1))

########################################################################################################################

# 2) THE QUADRATIC ASSIGNMENT PROCEDURE 

# 3.1. Multiple and logistic regression models (quick review)

# Example of a multiple regression (MR) test question: Are managers better connected than employees? 
attributes <- attributes[attributes$matrix == 1,]
print(attributes)
# (let's use betweenness centrality to check that)
attributes$btw_cent <- centr_betw(m1)$res
print(attributes)

mean(attributes[attributes$position == "manager",]$btw_cent)
mean(attributes[attributes$position == "employee",]$btw_cent)

model1 <- glm(data=attributes, # dataset
              formula=btw_cent ~ position, # formula
              family=gaussian(link='identity')) # family
summary(model1)

# Even after controlling for gender differences? (control)
model2 <- glm(data=attributes, 
              formula=btw_cent ~ position + gender,
              family=gaussian(link='identity')) 
summary(model2)

# Example of a logistic regression (LR) test
attributes$deg_5 <- 1*(degree(m1) >= 5) # nodes with 5 or above degree
print(attributes)

model3 <- glm(data=attributes,
              formula = deg_5 ~ position,
              family = binomial(link = 'logit')) # a binomial distribution, the logit
summary(model3)

# 3.2. LR-QAP (for non-dichotomous networks: MR-QAP)

# QAP enables running logistic regression tests when the data is not a vector of attributes (like the node's position
# or gender), but a matrix of relations. How does it do it? The LR-QAP algorithm proceeds in two steps. In step 1, it 
# performs a regular logistic regression across the cells of the response matrix and the predicting matrices. In step 2,
# it randomly permutes both rows and columns of the response matrix and recomputes the regression, storing resultant 
# values of all coefficients. This is repeated x times (often thousands) to estimate standard errors for the statistics 
# of interest. Then, for each coefficient, the programme counts the proportion of random permutations that yielded a 
# coefficient as extreme as the one computed in step 1, to test its significance. 

# For us, it is important to recall that:
# a) As the assumption of independent observations does not hold, to test the significance of the results we create a 
#    sampling distribution based on permutations of the existing data
# b) We can control for structural properties, say whether an actor is an isolate, has high degree, reciprocity, or
#    assortativity/homophily in certain dimensions (like gender)
# c) At a more practical level, to run a LR-QAP, we must transform all our variables to a matrix form

# Now we can ask questions like, do nodes with the same position connect with each other more often than they do with
# nodes of a different position? Or do differences in degree matters: nodes with high degree connect with other nodes 
# with high degree?
same_position <- 1*outer(attributes$position,attributes$position,FUN='==') # same position
diff_deg <- abs(outer(degree(m1),degree(m1),FUN='-')) # absolute differences in degree

matrix_1 # check that your matrices have no missing data
diag(matrix_1) <- 0 # remove all NAs for QAP

model4 <- netlogit(y=matrix_1, x=list(same_position,diff_deg),
                   nullhyp='qap', # By default, netlogit uses Dekker's semi-partialling routine (qapspp)
                   reps=100) # Number of permutations: increase to 1000 to reach a precision of 3 decimals (p-values)
summary(model4)

# EXERCISE: Run a (standard) MR to check whether female employees have higher degree than male employees in network 2. 
# Find out whether there exists gender homophily in network 3.

########################################################################################################################

# 3) MATHEMATICAL MODELS FOR NETWORK GRAPHS

rm(list=ls())

# 3.1) Classical Random Graph model (Erdos-Renyi)
set.seed(123)
rg <- sample_gnp(n=25, # nodes 
                 0.15, # prob. of a tie
                 directed=FALSE,loops=FALSE)

par(mfrow=c(1,2))
plot(rg,layout=layout_with_fr(rg))

# Let's check some of its properties:
edge_density(rg) # density 

is_connected(rg)
diameter(rg)
mean_distance(rg) # average path length

hist(degree(rg),xlab='Degree',ylab='Frequency') # Degree distribution
transitivity(rg) # very low compared to real-world networks

# 3.2) Small world (Caveman)
sworld_1 <- sample_smallworld(1, # number of draws
                              25, # nodes
                              3, # neighbourhood of size x
                              0) # rewiring probability [0,1]

par(mfrow=c(3,2))
plot(sworld_1,layout=layout_in_circle(sworld_1),main='rewiring prob = 0 (circle)') # without rewiring, we obtain a ring
plot(sworld_1,layout=layout_with_fr(sworld_1),main='rewiring prob = 0 (Fruchterman-Reingold)')

set.seed(123)
sworld_2 <- sample_smallworld(1, 
                              25,
                              3, 
                              0.05) # rewiring probability = 5%

plot(sworld_2,layout=layout_in_circle(sworld_2),main='rewiring prob = 0.05 (circle)') 
plot(sworld_2,layout=layout_with_fr(sworld_2),main='rewiring prob = 0.05 (Fruchterman-Reingold)') # still keep some of the ring-like structure

set.seed(123)
sworld_3 <- sample_smallworld(1, 
                              25, 
                              3, 
                              0.2) # rewiring probability = 20%

plot(sworld_3,layout=layout_in_circle(sworld_3),main='rewiring prob = .20 (circle)') 
plot(sworld_3,layout=layout_with_fr(sworld_3),main='rewiring prob = 0.20 (Fruchterman-Reingold)')

# Let's compare sworld_1 and sworld_3: average path and clustering
mean_distance(sworld_1)
transitivity(sworld_1)

mean_distance(sworld_3) # smaller average distance
transitivity(sworld_3) # lower transitivity

# To see this better, let us use a larger network (N = 100) 
results <- vector('list',length=length(seq(0,1,by=.01)))
names(results) <- seq(0,1,by=.01)
for(i in seq_along(results)){
  results[[i]] <- vector('list',length=100)
}

# A hundred samples per prob_value (0, 0.01, 0.02, ... , 1.00)
for(i in names(results)){
  for(j in 1:100){
    results[[i]][[j]] <- sample_smallworld(1,100,3,i) # 100 nodes, 3 neighbours
  }
}

res_clustering <- res_ave_path <- results

for(i in seq_along(results)){
  for(j in seq_along(results[[i]])){
    res_clustering[[i]][[j]] <- transitivity(results[[i]][[j]])
    res_ave_path[[i]][[j]] <- mean_distance(results[[i]][[j]])
  }
  res_clustering[[i]] <- mean(unlist(res_clustering[[i]])) # the mean value of 100 samples (clustering)
  res_ave_path[[i]] <- mean(unlist(res_ave_path[[i]])) # the mean value of 100 samples (average path)
}

res_clustering <- do.call('cbind',res_clustering)
res_ave_path <- do.call('cbind',res_ave_path)

results <- data.frame(p_value = colnames(res_clustering),
                      clustering = as.vector(res_clustering),
                      ave_path = as.vector(res_ave_path))
print(results)

results$p_value <- log10(as.numeric(results$p_value)) # the log transformation of p_value

par(mfrow=c(1,2))
plot(results$p_value,results$clustering,type='l',col='blue',xlab='Rewiring probability',ylab='Transitivity')
plot(results$p_value,results$ave_path,type='l',col='red',xlab='Rewiring probability',ylab='Mean average path')

# 3.3) Preferential attachment
set.seed(123)
pref_att <- sample_pa(100,
                      power=1,
                      directed=FALSE)

plot(pref_att)

edge_density(pref_att) # density 

is_connected(pref_att) # one component
diameter(pref_att) # large distances
mean_distance(pref_att) # average path length

hist(degree(pref_att),xlab='Degree',ylab='Frequency') # Degree distribution: few very popular nodes
transitivity(pref_att) # zero transitivity

# End of the script