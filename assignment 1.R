# SNA (Linkoping University)
# Assignment 1 (solutions) 
# Written by Jose Luis Estevez

# DATA LOADING
rm(list=ls())
attributes <- read.csv('attributes.csv') # attributes

# Here I am loading all 17 networks at once
name_ntw <- list.files()
(name_ntw <- name_ntw[startsWith(name_ntw,'ntw')])

networks <- list() # empty list
for(i in name_ntw){
  networks[[i]] <- read.csv(i,header=TRUE,row.names=1)
  colnames(networks[[i]]) <- rownames(networks[[i]]) # this is to correct the X added to the column names
  networks[[i]] <- as.matrix(networks[[i]]) # turn into a matrix object
}

########################################################################################################################

# DICHOTOMISATION
for(i in seq_along(networks)){
  # I chose to turn 2 into 1, all the rest to zero (observe that I am keeping the NAs)
  networks[[i]][!is.na(networks[[i]]) & networks[[i]] %in% c(-2,-1,0,1)] <- 0
  networks[[i]][!is.na(networks[[i]]) & networks[[i]] == 2] <- 1
  
  diag(networks[[i]]) <- 0 # don't forget to remove the diagonal
}

########################################################################################################################

# Matrices to igraph objects
library('igraph')
for(i in names(networks)){
  networks[[i]] <- graph_from_adjacency_matrix(networks[[i]],mode='directed')
  # This is to add attributes to the nodes of the igraph object
  V(networks[[i]])$female <- attributes[attributes$time == substr(names(networks[i]),6,6) & 
                                          attributes$class == substr(names(networks[i]),8,11),]$female
}

########################################################################################################################

# NETWORK DESCRIPTIVES

# This is function just to get descriptives for every network
describe_ntw <- function(igraph_obj){ # it requires a single argument (the network to be described in igraph format)
  # Create a data frame (this will be the output of the function)
  output <- data.frame(matrix(NA,nrow=11,ncol=1))
  rownames(output) <- c('nodes','ties','density','average out/indegree','SD outdegree','SD indegree','reciprocity',
                        'transitivity','isolates','degree assortativity','gender assortativity')
  colnames(output) <- 'Value'
  
  output[1,1] <- length(V(igraph_obj)) # number of nodes or vertices
  output[2,1] <- length(E(igraph_obj)) # number of ties or edges
  output[3,1] <- edge_density(igraph_obj)*100 # density (perc)
  output[4,1] <- mean(degree(igraph_obj,mode='out')) # average out or indegree
  output[5,1] <- sd(degree(igraph_obj,mode='out')) # SD outdegree
  output[6,1] <- sd(degree(igraph_obj,mode='in')) # SD indegree
  output[7,1] <- reciprocity(igraph_obj)*100 # reciprocity
  output[8,1] <- transitivity(igraph_obj)*100 # transitivity
  output[9,1] <- sum(degree(igraph_obj,mode='all')==0) # isolates
  output[10,1] <- assortativity.degree(igraph_obj) # degree assortativity
  output[11,1] <- assortativity.nominal(igraph_obj,as.factor(V(igraph_obj)$female)) # gender assortativity
  return(output)
}

# Apply the function to all networks
ntw_descriptives <- lapply(networks,describe_ntw)
ntw_descriptives <- do.call(cbind,ntw_descriptives) # turn list into a dataset
colnames(ntw_descriptives) <- substr(names(networks),5,11) # rename columns
round(ntw_descriptives,2) # If you use the same dichotomisation procedure, you can check your results here
write.table(round(ntw_descriptives,2),'network_descriptive.csv',sep=',',row.names=TRUE)

########################################################################################################################

# VISUALISATION OF DEGREES

# Here I am just using network 2200 (t1) as an example.
par(mfrow=c(2,1))
hist(degree(networks$ntw_t1_2200.csv,mode='out'),xlim=c(0,10), breaks=10,
     main="Outdegree distribution", xlab="Outdegree", col="tomato")
hist(degree(networks$ntw_t1_2200.csv,mode='in'),xlim=c(0,10), breaks=10,
     main="Indegree distribution", xlab="Indegree", col="royalblue")
par(mfrow=c(1,1))

########################################################################################################################

# STUDENTS DESCRIPTIVES 

# Creation of a dataframe per classroom and time
attributes_by_class <- list()
for(i in names(networks)){
  attributes_by_class[[i]] <- attributes[attributes$time == substr(names(networks[i]),6,6) & 
                                           attributes$class == substr(names(networks[i]),8,11),]
}

describe_nodes <- function(data){ 
  output <- data.frame(matrix(NA,nrow=8,ncol=1))
  rownames(output) <- c('Girls (%)','Roma (%)','Grammar (mean)','Grammar (SD)',
                        'Maths (mean)','Maths (SD)','Science (mean)','Science (SD)')
  colnames(output) <- 'Value'
  
  output[1,1] <- sum(data$female)/length(data$female)*100
  output[2,1] <- sum(data$roma,na.rm=TRUE)/length(data$roma)*100 # don't forget the argument na.rm=TRUE here
  output[3,1] <- mean(data$grades_grammar)
  output[4,1] <- sd(data$grades_grammar)
  output[5,1] <- mean(data$grades_maths)
  output[6,1] <- sd(data$grades_maths)
  output[7,1] <- mean(data$grades_sciences)
  output[8,1] <- sd(data$grades_sciences)
  return(output)
}

nodes_descriptives <- lapply(attributes_by_class,describe_nodes)
nodes_descriptives <- do.call(cbind,nodes_descriptives)
colnames(nodes_descriptives) <- substr(names(networks),5,11) # rename columns
round(nodes_descriptives,2)
write.table(round(nodes_descriptives,2),'nodes_descriptive.csv',sep=',',row.names=TRUE)

########################################################################################################################

# VISUALISATION OF NETWORK

# Here I am just including a couple of examples

jpeg(filename='networks1.jpeg',width=12,height=8,units='in',res=500)

par(mfrow=c(1,2))
set.seed(123)
layout1 <- layout_with_fr(networks$ntw_t1_2200.csv)
plot(networks$ntw_t1_2200.csv,
     vertex.color = ifelse(V(networks$ntw_t1_2200.csv)$female == 1, "magenta", "skyblue"),
     edge.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.3,
     # for visualisation purposes, you can use mathematical transfomations like square roots or logarithms
     vertex.size = 3 + sqrt(betweenness(networks$ntw_t1_2200.csv)), 
     vertex.label = "",
     layout = layout1,
     main = "Classroom 2200 (time 1)")
legend("topright",bty = "n",legend=c('girl','boy'),fill=c('magenta','skyblue'))

set.seed(123)
layout2 <- layout_with_kk(networks$ntw_t3_5100.csv)
plot(networks$ntw_t3_5100.csv,
     vertex.color = ifelse(V(networks$ntw_t3_5100.csv)$female == 1, "magenta", "skyblue"),
     edge.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.3,
     vertex.size = 3 + log(betweenness(networks$ntw_t3_5100.csv)+1),
     vertex.label = "",
     layout = layout2,
     main = "Classroom 5100 (time 3)")
legend("topright",bty = "n",legend=c('girl','boy'),fill=c('magenta','skyblue'))

dev.off()

########################################################################################################################

# FOR SUPPORTING THE THOUGHT EXERCISE

# REMOVE THE STUDENT WITH HIGHEST BETWEENNESS CENTRALITY
# Based on betweenness 
betweenness(networks$ntw_t1_2200.csv) 
max(betweenness(networks$ntw_t1_2200.csv)) # the actor with highest betweenness centrality
(remove_from_t1_2200 <- which(betweenness(networks$ntw_t1_2200.csv) == max(betweenness(networks$ntw_t1_2200.csv)))) 
# student 2203 (position 3)
x1 <- delete.vertices(networks$ntw_t1_2200.csv,remove_from_t1_2200) # Let's remove this person
layout1 <- layout1[-3,]

betweenness(networks$ntw_t3_5100.csv) 
max(betweenness(networks$ntw_t3_5100.csv)) 
(remove_from_t3_5100 <- which(betweenness(networks$ntw_t3_5100.csv) == max(betweenness(networks$ntw_t3_5100.csv)))) 
# Student 5171 (position 24)
x2 <- delete.vertices(networks$ntw_t3_5100.csv,remove_from_t3_5100)
layout2 <- layout2[-24,] 

# REPEAT THE PROCEDURE (SEQUENTIAL PROCEDURE)
betweenness(x1) 
max(betweenness(x1)) # the actor with highest betweenness centrality
(remove_from_x1 <- which(betweenness(x1) == max(betweenness(x1)))) 
# student 2203 (position 3)
x1 <- delete.vertices(x1,remove_from_x1) # Let's remove this person
layout1 <- layout1[-19,]

# Based on popularity
betweenness(x2) 
max(betweenness(x2)) 
(remove_from_t3_5100 <- which(betweenness(x2) == max(betweenness(x2)))) 
# Student 5171 (position 24)
x2 <- delete.vertices(x2,remove_from_t3_5100)
layout2 <- layout2[-20,] 

# REPEAT ONE THIRD TIME
betweenness(x1) 
max(betweenness(x1)) # the actor with highest betweenness centrality
(remove_from_x1 <- which(betweenness(x1) == max(betweenness(x1)))) 
# student 2203 (position 3)
x1 <- delete.vertices(x1,remove_from_x1) # Let's remove this person
layout1 <- layout1[-6,]

# Based on popularity
betweenness(x2) 
max(betweenness(x2)) 
(remove_from_t3_5100 <- which(betweenness(x2) == max(betweenness(x2)))) 
# Student 5171 (position 24)
x2 <- delete.vertices(x2,remove_from_t3_5100)
layout2 <- layout2[-11,] 

# New visualisation
jpeg(filename='networks2.jpeg',width=12,height=8,units='in',res=500)
par(mfrow=c(1,2))

plot(x1,
     vertex.color = ifelse(V(x1)$female == 1, "magenta", "skyblue"),
     edge.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.3,
     vertex.size = 3 + sqrt(betweenness(x1)), 
     vertex.label = "",
     layout = layout1,
     main = "Classroom 2200 (time 1)")
legend("topright",bty = "n",legend=c('girl','boy'),fill=c('magenta','skyblue'))

plot(x2,
     vertex.color = ifelse(V(x2)$female == 1, "magenta", "skyblue"),
     edge.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.3,
     vertex.size = 3 + log(betweenness(x2)+1),
     vertex.label = "",
     layout = layout2,
     main = "Classroom 5100 (time 3)")
legend("topright",bty = "n",legend=c('girl','boy'),fill=c('magenta','skyblue'))

dev.off()