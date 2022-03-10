########################################################################################################################

# SOCIAL NETWORK ANALYSIS
# Computational Social Science - Linkoping University
# 4 February 2022
# Jose Luis Estevez (jose.luis.estevez@liu.se or jose.luis.a.estevez.navarro@liu)
# The script reuses part of previous scripts written by Dorottya Kisfalusi and Andras Voros

########################################################################################################################
# PLEASE NOTE: All code and data used in the Kolaczyk-Cs?rdi book have been made available in the 
# R package sand, distributed through the CRAN archive.
# the codes are also available here: https://github.com/kolaczyk/sand

########################################################################################################################

# 0. PACKAGES FOR SOCIAL NETWORK ANALYSIS IN R

# There are plenty of packages tailored to perform different types of SNA in R. Among these, 
# - "sna" 
# - "igraph"
# - "network"
# - "ggraph"  for visualisations (uses ggplot2)
# - "statnet"  meta-package; herein you have "ergm", "tergm", etc.
# - "RSiena"  for Stochastic Actor-Oriented Models (SAOMs)
# - "intergraph"  for conversion between network data objects (igraph, network, matrix, edge list, etc.)
# - "influenceR"  for measures of structural importance (effective network size, constraint, bridging)
# - "egor"  for ego-nets (uses igraph objects)
# - "netseg"  for assortativity patterns, homophily, segregation, within-group mixing, etc.
# - "netropy" 
# - "networkdata"  a large sample of networks in igraph format (to install networkdata, see: https://github.com/schochastics/networkdata)
# - etc.

# If not yet, install "sna", "igraph", and "networkdata"
# install.packages("sna")
# install.packages("igraph")
# install.packages("remotes")
# remotes::install_github("schochastics/networkdata")

# Loading and detaching packages
library(sna) # beware of masked functions (att, order from package:base)
library(igraph)
# Sometimes convenient to detach packages
# detach(package:sna,unload = TRUE)

# Or call functions by specifying the package they are from
# sna::betweenness()
# igraph::betweenness()

########################################################################################################################

### 1. LOADING THE DATA 

# Network data can be stored in different ways (edges lists, matrices, igraph objects, etc) 
# A different types of files (.csv, .xlsx, .R, etc.)

starwars_iv <- networkdata::movie_656 # as igraph object
class(starwars_iv)
starwars_iv
V(starwars_iv)$name
E(starwars_iv)$weight
plot(starwars_iv)

thegodfather <- as_edgelist(networkdata::movie_304) # as edge list
thegodfather
TGF_ntw <- network(thegodfather,directed=FALSE) # turned into a network object (this is what ergm uses)
class(TGF_ntw)
TGF_ntw
network.vertex.names(TGF_ntw)
gplot(TGF_ntw,displaylabels=TRUE,label.cex=.5,usearrows = FALSE)

killbill <- as.matrix(get.adjacency(networkdata::movie_401)) # as matrix
View(killbill)
KB_igraph <- graph_from_adjacency_matrix(killbill,mode='undirected')
class(KB_igraph)
plot(KB_igraph)

# From now on, I will presume we depart from a matrix object (like in the assignment)

### The Glasgow dataset
# The dataset we use now is a subsample of the Teenage Friends and Lifestyle Study conducted in Scotland between 1995 
# and 1997. The dataset contains information on friendships, other dyadic variables (e.g. distance between homes), and 
# individual attributes (e.g. sex, drinking) of students from 3 waves. We use a subsample of a whole student cohort 
# comprising of 160 students from a single school.

# Further information about the data and other variables can be downloaded from the RSiena website: 
# http://www.stats.ox.ac.uk/~snijders/siena/Glasgow_data.htm

# first make sure to clear your workspace
rm(list=ls())

# check your working directory 
getwd()
# setwd("THE FOLDER WITH THE DATA FILES")
setwd("C:/Users/247379/OneDrive - MUNI/Plocha/SNA labs/week2lab")

# the data is in .RData format this time, which is R's own file type,
# let's load it
load("Glasgow-friendship.RData")
ls() 

# we are only going to use data from the first two waves, we can remove the last object
rm(friendship.3)

# now load the demographic and geopgraphic datasets
load("Glasgow-demographic.RData") # sex (1 boy, 2 girl), age
load("Glasgow-geographic.RData") # distance from school and others
ls()

# we remove what we do not need
rm(angle.1, angle.2, angle.3, dist.school, distance.3)
ls() 

### 2. INSPECTION OF THE DATA OBJECTS

# In all cases, before we start analyzing a new dataset, we need to learn about the data. This means two things in R:
# Step 1: learn about the objects (class, mode, no. of observations, value codes, etc.)
# Step 2: learn about the data (distributions, missings, etc.)

# the first step on the friendship networks
class(friendship.1)
mode(friendship.1)
dim(friendship.1)
rownames(friendship.1)
colnames(friendship.1)
View(friendship.1)
# are there the same students in the first two waves?
all.equal(rownames(friendship.1), rownames(friendship.2)) # note: the order matters for all.equal()

rownames(friendship.1) <- rownames(friendship.1)[order(rownames(friendship.1),decreasing = TRUE)]
rownames(friendship.1)
all.equal(rownames(friendship.1), rownames(friendship.2))
rownames(friendship.1) <- rownames(friendship.1)[order(rownames(friendship.1),decreasing = FALSE)]
all.equal(rownames(friendship.1), rownames(friendship.2))

# what kind of values are there in the data?
table(friendship.1,useNA="always")
table(friendship.2,useNA="always") # 0 = not friend; 2 = best friend; 1 = just a friend; 10 = structural missing

# let's have a look at the other variables
table(sex.F) # 2 - girl, 1 - boy
hist(age) # age in years
hist(distance.1) # distance between students' homes in km at t1
hist(distance.2) # the same at t2

########################################################################################################################

# 2. MATRIX BASIC TRANSFORMATION: DICHOTOMISATION (recoding), TRANSPOSING, SYMMETRIZATION

# 2.1. Dichotomisation
# Say we want to turn 1 and 2 into ties (1), and set the structural missing to non-ties (0) 
friendship.1[friendship.1 %in% c(1,2)] <- 1
friendship.1[friendship.1 == 10] <- 0
table(friendship.1,useNA="always")

friendship.2[friendship.2 %in% c(1,2)] <- 1
friendship.2[friendship.2 == 10] <- 0
table(friendship.2,useNA="always")

# Remove the diagonal just in case
diag(friendship.1) <- 0
diag(friendship.2) <- 0

# if you would like to recode missing data into 0:
friendship.1[is.na(friendship.1)] <- 0
friendship.2[is.na(friendship.2)] <- 0

table(friendship.1,useNA="always")
table(friendship.2,useNA="always")

# 2.2. Transposing
fr1_t <- t(friendship.1)

friendship.1['s005','s006']
friendship.1['s006','s005']

fr1_t['s005','s006']
fr1_t['s006','s005']

# 2.3. Symmetrization (turn a directed matrix into undirected)
fr1_weak <- symmetrize(friendship.1,rule='weak') # max-symmetrizing
fr1_strong <- symmetrize(friendship.1,rule='strong') # min-symmetrizing

network(friendship.1,directed=TRUE)
network(fr1_weak,directed=FALSE)
network(fr1_strong,directed=FALSE)

gplot(network(friendship.1,directed=TRUE))
gplot(network(fr1_weak,directed=FALSE),usearrows = FALSE)
gplot(network(fr1_strong,directed=FALSE),usearrows = FALSE)
# As you might have noticed, nodes change position from chart to chart, let's solve this with some plotting tricks 

########################################################################################################################

# 3. NETWORK VISUALIZATION IN sna

# let's describe the friendship networks at the two time points!
# from here on we will mostly use functions from the sna package

# how does the network look like - plot it!
?gplot
plot1 <- gplot(friendship.1)
plot2 <- gplot(friendship.2)

# Two things:
# 1. if you run gplot on the same data several times, it will be a little different
#    because the node placement algorithm is not exact but iterative
# 2. related to node placement, it is hard to compare the two plots - what can we do?

# what is saved in plot1 and plot2? you would never guess.
head(plot1) # head() displays the first part of an object
# so, what if...
gplot(friendship.1, coord=plot1) # now it will always look the same
# and
gplot(friendship.2, coord=plot1)
# it looks quite horrible, but here is a little trick:
plot12 <- gplot(friendship.1 + friendship.2)
# and now plot them both with the same node coordinates
gplot(friendship.1, coord=plot12)
gplot(friendship.2, coord=plot12)

# let's pimp our plots: they need a title, and perhaps boys and girls could get different colors
gplot(friendship.1, coord=plot12, main="Friendship network - wave 1", vertex.col=sex.F)
gplot(friendship.2, coord=plot12, main="Friendship network - wave 2", vertex.col=sex.F)

# Colors overlap between nodes and ties, let's improve this
gplot(friendship.1, coord=plot12, main="Friendship network - wave 1", vertex.col=sex.F,
      edge.col ='darkgray',vertex.border='black')
gplot(friendship.2, coord=plot12, main="Friendship network - wave 2", vertex.col=sex.F,
      edge.col ='darkgray',vertex.border='black')

### Graphical parameters:
# gplot, just like other plotting functions in R, uses a GRAPHICS DEVICE
# this has many features that can be changed by equally many GRAPHICAL PARAMETERS
# examples: how large should be the plotting region? what kind of font should be used for
# the title? how big should symbols be? (the default setup is roughly what you see now)
#
# type ?par whenever you want to find out about the plotting settings you can adjust
# the help page also suggests good and bad practices for restoring the default settings
# or use google to find out about what the "R community" suggests for dealing with par()
###

# the graphical paramters are set to the default when you restart R,
# so you need to change them before you run gplot, like this:
par(mfrow=c(1,2))
# now plot your graphs again (still looks bad, but we will do better later)
gplot(friendship.1, coord=plot12, main="Friendship network - wave 1", vertex.col=sex.F,
      edge.col ='darkgray',vertex.border='black')
gplot(friendship.2, coord=plot12, main="Friendship network - wave 2", vertex.col=sex.F,
      edge.col ='darkgray',vertex.border='black')
# after you are done, you need to reset the graphical parameter you changed
par(mfrow=c(1,1))

# Personally, for further customisation, I prefer to use 'igraph' or 'ggraph'

########################################################################################################################

### 4. NETWORK DESCRIPTIVES

# some basic functions in the sna package 

# a) the density of the networks
# Density is the proportion of existing ties on all ties.
sum(friendship.1) / (nrow(friendship.1) * (nrow(friendship.1) - 1))

gden(friendship.1)
gden(friendship.2)
# NAs are omitted by default - ?gden to see other options

# b) reciprocity

?grecip
grecip(friendship.1, measure="edgewise")
grecip(friendship.2, measure="edgewise")

# c) transitivity

# The transitivity index we use is the proportion of closed two-paths
#  a->b->c => a->c
# NAs should be handled!

?gtrans
gtrans(friendship.1, measure = "weak")
gtrans(friendship.2, measure = "weak")

# d) the degree distribution

# Our data comes in an adjacency matrix 
# Rows indicate senders, columns receivers
# The indegree centrality is thus simply the column sums
# R has a function in the base package to calculate this

colSums(friendship.1)

# The outdegree centrality is thus simply the row sums
# R has a function in the base package to calculate this

rowSums(friendship.1)

?sna::degree
(outdeg.1 <- sna::degree(friendship.1, cmode="outdegree")) # it returns a vector of degrees
indeg.1 <- sna::degree(friendship.1, cmode="indegree")
outdeg.2 <- sna::degree(friendship.2, cmode="outdegree")
indeg.2 <- sna::degree(friendship.2, cmode="indegree")

# Mean, standard deviation, minimum, maximum
mean(indeg.1)
mean(indeg.2)

sd(indeg.1)
sd(indeg.2)

# min and max
range(indeg.1)
range(indeg.2)

# which individuals have the lowest and highest indegrees?
which(indeg.1==12)
which(indeg.1==0)

# normally, we would visualize the distribution of degrees on histograms
?hist
hist(outdeg.1)
hist(indeg.1)
# the outdegree dispersion is low in this case, and the histogram looks bad
# we can fix this manually
hist(outdeg.1, breaks=7)
# to make the two distributions better comparable
hist(outdeg.1, xlim=c(0,12), ylim=c(0,50), breaks=7)
hist(indeg.1, xlim=c(0,12), ylim=c(0,50), breaks=13)
# now plot them both
par(mfrow=c(1,2))
hist(outdeg.1, xlim=c(0,12), ylim=c(0,50), breaks=7)
hist(indeg.1, xlim=c(0,12), ylim=c(0,50), breaks=13)
par(mfrow=c(1,1)) # remember to reset graphical parameters

# finally, let's make a relatively nice plot for in- and outdegrees from both waves
par(mfrow=c(2,2))
hist(outdeg.1, xlim=c(0,12), ylim=c(0,50), breaks=7,
     main="Outdegree distribution in wave 1", xlab="Outdegree", col="blue")
hist(indeg.1, xlim=c(0,12), ylim=c(0,50), breaks=13,
     main="Indegree distribution in wave 1", xlab="Indegree", col="red")
hist(outdeg.2, xlim=c(0,12), ylim=c(0,50), breaks=7,
     main="Outdegree distribution in wave 2", xlab="Outdegree", col="blue")
hist(indeg.2, xlim=c(0,12), ylim=c(0,50), breaks=13,
     main="Indegree distribution in wave 2", xlab="Indegree", col="red")
par(mfrow=c(1,1))

########################################################################################################################

# 5. THE igraph PACKAGE

# The igraph package has some functions with which you can calculate other descriptives and make nice figures. One of 
# those functions allow us to measure segregation

# there is a nice and easy-to-use segregation measure in igraph
?assortativity # there are three types of assortativity coefficients
# interpretation: positive: homophily , 0: random, negative: heterophily 

# assortativity: for continuous data
# assortativity.nominal: for nominal data
# assrotatibity.degree: fo network data

# try out what happens when you feed the network matric to the function
assortativity(friendship.1)
# igraph functions often work with special igraph objects only
# the package also comes with a function that is able to transform
# adjacency matrices to igraph objects
graph.1 <- graph.adjacency(friendship.1)
graph.2 <- graph.adjacency(friendship.2)
graph.1 

# age homophily
# unfortunately, the function has some issues with missings at the moment
# there is one missing age, which we now set to the median age
age[is.na(age)] <- median(age, na.rm=TRUE) 

# tendencies for friendship between students of similar/different age at the two waves
assortativity(graph.1, age)
assortativity(graph.2, age)

# second, we asses whether there is gender homophily in friendship choice
assortativity.nominal(graph.1, sex.F)
assortativity.nominal(graph.2, sex.F)

# "degree homophily"
# do students with similar number of friends tend to be friends with each other?
assortativity.degree(graph.1)
assortativity.degree(graph.2)
# this calls the degree function, calculates degrees, and uses them as node attributes
# strangely, you cannot specify in assortativity.degree whether you want to use
# in- or outdegree, it only works with total degree at the moment


### 6. NETWORK VISUALIZATION WITH THE igraph PACKAGE

# let's make nice plots of the two friendship graphs
# visualization is something that igraph is good at

# we already have the two networks as igraph objects
# first, set up a layout that is good for both networks (using the same trick as before)
graph.12 <- graph.adjacency(friendship.1 + friendship.2)
myLayout <- layout.fruchterman.reingold(graph.12)

# then plot the two networks with node colors corresponding to the sex of students
par(mfrow = c(1, 2))
plot(graph.1,
     vertex.color = ifelse(sex.F == 2, "red", "darkblue"),
     vertex.shape = ifelse(sex.F == 2, "square", "circle"),
     edge.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.1,
     vertex.size = 4,
     vertex.label = "",
     layout = myLayout,
     main = "Friendship network - wave 1")
plot(graph.2,
     vertex.color = ifelse(sex.F == 2, "red", "darkblue"),
     vertex.shape = ifelse(sex.F == 2, "square", "circle"),
     edge.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.1,
     vertex.size = 4,
     vertex.label = "",
     layout = myLayout,
     main = "Friendship network - wave 2")
par(mfrow=c(1,1))
# if the plots don't look good on your machine, play around a bit with the plotting arguments

# can we make the node sizes proportionate to the indegree of students?
par(mfrow = c(1, 2))
plot(graph.1,
     vertex.color = ifelse(sex.F == 2, "red", "darkblue"),
     vertex.shape = ifelse(sex.F == 2, "square", "circle"),
     edge.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.1,
     vertex.size = indeg.1,
     vertex.label = "",
     layout = myLayout,
     main = "Friendship network - wave 1")
plot(graph.2,
     vertex.color = ifelse(sex.F == 2, "red", "darkblue"),
     vertex.shape = ifelse(sex.F == 2, "square", "circle"),
     edge.color = "black",
     edge.width = 1,
     edge.arrow.size = 0.1,
     vertex.size = indeg.2,
     vertex.label = "",
     layout = myLayout,
     main = "Friendship network - wave 2")
par(mfrow=c(1,1))

# END OF THE SCRIPT