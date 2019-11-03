#################################################################
#
#   Script for running sampling bias and extinction metrics
#       Deep-time records of Mass Extinctions workshop
#           Nov 4, 2019
#
#################################################################

##### Step 0. libraries and code required to run this script #####

library(ichthyoliths) # for the rangechart function
library(extraDistr) # for the log-series distribution
library(viridis) # for the colors in the legend
library(vegan) # for rarefaction analysis
library(gtools) # for combining sample datasets function

source('sampling_bias_functions.R')
# note, if you don't install ichthyoliths, you can still get the rangechart function - just go to the sampling_bias_scripts.R script and un-comment it, then re-run the source command. It's function #5.

## All libraries except ichthyoliths are available on CRAN. to install, uncomment and run: 
# install.packages("extraDistr") # for the log-series distribution
# install.packages("viridis") #for the colors in the legend
# install.packages("vegan")
# install.packages("gtools")

## to install the ichthyoliths package, you need to be running devtools, if you don't have it already
# install.packages("devtools")
# devtools::install_github('esibert/ichthyoliths')




##### Step 1. create distributions #####
## dist types currently in this are: 
# normal (norm), parameters: mean, sd
# uniform (unif), parameters: min, max (only max really matters)
# lognormal (lnorm), parameters: meanlog, sdlog
# log-series (lgser), parameter: theta
# there are default parameters for each distribution built in, however you can change those specifically
#   in the function call
## hist = TRUE produces a histogram of species-abundances for the distribution so you can see what it is
## currently, the sample creates a distribution of 1 million (1e6 individuals) - you can change this by specifying the parameter n in any of the calls


species.dist <- select.species.distribution(dist.type = 'lgser', theta = 0.9, hist = TRUE)  #default theta=0.8
species.dist <- select.species.distribution(dist.type = 'lgser')
species.dist <- select.species.distribution(dist.type = 'norm') 

species.dist <- select.species.distribution(dist.type = 'unif')

species.dist <- select.species.distribution(dist.type = 'lnorm')

##### Step 2. sample the distribution #####
# time.bins is however many time bins you would like to include in your sample matrix
# sample.sizes is the number of individual fossils you're simulating sampling from the distribution


time.bins <- 50
sample.sizes <- 50
sample.matrix <- build.sample.matrix(species.dist = species.dist, time.slices = time.bins, sample.size = sample.sizes)

# you can hold sample size constant, as above, or create a variable sample size, using, for example, the code below. 
sample.sizes <- sample(x = 10:100, size = time.bins, replace = TRUE)


##### Step 3. plot samples in range chart - look at the raw data you've collected #####
# rangechart is a function currently available in the ichthyoliths package. 
# I have also included it in the sampling_bias_functions so you don't need to download 
# and install ichthyoliths if that's too much trouble. However, I strongly recommend
# installing ichthyoliths, because you'll get the help file for rangechart, which is super useful

rangechart(sample.matrix, reorder = 'lad.by.fad', count.breaks = c(0, 1, 2, 3, 4, 5, 6), 
           col.points = 'by.count', cols.vec ='viridis', 
           cex.points = 'by.count')
legend('topleft', horiz = TRUE, legend = c('1','2','3','4','5','6+'), col = viridis(6), 
       pt.cex = seq(1, 2, length.out = 6), pch = 16)
mtext(side = 1, text = 'species', line = 2.5)

##### Ecological Analyses ##### 
## Rarefaction curves
rarecurve(sample.matrix)

## shannon-weiner diversity index
diversity(sample.matrix)
plot(1:time.bins, diversity(sample.matrix))

## number of species (raw)
specnumber(sample.matrix)

##### Step 5. Foote rates #####

foote <- foote.metrics(sample.matrix)
foote.plot(foote)

##### Step 6. combine two distinctly different sample sets (e.g. an extinction!) across a horizon... #####

species.dist1 <- select.species.distribution(dist.type = 'lgser', theta = 0.8, hist = TRUE) # species 1:n
species.dist2 <- select.species.distribution(dist.type = 'norm') # species 1:n
species.dist3 <- select.species.distribution(dist.type = 'unif')

time.bins <- 20

#sample #1 - log-series distribution
sample.sizes <- sample(x = 10:100, size = time.bins, replace = TRUE)
sample.matrix1 <- build.sample.matrix(species.dist = species.dist1, time.slices = time.bins, sample.size = sample.sizes)

# sample #2 - normal distribution
sample.sizes <- sample(x = 10:100, size = time.bins, replace = TRUE)
sample.matrix2 <- build.sample.matrix(species.dist = species.dist2, time.slices = time.bins, sample.size = sample.sizes)

# combined matrix and analyses
sample.matrix <- combine.sample.matrices(sample.matrix2, sample.matrix1)
rangechart(sample.matrix, reorder = 'lad.by.fad', count.breaks = c(0, 1, 2, 3, 4, 5, 6), 
           col.points = 'by.count', cols.vec ='viridis', 
           cex.points = 'by.count')
rarecurve(sample.matrix)
foote <- foote.metrics(sample.matrix)
foote.plot(foote)




# sample #3 - log-series distribution, but scrambling the column names so it's different species being sampled
sample.sizes <- sample(x = 10:100, size = time.bins, replace = TRUE)
sample.matrix3 <- build.sample.matrix(species.dist = species.dist1, time.slices = time.bins, sample.size = sample.sizes)
colnames(sample.matrix3) <- sample(colnames(sample.matrix3))
colnames(sample.matrix3) <- as.numeric(colnames(sample.matrix3))+15


sample.matrix <- combine.sample.matrices(sample.matrix3, sample.matrix1)
rangechart(sample.matrix, reorder = 'lad.by.fad', count.breaks = c(0, 1, 2, 3, 4, 5, 6), 
           col.points = 'by.count', cols.vec ='viridis', 
           cex.points = 'by.count')
rarecurve(sample.matrix)
foote <- foote.metrics(sample.matrix)
foote.plot(foote)


##### 7. Rarefaction example - sampling the same distribution, but one with 10 samples, and one with 100 samples #####

species.dist <- select.species.distribution(dist.type = 'norm') # species 1:n
time.bins <- 5

sample.matrix1 <- build.sample.matrix(species.dist = species.dist, time.slices = time.bins, sample.size = 10)
sample.matrix2 <- build.sample.matrix(species.dist = species.dist, time.slices = time.bins, sample.size = 100)

sample.matrix <- combine.sample.matrices(sample.matrix2, sample.matrix1)
rangechart(sample.matrix, reorder = 'lad.by.fad', count.breaks = c(0, 1, 2, 3, 4, 5, 6), 
           col.points = 'by.count', cols.vec ='viridis', 
           cex.points = 'by.count')
rarecurve(sample.matrix)


##### 8. rarefaction example: sampling the same number of fossils, but completely different distributions
species.dist1 <- select.species.distribution(dist.type = 'norm', sd = 2) # species 1:n
species.dist2 <- select.species.distribution(dist.type = 'norm', sd = 10) # species 1:n
time.bins <- 5

sample.matrix1 <- build.sample.matrix(species.dist = species.dist1, time.slices = time.bins, sample.size = 100)
sample.matrix2 <- build.sample.matrix(species.dist = species.dist2, time.slices = time.bins, sample.size = 100)
sample.matrix <- combine.sample.matrices(sample.matrix2, sample.matrix1)
rangechart(sample.matrix, reorder = 'lad.by.fad', count.breaks = c(0, 1, 2, 3, 4, 5, 6), 
           col.points = 'by.count', cols.vec ='viridis', 
           cex.points = 'by.count')
rarecurve(sample.matrix)

# error bars
plot(specaccum(sample.matrix))
