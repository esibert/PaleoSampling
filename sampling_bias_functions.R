##################################################################
#
#       Functions for paleontological sampling exercies
#
#
###################################################################

#### 1. Create null distribution of species-abundances
# Choose a distribution, sample it, and make the lowest species number equal to 1. 
# dist types currently in this are: normal (norm), uniform (unif), lognormal (lnorm), log-series (lgser)
# there are default parameters for each distribution, however you can change those specifically
#   in the function call

select.species.distribution <- function(dist.type, n=1e6, histogram = TRUE, ...) {
    
    args <- list(...)
    
    # normal distribution
    if(dist.type == 'normal' | dist.type == 'norm') {
        if(is.null(args[['mean']])) {args$mean <- 20}
        if(is.null(args[['sd']])) {args$sd <- 8}
        species.dist <- as.integer(rnorm(n=n, mean = args$mean, sd = args$sd))
    }
    
    # uniform distribution
    if(dist.type == 'uniform' | dist.type == 'unif') {
        if(is.null(args[['min']])) {args$min = 1}
        if(is.null(args[['max']])) {args$max = 100}
        species.dist <- as.integer(runif(n=n, min = args$min, max = args$max))
    }
    
    # lognormal distribution
    if(dist.type == 'lognormal' | dist.type == 'lnorm') {
        if(is.null(args[['meanlog']])) {args$meanlog = 2}
        if(is.null(args[['sdlog']])) {args$sdlog = .5}
        species.dist <- as.integer(rlnorm(n=n, meanlog = args$meanlog, sdlog = args$sdlog))
    }
    
    # log-series distribution
    # From the literature, log-series distribution is the best fit for species abundance curves: 
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5183127/ 
    # "An extensive comparison of species-abundance distribution models", Baldridge et al 2016
    
    if(dist.type == 'log-series' | dist.type == 'lgser') {
        if(is.null(args[['theta']])) {args$theta = 0.8}
        species.dist <- rlgser(n=n, theta = args$theta) #note don't need to convert these to integers, already are...
    }
    
    
    ## Clean up species dist to have minimum species number = 1
    add.spp <- (-min(species.dist)) + 1
    species.dist <- species.dist + add.spp
    
    if(histogram == TRUE) {hist(species.dist)}
    
    return(species.dist)
}

#### 2. Sample the distribution
# Sample (without replacement) within that distribution at set time points to create "occurrance matrix" 
# note that can sample with constant sample size at each tiem bin or with variable sample size, in 
#    the form of a vector, length time.slices
build.sample.matrix <- function(species.dist = species.dist, time.slices = 10, sample.size = 50) {
    # create matrix of 0's
    sample.matrix <- matrix(data = 0, nrow = time.slices, ncol = max(species.dist), 
                            dimnames = list(c(1:time.slices),c(min(species.dist):max(species.dist))))
    
    # Account for variability in sample size input: 
    if(length(sample.size) == 1) {sample.size <- rep(sample.size, time.slices)}
    if(length(sample.size) != time.slices) stop ('variable sample size does not match time slices')
    
    # populate the matrix by sampling each time slice
    for(i in 1:time.slices) {
        fossil.sample <- sample(species.dist, sample.size[i])
        sample.matrix[i,] <- tabulate(fossil.sample, nbins = max(species.dist))
    }
    
    # remove the zero-columns
    sample.matrix <- sample.matrix[, colSums(sample.matrix != 0) > 0]
    
    # return the sample matrix for analysis
    return(sample.matrix)
}

## 2a - combine multiple distinct sample sets, e.g. ecological shift
combine.sample.matrices <- function(mat1, mat2) { #need to fix time slices later
    mat1 <- as.data.frame(mat1)
    mat2 <- as.data.frame(mat2)
    bind <- gtools::smartbind(list = list(mat1, mat2), fill = 0)
    rownames(bind) <- 1:length(bind[,1])
    bind <- as.matrix(bind)
    return(bind)
}


##### 3. sample.age.ranges
# creates LAD and FAD data frame for species sampled

sample.age.ranges <- function(sample.matrix) {
    depths <- rownames(sample.matrix)
    
    index.pos <- function(x) {
        return((1:length(x))[x > 0])
    }
    
    presence.index <- apply(sample.matrix, 2, index.pos)
    
    fads <- as.numeric(depths[as.numeric(lapply(presence.index, max))])
    lads <- as.numeric(depths[as.numeric(lapply(presence.index, min))])
    
    # give the species names appropriately
    
    ad <- data.frame(fads, lads)
    
    rownames(ad) <- colnames(sample.matrix)
    
    return(ad)
}


##### 4. Foote metrics calculations 


foote.metrics <- function(sample.matrix) {
    
    # calculate lad/fad datums
    ad <- sample.age.ranges(sample.matrix)
    time.slices <- as.numeric(rownames(sample.matrix))
    
    dt<-c()                           # change in time (interval)
    for(i in 1:(length(time.slices))-1) {
        dt[i]<-(time.slices[i+1]-time.slices[i]) }
    
    
    Nfl<-c()  # singletons (no boundary crossers)
    Nbl<-c()  # bottom crossers only
    Nft<-c()  # top crossers only
    Nbt<-c()  # both bottom and top crossers
    for (i in 1:length(time.slices)) {
        #subset of all species present or assumed present at the time point [i]
        ad.sub<-subset(ad, fads >= time.slices[i] & lads <= time.slices[i])
        #how many of each type of occurrance are present? 
        Nfl[i] <- length(subset(ad.sub, fads==lads)[,1])  # "singletons" (no boundary cross)
        Nbl[i] <- length(subset(ad.sub, fads>=time.slices[i] & lads==time.slices[i] & fads!=lads)[,1]) # bottom crossers only
        Nft[i] <- length(subset(ad.sub, fads==time.slices[i] & lads<=time.slices[i] & fads!=lads)[,1]) # top crossers only 
        Nbt[i] <- length(subset(ad.sub, fads> time.slices[i] & lads< time.slices[i])[,1]) # both top and bottom  crossers
    }
    
    ## calculate singletons, crossers present in each time bin - yay for for-loops?
    
    # b) calculate relevant metrics
    Ntot <- apply(cbind(Nfl, Nbl, Nft, Nbt), 1, sum) # total Diversity observed
    Nb   <- apply(cbind(Nbl, Nbt), 1, sum)           # All bottom boundary crossers
    Nt   <- apply(cbind(Nft, Nbt), 1, sum)           # all top boundary crossers
    No   <- apply(cbind(Nfl, Nft), 1, sum)           # number of originations
    Ne   <- apply(cbind(Nfl, Nbl), 1, sum)           # number of extinctions
    Ndiv <- apply(cbind(Nb, Nt), 1, function(x) sum(x)/2)     # estimated mean standing diversity
    pp   <- apply(cbind(Nbt, Nt, c(dt,0)), 1, function(x) {-log(x[1]/x[2]) / x[3]}) # per capita origination
    qq   <- apply(cbind(Nbt, Nb, c(0,dt)), 1, function(x) {-log(x[1]/x[2]) / x[3]}) # per capita extinction
    
    # for plotting origination/extinction estimates
    ylims <- na.exclude(c(pp, qq))
    ylims <- ylims[is.finite(ylims)]
    ylims <- c(0, max(ylims))
    
    foote.list <- list(Ntot = Ntot, Ndiv = Ndiv, pp = pp, qq = qq, ylims = ylims, time.slices = time.slices)
    
    return(foote.list)
}

foote.plot <- function(foote) {
    foote <- foote.metrics(sample.matrix)
    plot(foote$time.slices, foote$pp, col='blue', pch=16, type='o',  ylim = foote$ylims, 
         main='Foote 2000 boundary crosser extinction and origination', 
         xlab = 'ages', ylab = 'origination/extinction rate')  #AgeID.unique
    points(foote$time.slices, foote$qq, col='red', pch=16, type='o')
    legend ('topright', legend=c('origination', 'extinction'), col=c('blue', 'red'), lty=1, pch=16)
    
    # add diversity estimated
    par(new=T)
    plot(foote$time.slices, foote$Ndiv, axes=F, type = 'o', xlab = '', ylab = '')
    axis(4)
    
    # add diversity total
    points(foote$time.slices, foote$Ntot, pch=16)
}

##### 5. rangechart function from ichthyoliths #####
# rangechart <- function(counts, ages = NULL, taxa = NULL, tax.cat = NULL, reorder = NULL,
#                        normalize.counts = FALSE, count.breaks = NULL,
#                        cex.xaxis = 1, cex.yaxis = 1, yaxis.ticks = FALSE,
#                        llwd = 1, llcol = 'gray70', llty = 3,
#                        baselines = FALSE, blwd = 0.5, blcol = 'lightblue', blty = 3,
#                        col.points = 'gray70', cols.vec = NULL,
#                        pch.points = 16, pch.vec = NULL,
#                        cex.points = 1, largesize = 1,
#                        xaxis.labels = c('names', 'numeric', 'alphanum'), print.xaxis = FALSE, ...) {
#     
#     ##### set up the dataset #####
#     # turn a data.frame into a matrix if necessary.
#     if(class(counts) == 'data.frame') {
#         counts <- as.matrix(counts)
#     }
#     
#     
#     # Ages should be rownames of the counts table, and in increasing order
#     if(missing(ages)) {
#         ages <- as.numeric(rownames(counts))
#     }
#     else {
#         rownames(counts) <- ages
#     }
#     
#     # if the ages are not in increasing order, sort them and the counts table to be so
#     if(is.unsorted(ages) == TRUE) {
#         age.increasing <- sort(ages, index.return = TRUE)$ix
#         counts <- counts[age.increasing, ]
#         ages <- as.numeric(rownames(counts))
#     }
#     
#     # taxa should be column-names of the counts table, order doesn't matter.
#     if(missing(taxa)) {
#         taxa <- as.character(colnames(counts))
#     }
#     else {
#         colnames(counts) <- taxa
#     }
#     
#     original.taxa <- taxa #useful for matching tax-cat later too.
#     
#     
#     # clear NA values, if any, by replacing with zeros
#     if (sum(is.na(counts)) > 0) {
#         warning(paste(sum(is.na(counts)), "missing values in count matrix replaced with zeros"))
#         counts[is.na(counts)] <- 0
#     }
#     
#     # FAD: First (oldest) occurance datum calls the maximum index (mapped to the ages values)
#     # of a non-zero count value for each taxa column in the counts matrix
#     fad <- ages[apply(counts, 2, function(x) {max(which (x!=0))})]
#     
#     # LAD: Last (youngest) occurance datum calls the minimum index (mapped to the ages values)
#     # of a non-zero count value for each taxa column in the counts matrix
#     lad <- ages[apply(counts, 2, function(x) {min(which (x!=0))})]
#     
#     
#     ### Normalize the counts if they want to be normalized
#     if(normalize.counts == TRUE) {
#         norm.row<-function(row) {
#             row/sum(row)
#         }
#         norm.counts<-apply(counts, 1, norm.row)  #normalize the matrix by rows
#         counts<-t(norm.counts)   #for some reason I have to transpose the output back to the normal "counts" form
#         counts <- 100 * counts #make this a percentage instead of a decimal value.
#     }
#     
#     ### Group the counts into bins if you'd like them to be binned.
#     if(!is.null(count.breaks)) {
#         count.breaks <- c(count.breaks, max(counts)+1)
#         
#         for(i in 1:length(count.breaks)-1) {
#             counts[counts > count.breaks[i] & counts <= count.breaks[i+1] ] = i
#         }
#         
#     }
#     
#     ##### reorder counts #####
#     if(!is.null(reorder)) {
#         
#         # fad.by.lad (origination)
#         if(reorder == 'fad.by.lad') {
#             # First reorder the counts by LAD and recalculate FAD
#             reorder.vect <- sort(lad, decreasing = TRUE, index.return = TRUE)$ix #pulls index of order by fads
#             counts <- counts[, reorder.vect]
#             fad <- ages[apply(counts, 2, function(x) {max(which (x!=0))})]
#             
#             # Next, reorder the counts by FAD
#             reorder.vect <- sort(fad, decreasing = TRUE, index.return = TRUE)$ix #pulls index of order by fads
#             counts <- counts[, reorder.vect]
#             
#         }
#         
#         # lad.by.fad (extinction)
#         else if(reorder == 'lad.by.fad') {
#             # First reorder the counts by FAD and recalculate LAD
#             reorder.vect <- sort(fad, decreasing = TRUE, index.return = TRUE)$ix #pulls index of order by fads
#             counts <- counts[, reorder.vect]
#             lad <- ages[apply(counts, 2, function(x) {min(which (x!=0))})]
#             
#             # Next, reorder the counts by LAD
#             reorder.vect <- sort(lad, decreasing = TRUE, index.return = TRUE)$ix #pulls index of order by fads
#             counts <- counts[, reorder.vect]
#         }
#     }
#     
#     
#     ### re-generate 'taxa', 'ages', 'fad' and 'lad' and 'tax.cat' from the updated counts table
#     
#     taxa <- as.character(colnames(counts))
#     ages <- as.numeric(rownames(counts))
#     fad <- ages[apply(counts, 2, function(x) {max(which (x!=0))})]
#     lad <- ages[apply(counts, 2, function(x) {min(which (x!=0))})]
#     
#     if(!is.null(tax.cat)) { tax.cat <- tax.cat[match(taxa, original.taxa)] }
#     
#     ##### set up the graphical parameters #####
#     ## xaxis.labels
#     if(missing(xaxis.labels)) { xaxis.labels <- 'names' }
#     if(xaxis.labels == 'names') { xaxis.lab <- taxa }
#     if(xaxis.labels == 'numeric') { xaxis.lab <- 1:length(taxa) }
#     if(xaxis.labels == 'alphanum') { xaxis.lab <- match(taxa, original.taxa) }
#     
#     ## colors of points
#     if(missing(cols.vec)) {cols.vec = 'gray'}
#     
#     if(cols.vec[1] == 'rainbow') {
#         colors<-rainbow(max(counts), end=5/6)
#     }
#     else if (cols.vec[1] == 'viridis') {
#         colors <- viridis::viridis(max(counts))
#     }
#     else colors<-cols.vec
#     
#     ## sizes of points (and scale of the whole thing...)
#     sizes<-seq(0,(max(counts)-1), 1)
#     sizes<-((sizes/(max(sizes)/largesize)) + largesize)
#     
#     
#     ##### Actually make the plot #####
#     ### make blank plot with appropriate dimensions, suppress x- and y- axes
#     plot(1:ncol(counts), ylim = c(max(ages), min(ages)),
#          type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "Age (Ma)", ...)
#     
#     ### add segments to the plot
#     segments(1:ncol(counts), lad, 1:ncol(counts), fad,
#              lwd = llwd, col = llcol, lty = llty, ...)
#     
#     if (baselines == TRUE) {
#         segments(1:ncol(counts), fad, 1:ncol(counts), rep(par()$usr[3], ncol(counts)),
#                  col = blcol, lty = blty, lwd = blwd, ...)
#     }
#     
#     ### Add points to the plot
#     for (i in 1:ncol(counts)) {
#         num.val<-c(counts[,i][counts[,i]>0])
#         
#         # point characters
#         if (pch.points == 'by.count') {
#             pts.pch<-pch.vec[num.val]
#         }
#         else if (pch.points == 'by.category') {
#             if(class(tax.cat)=='factor') {tax.cat <- as.numeric(tax.cat)}
#             pts.pch <- pch.vec[tax.cat[i]]
#         }
#         else { pts.pch<-pch.points }
#         
#         # point colors
#         if (col.points == 'by.count') {
#             pts.cols<-c(colors[num.val])
#         }
#         else if (col.points == 'by.category') {
#             pts.cols <- colors[tax.cat[i]]
#         }
#         else { pts.cols<-col.points }
#         
#         # point size
#         if (cex.points == 'by.count') {
#             pts.cex<-c(sizes[num.val])
#         }
#         else if(cex.points == 'by.category') {
#             pts.cex <- sizes[tax.cat[i]]
#         }
#         else { pts.cex<-cex.points }
#         
#         # point y-values
#         plocs <- ages[(counts > 0)[, i]]
#         
#         #actually add the points
#         points(rep(i, length(plocs)), plocs, cex=pts.cex, pch=pts.pch, col=pts.cols,
#                ...)
#     }
#     
#     ### add axes
#     axis(1, at = 1:ncol(counts), cex.axis = cex.xaxis, labels = xaxis.lab,
#          las = 3)
#     axis(2, las = 1, cex.axis = cex.yaxis)
#     if(yaxis.ticks == TRUE) {axis(2, at = ages, labels = FALSE, tck = -0.01)}
#     
#     
#     ##### print taxa in list #####
#     # this has to be the last thing that the function does, because R stops after a return value
#     if(print.xaxis == TRUE) {
#         if(xaxis.labels == 'alphanum') {
#             return(original.taxa)
#         }
#         else {
#             return(taxa)
#         }
#     }
#     
# }
