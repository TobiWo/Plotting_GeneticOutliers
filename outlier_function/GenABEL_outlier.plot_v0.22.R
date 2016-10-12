#' Description
#' This function can be used for every cmd- or principal component matrix.
#' However, it was designed as a supporting function for the GWAS-workflow with GenABEL.
#' It is important to have the sample names or ID names within the rownames.
#' The function plots different mds plots to verify genetic outliers obtained in a PC1 vs. PC2 plot
#' Further it returns a vector of IDs which represent the genetically homogen population for a GWAS
#' The input parameters are as follows:
#' @param x         a matrix of numeric values (for the GenABEL-vignette it is a object returned by cmdscale)
#' @param clusters  the number of clusters obtained from a initial MDS-plot (including the main cluster)
#'                  the function will mark the clusters in every single plot based on the outlier-clusters from the PC1 vs. PC2 plot
#' @param n.pc       the number of pcas you want to plot
#'                  e.g. if you choose 4 pcas the function will return 6 plots which are namely:
#'                  PC1 vs. PC2; PC1 vs. PC3; PC1 vs. PC4; PC2 vs. PC3; PC2 vs. PC4; PC3 vs. PC4 
#'                  Note: If the user want to plot the combinations of 5 or more PCAs, the labels of the x- and y-axis are
#'                        disabled because of space issues. The order of plots is always from left to right and top to bottom.
#' @param plot.path under default conditions all plots will be merged to one single plot
#'                  e.g. if you have 6 pcas to display you will receive a plot with 3x2 PC-plots
#'                  if you set a path, there will be no merging, instead every single plot will be saved as pdf to the supplied path
#'                  Note: It is recommended to store the plots instead of merging while plotting more than 7 PCAs
#' @return case1    multiple PC-plots
#'         case2    multiple PC-plots and IDs of the main-cluster which should be kept in the analysis   


plot.GeneticOutlier <- function(x, clusters, n.pc, return.main=FALSE, plot.path=NULL) {
  
  # Initial check, whether enough PCs has been computed
  if(n.pc>ncol(GWAS_test.mds)){
    stop("Dataset is smaller than supplied number of n.pcs!")
  }

  # saving the original plot parameters for margins and outer margins
  org_mar <- c(5.1,4.1,4.1,2.1)
  org_oma <- c(2,2,2,2)
  
  # save x to object data --> because later the function uses x for coordinates and therefore x will be overwritten
  data <- x
  
  # chech whether plot.path is null or not and create new boolean which makes it easier for the following if-statements
  if(is.null(plot.path)){
    plot.path_bool <- FALSE
  } else if (!is.null(plot.path)){
    plot.path_bool <- TRUE
  }
  
  # build the cluster object with kmeans() based on 10000 random checks
  # assign clusters automatically to appropriate names
  # clusters are now visible within the function
  km <- kmeans(data, centers = clusters, nstart = 10000)
  for (i in 1:clusters) {
    assign(x = paste("cl", i, sep = ""), value = names(which(km$cluster==i)))
  }
  
  # the following two steps are necessary because the main cluster (cluster with the most samples) is not necessarily the last one
  # therefore the number of the main cluster is important for later plotting
  # creates a vector with the number of IDs in each cluster
  # "maximum" points to the position of main cluster which should include most samples
  laenge = NULL
  for(i in 1:clusters){
    l = length(x = get(paste("cl", i, sep = "")))
    laenge = c(laenge, l)
    maximum = which(laenge==max(laenge))
  }
  
  # if maximum equals the number of clusters you provide to the function, nothing particular will happen
  # otherwise the main cluster and the last cluster are changed with each other
  # e.g. main cluster == 3, provided cluster number and therfore the last cluster == 5 --> ...
  # ...the code changes main cluster to 5 and the last cluster (here 5) to 3 using temporary variables
  if(maximum==clusters){
    maximum==maximum
  }else {
    temp1 = get(x = paste("cl", maximum, sep = ""))
    temp2 = get(x = paste("cl", clusters, sep = ""))
    assign(x = paste("cl", clusters, sep = ""), value = temp1)
    assign(x = paste("cl", maximum, sep = ""), value = temp2)
  }

  # sets the plotting options (number of plots == n)
  # the automatic scaling of the number of plots, combined to one plot or separately saved in a folder, is based on combinatorics and particular on:
  # n.pc!/(n.pc-2)!*2! --> n.pc is the argument within the function which represent the number of principal components to plot 
  # the two comes from the x- and y-axis (2 axes)
  #
  # if 1>=n<=3, you will receive a n(row)*1(column) plot
  # else it checks for two conditions:
  # first: (round up of first factorial expression / 2) / (second factorial expression / 2) == 1 --> ...
  # ... if the expression == 1, the number of rows is based on the division: factorial expression / 2 and ...
  # ...subsequently the number of cols is based on the division: factorial expression / rows 
  #
  # second: (round up of first factorial expression / 2) / (second factorial expression / 2) !== 1 --> ...
  # ... if the expression != 1, the number of rows is based on the division: factorial expression / 3 and ...
  # ...subsequently the number of cols is based on the division: factorial expression / rows
  #
  # further, because of space reasons, for more than 4 plots the margins and outer margins are changed ...
  # ... this will eliminate the axes-labels
  # I recommend to separately save the plots if one want to plot more than 5 PCs
  
  if(n.pc<=3) {
    par(mfrow = c(n.pc,1))
  }else {
    if(ceiling((factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))/2)/((factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))/2)==1){
      rows = (factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))/2
      cols = (factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))/rows
      if(n.pc>4){
        par(oma=c(0.1,0.1,0.1,0.1), mar = c(1.5, 3.1, 1, 0.4))
      }else {
        par(mar=c(4.1,4.1,1.0,2.1))
      }
      par(mfrow = c(rows, cols))
    }else if (ceiling((factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))/2)/((factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))/2)!=1){
      rows = ceiling((factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))/3)
      cols = ceiling((factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))/rows)
      if(n.pc>4){
        # since the above "else if" is only true for more than 4 pcs, no further "else" is needed
        par(oma=c(0.1,0.1,0.1,0.1), mar = c(1.5, 3.1, 1, 0.4))
      }
      par(mfrow = c(rows, cols))
    }
  }
  
  # printing small messages depending whether a path for saving was supplied
  # if one was supplied, every single plot will be saved in the respective folder
  if(!plot.path_bool){
    cat("Merging plots...")
  }else{
    cat("Saving single plots...")
    par(mfrow = c(1,1))
  }
  
  # the main plotting functionality
  pc1 = 1
  pc2 = 2
  for(i in 1:(factorial(n.pc)/(factorial(n.pc-2)*factorial(2)))){
    
    # if a path was supplied every single plot will be stored as pdf
    if(plot.path_bool) {
      pdf(file = paste(plot.path, "/", "Genetic_outlier_test_PC", pc1, "_PC", pc2, "_" ,Sys.Date() ,".pdf",sep = ""))
    }
    
    # the plotting function of two PCs, the size of the plot is adapted since the naming of a cluster could lie out of the plotting area  
    plot(data[,c(pc1,pc2)], xlim = c(min(data[,pc1])+(0.2*min(data[,pc1])),max(data[,pc1])+(0.2*max(data[,pc1]))),
         ylim = c(min(data[,pc2])+(0.2*min(data[,pc2])), max(data[,pc2])+(0.2*max(data[,pc2]))), xlab = paste("PC", pc1, sep = ""), 
         ylab = paste("PC", pc2, sep = ""))
    
    # this for-loop marks the outlier samples as red dots
    # since it falls back to the above calculated clusters, the plotted outliers are always the same
    for(j in 1:length(km$cluster)){
      for(b in 1:(clusters-1)) {
        id <- names(km$cluster[j])
        cl = paste("cl", b, sep = "")
        if(id%in%get(cl)) {
          x <- data[which(rownames(data)==id),pc1]
          y <- data[which(rownames(data)==id),pc2]
          points(x = x, y = y, pch = 16, col = "red")
        }
      }
    }
    rm(x, y, j, id, cl, b)
    
    # adding the cluster name to the plot, based on the outlier which is in the middle of the cluster (middle regarding the ID-vector 'cl')
    for(i in 1:(clusters-1)) {
      cl = get(x = paste("cl", i, sep = ""))
      id = cl[round(length(cl)/2)]
      x <- data[which(rownames(data)==id),pc1]
      y <- data[which(rownames(data)==id),pc2]
      text(x = x, y = y+(y*0.15), labels = paste("cluster", i, sep = ""), col = "red")
    }
    
    # if a path was supplied, the actual plot is only saved if dev.off() was called
    if(plot.path_bool){
      dev.off()  
    }
    
    # increasing PC2
    pc2 = pc2+1
    
    # if PC2 getting bigger than the number of PCs supplied (n.pc), PC1 is increased.
    # PC2 is as well increased, based on the new defined PC1 
    if(pc2 > n.pc) {
      pc1 <- pc1+1
      pc2 <- pc1+1
    }
  }

  # sets the graphical parameters mar, oma and mfrow back to the default settings 
  par(mar = org_mar, oma = org_oma ,mfrow = c(1,1))
  
  # if return.main is true, the IDs which represent the main cluster are send back from the function
  # one can directly use this vector to eliminate the outliers from the data-set
  if(!return.main){
    cat("\n")
    return(NULL)
  }else if(return.main){
    return(get(x = paste("cl", clusters, sep = "")))
  }
}









