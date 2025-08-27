#Load Necessary packages
install.packages('mapview')
install.packages('fda')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('RColorBrewer')
install.packages('cluster')
library(mapview)
library(fda)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(cluster)

#Reading the dataset
caspian<-read.csv("/Users/abhinavsingh/Downloads/caspian.csv")
View(caspian)
#Dimension of the dataset
dim(caspian)

####caspian####DATA TRANSFORMATION AND EXPLORATORY ANALYSIS####

#Check for missing values
sum(is.na(caspian))
sapply(caspian, function(x) sum(is.na(x)))

#subset of Caspian
set.seed(30)
caspianSub<- caspian %>% sample_n(50)

#map
mapview(caspianSub, xcol='x', ycol='y', crs=4269, grid=FALSE, col.regions='lightblue')
View(caspianSub[,4:408])

#Transforming the data set
tempt<- t(caspianSub[,4:408])
View(tempt)

#Yearly variation
daytime <- seq(from = 1995.5, to = 2012.30, length.out = 405)
daytime
par(mfrow=c(1,1))
matplot(daytime,tempt,type='l',xlab='Year', ylab='Temperature(Kelvin)',main="LSWT Time Series at Multiple Grid Locations (1995–2012)",axes=F)
axis(2)
axis(side=1, at=seq(1995,2012))

#Mean line
mean_line <- rowMeans(tempt)

# Mean plot on LSWT 
matplot(daytime,tempt,type='l',xlab='Year', ylab='Temperature(Kelvin)',main="Mean Lake Surface Water Temperature",col='lightblue' ,axes=F)
axis(2)
axis(side=1, at=seq(1995,2012))
lines(daytime, mean_line,col = "red", lwd = 2, type='l')
#line for interpretative purpose
abline(h=279.15, col='blue', lty='dashed')

#Seasonal variation(yearly variation) 1996
season<-seq(from=1, to=365, length.out=24)
season
#subset of 1996 from tempt
tempt[15:38,]
matplot(season,tempt[15:38,],xlab='Days', ylab='Temperature(Kelvin)',main="Seasonal variation year 1996",type='l',axes=F)
axis(2)
axis(side=1, at=seq(1,365, by=5))
#horizontal line to compare seasonal variation
abline(h=295, col='red', lty='dashed')


#Seasonal variation(yearly variation) 2005
matplot(season,tempt[231:254,],xlab='Days', ylab='Temperature(Kelvin)',main="Seasonal variation year 2005",type='l',axes=F)
axis(2)
axis(side=1, at=seq(1,365, by=5))
abline(h=295, col='red', lty='dashed')


#Seasonal variation(yearly variation) 2009
tempt[327:350,1]
matplot(season,tempt[327:350,],xlab='Days', ylab='Temperature(Kelvin)',main="Seasonal variation year 2009",type='l',axes=F)
axis(2)
axis(side=1, at=seq(1,365, by=5))
abline(h=295, col='red', lty='dashed')

#Seasonal variation(yearly variation) 2010
tempt[351:374,1]
matplot(season,tempt[351:374,],xlab='Days', ylab='Temperature(Kelvin)',main="Seasonal variation year 2010",type='l',axes=F)
axis(2)
axis(side=1, at=seq(1,365, by=5))
abline(h=295, col='red', lty='dashed')


####FORMAL ANALYSIS#####


#Fit without smoothing
#For Saturated B-spline as many breaks as datapoints
num_breaks<- length(season)-4+2
num_breaks
breaks <- seq(1, 365, length.out = num_breaks)
breaks
dayrng = range(season)

#APPROACH-1
#Bbasis fit by specifying number of basis

#Smaller number of basis for representation
smallerBasis = create.bspline.basis(dayrng,nbasis=5,norder=2)
#Increasing number of basis for representation
IncreasedBasis = create.bspline.basis(dayrng,nbasis=8,norder=4)

#####Saturated fit for our research#####
bbasis = create.bspline.basis(dayrng,nbasis = 22,norder=4)
#number of basis applied can be looked using:
bbasis$nbasis
#Knots applied at positions can be looked using:
bbasis$params

#plot saturated basis
plot(bbasis)
#plot w smaller basis
plot(smallerBasis, main='B-spline fit with 5 basis and degree 1')
#plot w cubic poly and increased number of basis
plot(IncreasedBasis, main='B-spline fit with 8 basis and degree 3')


# Smoother fit for one pixel for smaller number of basis (Season 1996)
tempSmooth0 = smooth.basis(season,tempt[15:38,1],smallerBasis)
smooth_vals <- eval.fd(season, tempSmooth0$fd)

y_all <- c(tempt[15:38, 1], smooth_vals)
plot(season, tempt[15:38, 1], 
     ylim = range(y_all), xlab='Days', ylab='Temperature(Kelvin)',main="B-spline fit with 5 basis") 
lines(season, smooth_vals, col = "lightblue", lwd = 2)


# Smoother fit for one pixel for Increased number of basis (Season 1996)
tempSmooth00 = smooth.basis(season,tempt[15:38,1],IncreasedBasis)
smooth_vals <- eval.fd(season, tempSmooth00$fd)

y_all <- c(tempt[15:38, 1], smooth_vals)
plot(season, tempt[15:38, 1], 
     ylim = range(y_all), xlab='Days', ylab='Temperature(Kelvin)',main="Cubic spline fit, Increased number of  basis") 
lines(season, smooth_vals, col = "lightblue", lwd = 2)


# Cubic spline fit for one pixel with saturated basis (Season 1996)
tempSmooth000 = smooth.basis(season,tempt[15:38,1],bbasis)
smooth_vals <- eval.fd(season, tempSmooth000$fd)

y_all <- c(tempt[15:38, 1], smooth_vals)
plot(season, tempt[15:38, 1], 
     ylim = range(y_all), xlab='Days', ylab='Temperature(Kelvin)',main="Cubic spline fit for one pixel with saturated basis")  
lines(season, smooth_vals, col = "lightblue", lwd = 2)

########################################
########################################
###MAIN FOR FURTHER FORMAL ANALYSIS#####
########################################
########################################
###Smoother for all pixels chosen for research (Season 1996)###
tempSmooth1 = smooth.basis(season,tempt[15:38,],bbasis)
# Overlay smooth curve
plot(tempSmooth1, xlab='Days', ylab='Temperature(Kelvin)',main="B-spline fit w/o penalty(Year-1996)" )

#Smoothing
#Optimal value for lambda(Penalty coefficient)
set.seed(30)
lambdas = seq(70,120)    #'  lambdas to look over
lambdas

mean.gcv = rep(0,length(lambdas))
curv.Lfd = int2Lfd(2)

for(ilam in 1:length(lambdas)){
  #'  Set lambda
  curv.fdPari = fdPar(bbasis,curv.Lfd,ilam)
  curv.fdPari$lambda = lambdas[ilam]
  #'  Smooth
  tempSmoothi = smooth.basis(season,tempt[15:38,],curv.fdPari)
  #'  Record average gcv
  mean.gcv[ilam] = mean(tempSmoothi$gcv)
}
#'  We can plot what we have
par(mfrow=c(1,1))
plot(lambdas,mean.gcv,type='b',log='x', xlab='Lambda', ylab='GCV',main="Genralised Cross Validation")

#'  Lets select the lowest of these and smooth

best = which.min(mean.gcv)
best
lambdabest = lambdas[best]
lambdabest
abline(v=lambdabest, lty='dashed', col='red')

curv.fdPar = fdPar(bbasis,curv.Lfd, lambdabest)
curv.fdPar$lambda
#Final smoothing with optimal lambda value
tempSmooth1 = smooth.basis(season,tempt[15:38,],curv.fdPar)
plot(tempSmooth1$fd, xlab='Days', ylab='Temperature(Kelvin)',main="B-spline fit with smoothing penalty(Year-1996)")

#Visual of the Mean on the smoother
plot(tempSmooth1$fd,xlab='Days', ylab='Temperature(Kelvin)',main='Mean smoothed curve in red', col=4)
lines(mean.fd(tempSmooth1$fd), col=2, lwd=4)


#Residual check, not necessarily required checking just for sanity/good practice
y_hat <- eval.fd(season, tempSmooth1$fd)
# Pointwise residuals (same shape)
res_mat <- tempt[15:38,] - y_hat  
#Spaghetti Plot shows random noise 
matplot(season, res_mat, type = "l", lty = 1, col = "grey70",
        xlab = "Time(Days)", ylab = "Residual", main='Raw Residuals')
abline(h = 0, lty = 2)
##(random spikes arouond zero -no systematic structure i.e not oversmoothing)


#Contour plot 
temp.cor = cor.fd(season,tempSmooth1$fd)
filled.contour(season,season,temp.cor,xlab='Time(Days)',ylab='Time(Days)', main="Functional correlation surface plot")


###Functional PCA###
tempfd<- tempSmooth1$fd
tempPCA = pca.fd(tempfd,nharm=5)

#Mathematical check
cumsum(tempPCA$varprop)

#Graphical Check
par(mfrow=c(1,1))
plot(cumsum(tempPCA$varprop),type='b', xlab="Number of components", 
     ylab="variation explained", main="Proportion of variation explained by individual components")
abline(h=0.95, col='red', lty='dashed')

#looking at PCA
plot(tempPCA$harmonics[1:2], main="First 2 principal components")

par(mfrow=c(1,2))
plot(tempPCA,harm=1:2)

#Storing score in a new variable
scores<-tempPCA$scores


### K-Means Clustering Algorithm ###

#Elbow plot for optimal K
set.seed(123)

# Create vector to store WCSS for each k
wcss <- numeric()

# Try k from 1 to 10
for (k in 1:10) {
  km <- kmeans(scores, centers = k, nstart = 125)
  wcss[k] <- km$tot.withinss
}

# Plot the elbow plot
par(mfrow=c(1,1))
plot(1:10, wcss, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow Method for Choosing Optimal k")

##Cluster with 3 centers
set.seed(123)
clusters <- kmeans(scores, centers = 3, nstart = 25)$cluster

#visual rep
caspianSub$cluster <- as.factor(clusters)
clusters <- factor(caspianSub$cluster)

# Color pallete for 5 clusters
cluster_colors <-brewer.pal(3, "Set2")

# Map colors to cluster levels
names(cluster_colors) <- levels(clusters)

#Mapping on spatial grid - interpretative view
mapview(caspianSub, xcol='x', ycol='y',zcol="cluster",col.regions=cluster_colors[as.character(clusters)], crs=4269, grid=FALSE)

##Final plot ( time series) (colors match accurately with mapview above)
plot(tempSmooth1$fd, col=cluster_colors[as.character(clusters)], axes=T, xlab = 'Days', ylab = "Temperature(Kelvin)")

# manual legend labels, will produce different results when run again
legend("topleft",                  
       legend = c("South", "North", "Central"),     
       col = cluster_colors,        # same colors as plot
       lty = 1,                     # line type
       cex = 0.7, lwd = 6)                   # size


##Cluster with 5 centers
set.seed(123)
clusters <- kmeans(scores, centers = 5, nstart = 25)$cluster

#visual rep
caspianSub$cluster <- as.factor(clusters)
clusters <- factor(caspianSub$cluster)

# Color pallete for 5 clusters
cluster_colors <-brewer.pal(5, "Set2")

# Map colors to cluster levels
names(cluster_colors) <- levels(clusters)

#Mapping on spatial grid - interpretative view
mapview(caspianSub, xcol='x', ycol='y',zcol="cluster",col.regions=cluster_colors[as.character(clusters)], crs=4269, grid=FALSE)

##Final plot ( time series)
plot(tempSmooth1$fd, col=cluster_colors[as.character(clusters)], axes=T, xlab = 'Days', ylab = "Temperature(Kelvin)")
legend("topleft",                  # position of the legend
       legend = c("South/Southwest", "North", "Central/Central East", "Central/CentralWest", "South/SouthEast"),     # cluster labels
       col = cluster_colors,        # same colors as plot
       lty = 1,                     # line type
       cex = 0.7, lwd = 6)                   # size


####APPENDIX#####

##Seasonal Mann-Kendall test

# Example for first pixel column
pixel_series <- tempt[, 14]  # change to your actual pixel column

# Create ts object with correct seasonal frequency
lswt_ts <- ts(pixel_series, start = c(1995, 1), frequency = 24) 

# Run Seasonal Mann–Kendall test
res <- smk.test(lswt_ts, alternative = "two.sided", continuity = TRUE)
summary(res)

# Function to run Seasonal Mann-Kendall test for one pixel
run_smk <- function(pixel_series) {
  ts_data <- ts(pixel_series, start = c(1995, 1), frequency = 24)
  
  res <- tryCatch({
    smk.test(ts_data, alternative = "two.sided", continuity = TRUE)
  }, error = function(e) return(NA))  # return NA if error
  
  return(res)
}

# Apply to all pixel columns
smk_results <- apply(tempt, 2,  run_smk)

p_values <- sapply(smk_results, function(r) ifelse(is.list(r), r$p.value, NA))
tau_values <- sapply(smk_results, function(r) ifelse(is.list(r), r$tau, NA))

# Combine into a results data frame
smk_summary <- data.frame(
  Pixel = caspianSub[,1],
  Tau = tau_values,
  P_value = p_values
)

sum(smk_summary$P_value<0.05)
##17 locations came out to be significant

