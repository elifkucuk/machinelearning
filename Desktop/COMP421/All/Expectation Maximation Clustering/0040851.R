library(MASS)
library(mvtnorm)

set.seed(521)

##helper methods
my_plot<-function(X, centroids, assignments){
  if(is.null(assignments)){
    plot(X[,1], X[,2], type = "p", pch = 19, col = "lightgray", las = 1,
         xlim = c(-6, 6), ylim = c(-6, 6),
         xlab = "x1", ylab = "x2")
  }

  else{
    renk<-c("red","green","blue","pink","yellow")
    plot(X[,1], X[,2], type = "p", pch = 19, col =renk[assignments], las = 1,
         xlim = c(-6, 6), ylim = c(-6, 6),
         xlab = "x1", ylab = "x2")
  }
  if(!is.null(centroids)){
    points(centroids[,1], centroids[,2], col = "black", pch = 19, cex = 1) 
  }
}

add <- function(x) Reduce("+", x)



## Question 1
X1 <- mvrnorm(50, c(+2.5, +2.5), matrix(c(.8, -0.6, -0.6, .8), 2, 2))
X2 <- mvrnorm(50, c(-2.5, +2.5), matrix(c(.8, 0.6, 0.6, .8), 2, 2))
X3 <- mvrnorm(50, c(-2.5, -2.5), matrix(c(.8, -0.6, -0.6, .8), 2, 2))
X4 <- mvrnorm(50, c(2.5, -2.5), matrix(c(.8, 0.6, 0.6, .8), 2, 2))
X5 <- mvrnorm(100, c(0, 0), matrix(c(1.6, 0.0, 0.0, 1.6), 2, 2))
X<-rbind(X1,X2,X3,X4,X5)

my_plot(X, NULL, NULL)

#Question 2


K<-5
N<-length(X[,1])

#In initialization we take choose random points.

updateAssignments<-function(){
  
  #dist method evaluates distance for each row to other rows and returns a matrix where each row holds
  # distances to all other rows. So 1 row has 305 column.  
  D <- as.matrix(dist(rbind(centroids, X), method = "euclidean"))
  #We want to keep only distances to dectroid. So We take the first 5 rows and from column 6 to the end of the columns 305 in our case
  D <- D[1:nrow(centroids), (nrow(centroids) + 1):(nrow(centroids) + nrow(X))]
  #For each column we find the least row value which will point to centroid that is closest to the point indicated by column index. 
  assignments <<- sapply(1:ncol(D), function(i) {which.min(D[,i])})
  my_plot(X,centroids, assignments)
  
}

updateCentroid_K_Means<-function(){
  
  centroids<<-t(sapply(seq(1:K), function(k){colMeans(X[assignments==k,])}))
  my_plot(X,centroids, assignments)
}

centroids<-X[sample(seq(1:N),5),]
my_plot(X,centroids, NULL)
updateAssignments()
updateCentroid_K_Means()
updateAssignments()

# Question 3

#Initialization for S and P using K-means
S<-t(sapply(seq(1:K), function(k){cov(X[assignments==k,])}))

P<-sapply(seq(1:K), function(k){sum(assignments==k)/N})


#Question 4


E_Step<-function(){
  
  h_numerator<-list()
  for(k in seq(1:K)){
   h_numerator[[k]]<-sapply(seq(1:N), function(i) {P[k]*(det(matrix(S[k,],2,2))^(-1/2))*exp(-(1/2)*(t(X[i,]-centroids[k,])%*%chol2inv(chol(matrix(S[k,],2,2)))%*%(X[i,]-centroids[k,])))})
  }

  h_denominator<-add(h_numerator)

  h<<-sapply(seq(1:K), function(k) {h_numerator[[k]]/h_denominator})
  
  #assignments<<-sapply(seq(1:N), function(i) {which.max(h[i,])})
  
  #P<<-sapply(seq(1:K), function(k){sum(assignments==k)/N})
  P<<-sapply(seq(1:K), function(k){sum(h[,k])/sum(h)})
}


M_Step<-function(h){
  centroids<<-t(sapply(seq(1:K), function(k){(h[,k]%*%X)/sum(h[,k])}))
  S<<-t(sapply(seq(1:K), function(k){colSums(t(sapply(seq(1:N), function(i){h[i,k]*(X[i,]-centroids[k,])%*%t(X[i,]-centroids[k,])})))/sum(h[,k])}))
}


runtime<-100


for(l in seq(1:runtime)){
  E_Step()
  
  M_Step(h)
  assignments<<-sapply(seq(1:N), function(i) {which.max(h[i,])})
  my_plot(X,centroids,assignments)
  
}



#Question 5


#Following soltion worked better for me (faster) then finding multivariate normal distr. on x1 and x2 grid.
#This was taken from web source https://stats.stackexchange.com/questions/9898/how-to-plot-an-ellipse-from-eigenvalues-and-eigenvectors-in-r

 for(k in seq(1:K)){
 mu    <- centroids[k,]                               
 S_i      <- matrix(S[k,],2,2)                        
 decom     <- chol(S_i)                               
 t <- seq(0, 2*pi, length.out=200)                    
 elips    <- 1 * cbind(cos(t), sin(t)) %*% decom      
 center_elips <- sweep(elips, 2, mu, "+")             
 lines(center_elips, type="l", lwd=2, asp=1)          
 }


m1<-c(+2.5, +2.5)
S1<-c(.8, -0.6, -0.6, .8)
m2<-c(-2.5, +2.5)
S2<-c(.8, 0.6, 0.6, .8)
m3<-c(-2.5, -2.5)
S3<-c(.8, -0.6, -0.6, .8)
m4<- c(2.5, -2.5)
S4<-c(.8, 0.6, 0.6, .8)
m5<-c(0, 0)
S5<-c(1.6, 0.0, 0.0, 1.6)


m<-rbind(m1,m2,m3,m4,m5)

S_real<-rbind(S1,S2,S3,S4,S5)
for(k in seq(1:K)){
  
  mu    <- m[k,]                               # data centroid -> colMeans(dataMatrix)
  S_i      <- matrix(S_real[k,],2,2) # covariance matrix -> cov(dataMatrix)
  decom     <- chol(S_i)                               # Cholesky decomposition
  t <- seq(0, 2*pi, length.out=200)          # angles for ellipse
  elips    <- 1 * cbind(cos(t), sin(t)) %*% decom  # ellipse scaled with factor 1
  center_elips <- sweep(elips, 2, mu, "+")               # center ellipse to the data centroid
  lines(center_elips, type="l", lty=2, lwd=2, asp=1)            # plot ellipse
}

