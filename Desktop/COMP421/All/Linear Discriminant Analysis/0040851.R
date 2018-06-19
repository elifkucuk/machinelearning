#read data in memory
training_digits <- read.csv("hw06_mnist_training_digits.csv", header = FALSE)
training_labels <- read.csv("hw06_mnist_training_labels.csv", header = FALSE)


# get X and y values for training digits
X_train <- as.matrix(training_digits)/255
y_train <- training_labels[,1]

#read data in memory
test_digits <- read.csv("hw06_mnist_test_digits.csv", header = FALSE)
test_labels <- read.csv("hw06_mnist_test_labels.csv", header = FALSE)


# get X and y values for training digits
X_test <- as.matrix(test_digits)/255
y_test <- test_labels[,1]


### Question 2 ####
linear_disc_analysis<-function(X, y, R){
  # get number of overall samples number of features, number of class and number of samples for each class
  N <- length(y)
  D <- ncol(X)
  K<-length(unique(y))
  N_i<-sapply(seq(1:K), function(k){sum(y==k)})
  

  m_before<-c()
  for(k in seq(1:K)){
    m<-colMeans(X[which(y==k),])
    m_before<-rbind(m_before,m)
  }
  
  #Helper add method is defined to help add the list members.
  add <- function(x) Reduce("+", x)
  
  
  #Within class scatter matrix is a list S.
  S<-list()
  for(k in seq(1:K)){
    l<-lapply(which(y==k), function(i){(X[i,]-m_before[k,])%*%t(X[i,]-m_before[k,])})
    S[[k]]<-add(l)
  }
  
  #Total within class scatter is S_w
  S_w<-add(S)
  
  #Overall mean is mean of all class
  overall_mean<-colMeans(m_before)
  
  
  #Between class scatter matrix is S_b
  S_bs<-lapply( seq(1:K), function(i) {N_i[i]*((m_before[i,]-overall_mean)%*%t(m_before[i,]-overall_mean))})
  S_b<-add(S_bs)
  
  # W is equal to largest eigenvector of (Sw)^-1*Sb 
  
  #To take the inverse of S_w we add  diagonal matric with value 1e-10
  A<-diag(rep(1e-10,784))
  S_w<-S_w+A
  S_w_in<-chol2inv(chol(S_w))
  
  #The vector we will take the eigenvector of
  val<-S_w_in%*%S_b 
  
  # calculate the eigenvalues and eigenvectors and take the largest eigenvector as W
  decomposition <- eigen(val, symmetric = TRUE)
  decomposition$values
  W<-decomposition$vectors[,1:R]
  
  # calculate two-dimensional projections for training set
  return (W)
  
}

### Question 3  ####

#Calculate Z values for training data set
W<-linear_disc_analysis(X_train,y_train, 2)
Z_train <- (X_train) %*% W


point_colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6")
# plot two-dimensional projections for training set
plot(Z_train[,1], Z_train[,2], type = "p", pch = 19, col = point_colors[y_train], cex = 0,
     xlab = "PC1", ylab = "PC2", las = 1)
text(Z_train[,1], Z_train[,2], labels = y_train %% 10, col = point_colors[y_train])

#Calculate Z values for test 
Z_test <- (X_test) %*% W
# plot two-dimensional projections for test set
plot(Z_test[,1], Z_test[,2], type = "p", pch = 19, col = point_colors[y_test], cex = 0,
     xlab = "PC1", ylab = "PC2", las = 1)
text(Z_test[,1], Z_test[,2], labels = y_test %% 10, col = point_colors[y_test])




#### Question 4 ####
euclidian<-function(x,y){
  D<-length(x)
  squared_dis<-sapply(seq(1:D), function(d){(x[d]-y[,d])^2})
  dist<-sqrt(rowSums(squared_dis))
  return (dist)  
}

# knn specifies the number of neighbors
knn <- 5

#N is the number of data points, K is the number of classes 
N <- length(y_test)
K<-length(unique(y_train))

#find_KNN function is defined to find the prediction for a point using the k-nearest neighbor with euclidian distance
find_KNN<-function(x, X, Y){
  index_of_nearest_k<-order(euclidian(x, X), decreasing = FALSE)[1:(knn)]
  p_head<-sapply(seq(1:K), function (k){(sum(Y[index_of_nearest_k]==k))/knn})
  return (which.max(p_head))
}

# for r from 1 to 9 we learn W using X_train and Y_train and predict for X_test. 
#And calculate accuracy. 
R<-seq(1:9)
accuracy<-c()
for(r in R ){
  W<-linear_disc_analysis(X_train,y_train, r)
  Z_test <- (X_test) %*% W
  Z_train<-(X_train) %*% W
  p_head <- apply(Z_test, 1, find_KNN, Z_train, y_train)
  accuracy<-cbind(accuracy, sum(p_head==y_test)/N)
}
plot(R, accuracy*100, type="o", pch=19, ylab="Classification Accuracy %", xlab="R")
