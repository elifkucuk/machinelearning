# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

library(MASS)

#set.seed(421)


# Parameters for the random data points to be generated.
class_size1<-100
class_size2<-100
class_size3<-100

class_means1<-c(0.0,1.5)
class_means2<-c(-2.5,-3.0)
class_means3<-c(2.5,-3.0)

class_cov1<-matrix(c(1,0.2,0.2,3.2),2,2)
class_cov2<-matrix(c(1.6,-0.8,-0.8,1.0),2,2)
class_cov3<-matrix(c(1.6,0.8,0.8,1.0),2,2)

#Generation of the random data points from three bivariate Gaussian densities.
points1 <- mvrnorm(n = class_size1, mu = class_means1, Sigma = class_cov1)

points2 <- mvrnorm(n = class_size2, mu = class_means2, Sigma = class_cov2)

points3 <- mvrnorm(n = class_size3, mu = class_means3, Sigma = class_cov3)


#y is the sample size*class size matrix. Each column belongs to one class type with 1 and 0 binary value
#indicating the class of the data point.
y<-rbind(cbind(rep(1, class_size1), rep(0, class_size2), rep(0, class_size3)),
         cbind(rep(0, class_size1), rep(1, class_size2), rep(0, class_size3)),
         cbind(rep(0, class_size1), rep(0, class_size2), rep(1, class_size3)))

X<-rbind(points1, points2, points3)

z<- cbind(X,y)

print(z)

#Generated sample data is plotted.
plot(points1[,1], points1[,2], type = "p", pch = 19, col = "red", las = 1,
     xlim = c(-6, 6), ylim = c(-6, 6),
     xlab = "x1", ylab = "x2")
points(points2[,1], points2[,2], type = "p", pch = 19, col = "green")
points(points3[,1], points3[,2], type = "p", pch = 19, col = "blue")


#Eta is the step size, epsilon is the stopping condition.
eta <- 0.01
epsilon <- 1e-3

#Initializing the weights for each class 
w1 <- runif(ncol(X), min = -0.01, max = 0.01)
w10 <- runif(1, min = -0.01, max = 0.01)

w2 <- runif(ncol(X), min = -0.01, max = 0.01)
w20 <- runif(1, min = -0.01, max = 0.01)

w3 <- runif(ncol(X), min = -0.01, max = 0.01)
w30 <- runif(1, min = -0.01, max = 0.01)

w<-cbind(w1, w2, w3)
w0<-cbind(t(w10), t(w20), t(w30))


#gradient_w function takes the exact y, estimated y and all the sample data
#and it calculates where to move in order to decrease the objective function.
gradient_w<-function(y, y_pred, X){

  a<-(-colSums(matrix(y[,1] - y_pred[,1], nrow = nrow(X), ncol = ncol(X), byrow = FALSE) * X))
  b<-(-colSums(matrix(y[,2] - y_pred[,2], nrow = nrow(X), ncol = ncol(X), byrow = FALSE) * X))
  c<-(-colSums(matrix(y[,3] - y_pred[,3], nrow = nrow(X), ncol = ncol(X), byrow = FALSE) * X))
  return (cbind(a,b,c))
}

gradient_w0 <- function(y, y_pred) {
  a<-(-sum(y[,1] - y_pred[,1]))
  b<-(-sum(y[,2] - y_pred[,2]))
  c<-(-sum(y[,3] - y_pred[,3]))
  return (cbind(t(a), t(b), t(c)))
}



#softmax function takes the weight matrix and sample data
# to estimate y predicted.
softmax <- function(w, w0, x) {
  i1<-exp(w[,1]%*%x+w0[,1]) / (1 + sum( exp(w[,1]%*%x+w0[,1]), exp(w[,2]%*%x+w0[,2]), exp(w[,3]%*%x+w0[,3])))
  i2<-(exp(w[,2]%*%x+w0[,2]) / (1 + sum( exp(w[,1]%*%x+w0[,1]), exp(w[,2]%*%x+w0[,2]), exp(w[,3]%*%x+w0[,3]) )))
  i3<-(exp(w[,3]%*%x+w0[,3]) / (1 + sum( exp(w[,1]%*%x+w0[,1]), exp(w[,2]%*%x+w0[,2]), exp(w[,3]%*%x+w0[,3]))))

  return (cbind(i1,i2,i3))
}



#initialize y_pred matrix to store the predicted y values.
y_pred<-matrix(0,300,3)

#objective_values will store the each objective value estimated in each iteration. 
objective_values <- c()
j=0

#keep initializing until keepGoing is less than epsilon.
keepGoing<-c(1,1,1)

while (keepGoing>epsilon) {

  #using softmax predict all the y values by using the weights estimated before.
  t<-1
  while(t<301){
    y_pred[t,] <- softmax(w, w0, X[t,])
    t<-t+1
  }
  w_pre <- w
  w0_pre <- w0

  #store objective values in terms of sum(y*log(y_pred)+(1-y)*log(1-y_pred)) y and y_pred are either 1 or 0.
  objective_values <- c(objective_values, -sum(y * log(y_pred + 1e-100) + (1 - y) * log(1 - y_pred + 1e-100)))

  #estimate new w and w0 using the gradient functions defined.
  w <- w - eta * gradient_w(y, y_pred, X)
  w0 <- w0 - eta * gradient_w0(y, y_pred)

  #check the continuation condition keepgoing. The condition for stopping is convergence.
  keepGoing<-sqrt((w0 - w0_pre)^2 + sum((w - w_pre)^2))

  j<-j+1
}


#I round y_pred to 1 or 0 to make data more clean. 
t<-1
while(t<301){
  k<-1
  while(k<4){
  y_pred[t,k] <-round(y_pred[t,k])
  k<-k+1
  }
  t<-t+1
}


#Plot objective function
plot(1:j, objective_values,
     type = "l", lwd = 2, las = 1,
     xlab = "Iteration", ylab = "Error")



confusion_matrix<-matrix(0,3,3)


y_pred_subs<-matrix(1,300,1)
y_pred_subs[y_pred[,2]==1]<-2
y_pred_subs[y_pred[,3]==1]<-3

y_subs<-matrix(1,300,1)
y_subs[y[,2]==1]<-2
y_subs[y[,3]==1]<-3

confusion_matrix<-table(y_subs,y_pred_subs)

print(confusion_matrix)

softmax_plot<- function(x1, x2){
  x<-matrix(c(x1,x2),2,1)
  l<-matrix(0,1,3)
  l<-softmax(w,w0,x)
  j<-1
  while(j<4){
    l[j] <-round(l[j])
    if(l[j]==1){
      return (j)
    }
    j<-j+1
  }
  return (0)
}


#In order to plot the data we first define the interval int for x1 and x2 axis
int<-seq(from=-6, to=6, by=0.06)

#To be able to find our way in plot, for eac point we mark corresponding
#x1 grid points and x2 grid points in the rectangular area.
x1_grid <- matrix(int, nrow = length(int), ncol = length(int), byrow = FALSE)
x2_grid <- matrix(int, nrow = length(int), ncol = length(int), byrow = TRUE)



#Using the function f defined we estimate and store the class_values for each data point in 6 by 6 plot area.
class_values <- matrix(mapply(softmax_plot, x1_grid, x2_grid), nrow(x2_grid), ncol(x2_grid))

plot(points1[,1], points1[,2], type = "p", pch = 19, col = "red", las = 1,
     xlim = c(-6, 6), ylim = c(-6, 6),
     xlab = "x1", ylab = "x2")
points(points2[,1], points2[,2], type = "p", pch = 19, col = "green")
points(points3[,1], points3[,2], type = "p", pch = 19, col = "blue")

#Draws  black circle around the data points in which y_predicted is not same with actual y values.
points(X[y_pred_subs != y_subs, 1], X[y_pred_subs != y_subs, 2], cex = 1.5, lwd = 2)

#For each grid point take the points that have class_value 1,2 and 3 and paint them in corresponding color, red, green and blue.
points(x1_grid[class_values==1], x2_grid[class_values==1 ], col = rgb(red = 1, green = 0, blue = 0, alpha = 0.01), pch = 16)
points(x1_grid[class_values==2], x2_grid[class_values==2 ], col = rgb(red = 0, green = 1, blue =0 , alpha = 0.01), pch = 16)
points(x1_grid[class_values==3], x2_grid[class_values==3 ], col = rgb(red = 0, green = 0, blue = 1, alpha = 0.01), pch = 16)

contour(int, int, class_values, levels = c(0), add = TRUE, lwd = 2, drawlabels = FALSE)


