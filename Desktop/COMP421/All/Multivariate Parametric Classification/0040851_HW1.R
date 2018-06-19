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

#Collecting data points from three bivariate Gaussian densities in one vector.
X<-rbind(points1,points2, points3)

#Y values (Class type) for the corresponding data points.
y<-c(rep(1,class_size1), rep(2, class_size2), rep(3, class_size3))

z<-cbind(X,y)

write.csv(x=cbind(X,y),file = "hw1_data_set.csv")



sample_mean1<-c(mean(points1[,1]),mean(points1[,2]))
sample_mean2<-c(mean(points2[,1]), mean(points2[,2]))
sample_mean3<-c(mean(points3[,1]), mean(points3[,2]))

#Sample mean matrix
sample_means<-matrix(c(sample_mean1, sample_mean2, sample_mean3),2,3)

 sample_cov1<-matrix(c(cov(points1[,1], points1[,1]), cov(points1[,2],points1[,1]),
                        cov(points1[,1], points1[,2]), cov(points1[,2],points1[,2])),2,2)

 sample_cov2<-matrix(c(cov(points2[,1], points2[,1]), cov(points2[,2],points2[,1]),
                      cov(points2[,1], points2[,2]), cov(points2[,2],points2[,2])),2,2)

 sample_cov3<-matrix(c(cov(points3[,1], points3[,1]), cov(points3[,2],points3[,1]),
                       cov(points3[,1], points3[,2]), cov(points3[,2],points3[,2])),2,2)

 #Sample covariance matrix
sample_covariance<-array(c(cov(points1[,1], points1[,1]), cov(points1[,2],points1[,1]),
                    cov(points1[,1], points1[,2]), cov(points1[,2],points1[,2]),
                    cov(points2[,1], points2[,1]), cov(points2[,2],points2[,1]),
                    cov(points2[,1], points2[,2]), cov(points2[,2],points2[,2]),
                    cov(points3[,1], points3[,1]), cov(points3[,2],points3[,1]),
                    cov(points3[,1], points3[,2]), cov(points3[,2],points3[,2])),
                  dim=c(2,2,3), dimnames = NULL)

#Class priors
class_priors<-c((class_size1/(class_size1+class_size2+class_size3)),
                (class_size2/(class_size1+class_size2+class_size3)),
                (class_size3/(class_size1+class_size2+class_size3)))

print(sample_means)

print(sample_covariance)

print(class_priors)


W1<-(-1/2)*chol2inv(chol(sample_cov1))
w1<-chol2inv(chol(sample_cov1))%*%sample_mean1
w10<-(-1/2)*t(sample_mean1)%*%chol2inv(chol(sample_cov1))%*%sample_mean1-(1/2)*log(det(sample_cov1))+log(class_priors[1])

#Score function for class 1.
g1<-function(x){
(t(x)%*%W1)%*%x + t(w1)%*%x + w10
}

W2<-(-1/2)*chol2inv(chol(sample_cov2))
w2<-chol2inv(chol(sample_cov2))%*%sample_mean2
w20<-(-1/2)*t(sample_mean2)%*%chol2inv(chol(sample_cov2))%*%sample_mean2-(1/2)*log(det(sample_cov2))+log(class_priors[2])

#Score function for class 2
g2<-function(x){
  (t(x)%*%W2)%*%x + t(w2)%*%x + w20
}


W3<-(-1/2)*chol2inv(chol(sample_cov3))
w3<-chol2inv(chol(sample_cov3))%*%sample_mean3
w30<-(-1/2)*t(sample_mean3)%*%chol2inv(chol(sample_cov3))%*%sample_mean3-(1/2)*log(det(sample_cov3))+log(class_priors[3])

#Score function for class 3
g3<-function(x){
  (t(x)%*%W3)%*%x + t(w3)%*%x + w30
}



i<-1

#Estimating the y_predicted using the score functions for class 1,2 and 3. Data point belongs to the class
#with highest score function.
y_predicted<-matrix(cbind(X,0),300,3)
trans_X<-t(X)
while(i<301){

  bir<-g1(trans_X[,i])
  iki<-g2(trans_X[,i])
  uc<-g3(trans_X[,i])
  if(bir>iki && bir>uc){
    y_predicted[i,3]<-1
  }
  else if(iki>bir && iki>uc){
    y_predicted[i,3]<-2
  }
  else{
    y_predicted[i,3]<-3
  }
  i<-i+1
}


y_predicted_list<-y_predicted[,3]

j<-1

#Confusion matrix
confusion_matrix<-matrix(0,3,3)
while(j<301){
  confusion_matrix[y_predicted[j,3], z[j,3]]<-confusion_matrix[y_predicted[j,3], z[j,3]]+1
  j<-j+1
}


#In order to plot the data we first define the interval int for x1 and x2 axis
int<-seq(from=-6, to=6, by=0.06)

#To be able to find our way in plot, for eac point we mark corresponding
#x1 grid points and x2 grid points in the rectangular area.
x1_grid <- matrix(int, nrow = length(int), ncol = length(int), byrow = FALSE)
x2_grid <- matrix(int, nrow = length(int), ncol = length(int), byrow = TRUE)

#For each point in the rectangular plot area using the score functions defined for
#class 1,2,3 we define a function that estimates the class that point belongs to.

class_values <- matrix(0, c(length(int), length(int)))
f <- function(x1, x2) {
  x<-matrix(c(x1,x2),2,1)
  bir<-g1(x)
  iki<-g2(x)
  uc<-g3(x)
  if(bir>iki && bir>uc){
    1
  }
  else if(iki>bir && iki>uc){
    2
  }
  else{
    3
  }
}

#Using the function f defined we estimate and store the class_values for each data point in 6 by 6 plot area.
class_values <- matrix(mapply(f, x1_grid, x2_grid), nrow(x2_grid), ncol(x2_grid))

#Plot the data points belonging to class 1 (red) , 2 (green) and 3 (blue)
plot(points1[,1], points1[,2], type = "p", pch = 19, col = "red", las = 1,
     xlim = c(-6, 6), ylim = c(-6, 6),
     xlab = "x1", ylab = "x2")
points(points2[,1], points2[,2], type = "p", pch = 19, col = "green")
points(points3[,1], points3[,2], type = "p", pch = 19, col = "blue")

#Draws  black circle around the data points in which y_predicted is not same with actual y values.
points(X[y_predicted_list != y, 1], X[y_predicted_list != y, 2], cex = 1.5, lwd = 2)

#For each grid point take the points that have class_value 1,2 and 3 and paint them in corresponding color, red, green and blue.
points(x1_grid[class_values==1], x2_grid[class_values==1 ], col = rgb(red = 1, green = 0, blue = 0, alpha = 0.01), pch = 16)
points(x1_grid[class_values==2], x2_grid[class_values==2 ], col = rgb(red = 0, green = 1, blue =0 , alpha = 0.01), pch = 16)
points(x1_grid[class_values==3], x2_grid[class_values==3 ], col = rgb(red = 0, green = 0, blue = 1, alpha = 0.01), pch = 16)

contour(int, int, class_values, levels = c(0), add = TRUE, lwd = 2, drawlabels = FALSE)

