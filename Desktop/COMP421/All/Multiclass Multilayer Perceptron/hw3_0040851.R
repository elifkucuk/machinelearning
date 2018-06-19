
library(MASS)

set.seed(221)


# Parameters for the random data points to be generated.
class_size<-50



class_means<-matrix(c(2.0,2.0,-4.0,-4.0,-2.0,2.0,4.0,-4.0,
                      -2.0,-2.0,4.0,4.0,2.0,-2.0,-4.0,4.0),2,8)


class_covs<-array(c(0.8,-0.6,-0.6,0.8,
                     0.4, 0.0,0.0,0.4,
                     0.8,0.6,0.6,0.8,
                     0.4,0.0,0.0,0.4,
                     0.8,-0.6,-0.6,0.8,
                     0.4,0.0,0.0,0.4,
                     0.8,0.6,0.6,0.8,
                     0.4,0.0,0.0,0.4),dim=c(2,2,8), dimnames = NULL)


#Generation of the random data points from three bivariate Gaussian densities.
pointList<-list()

for(i in seq(1:8)){
  pointList[[length(pointList)+1]] <- mvrnorm(n = class_size, mu = class_means[,i], Sigma = class_covs[,,i])

}

X<-rbind(pointList[[1]], pointList[[2]], pointList[[3]], pointList[[4]],
         pointList[[5]], pointList[[6]], pointList[[7]], pointList[[8]])

#y_binary is the N*Class_type_Size matrix. Each column belongs to one class type with 1 and 0 binary value
#indicating the class of the data point.
y_binary<-rbind(cbind(rep(1, (class_size*2)), rep(0, (class_size*2)), rep(0, (class_size*2)), rep(0, (class_size*2))),
         cbind(rep(0, (class_size*2)), rep(1, (class_size*2)), rep(0, (class_size*2)), rep(0, (class_size*2))),
         cbind(rep(0, (class_size*2)), rep(0, (class_size*2)), rep(1, (class_size*2)), rep(0, (class_size*2))),
         cbind(rep(0, (class_size*2)), rep(0, (class_size*2)), rep(0, (class_size*2)), rep(1, (class_size*2)))
          )
#non_binary representation of classes for each data point
y<-c()
for(i in 1:4){
  y<-rbind(y, matrix(rep(i, (class_size*2)), class_size*2, 1))
}


Data<-cbind(X,y)


#Parameters defined for multiclass multilayer perceptron algorithm
eta<-0.1
epsilon<-1e-3
H<-20
max_iteration<-200

#Parameters defined for the problem
K<-length(y_binary[1,])
N<-length(y[,1])
D<-length(X[1,])


#Sigmoid function is defined. It will be used to create Z values.
sigmoid<-function(a){
  return (1/(1+exp(a)))
}

#Softmax function is defined. It will be used to transform values between 1 and 0
#for Y_predicted values.
softmax<-function(o){
  exp_o<-sapply(o, function(x) exp(x))
  ret<-sapply(exp_o, function(x) (x/(sum(exp_o))) )
  return(ret)
}

#softmax_All  function applies softmax function for all rows
#in one matrix.
softmax_ALL<-function(o){
  ret<-apply(o, 1, softmax)
  return (t(ret))
}

#gradient_v estimates the change in v for given 1 row of y_truth,
#y_predicted and z
gradient_v<- function(y_truth, y_predicted, z ){
  grad<-eta *  c(1, z) %*% t((y_truth - y_predicted))
  return (grad)
}

#gradient_w estimates the change in w for given 1 row of y_truth,
#y_predicted, z and x.
gradient_w<-function(y_truth, y_predicted, z, x){
  grad<-matrix(0, D+1, H)
  for (h in 1:H) {
    grad[,h]<- eta * as.vector((t(v[h,])%*% (y_truth - y_predicted)))*as.vector((z[h] %*% t((1 - z[h])))) * c(1, x)
  }
  return (grad)
}

#objective_val estimates the error of given y_truth table and y_pred table.
objective_val<-function(y_truth, y_pred){
  ret<- y_pred * log(y_truth + 1e-100)
  return (-sum(colSums(ret)))
}


#I round y_pred to 1 or 0 to make data more clean.
round_y<-function(y_pred){
  for(t in 1:N){
    k<-1
    for(k in 1:K){
      y_pred[t,k] <-round(y_pred[t,k])
    }
  }
  return (y_pred)
}




#Starting point for W, v, Z

#W is  matrix of dimension (D+1)*H
W <- matrix(runif((D + 1) * H, min = -0.01, max = 0.01), D + 1, H)
#v is  matrix of dimension (H+1)*K
v <- matrix(runif(K*(H + 1), min = -0.01, max = 0.01), H+1, K)

#Z is  matrix of dimension (N)*H
Z <- sigmoid(cbind(1, X) %*% W)

#y_pred_binary gives the predicted value for classes associated with each input.
#y_pred_binary is  matrix of dimension (N)*K
y_pred_binary <- softmax_ALL((cbind(1, Z) %*% v))

#We store each objective value in an array called objective_values.
objective_values <- c(objective_val(y_binary, y_pred_binary))
iteration<-1


# learn W and v using gradient descent and online learning
while (1) {
  for (u in sample(N)) {
    # calculate hidden nodes
    Z[u,] <- sigmoid(c(1, X[u,]) %*% W)
    # calculate output node
    y_pred_binary[u,] <- softmax(c(1, Z[u,]) %*% v)
    #Update the weights associated with hidden and output nodes.
    v<- v + gradient_v(y_binary[u,], y_pred_binary[u,], Z[u,])
    W<-W+gradient_w(y_binary[u,], y_pred_binary[u,], Z[u,], X[u,])
  }

  Z <- sigmoid(cbind(1, X) %*% W)
  y_pred_binary <- softmax_ALL(cbind(1, Z) %*% v)
  objective_values <- c(objective_values, objective_val(y_binary, y_pred_binary))

  if (abs(objective_values[iteration + 1] - objective_values[iteration]) < epsilon | iteration >= max_iteration) {
    break
  }
  iteration <- iteration + 1
}

#Plot the objective function by each iteration.
plot(1:(iteration + 1), objective_values,
     type = "l", lwd = 2, las = 1,
     xlab = "Iteration", ylab = "Error")


#round to 1 or 0 for each y_pred_binary element.
y_pred_binary<-round_y(y_pred_binary)

#y_predicted is a column vector with N rows. It gives the class type values for each input in 1,2,3,4 fashion.
y_predicted<-c(1:400)
y_predicted[which(y_pred_binary[,1]==1)]<-1
y_predicted[which(y_pred_binary[,2]==1)]<-2
y_predicted[which(y_pred_binary[,3]==1)]<-3
y_predicted[which(y_pred_binary[,4]==1)]<-4


#Confusion table is created using table function and column vectors y_predicted and y.
confusion_table<-table(y_predicted,y)
print(confusion_table)

#Predict_plot function takes one input data anc calculates the predicted y
predict_plot<- function(x1, x2){
  z_plot<-sigmoid(c(1,x1,x2)%*% W)
  y_pred_plot<-softmax(cbind(1,z_plot) %*%v)
  z<-max.col(t(y_pred_plot))
  return (z)
}

# Define intervals for x1 and x2 for plot
x1_interval <- seq(from = -6, to = +6, by = 0.06)
x2_interval <- seq(from = -6, to = +6, by = 0.06)
x1_grid <- matrix(x1_interval, nrow = length(x1_interval), ncol = length(x1_interval), byrow = FALSE)
x2_grid <- matrix(x2_interval, nrow = length(x2_interval), ncol = length(x2_interval), byrow = TRUE)


#Estimated class_values for associated plot points are calculated.
class_value <- matrix(mapply(predict_plot, x1_grid, x2_grid), nrow(x2_grid), ncol(x2_grid))



#Plot the sample points
plot(X[y == 1, 1], X[y == 1, 2], type = "p", pch = 19, col = "red",
     xlim = c(-6, +6),
     ylim = c(-6, +6),
     xlab = "x1", ylab = "x2", las = 1)
points(X[y == 2, 1], X[y == 2, 2], type = "p", pch = 19, col = "green")
points(X[y == 3, 1], X[y == 3, 2], type = "p", pch = 19, col = "blue")
points(X[y == 4, 1], X[y == 4, 2], type = "p", pch = 19, col = "magenta")

#Circle the sample points that are predicted wrong.
points(X[y_predicted != y, 1], X[y_predicted != y, 2], cex = 1.5, lwd = 2)

#Draw the decision boundaries associated with each class.
points(x1_grid[class_value==1 ], x2_grid[class_value==1], col = rgb(red = 1, green = 0, blue = 0, alpha = 0.01), pch = 16)
points(x1_grid[class_value==2], x2_grid[class_value==2], col = rgb(red = 0, green = 1, blue =0 , alpha = 0.01), pch = 16)
points(x1_grid[class_value==3 ], x2_grid[class_value==3], col = rgb(red = 0, green = 0, blue = 1, alpha = 0.01), pch = 16)
points(x1_grid[class_value==4], x2_grid[class_value ==4], col = rgb(red =1, green =0, blue =1 , alpha = 0.01), pch = 16)

#Countour the decisicon boundaries.
contour(x1_interval, x2_interval, class_value, levels = c(1,2,3,4), add = TRUE, lwd = 2, drawlabels = FALSE)


