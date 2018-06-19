# read data into memory
data_set <- read.csv("hw04_data_set.csv")

N<-133
N_train<-100
N_test<-33
h<-3

#choose random points between 1 and 133
indices<-c(1:N)
train_indices<-sample(c(1:N), N_train)
test_indices<-setdiff(indices,train_indices)

# seperate training and test data
data_train <-  data_set[train_indices,]
data_test <-  data_set[test_indices,]


# get number of classes and number of samples
#choose the origin as the floor function of minimum x value and bindwidth as 3
minimum_value<-floor(min(data_train[,1]))
maximum_value<-ceiling(max(data_train[,1]))


#Data interval for plot
data_interval <- seq(from = minimum_value, to = maximum_value, by = 0.1)

#Question 2
#Firstly define the bind_intervals
bind_interval <- cbind(seq(from = minimum_value, to = maximum_value-h, by = h),
                       seq(from=minimum_value+h, to=maximum_value, by=h))

#Define a function that corresponds to g(x) for regressogram
reg_est<-function(data__){
num<-sum(sapply(1:length(bind_interval[,1]), function(x) {sum((data__>=bind_interval[x,1] & data__<bind_interval[x,2])*(data_train[,1]>=bind_interval[x,1] & data_train[,1]<bind_interval[x,2]) * 
                                                            data_train[,2])}))
den<-sum(sapply(1:length(bind_interval[,1]), function(x) {sum((data__>=bind_interval[x,1] & data__<bind_interval[x,2])*(data_train[,1]>=bind_interval[x,1] & data_train[,1]<bind_interval[x,2]))}))
return (num/den)
}

##Draw the regrossogram on the defined data_interval

p_head_plot<-sapply(data_interval, reg_est)

plot(data_train[,1],  data_train[,2], type = "p", pch = 20, col = "RED",
     ylab = "density", xlab = "x")
points(data_test[,1], data_test[,2], type = "p", pch = 20, col = "BLUE" )
lines(data_interval, p_head_plot, type = "l", lwd = 2, col = "black")

#Question 3

#Defined a function to calculate mean squared error. I use this function in Question 3-5-7
find_mse<-function(p_h, s){
  error <- mean((data_test[,2]-p_h)^2)
  RMSE<-sqrt(sum(error))
  output_question5<-sprintf("%s=> RMSE is %f when h is %d",s, RMSE, h )
  output_question5
}

p_head<-sapply(data_test[,1], reg_est)

find_mse(p_head, "Regressogram")

#Question 4

# w function defined. Bin is symmetric around x 
running_mean_weight<-function(u){
  value<-0
  if (abs(u)<1){
    value<-1
  }
  return (value)
}

running_mean_smoother<-function(point, train_data){
  u<-(point-train_data[,1])/h
  up<-sum(sapply(u, running_mean_weight)*train_data[,2])
  down<-sum(sapply(u, running_mean_weight))
  value<-(up/down)
  return (value)
}

p_head<- sapply(data_interval, running_mean_smoother, data_train)

plot(data_train[,1], data_train[,2], type = "p", pch = 20, col = "RED",
     ylab = "density", xlab = "x")
points(data_test[,1], data_test[,2], type = "p", pch = 20, col = "BLUE" )
lines(data_interval, p_head, type = "l", lwd = 2, col = "black")

#Question 5
p_head<-learn_running_mean_smoother_fun(data_test[,1], data_train)

find_mse(p_head, "Running Mean Smoother")

#Question 6

Kernel<-function(x){
value<-(1/sqrt(2*pi)) * exp(-x^2/2)
}


# kernel estimator

kernel_est<-function(data__){
  numerator <- sapply(data__, function(x) {sum(1 / sqrt(2 * pi) * exp(-0.5 * (x - data_train[,1])^2 ) * data_train[,2])}) 
  denominator <- sapply(data__, function(x) {sum(1 / sqrt(2 * pi) * exp(-0.5 * (x - data_train[,1])^2 ))}) 
  return (numerator/denominator)
}

p_head<-kernel_est(data_interval)
plot(data_train[,1],  data_train[,2], type = "p", pch = 20, col = "RED",
     ylab = "density", xlab = "x")
points(data_test[,1], data_test[,2], type = "p", pch = 20, col = "BLUE" )
lines(data_interval, p_head, type = "l", lwd = 2, col = "black")

#Question 7
find_mse(kernel_est(data_test[,1]), "Kernel Estimator")




