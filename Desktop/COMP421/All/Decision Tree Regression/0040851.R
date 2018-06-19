#read data given for hw 5
data_set<-read.csv("hw05_data_set.csv")

set.seed(521)
N_training<-100
N_test<-33

#Question 1
# separate the train indices from test
indice_train<-sample(seq(1:length(data_set[,1])), N_training)
indice_test<-setdiff(seq(1:length(data_set[,1])), indice_train)
data_train<-data_set[indice_train,]
data_test<-data_set[indice_test,]

#Question 2

#g function evaluates the estimated value in a node
g<-function(indices_data){
  return (sum(data_train[indices_data,2])/length(indices_data))
}

#E_m function evaluates mean squared error. This defines the goodness of split.
E_m<-function(indices_data){
  return ((1/length(indices_data))*sum((data_train[indices_data,2]-g(indices_data))^2))
}

#Data structure defined for the decision tree.
need_splits<-c()
node_indices<-list()
is_terminal<-c()
node_splits<-c()
node_terms<-c()
split_scores<-c()
best_split<-c()
best_score<-c()


#Initialize the data structure. 
need_splits<-c(TRUE)
is_terminal<-c(FALSE)
node_indices<-list(1:length(data_train[,1]))
node_terms<-c()
P<-10

#an algorithm to learn a decision tree. 
while(1){
  split_nodes<-which(need_splits)
  if(length(split_nodes)==0){
    break
  }
  else{
    for(node_split in split_nodes){
      data_indices<-node_indices[[node_split]]
      #node_terms[node_split]<-mean(data_train[data_indices,2])
      node_terms[node_split]<-g(data_indices)
      need_splits[node_split]<-FALSE
      if(length(unique(data_train[data_indices,2]))<=P){
        is_terminal[node_split]<-TRUE
      }
      else if(length(unique(data_train[data_indices,1]))==1){
        is_terminal[node_split]<-TRUE
      }
      else{
        is_terminal[node_split]<-FALSE
        unique_values<-sort(unique(data_train[data_indices,1]))
        split_points<-(unique_values[-1]+unique_values[-(length(unique_values))])/2
        split_scores<-c()
        best_split<-c()
        best_score<-c()
        for(s in seq(1:length(split_points))){
          left_indices<-data_indices[which(data_train[data_indices,1]<split_points[s])]
          right_indices<-data_indices[which(data_train[data_indices,1]>split_points[s])]
          split_scores[s]<-sum((length(left_indices)/length(data_indices))*E_m(left_indices),
                               (length(right_indices)/length(data_indices))*E_m(right_indices))
        }
        s<-sprintf("split_scores: %s", paste(data_indices, collapse=","))
        #print(s)
        best_score<-min(split_scores)
        best_split<-split_points[which.min(split_scores)]
        
        node_splits[node_split]<-best_split
        
        #create left indices
        left_indices<-data_indices[which(data_train[data_indices,1]<best_split)]
        node_indices[[2*node_split]]<-left_indices
        need_splits[2*node_split]<-TRUE
        is_terminal[2*node_split]<-FALSE
        
        right_indices<-data_indices[which(data_train[data_indices,1]>best_split)]
        node_indices[[(2*node_split)+1]]<-right_indices
        need_splits[(2*node_split)+1]<-TRUE
        is_terminal[(2*node_split)+1]<-FALSE
      }
    }
  }
}

#In order to print the rules defined by the algorithm.
terminal_nodes <- which(is_terminal)
for (terminal_node in terminal_nodes) {
  index <- terminal_node
  rules <- c()
  while (index > 1) {
    parent <- floor(index / 2)
    if (index %% 2 == 0) {
      # if node is left child of its parent
      rules <- c(sprintf("x < %g", node_splits[parent]), rules)
    } else {
      # if node is right child of its parent
      rules <- c(sprintf("x >= %g", node_splits[parent]), rules)
    }
    index <- parent
  }
  rules<-sprintf("{%s} => [%s]", paste0(rules, collapse = " AND "), paste0(node_terms[terminal_node], collapse = "-"))
}





#Question 3

#Decision tree function given the P value learns a decision tree and evaluates y_predicted for the data__
#This function is same as the algorithm above.  
decision_tree<-function(P, data__){
  #Data Structure Needed
  
  #Initialize data structure with root
  need_splits<-c(TRUE)
  is_terminal<-c(FALSE)
  node_indices<-list(1:length(data_train[,1]))
  node_terms<-c()
  
  while(1){
    split_nodes<-which(need_splits)
    if(length(split_nodes)==0){
      break
    }
    else{
      for(node_split in split_nodes){
        data_indices<-node_indices[[node_split]]
        #node_terms[node_split]<-mean(data_train[data_indices,2])
        node_terms[node_split]<-g(data_indices)
        need_splits[node_split]<-FALSE
        if(length(unique(data_train[data_indices,2]))<=P){
          is_terminal[node_split]<-TRUE
        }
        else if(length(unique(data_train[data_indices,1]))==1){
          is_terminal[node_split]<-TRUE
        }
        else{
          is_terminal[node_split]<-FALSE
          unique_values<-sort(unique(data_train[data_indices,1]))
          split_points<-(unique_values[-1]+unique_values[-(length(unique_values))])/2
          split_scores<-c()
          best_split<-c()
          best_score<-c()
          for(s in seq(1:length(split_points))){
            left_indices<-data_indices[which(data_train[data_indices,1]<split_points[s])]
            right_indices<-data_indices[which(data_train[data_indices,1]>split_points[s])]
            split_scores[s]<-sum((length(left_indices)/length(data_indices))*E_m(left_indices),
                                 (length(right_indices)/length(data_indices))*E_m(right_indices))
          }
          s<-sprintf("split_scores: %s", paste(data_indices, collapse=","))
          #print(s)
          best_score<-min(split_scores)
          best_split<-split_points[which.min(split_scores)]
          
          node_splits[node_split]<-best_split
          
          #create left indices
          left_indices<-data_indices[which(data_train[data_indices,1]<best_split)]
          node_indices[[2*node_split]]<-left_indices
          need_splits[2*node_split]<-TRUE
          is_terminal[2*node_split]<-FALSE
          
          right_indices<-data_indices[which(data_train[data_indices,1]>best_split)]
          node_indices[[(2*node_split)+1]]<-right_indices
          need_splits[(2*node_split)+1]<-TRUE
          is_terminal[(2*node_split)+1]<-FALSE
        }
      }
    }
  }
  y_predicted <- rep(0, length(data__))
  for (i in seq(1:length(data__))) {
    index <- 1
    while (1) {
      if (is_terminal[index] == TRUE) {
        y_predicted[i] <- node_terms[index]
        break
      } else {
        if (data__[i] < node_splits[index]) {
          index <- index * 2
        } else {
          index <- index * 2 + 1
        }
      }
    }
  }
  terminal_nodes <- which(is_terminal)
  for (terminal_node in terminal_nodes) {
    index <- terminal_node
    rules <- c()
    while (index > 1) {
      parent <- floor(index / 2)
      if (index %% 2 == 0) {
        # if node is left child of its parent
        rules <- c(sprintf("x < %g", node_splits[parent]), rules)
      } else {
        # if node is right child of its parent
        rules <- c(sprintf("x >= %g", node_splits[parent]), rules)
      }
      index <- parent
    }
    rules<-sprintf("{%s} => [%s]", paste0(rules, collapse = " AND "), paste0(node_terms[terminal_node], collapse = "-"))
  }
  
  return (y_predicted)
}


#Question 3
data_interval <- seq(from = min(data_train[,1]), to = max(data_train[,1]), by = 0.1)
p_head<-decision_tree(10,data_interval)
plot(data_train[,1], data_train[,2], type = "p", pch = 20, col = "BLUE",
     ylab = "y", xlab = "x")
points(data_test[,1], data_test[,2], type = "p", pch = 20, col = "RED" )
lines(data_interval, p_head, type = "l", lwd = 2, col = "black")


#Question 4

#RMSE function given the y_truth and p_head evaluates the mean squared error
P<-10
RMSE<-function(y_truth, p_h){
  val<-(sqrt(sum((y_truth-p_h)^2)/length(y_truth)))
  s<-sprintf("RMSE is %f when P is %d", val, P)
  print(s)
  return (val)
}

p_head<-decision_tree(P, data_test[,1])
error<-RMSE(data_test[,2],p_head)


#Question 5
#Plotting the RMSE for decision tree wrt P in range 1 and 20.
Ps<-seq(from=1 , to=20, by=1)
error<-c()
for(p in Ps){
  P<-p
  p_head<-decision_tree(P, data_test[,1])
  error<-c(error, RMSE(data_test[,2],p_head))
}
plot(Ps,error , type="l", xlab="P", ylab="RMSE", main="RMSE wrt P")
