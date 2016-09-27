### Store functions that are used for gradient descent:
### peak, grx, eval.fik, eval.fik.m, list.prod
#'List product
#'
#' Returns the product of all elements in a list, e.g. a vector if the elements of the list are vectors, etc.
#' @param L list used
#' @param col If it is a list of matrices, and only one column should be used, name of the column.
#' @keywords List handling
#' #Write example for list_prod
list_prod<-function(L,col=NULL){
  if(is.null(col)){
    if(length(L)>1){
      result<-L[[1]]
      for(i in 2:length(L)){
        result<-result*L[[i]]
      }
    }
    else{
      return(L[[1]])
    }
  }
  else{
    if(length(L)>1){
      result<-L[[1]][,col]
      for(i in 2:length(L)){
        result<-result*L[[i]][,col]
      }
    }
    else{
      return(L[[1]][,col])
    }
  }
  return(result)
}

fik.from.al<-function(al,id=NULL,keep.all.poss = TRUE){
  
  spare<-list_prod(al)
  #tab<-table(u)
  if(keep.all.poss){
    return(spare)
  }
  else{
    fik<-matrix(nrow = length(unique(id)),ncol = ncol(al[[1]]))
    u<-unique(id)
    # fik<- as.data.frame(cbind(spare,id =id)) %>% group_by(id) %>% summarise_all(funs(sum)) # takes longer on matrices of 100 rows
    for(i in 1:length(u)){
      if(sum(id==u[i])>1){ ##more than one possibility for a mutation
        fik[i,]<-apply(X = spare[id==u[i],],MARGIN = 2, sum) ##normalize by sum of possibilities...
      }
      else{ ## only one possibility for a mutation
        fik[i,]<-spare[id==u[i],]
      }
    }
    
  }
  fik[fik==0]<-.Machine$double.xmin ## replace by machine limit to avoid the log(0) issue
  #print(cbind(fik2,fik))
  #return(as.matrix(fik[,-1]))
  return(as.matrix(fik))
}


eval.fik<-function(Schrod,centers,weights,keep.all.poss=TRUE,alpha,adj.factor,epsilon,log = FALSE){
  al<-list()
  if(is.list(centers)){
    centers<-unlist(centers)
  }
  idx<-0
  if(log){
    al<-matrix(data = 0,nrow=nrow(Schrod[[1]]),ncol=length(weights))
  }
  else{
    al<-matrix(data = 1,nrow=nrow(Schrod[[1]]),ncol=length(weights))
  }
  for(i in 1:length(Schrod)){ ## i is a sample
    Alt<-Schrod[[i]]$Alt
    Depth<-Schrod[[i]]$Depth
    adj<-adj.factor[,i]
    for(k in 1:length(weights)){ ## k is a clone
      idx<-idx+1
      pro<-centers[idx]*adj
      test<-pro <=1 & pro >=0
      #pro_0<-pro
      #pro_0[pro>1 | pro<0]<-0
      if(log){
        al[test,k]<-al[test,k]+dbinom(x =Alt[test] ,size = Depth[test],prob = pro[test],log = TRUE)
      }
      else{
        al[test,k]<-al[test,k]*dbinom(x =Alt[test] ,size = Depth[test],prob = pro[test],log = FALSE)
        
      }
      al[!test,k]<-NA
      
    }
  }
  
  al
}

#' Eval probability for M step
#' Computes the log directly as log density is faster to compute
#' 
#' @param Schrod The shcrodinger list of matrices
#' @param centers centers of the clusters
#' @param weights weight of each cluster
#' @param adj.factor The adjusting factor, taking into account contamination, copy number, number of copies
#' @param epsilon Stop value: maximal admitted value of the difference in cluster position and weights
#' @param log Should it compute the log distribution (TRUE) or probability (FALSE)
#' @param alpha Corrective value (if some variants are more likely than others) 
#' between two optimization steps. If NULL, will take 1/(median depth).
eval.fik.m<-function(Schrod,centers,weights,alpha,adj.factor,epsilon,log = TRUE){
  spare<-eval.fik(Schrod = Schrod,
                  centers=centers,
                  weights = weights,
                  alpha = alpha, 
                  adj.factor= adj.factor,
                  epsilon= epsilon,
                  log = log)
  spare
}

#' Computes gradient of function
#'
#' @param fik Evaluation of fik for previous iteration
#' @param adj.factor Factor to compute the probability: makes transition between the cellularity of the clone and the frequency observed
#' @param centers vector with cellularity of each clone (numeric vector, ordered by samples)
#' @param Alt Matrix with number of draws in rows for a mutation/possibility, and samples in columns
#' @param Depth Matrix with number of not draws (Depth - Alt) in rows for a mutation/possibility, and samples in columns
grbase<-compiler::cmpfun(function(fik,adj.factor,centers,Alt,Depth){
  result<-numeric(length = length(centers))
  centers.per.sample<-length(centers)/ncol(adj.factor)
  ## fik has the Schrod possibilities in rows and clones in cols
  ## adj.factors has mutations in row and samples in cols
  index<-0
  for(s in 1:ncol(adj.factor)){
    for(i in 1:centers.per.sample){
      index<-index+1
      pro<-adj.factor[,s]*centers[index]
      # test<-pro >= 0 & pro <= 1
      # 
      # spare<-sum(fik[test,i]*
      #     {Alt[test,s]/centers[index]} - {adj.factor[test,s]*( Depth[test,s] - Alt[test,s])}/
      #     {1-adj.factor[test,s]*centers[index]}
      # )
      
      
      spare<--sum(
        fik[,i]*{
          
          {Alt[,s]/centers[index]} - {adj.factor[,s]*( Depth[,s] - Alt[,s])}/
          {1-adj.factor[,s]*centers[index]}
        })
      if(is.infinite(spare)){
        result[index]<-sign(spare)*.Machine$double.xmax
      }
      else{
        result[index]<-spare
      }
    }
  }
  result
},options = list(optimize = 3)
)