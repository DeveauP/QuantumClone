### Store functions that are used for gradient descent:
### peak, grx, eval.fik, eval.fik.m, list.prod

#'List product
#'
#' Returns the product of all elements in a list, e.g. a vector if the elements of the list are vectors, etc.
#' @param L list used
#' @param col If it is a list of matrices, and only one column should be used, name of the column.
#' @keywords List handling
#' @examples 
#' list_prod(list(matrix(1:4,nrow = 2),matrix(1:4,nrow = 2)))
#' @export
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
  test<-is.infinite(spare)
  if(sum(test)){
    spare[test]<--log(.Machine$double.xmax)
  }
  spare
}

#' Computes gradient of function
#'
#' @param fik Evaluation of fik for previous iteration
#' @param adj.factor Factor to compute the probability: makes transition between the cellularity of the clone and the frequency observed
#' @param centers vector with cellularity of each clone (numeric vector, ordered by samples)
#' @param Alt Matrix with number of draws in rows for a mutation/possibility, and samples in columns
#' @param Depth Matrix with number of not draws (Depth - Alt) in rows for a mutation/possibility, and samples in columns
#' @examples 
#' fik<-matrix(c(1,0,0,1),nrow = 2)
#' adj.factor<-matrix(1/2,nrow =2 ,ncol =1)
#' centers<-c(0.25,0.75)
#' Alt<-c(125,375)
#' Depth<-c(1000,1000)
#' grbase(fik,adj.factor,centers,Alt,Depth)
#' @export
grbase<-compiler::cmpfun(function(fik,adj.factor,centers,Alt,Depth){
  result<-numeric(length = length(centers))
  ## fik has the Schrod possibilities in rows and clones in cols
  ## adj.factors has mutations in row and samples in cols
  if(is.matrix(adj.factor) && ncol(adj.factor)>1){
    ###
    # MULTISAMPLE CASE
    ###
    index<-0
    centers.per.sample<-length(centers)/ncol(adj.factor)
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
        
        
        spare<-fik[,i]*{
          {Alt[,s]/centers[index]} - {adj.factor[,s]*( Depth[,s] - Alt[,s])}/
          {1-adj.factor[,s]*centers[index]}
        }
        spare[is.infinite(spare)]<-sign(is.infinite(spare))*log(.Machine$double.xmax)
        result[index]<--sum(spare)
      }
    }
  }else{
    ###
    # SINGLE SAMPLE
    ###
    centers.per.sample<-length(centers)
    
    for(i in 1:centers.per.sample){
      
      pro<-adj.factor*centers[i]
      spare<-fik[,i]*{
        {Alt/centers[i]} - {adj.factor*( Depth - Alt)}/
        {1-adj.factor*centers[i]}
      }
      spare[is.infinite(spare)]<-sign(is.infinite(spare))*log(.Machine$double.xmax)
      spare<--sum(spare)
      if(is.infinite(spare)){
        result[i]<-sign(spare)*log(.Machine$double.xmax)
      }
      else{
        result[i]<-spare
      }
    }
  }
  result
},options = list(optimize = 3)
)