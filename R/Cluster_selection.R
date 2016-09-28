# Cluster selection: With BIC, AIC, modified BIC or with variance

#' Bayesian Information Criterion
#'
#' Computes BIC from a list of outputs of EM algorithm, then returns the position with minimal BIC
#' @param EM_out_list list of outputs from EM.algo or FullEM
#' @param model.selection The function to minimize for the model selection: can be "AIC", 
#' "BIC", or numeric. In numeric, the BIC function is modified. If variance: returns max(abs(1 - Var(cluster)/expected(Var)))
#' @keywords EM clustering number
BIC_criterion<-function(EM_out_list,model.selection){
  ### Criterion should be minimized
  # Here we assimilate EM.output$val to -ln(L) where L is the likelihood of the model
  # BIC is written -2*ln(L)+k*ln(k)
  # Generalized BIC is written -2*ln(L)+q * k*ln(k)
  # AIC is written 2*k - 2*ln(L)
  
  if(is.numeric(model.selection)){
    ### Modified BIC to relax or add constraints on model selection
    # if q > 1 adding explicative variables should explain observed values better => control overfitting
    # if q = 1 BIC
    # if q < 1 adding explicative variables is less costly
    
    
    Bic<-numeric()
    if(length(EM_out_list)==0){
      return(0)
    }
    Mut_num<-nrow(EM_out_list[[1]]$EM.output$fik)
    for(i in 1:length(EM_out_list)){
      k<-length(EM_out_list[[i]]$EM.output$centers[[1]])
      Bic[i]<-2*EM_out_list[[i]]$EM.output$val+model.selection * k *log(Mut_num)
    }
    W<-which.min(Bic)
    L<-0
    return(Bic)
  }
  else if(model.selection == "BIC"){
    Bic<-numeric()
    if(length(EM_out_list)==0){
      return(0)
    }
    Mut_num<-nrow(EM_out_list[[1]]$EM.output$fik)
    
    for(i in 1:length(EM_out_list)){
      k<-length(EM_out_list[[i]]$EM.output$centers[[1]])
      Bic[i]<-2*EM_out_list[[i]]$EM.output$val+k*log(Mut_num)
    }
    return(Bic)
  }
  else if(model.selection == "AIC"){
    Aic<-numeric()
    if(length(EM_out_list)==0){
      return(0)
    }
    for(i in 1:length(EM_out_list)){
      Aic[i]<-2*EM_out_list[[i]]$EM.output$val+2*length(EM_out_list[[i]]$EM.output$centers[[1]])
    }
    return(Aic)
    
  }

}

#' Hard clustering based on EM output
#'
#' Attributes a mutation to its most likely clone based on the output of the EM algorithm
#' @param EM_out Output from EM.algo or FullEM
#' @keywords EM Hard clustering
hard.clustering<-function(EM_out){
  clust<-apply(X = EM_out$fik,MARGIN = 1,FUN = function(z) {
    if(sum(z==max(z))>1){ ### Look for the multiple clones, and attribute with probability proportional to the weight
      if(max(z)>0){
        pos<-which(z==max(z))
        prob<-EM_out$weights[pos]/(sum(EM_out$weights[pos]))
        #sample(x = pos, size = 1, prob = prob))
        return(pos[which.max(prob)])
      }
      else{ ### all possibilities have 0 probability, so choose one randomly
        return(sample(1:length(z),size = 1))
      }
    }
    else{ # only one clone has maximal probability
      return(which.max(z))
    }
  })
  return(clust)
}