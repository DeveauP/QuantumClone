#' Data generation
#'
#' Creates a phylogenetic tree on simple assumptions:
#' 1. There is a common ancestor to all clones
#' 2. Each clone creates a partition of the space. One node will carry mutations, the other will not. 
#' A node that is not mutated can be partitioned.
#' 3. Nodes that have less than 2% difference in all samples to another existing clone will be discarded or
#' less than 1% cellularity in all samples.
#' 4. Nodes that have less than 10% cellularity in ALL samples cannot be split further.
#' WARNING: Tree_generation recreates a tree from data while phylo_tree_generation randomly creates a phylogeny
#' @param number_of_clones The wanted number of observable clones (meaning bearing at least 1 mutation)
#' @param number_of_samples The number of samples on which the data should be simulated
#' @keywords Data creation phylogeny
phylo_tree_generation<-function(number_of_clones,number_of_samples){
  
  if(number_of_clones<2){
    warning("Cannot create a phylogenetic tree with less than 2 clones")
    return(NA)
  }
  ### Initialize  
  First<-sample(x = 0:99,replace = T,size = number_of_samples)
  
  if(!check_leaf(First, matrix(rep(100,times = number_of_samples),
                               ncol = number_of_samples )
  )
  ){
    while(!check_leaf(First, matrix(rep(100,times = number_of_samples),
                                    ncol = number_of_samples )
    )){
      First<-sample(x = 0:99,replace = T,size = number_of_samples)
    }
  }
  
  Proportions<-matrix(
    c(rep(100,times = number_of_samples),
      First,
      100 - First),
    byrow = TRUE,
    ncol = number_of_samples,nrow = 3)
  Localization<-c('Ancestral',"L","R")
  can_split<-c(FALSE,
               check_split(Proportions[2,]),
               check_split(Proportions[3,]))
  
  
  probs<-update_probs(Proportions,can_split)
  j<-2
  mutated<-c(TRUE, TRUE,FALSE)
  
  while(j < number_of_clones){
    if(length(probs)==1){
      to_split<-which(can_split)
    }
    else{
      to_split<-sample(which(can_split),size = 1,prob = probs)
    }

    new_leaf<-sapply(X = Proportions[to_split,],FUN = function(p){
      sample(0:p,size = 1)
    }
    )
    opposite<-Proportions[to_split,]-new_leaf
    if(check_leaf(new_leaf,Proportions[mutated,]) && 
       check_leaf(opposite,Proportions[mutated,])){
      j<-j+1
      Proportions<-rbind(Proportions, 
                         matrix(c(new_leaf,opposite),
                                nrow = 2, byrow = TRUE )
      )
      can_split[to_split]<-FALSE
      can_split<-c(can_split,check_split(Proportions),check_split(opposite))
      probs<-update_probs(Proportions,can_split)
      Localization<-c(Localization,paste0(Localization[to_split],c("L","R")))
      mutated<-c(mutated,TRUE,FALSE)
    }
    
  }
  points<-data.frame(mutated,Localization,Proportions)
  colnames(points)<-c("mutated","Localization",1:number_of_samples)
  
  points
  
}

#' Check created leaf
#' 
#' Checks that created leaf has cellularity >1% in at least a sample, that it has at least 2% difference with another
#' existing clone in at least one sample, and that cellularities in all samples are greater than or equal to 0.
#' @param new_leaf A numeric vector to be added
#' @param Proportions_mutated Matrix with samples in columns, clones (carrying mutations) in rows.
check_leaf<-function(new_leaf,Proportions_mutated){
  if(sum(new_leaf<=1)==length(new_leaf) | sum(new_leaf<0)>0 ){ 
    # The cellularity of a clone is less than 1% in ALL samples
    # OR
    # The cellularity is lower than 0 (not biological)
    return(FALSE)
  }
  if(is.matrix(Proportions_mutated)){
    check<-apply(Proportions_mutated,
                 MARGIN = 1,
                 FUN = function(row){
                   max(abs(row-new_leaf))
                 })
  }
  else{
    check<-abs(Proportions_mutated-new_leaf)
  }
  if(sum(check<=2)){ ### There is a less than 2% difference between new clone and another clone in ALL samples
    return(FALSE)
  }
  TRUE
}

#' Check 
#' 
#' Check if node can be split further, i.e. cellularity > 10% in at least 1 sample
#' @param leaf Numeric vector
check_split<-function(leaf){
  ### Check that you still have room for splitting (at least on sample with more than 10% cellularity)
  if(max(leaf)<10){
    return(FALSE)
  }
  return(TRUE)
}

#' Update proportions
#' 
#' Creates vector with probability to be sampled
#' @param Proportions numeric matrix
#' @param can_split logical vector
update_probs<-function(Proportions,can_split){
  w<-which(can_split)
  probs<-numeric(length(w))
  for(i in 1:length(w)){
    if(is.matrix(Proportions)){
      probs[i]<-sum(Proportions[w[i],])
    }
    else{
      probs[i]<-Proportions[w[i]]
    }
  }
  probs/sum(probs)
}