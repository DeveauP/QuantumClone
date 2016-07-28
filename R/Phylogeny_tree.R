#' Phylogenetic tree
#'
#' Generates a list of possible trees based on the cellularity of each clone, and the spatial and temporal distribution of the samples.
#' Assumption is made the different clones are on different lines of the matrix
#' @param Clone_cellularities A dataframe with cellularities (ranging from 0 to 1) of each clone (rows) in each sample (columns)
#' @param timepoints A numeric vector giving the spatial and/or temporal distribution of the samples
#' @export
#' @keywords Clonal inference phylogeny
# Clone_cell<-cbind(QuantumClone::QC_output$EM.output$centers[[1]],QuantumClone::QC_output$EM.output$centers[[2]])
# print("Using clone cellularities:")
# print(Clone_cell)
# Tree_generation(Clone_cell)
Tree_generation<-function(Clone_cellularities,timepoints=NULL){
  if(is.data.frame(Clone_cellularities)){
    Clone_cellularities<-as.matrix(Clone_cellularities)
  }
  #Inclusion detection
  if(is.null(timepoints)){
    timepoints<-rep(1,times=ncol(Clone_cellularities))
  }
  else if(length(timepoints)!=ncol(Clone_cellularities)){
    warning('Temporal vector and number of samples are of different sizes')
  }
  nr<-nrow(Clone_cellularities)
  nc<-ncol(Clone_cellularities)
  
  ####Creating network matrix (inclusion), 1 at line i col j means i is included in j
  M<-matrix(0,nr,nr)
  for(i in 1:nr){
    for(j in 1:nr){
      if(i!=j & is_included(Clone_cellularities[i,],Clone_cellularities[j,])){
        M[j,i]<-1
      }
    }
  }
  ###Re-ordering so that you go from root to leaves
  
  
  S<-apply(M,1,sum)
  
  Clone_cellularities<-Clone_cellularities[order(S,decreasing = F),]
  S<-S[order(S,decreasing = F)]
  
  connexion_list<-matrix(0,nrow=(2*nr-1),ncol = (2*nr-1))
  
  if(is.matrix(Clone_cellularities)){
    ###adding non mutated clones
    connexion_list<-cbind(connexion_list,rbind(Clone_cellularities,matrix(0,nrow=nr-1,ncol=nc)))
  }
  else{
    connexion_list<-cbind(connexion_list,c(Clone_cellularities,rep(0,times=nc-1)))
  }
  start<-min(which(S>0))
  for(i in (start-1):nr){
    remove<-0
    if(i==(start-1)){
      connexion_list<-list(list(connexion_list,1))
    }
    else{
      for(k in 1:length(connexion_list)){
        if(is.matrix(Clone_cellularities)){
          t<-add_leaf_list(leaf = Clone_cellularities[i,],connexion_list = connexion_list[[k]],
                           timepoints = timepoints,d = nr,selector_position = i)
        }
        else{
          t<-add_leaf_list(leaf = Clone_cellularities[i],connexion_list = connexion_list[[k]],
                           timepoints = timepoints,d=nr,selector_position = i)
        }
        if(!is.list(t)){
          if(remove==0){
            remove<-k
          }
          else{
            remove<-c(remove,k)
          }
        }
        
        else if(length(t)==2 & !is.list(t[[1]])){
          connexion_list[[k]]<-t
        }
        else{
          connexion_list[[k]]<-t[[1]]
          connexion_list<-c(connexion_list,t[-1])
        }
      }
    }
    if(sum(remove)>0){
      connexion_list<-connexion_list[-remove]
    }
  }
  return(connexion_list)
}

#' Phylogenetic tree leaf
#'
#' Adds a leaf to an already built tree. Output is a list of all possibilities.
#' @param leaf A vector of cellularities (ranging from 0 to 1)
#' @param connexion_list List containing 1. An interaction matrix concatenated with the cellularity of each cluster (one line per cluster)
#' @param timepoints A numeric vector giving the spatial and/or temporal distribution of the samples
#' @param d The initial number of clusters
#' @param selector_position The row of the studied leaf in the data frame.
#' @keywords Clonal inference phylogeny

add_leaf_list<-function(leaf,connexion_list,timepoints,d,selector_position){
  
  Exclude<-apply(X = connexion_list[[1]][,1:(2*d-1)],MARGIN = 1,FUN = sum)==2
  
  if(ncol(connexion_list[[1]])-2*d>0){
    Inclusion<-apply(X = connexion_list[[1]][,(2*d):(ncol(connexion_list[[1]]))],MARGIN = 1,FUN = is_included,leaf)
  }
  else{
    Inclusion<-sapply(X = connexion_list[[1]][,(2*d)],FUN = is_included,leaf)
  }
  Inclusion[selector_position]<-F
  
  S<-sum(Inclusion & !Exclude)
  if(S==0){
    return(NA)
  }
  else if(S>=1){
    w<-which(Inclusion & !Exclude)
    result<-list()
    t<-min(timepoints[leaf>0])
    if(ncol(connexion_list[[1]])-2*d>0){
      size_at_t<-apply(X = connexion_list[[1]][,(2*d):(ncol(connexion_list[[1]]))],MARGIN = 1,FUN = sum)
    }
    else{
      size_at_t<-sapply(X = connexion_list[[1]][,(2*d)],FUN = sum)
    }
    prob_at_t<-size_at_t/(sum(size_at_t[w]))
    if(length(w)==1){
      spare<-connexion_list[[1]]
      spare[w,c(selector_position,selector_position+d-1)]<-1
      spare[selector_position+d-1,(2*d):(ncol(connexion_list[[1]]))]<-spare[w,(2*d):(ncol(connexion_list[[1]]))]-leaf
      prob<-prob_at_t[w]*connexion_list[[2]]
      result<-list(list(spare,prob))
    }
    else{
      for(i in w){
        spare<-connexion_list[[1]]
        spare[i,c(selector_position,selector_position+d-1)]<-1
        spare[selector_position+d-1,(2*d):(ncol(connexion_list[[1]]))]<-spare[i,(2*d):(ncol(connexion_list[[1]]))]-leaf
        prob<-prob_at_t[i]*connexion_list[[2]]
        result<-c(result,list(list(spare,prob)))
      }
    }
  }
  return(result)
}

#' Length
#'
#' Computes the length from the clone on the n-th row of the matrix, to the most ancestral clone
#' @param matrix The interaction matrix of the tree (1 on the i-th row j-th column means "clone j is the progeny of clone i")
#' @param n Index of the clone in the matrix
#' @keywords Clonal inference phylogeny

longueur<-function(matrix,n){
  if(sum(matrix[,n])==0){
    return(0)
  }
  else{
    return(longueur(matrix,which(matrix[,n]==1))+1)
  }
}



#' Graphic position
#'
#' Computes the position of a node on the graph, based on the interaction matrix.
#' @param matrix The interaction matrix of the tree (1 on the i-th row j-th column means "clone j is the progeny of clone i")
#' @param d Initial number of clones
#' @param n Index of the clone of interest in the matrix
#' @keywords Clonal inference phylogeny

find_x_position<-function(matrix,n,d){
  if(sum(matrix[,n])==0){
    return(0)
  }
  else if(n<=d){
    return(find_x_position(matrix,which(matrix[,n]==1),d)-1/2**longueur(matrix,n))
  }
  else{
    return(find_x_position(matrix,which(matrix[,n]==1),d)+1/2**longueur(matrix,n))
  }
}

#' Group theory
#'
#' Clone2 is included in Clone1 if all values of Clone2 are lower or equal to the ones in Clone1 at the same position. Returns TRUE is Clone2 is included in Clone1.
#' @param Clone1 Numeric vector, representing the cellularity of Clone1 in different samples
#' @param Clone2 Numeric vector, representing the cellularity of Clone2 in different samples
#' @keywords Clonal inference phylogeny

is_included<-function(Clone1,Clone2){#Returns True if Clone2 is included in Clone1
  for(i in 1:length(Clone1)){
    if(Clone1[i]<Clone2[i]){
      return(F)
    }
  }
  return(T)
}