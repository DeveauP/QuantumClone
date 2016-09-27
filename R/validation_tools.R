#' Normalized Mutual Information
#' 
#' Compute the normalized mutual information to assess clustering quality
#' @param  QC_out output from QuantumClone clustering 
#' @examples 
#' Compute_NMI(QC_output)
#' @export
Compute_NMI<-function(QC_out){
  # Probabilities
  cluster<-as.numeric(as.character(QC_out$cluster))
  if(length(unique(cluster))!=max(cluster)){### If a cluster is unused
    spare<-cluster
    for(i in 2:max(cluster)){ # take k for each unused value (smaller than i) of a cluster
      spare[spare==i]<-i-sum(!(1:i %in% cluster[cluster<=i]))
    }
    cluster<-spare
  }
  P_cluster<-table(cluster)/length(cluster)
  P_clone<-table(QC_out$filtered.data[[1]]$Chr)/length(QC_out$filtered.data[[1]]$Chr)
  
  # Information entropy
  H_clone<--sum(P_clone*log(P_clone))
  H_cluster<--sum(P_cluster*log(P_cluster))
  
  A<-aggregate(rep(1, length(cluster)), 
               by = list(x= cluster,
                         y=QC_out$filtered.data[[1]]$Chr ),
               sum)
  
  L<-log(A[,3]/(length(cluster)*P_cluster[A[,1]]*P_clone[A[,2]]))
  NMI<-2*sum(A[,3]/length(cluster)*L)/(H_clone+H_cluster)
  return(NMI)
}

#' NMI
#' 
#' Computes the NMI based on the clustering
#' @param cut_tree a numeric vector of cluster selection
#' @param chr the ground truth for clusters
#' @return numeric value of NMI (between 0 and 1)
#' @examples
#' set.seed(123)
#' #1: Cluster data
#' FQC<-FlashQC(QuantumClone::Input_Example,conta = c(0,0),Nclus = 2:10)
#' 
#' #2: Compute NMI
#' NMI_cutree(FQC$cluster,chr = QuantumClone::Input_Example[[1]]$Chr)
#' @export
NMI_cutree<-function(cut_tree,chr){
  # Probabilities
  cluster<-as.numeric(as.character(cut_tree))
  clones<-chr
  
  P_cluster<-table(cluster)/length(cluster)
  P_clone<-table(clones)/length(clones)
  
  # Information entropy
  H_clone<--sum(P_clone*log(P_clone))
  H_cluster<--sum(P_cluster*log(P_cluster))
  
  A<-aggregate(rep(1, length(cluster)), 
               by = list(x= cluster,
                         y=clones ),
               sum)
  
  L<-log(A[,3]/(length(cluster)*P_cluster[A[,1]]*P_clone[A[,2]]))
  NMI<-2*sum(A[,3]/length(cluster)*L)/(H_clone+H_cluster)
  return(NMI)
  
}

#' Precision
#' 
#' Computes the precision based on the clustering
#' @param hx a numeric vector of cluster selection
#' @param Truth the ground truth for clusters
#' @return \describe{
#'  \item{TP}{ The number of true positive links}
#'  \item{TN}{ The number of true negative links}
#'  \item{FP}{ The number of false positive links}
#'  \item{FN}{ The number of false negative links}
#'  \item{Pr}{ The precision, defined by \eqn{Pr = \frac{TP}{TP+FP}}}
#'  \item{R}{The recall, defined by \eqn{R = \frac{TP}{TP+FN}}}
#'  \item{F1}{The F1 index, defined by \eqn{F1 = \frac{2\times P \times R}{P + R}}}
#'  \item{RI}{Rand Index, defined by \eqn{RI = \frac{TP+TN}{TP+TN+FP+FN}}}
#'  \item{validat}{Is positives + negatives equal to total number of links - returns absolute difference if false }
#'  }
#' @examples
#' set.seed(123)
#' #1: Cluster data
#' FQC<-FlashQC(QuantumClone::Input_Example,conta = c(0,0),Nclus = 2:10)
#' 
#' #2: Compute NMI
#' Precision_Recall(hx = FQC$cluster,Truth = QuantumClone::Input_Example[[1]]$Chr)
#' 
#' ### From Stanford NLP example:
#' cluster<-c(rep(1,6),rep(2,6),rep(3,5))
#' truth<-c(rep(1,5),2,
#'          1,rep(2,4),3,
#'          rep(1,2),rep(3,3))
#' Precision_Recall(cluster,truth)
#' @export
Precision_Recall<-function(hx,Truth){
  tot<-choose(n=length(hx),k=2) ### Total number of edges
  
  Tabhx<-table(as.character(hx)) 
  positives<-sum(choose(n = Tabhx,2)) ## Number of positive edges (TP + FP)
  
  Insersect_sizes<-aggregate(rep(1, length(hx)), 
                             by = list(x= hx,
                                       y=Truth ),
                             sum)
  TP<-sum(choose(n=Insersect_sizes[,3],k=2)) 
  
  Combinations<-expand.grid(1:(length(hx)-1),2:length(hx))
  Combinations<-Combinations[Combinations[,1]<Combinations[,2],]
  
  test<-hx[Combinations[,1]] != hx[Combinations[,2]]
  
  negatives<-sum(test) ### Number of negative edges (TN + FN)
  
  if(negatives + positives == tot){
    validat<-TRUE
  }else{
    validat<-abs(negatives + positives - tot)
  }
  
  TN<-sum( test & Truth[Combinations[,1]] != Truth[Combinations[,2]])
  FN<- negatives - TN
  FP<- positives - TP
  Pr<-TP/(TP+FP)
  R<-TP/(TP+FN)
  F1<-2*Pr*R/(Pr+R)
  RI<-(TP+TN)/tot
  
  return(list(TP = TP, TN = TN,FP = FP, FN = FN, 
         Pr = Pr, R = R, F1 = F1, RI =RI,
         validat = validat,
         nclus = length(unique(hx)))
  )
}
