
#' Data generation
#'
#' Creates plausible data as would be observed by genome sequencing
#' @param number_of_clones The wanted number of observable clones (meaning bearing at least 1 mutation)
#' @param number_of_mutations The total observed number of mutations (across all clones)
#' @param number_of_samples The number of samples on which the data should be simulated. Default is 2.
#' @param depth The depth of sequencing (does not account for contamination). Default is 100x
#' @param ploidy The general ploidy of the tumor. Default is 2. If "disomic" : only AB regions will be generated.
#' @param Random_clones Should the number of clones be generated randomly (sample(1:10))
#' @param contamination A numeric vector indicating the fraction of normal cells in each sample.
#' @param Subclonal.CNA.fraction Cell fraction of the subclone that has subclonal CNA
#' @keywords Data generation phylogeny
#' @export
#' @return list of dataframes containing observations of a tumor (1 dataframe / sample)
#' @examples
#' print("Generate small set of mutations from 2 differents clones...")
#' print("...in 1 sample, contaminated at 10% by normal cells")
#'  
#' QuantumCat(number_of_clones=2,number_of_mutations=50,number_of_samples=1,contamination=0.1)
QuantumCat<-function(number_of_clones,number_of_mutations,ploidy=2,depth=100,number_of_samples=2,Random_clones=FALSE,contamination=NULL,
                     Subclonal.CNA.fraction=NULL){
  if(Random_clones){
    number_of_clones<-sample(1:10)
    print(paste("Number of clones",number_of_clones))
  }
  #Tree generation
  if(number_of_clones<=0){
    warning("Invalid number of clones")
    return(NA)
  }
  else{#creates phylogenetic tree of clones
    Tree<-phylo_tree_generation(number_of_clones = number_of_clones,number_of_samples)
  }
  #print(Tree)
  Cellularities<-cbind(as.matrix(Tree[Tree$mutated,3:ncol(Tree)]))
  Clonal_attribution<-(rep(1,times = number_of_mutations))
  while(length(unique(Clonal_attribution))<number_of_clones){
    Clonal_attribution<-sample(x=(1:number_of_clones),size=number_of_mutations,replace = T)
  }
  if(is.null(contamination)){
    conta<-rep(0,times = number_of_samples)
  }
  else{
    conta<-contamination
  }
  A<-list()
  B<-list()
  number_of_copies<-list()
  Genotype<-list()
  i<-0
  test<-T
  if(ploidy=="disomic"){
    for(i in 1:number_of_samples){
      A[[i]]<-rep(1,times = number_of_mutations)
      B[[i]]<-rep(1,times = number_of_mutations)
      Genotype[[i]]<-rep("AB",times = number_of_mutations)
      number_of_copies[[i]]<-rep(1,times = number_of_mutations)
    }
  }
  else if(class(ploidy)=="character"){
    if(length(ploidy)==1){
      for(i in 1:number_of_samples){
        A[[i]]<-rep(strcount(x = ploidy,pattern = 'A'),times = number_of_mutations)
        B[[i]]<-rep(strcount(x = ploidy,pattern = 'B'),times = number_of_mutations)
        Genotype[[i]]<-rep(ploidy,times = number_of_mutations)
        number_of_copies[[i]]<-sample(x = 1:(A[[i]][1]),size = number_of_mutations,replace = TRUE)
      }
    }
    else{
      for(i in 1:number_of_samples){
        A[[i]]<-rep(strcount(x = ploidy[i],pattern = 'A'),times = number_of_mutations)
        B[[i]]<-rep(strcount(x = ploidy[i],pattern = 'B'),times = number_of_mutations)
        Genotype[[i]]<-rep(ploidy[i],times = number_of_mutations)
        number_of_copies[[i]]<-sample(x = 1:(A[[i]][1]),size = number_of_mutations,replace = TRUE)
      }
    }
  }
  else{
    for(n in 1:number_of_mutations){
      c<-rpois(number_of_samples,ploidy/2)
      while(sum(c==0)>0){
        for(s in which(c==0)){
          c[s]<-rpois(1,ploidy/2)
        }
      }
      for(i in 1:number_of_samples){
        if(length(A)<number_of_samples){
          A[[i]]<-c[i]
          B[[i]]<-sample(0:c[i],size=1)
          Genotype[[i]]<-paste(paste(rep('A',times=A[[i]]),collapse=''),paste(rep('B',times=B[[i]]),collapse=''),sep='')
        }
        else{
          A[[i]]<-c(A[[i]],c[i])
          B[[i]]<-c(B[[i]],sample(0:c[i],size=1))
          Genotype[[i]]<-c(Genotype[[i]],paste(paste(rep('A',times=tail(A[[i]],1)),collapse=''),paste(rep('B',times=tail(B[[i]],1)),collapse=""),sep=''))
        }
      }
    }
    
    for(j in 1:number_of_samples){
      number_of_copies[[j]]<-sample(x=1:A[[j]][1],size=1)
      
      for(i in 2:number_of_mutations){
        number_of_copies[[j]]<-c(number_of_copies[[j]],sample(x=1:A[[j]][i],size=1))
      }
    }
  }
  if(number_of_samples>1){
    Cell<-Cellularities[Clonal_attribution[1:number_of_mutations],]
  }
  else{
    Cell<-Cellularities[Clonal_attribution[1:number_of_mutations]]
  }
  recap<-list()
  if(number_of_samples>1){
    for(i in 1:number_of_samples){
        #VAF = Cell * Ncopies/(Ncancer+ {c*Nnormal/(1-c)})
      frequency<-Cell[,i]*number_of_copies[[i]]/{A[[i]]+B[[i]]+conta[i]/(1-conta[i])*2}
      SNP_depth<-rnbinom(n=number_of_mutations,size = 4.331601,mu = depth)
      Alt_depth<-numeric()
      for(j in 1:number_of_mutations){
        Alt_depth[j]<-rbinom(n = 1,size=SNP_depth[j],prob= frequency[j]/100)
      }
      recap[[i]]<-data.frame(Clonal_attribution,1:length(Clonal_attribution),Genotype[[i]],Cell[,i],number_of_copies[[i]],frequency,SNP_depth,Alt_depth)
      colnames(recap[[i]])<-c("Chr","Start",'Genotype','Cellularit','number_of_copies',"Frequency",'Depth',"Alt")
    }
  }
  else if(number_of_samples==1){
    frequency<-Cell*number_of_copies[[1]]/{A[[1]]+B[[1]]+conta/(1-conta)*2}
    SNP_depth<-rnbinom(n=number_of_mutations,size = 4.331601,mu = depth)
    Alt_depth<-numeric()
    for(j in 1:number_of_mutations){
      Alt_depth[j]<-rbinom(n = 1,size=SNP_depth[j],prob= frequency[j]/100)
    }
    Observed_freq<-Alt_depth/SNP_depth*100
    recap[[1]]<-data.frame(Clonal_attribution,1:length(Clonal_attribution),Genotype[[1]],Cell,number_of_copies[[1]],frequency,SNP_depth,Alt_depth)
    colnames(recap[[1]])<-c("Chr","Start",'Genotype','Cellularit','number_of_copies',"Frequency",'Depth',"Alt")
  }
  return(recap)
}