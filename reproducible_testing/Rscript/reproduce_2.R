### Source code to reproduce the filtering pipeline

strcount <- function(x, pattern='', split=''){
  unlist(lapply(strsplit(x, split),function(z) na.omit(length(grep(pattern, z)))))
}

### Create stringent dataset

QuantumCat_stringent<-function(number_of_clones = 6,number_of_mutations = 200,
                               ploidy = "AB",depth = 100,
                               contamination = c(0.3,0.4),min_depth = 50){
  if(depth<min_depth){
    error("The depth is lower than the minimal depth")
  }
  variants<-QuantumCat(number_of_clones = number_of_clones,number_of_mutations = number_of_mutations,
                       ploidy = ploidy,depth = depth,number_of_samples= 2,
                       contamination = contamination)
  
  test<-{variants[[1]]$Depth<min_depth} | {variants[[2]]$Depth < min_depth}
  
  while(sum(test)>0){ ### Check that all variants are above 50X depth
    for(i in 1:2){
      sparetest<-variants[[i]]$Depth<min_depth
      variants[[i]]$Depth[sparetest]<-rnbinom(n=sum(sparetest),size = 4.331601,mu = depth)
      variants[[i]]$Alt[sparetest]<-rbinom(n = sum(sparetest),
                                           size=variants[[i]]$Depth[sparetest],
                                           prob= variants[[i]]$Freq[sparetest]/100
      )
      
    }
    test<-{variants[[1]]$Depth<min_depth} | {variants[[2]]$Depth < min_depth}
  }
  return(variants)
}

QuantumCat_permissive<-function(fromQuantumCat,number_of_mutations = 200,
                                ploidy = "AB",depth = 100,
                                contamination = c(0.3,0.4),max_depth = 50, min_depth = 30){
  if(depth<min_depth | max_depth < min_depth){
    error("Depth error")
  }
  clones<-1:max(fromQuantumCat[[1]]$Chr) # retrieve clones
  # Now retrieve cellularity of each clone:
  cells<-list()
  for(i in 1:2){
    cells[[i]]<-sapply(clones,function(z){
      fromQuantumCat[[1]]$Cellularit[fromQuantumCat[[1]]$Chr==z]
    }
    )
  }
  
  
  
  test<-{variants[[1]]$Depth<min_depth} | {variants[[2]]$Depth < min_depth} | 
  {variants[[1]]$Depth>max_depth} | {variants[[2]]$Depth > max_depth}
  while(sum(test)>0){ ### Check that all variants are above 50X depth
    for(i in 1:2){
      sparetest<-variants[[i]]$Depth<min_depth
      variants[[i]]$Depth[sparetest]<-rnbinom(n=sum(sparetest),size = 4.331601,mu = depth)
      variants[[i]]$Alt[sparetest]<-rbinom(n = sum(sparetest),
                                           size=variants[[i]]$Depth[sparetest],
                                           prob= variants[[i]]$Freq[sparetest]/100
      )
      
    }
    test<-{variants[[1]]$Depth<min_depth} | {variants[[2]]$Depth < min_depth} | 
    {variants[[1]]$Depth>max_depth} | {variants[[2]]$Depth > max_depth}
    
  }
  return(variants)
}