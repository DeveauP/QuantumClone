reproduce<-function(iter = 9){
  set.seed(345)
  for(z in 1:iter){
    #CREATING DATA
    toy.data<-QuantumCat(number_of_clones = 6,number_of_mutations = 200,
                         ploidy = "AB",depth = 100,number_of_samples = 2,
                         contamination = c(0.3,0.4))
    keep<- sapply(X = 1:nrow(toy.data[[1]]),
                  FUN = function (z){
                    toy.data[[1]][z,"Depth"]>50 && toy.data[[2]][z,"Depth"]>50
                  })
    
    filtered.data<-lapply(X = toy.data, function(df){
      df[keep,]
    })
    noisy_data<-toy.data
    
    Chr<-sample(1:6,size = 200,replace = TRUE)
    Start <- (nrow(noisy_data[[1]])+1):(nrow(noisy_data[[1]])+200)
    purity<-c(0.7,0.6)
    for(j in 1:length(toy.data)){
      Cells<-sapply(1:6, function(z) unique(noisy_data[[j]]$Cellularit[noisy_data[[j]]$Chr==z]))
      
      Genotype <-rep(c("AAB","AABB"),times = c(150,50))
      NCh<-rep(c(3,4),times = c(150,50))
      
      number_of_copies<-c(sample(x = 1:2,size = 150,replace = TRUE,prob = c(3/4,1/4)),
                          sample(x = 1:2,size = 50,replace = TRUE)
      )
      Cellularit<-sapply(Chr,function(n) Cells[n])
      Frequency <- as.numeric(as.character(Cellularit))*purity[j]*number_of_copies/NCh
      Depth<-rnbinom(n=200,size = 4.331601,mu = 100)
      Alt<-rbinom(n=200,size = Depth,prob = Frequency/100)
      noisy_data[[j]]<-rbind(noisy_data[[j]],
                             data.frame(Chr = Chr,Start = Start,Genotype = Genotype,
                                        Cellularit = Cellularit,
                                        number_of_copies = number_of_copies,
                                        Frequency = Frequency, Depth = Depth, Alt = Alt))
      
      
    }
    drivers<-c(sample(x = which(keep),size = round(ndrivers/4),replace = FALSE),
               sample(x = (1:nrow(noisy_data[[1]]))[!keep],size= 3*round(ndrivers/4),replace = FALSE)
    )
    
    log_keep<-rep(FALSE,times = nrow(noisy_data[[1]]))
    log_keep[keep]<-TRUE
    
    log_drivers<-rep(FALSE,times = nrow(noisy_data[[1]]))
    log_drivers[drivers]<-TRUE
    
    
    drivers_data<-lapply(noisy_data,function(df) df[drivers,])
    drivers_id<-drivers_data[[1]]$Start
    
    id_drivers_less_stringent<-noisy_data[[1]]$Start[log_drivers & !log_keep]

    ### Clustering
    Paper_pipeline<-paper_pipeline(cleaned.data = filtered.data,drivers = drivers_data)
    All<-All_clustering(noisy = noisy_data)
    Filter_drivers<-filter_drivers(cleaned.data = filtered.data,
                                   noisy_data,
                                   id_drivers_less_stringent)
    
    
    ### Post processing
    if(nrow(Paper_pipeline$Driver_clusts)> nrow(drivers_data[[1]])){
      # some variants have multiple possibilities
      #colnames(Paper_pipeline$Driver_clusts)<-1:ncol(Paper_pipeline$Driver_clusts)
      Paper_pipeline$Driver_clusts<-cbind(Paper_pipeline$Driver_clusts,
                                          rep(drivers_data[[1]]$Start,
                                              times = strcount(x = as.character(drivers_data[[1]]$Genotype),
                                                               pattern  = "A")*strcount(x = as.character(drivers_data[[2]]$Genotype),
                                                                                        pattern  = "A")
                                          )
      )
      colnames( Paper_pipeline$Driver_clusts)<-c(1:(ncol(Paper_pipeline$Driver_clusts)-1),"Start")
      uStart<-unique(Paper_pipeline$Driver_clusts[,"Start"])
      cluster<-numeric(length = length(uStart))
      index<-0
      adjCol<-ncol(Paper_pipeline$Driver_clusts) -1 
      for(i in uStart){
        index<-index+1
        test<-Paper_pipeline$Driver_clusts[,"Start"]==i
        if(sum(test)>1){
          cluster[index]<-which.max(apply(X = Paper_pipeline$Driver_clusts[test,1:adjCol],
                                          MARGIN = 2,
                                          FUN = max)
          )
        }else{
          cluster[index]<-which.max(Paper_pipeline$Driver_clusts[test,1:adjCol])
        }
      }
      cluster[cluster==0]<-adjCol
      Paper_pipeline$probs<-Paper_pipeline$Driver_clusts
      Paper_pipeline$Driver_clusts<-cluster
    }else{
      Paper_pipeline$probs<-Paper_pipeline$Driver_clusts
      Paper_pipeline$cluster<-apply(X = Paper_pipeline$probs, MARGIN = 1,
                                    FUN = which.max)
    }
    ### Computing results 
    NMI<-c(Compute_NMI(Paper_pipeline$Clusters),
           Compute_NMI(All),
           Compute_NMI(Filter_drivers)
    )
    
    MaxDistErrClus<-c(MaxDistanceCluster(centers = Paper_pipeline$Clusters$EM.output$centers,
                                         filtered.data = Paper_pipeline$Clusters$filtered.data),
                      MaxDistanceCluster(centers = All$EM.output$centers,filtered.data = All$filtered.data),
                      MaxDistanceCluster(centers = Filter_drivers$EM.output$centers,filtered.data = Filter_drivers$filtered.data)
    )
    
    nclusters<-c(max(Paper_pipeline$Clusters$cluster),
                 max(All$cluster),
                 max(Filter_drivers$cluster))
    All_ord<-order(as.numeric(as.character(
      Filter_drivers$filtered.data[[1]]$Start[All$filtered.data[[1]]$Start %in% drivers_data[[1]]$Start]
    )))
    
    ### Starting error distance computation (on drivers)
    
    a<-sqrt((Paper_pipeline$Clusters$EM.output$centers[[1]][Paper_pipeline$Driver_clusts] - drivers_data[[1]]$Cellularit/100)**2
            + (Paper_pipeline$Clusters$EM.output$centers[[2]][Paper_pipeline$Driver_clusts] - drivers_data[[2]]$Cellularit/100)**2
    )
    
    reordered_log_drivers<-All$filtered.data[[1]]$Start %in% drivers_id
    b<-sqrt((All$EM.output$center[[1]][All$cluster[reordered_log_drivers]] - All$filtered.data[[1]]$Cellularit[reordered_log_drivers]/100 )**2+
              (All$EM.output$center[[2]][All$cluster[reordered_log_drivers]] - All$filtered.data[[2]]$Cellularit[reordered_log_drivers]/100 )**2)
    
    reordered_log_drivers<-Filter_drivers$filtered.data[[1]]$Start %in% drivers_id
    c<-sqrt((Filter_drivers$EM.output$center[[1]][Filter_drivers$cluster[reordered_log_drivers]] - Filter_drivers$filtered.data[[1]]$Cellularit[reordered_log_drivers]/100 )**2+
              (Filter_drivers$EM.output$center[[2]][Filter_drivers$cluster[reordered_log_drivers]] - Filter_drivers$filtered.data[[2]]$Cellularit[reordered_log_drivers]/100 )**2)
    
    
    MaxDrivEr <-c(max(a),
                  max(b),
                  max(c)
    )
    
    
    c<-mean(c)
    MeanDrivEr <-c(mean(a),
                   mean(b),
                   mean(c)
    )
    
    Summ<-data.frame(NMI = NMI,
                     Maximal.Error.Cluster.Distance = MaxDistErrClus,
                     Number.Clusters.Found = nclusters,
                     Maximal.Error.Driver.Distance = MaxDrivEr,
                     Mean.Error.Driver.Distance = MeanDrivEr,
                     id = c("Paper_pipeline","All","Filter_drivers")
    )
    if(z == 1){
      result<-Summ
    }
    else{
      result<-rbind(result,Summ)
      
    }
  }
  return(result)
}