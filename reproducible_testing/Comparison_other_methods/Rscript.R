if(!require("QuantumClone")){
	if(!require(devtools)){
		install.packages("devtools")
		library("devtools")
	}
	install_github("DeveauP/QuantumClone")
}
if(!require("sciClone")){
	if(!require(devtools)){
		install.packages("devtools")
		library("devtools")
	}
	install_github("genome/sciclone")
	
}
arg <- commandArgs()
print(arg)
range_clones<-eval(parse(text=arg[2])) #number of clones in the samples
range_mutations<-eval(parse(text=arg[3])) # number of mutations


range_ploidy<-as.character(arg[4])
print(range_ploidy)

range_samples<-eval(parse(text=arg[5]))
range_depth<-eval(parse(text=arg[6]))
range_contamination<-eval(parse(text=arg[7]))
number_of_tests<-eval(parse(text=arg[8]))


cnum<-function(x){as.numeric(as.character(x))}
################################
###########Start comparison#####
################################
###########Define functions#####
################################
set.seed(123)
strcount<-function(x,pattern='',split=''){
	unlist(lapply(strsplit(x,split),function(z) na.omit(length(grep(pattern,z)))))
}


Compare_to_sciClone<-function(number_of_tests,contamination = 0,number_of_clones = 4,number_of_mutations = 100,number_of_samples = 2,depth = 100,ploidy = "AB"){
  whole.data<-list()
  result<-data.frame()
  maxit<-3
  corr_BIC<-rep(2:10,each = maxit)*log(number_of_mutations)
  BIC_seq<-seq(from = -2, to = 10, by =0.5) 
  for(i in 1:number_of_tests){
    
    if(!grepl(pattern="ploidy",x = ploidy)){ # if tests are not on ploidy run sciClone
	start.data<-QuantumCat(number_of_clones = number_of_clones,number_of_mutations = number_of_mutations,
	ploidy = ploidy,depth = depth,number_of_samples = number_of_samples,contamination = rep(contamination, times = number_of_samples))
    		whole.data<-c(whole.data,start.data)
    	sci.data<-Format_to_sciClone(start.data,contamination = contamination)
    	ptm<-proc.time()
    	sciCloned<-sciClone(vafs = sci.data,sampleNames = as.character(1:number_of_samples),maximumClusters = 10,minimumDepth = 0)
    	sci.elapsed<-sum((proc.time()-ptm)[1:3])
    	ptm<-proc.time()
    	QCloned<-One_step_clustering(SNV_list = start.data,contamination = rep(contamination,times = number_of_samples),
		maxit = maxit,preclustering = TRUE,save_plot = FALSE,
		nclone_range=2:10,ncores = 4, keep.all.models = TRUE)
    	QC.elapsed<-sum((proc.time()-ptm)[1:3])

    }
    else{ # run QC with single possibility
	range_ploidy<-as.numeric(as.character(strsplit(x = ploidy, split = ":",fixed = TRUE)[[1]][2]))
	if(is.na(range_ploidy)) range_ploidy<-as.character(strsplit(x = ploidy, split = ":",fixed = TRUE)[[1]][2])
	
	start.data<-QuantumCat(number_of_clones = number_of_clones,number_of_mutations = number_of_mutations,
	ploidy = range_ploidy,depth = depth,number_of_samples = number_of_samples,contamination = rep(contamination, times = number_of_samples))
    		whole.data<-c(whole.data,start.data)
	ptm<-proc.time()

	QCloned.single<-One_step_clustering(SNV_list = start.data,contamination = rep(contamination,times = number_of_samples),
		maxit = maxit,preclustering = TRUE,save_plot = FALSE,
		nclone_range=2:10,ncores = 4,force.single.copy = TRUE,keep.all.models = TRUE)
    	QCloned.single.elapsed<-sum((proc.time()-ptm)[1:3])
   	ptm<-proc.time()
    	QCloned<-One_step_clustering(SNV_list = start.data,contamination = rep(contamination,times = number_of_samples),
		maxit = maxit,preclustering = TRUE,save_plot = FALSE,
		nclone_range=2:10,ncores = 4,keep.all.models = TRUE)
    	QC.elapsed<-sum((proc.time()-ptm)[1:3])

    }
    ptm<-proc.time()
    
    df_fpc<-matrix(0,number_of_mutations,number_of_samples)
    for(i in 1:length(start.data)){
	df_fpc[,i]<-start.data[[i]]$Alt/(start.data[[i]]$Depth)
    }
    ptm<-proc.time()
    fpc.res<-pamk(data=df_fpc,krange = 2:10)
    fpc.elapsed<-sum((proc.time()-ptm)[1:3])

    #### NMI calculation for sciClone ####
    if(!grepl(pattern="ploidy",x = ploidy)){
    	clus<-sciCloned@vafs.merged$cluster
    	P_cluster<-table(clus[clus>0])/length(clus[clus>0])
    	P_clone<-table(sciCloned@vafs.merged$chr[clus>0])/length(sciCloned@vafs.merged$chr[clus>0])
    	H_clone<--sum(P_clone*log(P_clone))
    	H_cluster<--sum(P_cluster*log(P_cluster))
    	A<-aggregate(rep(1, times = length(clus[clus>0])), by = list(x=clus[clus>0],y=sciCloned@vafs.merged$chr[clus>0] ), sum)
    	L<-log(A[,3]/(length(clus[clus>0])*P_cluster[A[,1]]*P_clone[A[,2]]))
    	NMI_sci<-2*sum(A[,3]/length(clus[clus>0])*L)/(H_clone+H_cluster)
    }
    else{

	NMI_QC.single<-numeric(length=length(BIC_seq))
    	Crits<-numeric(length = length(QCloned.single))
    	for(model in 1:length(QCloned.single)){
    	  	Crits<-QCloned.single[[model]]$Crit
    	}
    	mod_index<-1
   	for( k in BIC_seq){
		w <- which.min(Crits + (k-1)*corr_BIC)
		if(length(w)){
			NMI_QC.single[mod_index]<-Compute_NMI(QCloned.single[[w]])
		}
		else{
			NMI_QC.single[mod_index]<-NA
		}
		mod_index<-mod_index+1
    	}
    	names(NMI_QC.single)<-paste0("QC.BIC.single",BIC_seq)
    }
    ### NMI calculation for QClone, with a range of BIC ####
    NMI_QC<-numeric(length=length(BIC_seq))
    Crits<-numeric(length = length(QCloned))
    nclust<-numeric(length = length(BIC_seq))
    for(model in 1:length(QCloned)){
      Crits[model]<-QCloned[[model]]$Crit
    }
    mod_index<-1
    for( k in BIC_seq){
	w <- which.min(Crits + (k-1)*corr_BIC)
	if(length(w) && length(QCloned[[w]]$cluster)){
		NMI_QC[mod_index]<-Compute_NMI(QCloned[[w]])
		nclust[mod_index]<-length(QCloned[[w]]$EM.output$weights)
	}
	else{
		NMI_QC.single[mod_index]<-NA
		nclust[mod_index]<-NA
	}
	mod_index<-mod_index+1
    }
    names(NMI_QC)<-paste0("QC.BIC.",BIC_seq)
    names(nclust)<-paste0("nclust.BIC.",BIC_seq)


    ### NMI calculation for fpc ####
    P_cluster<-table(fpc.res$pamobject$cluster)/length(fpc.res$pamobject$cluster)
    P_clone<-table(start.data[[1]]$Chr)/length(start.data[[1]]$Chr)
    H_clone<--sum(P_clone*log(P_clone))
    H_cluster<--sum(P_cluster*log(P_cluster))
    A<-aggregate(rep(1, times = length(fpc.res$pamobject$cluster)), by = list(x= fpc.res$pamobject$cluster,y=start.data[[1]]$Chr ), sum)
    L<-log(A[,3]/(length(fpc.res$pamobject$cluster)*P_cluster[A[,1]]*P_clone[A[,2]]))
    NMI_fpc<-2*sum(A[,3]/length(fpc.res$pamobject$cluster)*L)/(H_clone+H_cluster)
	
    ### binding results
   if(!grepl(pattern="ploidy",x = ploidy)){
   spare<-data.frame(Id= i,t(nclust),
		NMI_sci = NMI_sci,
		t(NMI_QC),
		NMI_kmedoid = NMI_fpc,
 		Time.sci = sci.elapsed,
		Time.QC = QC.elapsed,
		Time.kmedoid = fpc.elapsed,
		Mutations_left.sci = length(clus[clus>0]),
		Mutations_left.QC = length(QCloned$cluster)
		)
   #colnames(spare)<-c("Id","NMI_sci","NMI_QC","NMI_kmedoid","Time(sci)","Time(QC)","Time(kmedoid)","Mutations_left(sci)","Mutations_left_QC")

   }

   else{
	spare<-data.frame(Id = i,t(nclust),
				t(NMI_QC),
				t(NMI_QC.single),
				NMI_kmedoid= NMI_fpc,
				Time.QC = QC.elapsed,
				Time.QC.single = QCloned.single.elapsed,
				Time.kmedoid = fpc.elapsed,
		Mutations_left_QC =length(QCloned$cluster)
	)
   	#colnames(spare)<-c("Id","NMI_sci","NMI_QC","NMI_kmedoid","Time(sci)","Time(QC)","Time(kmedoid)","Mutations_left(sci)","Mutations_left_QC")

   }
  result<-rbind(result,spare)
  }
  write.table(x=result,file="results_sciClone_QCclone.txt",quote=FALSE,sep="\t",row.names=FALSE)
  return(whole.data)
}
create_for_pyClone<-function(QCat){
  result<-list()
  # mutation_id ref_counts var_count normal_cn minor_cn major_cn
  for(i in 1:length(QCat)){
    print(head(QCat[[i]]))
    mutation_id<-apply(X = QCat[[1]][,c("Chr","Start")],MARGIN = 1,FUN = function(z) paste("chr",z[1],":",z[2],":sample",sep = ''))
    ref_counts<-QCat[[i]]$Depth-QCat[[i]]$Alt
    var_counts<-QCat[[i]]$Alt
    minor_cn<-sapply(X = as.character(QCat[[i]]$Genotype),FUN = function(z) strcount(x=z,pattern = "B"))
    major_cn<-sapply(X = as.character(QCat[[i]]$Genotype),FUN = function(z) strcount(x=z,pattern = "A"))
    normal_cn<-rep(2,times=dim(QCat[[1]])[1])
	print(head(cbind(mutation_id,ref_counts,var_counts,normal_cn,minor_cn,major_cn)))
	print(dim(cbind(mutation_id,ref_counts,var_counts,normal_cn,minor_cn,major_cn))) 
   result[[i]]<-cbind(mutation_id,ref_counts,var_counts,normal_cn,minor_cn,major_cn) 
 }
 #print("result")
 #print(head(result[[1]]))
 #print(head(result[[2]]))
 return(result)
}
Format_to_sciClone<-function(QuantumCat_out,contamination = 0){
  test<-rep(TRUE,times = dim(QuantumCat_out[[1]])[1])
  result<-list()
  for(i in 1:length(QuantumCat_out)){
    test<-test & QuantumCat_out[[i]]$Genotype=="AB"
  }
  for(i in 1:length(QuantumCat_out)){
    QCf<-QuantumCat_out[[i]][test,]
    
    result[[i]]<-data.frame(cbind(QCf$Chr,QCf$Start,
                             QCf$Depth-QCf$Alt,QCf$Alt,
                            (QCf$Alt/QCf$Depth)*100 /(1-contamination)))
  }
  return(result)
}
Filter_on_AB<-function(QuantumCat_out){
  test<-rep(TRUE,times = dim(QuantumCat_out[[1]])[1])
  result<-list()
  for(i in 1:length(QuantumCat_out)){
    test<-test & QuantumCat_out[[i]]$Genotype=="AB"
  }
  for(i in 1:length(QuantumCat_out)){
    result[[i]]<-QuantumCat_out[[i]][test,]
  }
  return(result)
}
###################################
#######Use functions###############
###################################
print("Starting with sciClone and QClone")
valid.data<-Compare_to_sciClone(number_of_tests,contamination = range_contamination,
                                number_of_clones = range_clones,number_of_mutations = range_mutations,
                                number_of_samples = range_samples,depth = range_depth,
                                ploidy = range_ploidy)

for(i in 1:length(valid.data)){
	if(range_samples>1){
		if((i%%range_samples)==1){
			#print('Condition OK')
			dir.create(paste("test",i%/%range_samples+1,sep=""),showWarnings=FALSE)	
	   		spare<-create_for_pyClone(valid.data[i:(i+range_samples-1)])
   			#print(length(spare))
   			for(j in 1:length(spare)){
   				write.table(x=spare[[j]],file=paste("./test",i%/%range_samples+1,"/sample",j,".tsv",sep=""),quote=FALSE,row.names=F,sep="\t")
				
			}
			save(valid.data, file = paste("./test",i%/%range_samples+1,"/NMI_sci_QC.Rda",sep=""))
		}   
  }
	else if(range_samples==1){
		dir.create(paste("test",i,sep=""),showWarnings=FALSE)	
	   	spare<-create_for_pyClone(valid.data[i:(i+range_samples-1)])
   		#print(length(spare))
   		for(j in 1:length(spare)){
   			write.table(x=spare[[j]],file=paste("./test",i,"/sample",1,".tsv",sep=""),quote=FALSE,row.names=F,sep="\t")
		}
			save(valid.data, file = paste("./test",i,"/NMI_sci_QC.Rda",sep=""))
	}
}
