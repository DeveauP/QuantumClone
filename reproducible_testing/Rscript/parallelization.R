# Testing the parallelization procedure

if(!require("QuantumClone")){
	install.packages("QuantumClone")
}
if(!require("microbenchmark")){
	install.packages("microbenchmark")
}

set.seed(123)
Start.data<-QuantumCat(4,100,"AB",contamination = c(0.3,0.3))

parallel<-function(n){
	### This function optimizes clones in a parallel way
	One_step_clustering(SNV_list = Start.data,
                        contamination =  c(0.3,0.3),
                        simulated = T,maxit = 1,
                        optim = "default",ncores = n,
                        save_plot =FALSE)
						}
parallel2<-function(n,m){
	### This function treats clones in a parallel way except if the number of trials per clone is higher than the number of clones to test
	One_step_clustering(SNV_list = Start.data,
                        contamination =  c(0.3,0.3),
                        simulated = T,maxit = m,
                        optim = "default",ncores = n,
                        save_plot =FALSE)
						}

mb<-microbenchmark(parallel(1),
					parallel(2),
					parallel(3),
					parallel(4),
					times = 5)

autoplot(mb)

mb2<-microbenchmark(
					parallel2(n = 4, k = 1),
					parallel2(n = 4, k = 5),
					parallel2(n = 4, k = 6),
					times = 1
				)
				
autoplot(mb2)