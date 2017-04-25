#' Plot cellularity
#'
#' 2D plot of cellularity based on the output of the EM
#' @param lis Output from Return_one_cell_by_mut, list of cellularities (one list-element per sample)
#' @param Sample_names Name of the samples.
#' @param output_dir Directory in which to save plots
#' @keywords Clonal inference plot
plot_cell_from_Return_out<-function(lis,Sample_names,output_dir=NULL){
    if(length(lis)>1){
        U<-expand.grid(1:length(lis),1:length(lis))
        U<-U[U[,1]<U[,2],]
        #    Sample_names<-lapply(lis,FUN = function(z) z[1,1])
        for(i in 1:nrow(U)){
            d<-ggplot2::qplot(x = lis[[U[i,1]]][,'Cellularity'], y = lis[[U[i,2]]][,'Cellularity'],asp = 1,
                              xlab=paste('Cellular prevalence',Sample_names[[U[i,1]]]),ylab=paste('Cellular prevalence',Sample_names[[U[i,2]]]),
                              main=paste('Cellular prevalence of all possibilities'))+ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+ggplot2::theme_bw()
            if(is.null(output_dir)){
                ggplot2::ggsave(filename = paste(Sample_names[[1]],'/', 'Cellularity', Sample_names[[U[i,1]]],"_",Sample_names[[U[i,2]]],'.pdf',sep=''),plot = d,width = 6.04,height = 6.04)
            }
            else{
                ggplot2::ggsave(filename = paste(output_dir,'/', 'Cellularity', Sample_names[[U[i,1]]],"_",Sample_names[[U[i,2]]],'.pdf',sep=''),plot = d,width = 6.04,height = 6.04)
                
            }
        }
    }
    else{
        d<-ggplot2::qplot(x = lis[[1]][,'Cellularity'], y = jitter(rep(0.5,times=length(lis[[1]][,'Cellularity'])),factor = 10),asp = 1,
                          xlab=paste('Cellular prevalence'),ylab='',
                          main=paste(Sample_names,'Cellular prevalence of all possibilities'))
        d<-d+ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+ggplot2::theme_bw()+ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                                                                                                  axis.ticks.y=ggplot2::element_blank(),
                                                                                                  axis.text.y = ggplot2::element_blank())
        if(is.null(output_dir)){
            ggplot2::ggsave(filename = paste(Sample_names[1],'/', 'Cellularity', Sample_names,'.pdf',sep=''),plot = d,width = 6.04,height = 6.04)
        }
        else{
            ggplot2::ggsave(filename = paste(output_dir,'/', 'Cellularity', Sample_names,'.pdf',sep=''),plot = d,width = 6.04,height = 6.04)
            
        }
    }
}

#' Plots
#'
#' Creates density plot when only one sample is given
#' @param EM_out output from the EM algorithm
#' @param contamination Numeric vector giving the proportion of normal cells in each samples
#' @keywords Clonal inference phylogeny

One_D_plot<-function(EM_out,contamination){
    theta=seq(from = 0,to = 1,by = 0.0001)
    p<-outer(theta,EM_out$filtered.data[[1]]$NC*(1-contamination)/ EM_out$filtered.data[[1]]$NCh)
    P<-choose(EM_out$filtered.data[[1]]$Depth,EM_out$filtered.data[[1]]$Alt)*p**(EM_out$filtered.data[[1]]$Alt)*(1-p)**(EM_out$filtered.data[[1]]$Depth-EM_out$filtered.data[[1]]$Alt)
    y<-apply(X = P,MARGIN = 1,FUN = sum)
    y<-y/sum(y)
    r<-ggplot2::qplot(x=theta,y=y,geom = "line",main="Density of probability of presence of a clone",xlab="Cell fraction",ylab="density")+ggplot2::theme_bw()
    return(r)
}

#'Plot with margin densities
#'
#'Adapted from http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
#'Uses gridExtra package
#' @param QClone_Output Output from QuantumClone algorithm
#' @keywords Plot Densities
#' @export plot_with_margins_densities
#' @examples
#' require(ggplot2)
#' require(gridExtra)
#' message("Using preclustered data:")
#' QC_out<-QuantumClone::QC_output
#' plot_with_margins_densities(QC_out)
#' @importFrom gridExtra grid.arrange
#' 
#' 
plot_with_margins_densities<-function(QClone_Output){
    if(length(QClone_Output$filtered.data)!=2){
        stop("This function can only take 2 samples at a time.")
    }
    sq<-floor(sqrt(max(QClone_Output$cluster)))+1
    
    
    main<-ggplot2::qplot(x=QClone_Output$filtered.data[[1]]$Cellularity,y=QClone_Output$filtered.data[[2]]$Cellularity,color=as.character(QClone_Output$cluster),
                         xlab="Cellularity diag",ylab="Cellulariy relapse",xlim=c(0,1),ylim=c(0,1))+ggplot2::theme_bw()+ggplot2::scale_color_discrete(guide = ggplot2::guide_legend(title="Cluster",ncol=sq))
    
    
    top<-ggplot2::ggplot(QClone_Output$filtered.data[[1]],
                         ggplot2::aes_string("Cellularity"))+ggplot2::geom_density(alpha=.5)+ggplot2::theme_bw()+ggplot2::theme(legend.position="none",
                                                                                                                                axis.title.x=ggplot2::element_blank())
    
    right<-ggplot2::ggplot(QClone_Output$filtered.data[[2]],
                           ggplot2::aes_string("Cellularity"))+ggplot2::geom_density(alpha=.5)+ggplot2::coord_flip()+ggplot2::theme_bw()+ggplot2::theme(legend.position="none",
                                                                                                                                                        axis.title.y=ggplot2::element_blank())
    
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(main)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    
    return(gridExtra::grid.arrange(top,
                                   legend,
                                   main+ggplot2::theme(legend.position="none"),
                                   right,
                                   ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4)))
    
}
#' Plot QC_output
#'
#' This function was implemented to re-plot easily the diagrams of clonality for changes/enhancement.
#' Returns a ggplot object
#' Uses ggplot2 package
#' @param QClone_Output Output from QuantumClone algorithm
#' @param simulated Was the data generated by QuantumCat?
#' @param sample_selected : number of the sample to be considered for plot (can be 1 or 2 samples)
#' @param Sample_names : character vector of the names of each sample (in the same order as the data)
#' @keywords Plot 
#' @importFrom ggplot2 aes_string ggplot theme_bw coord_cartesian xlab ylab
#' @export plot_QC_out
#' @examples
#' require(ggplot2)
#' message("Using preclustered data:")
#' QC_out<-QuantumClone::QC_output
#' plot_QC_out(QC_out,Sample_names = c("Diagnosis","Relapse"))
#' 
plot_QC_out<-function(QClone_Output,Sample_names=NULL, simulated = FALSE,sample_selected = 1:2){
    if(is.null(names(QClone_Output)) && sum(grepl(pattern = "Crit",x = names(QClone_Output[[1]])))>0){
        #### All models are kept
        if(is.null(Sample_names)){
            Sample_names<-unlist(lapply(X = QClone_Output[[1]]$filtered.data,FUN = function(df){
                df[1,1]
            }))
            
        }
        plot_df<-data.frame()
        for(i in 1:length(QClone_Output)){
            plot_df<-rbind(plot_df,
                           cbind(x = QClone_Output[[i]]$filtered.data[[1]]$Cellularity,
                                 y= QClone_Output[[i]]$filtered.data[[2]]$Cellularity,
                                 clone = QClone_Output[[i]]$cluster,
                                 Crit = QClone_Output[[i]]$Crit)
            )
        }
        
        plot_df$clone<-as.factor(plot_df$clone)
        
        result<-ggplot2::ggplot(data=plot_df,
                                ggplot2::aes_string(x = "x", y="y",
                                                    colour = "clone" )
        )+ggplot2::geom_point()+
            ggplot2::scale_colour_discrete(name='Clone')+
            ggplot2::theme_bw()+
            ggplot2::facet_wrap(facets = "Crit",
                                ncol = floor(sqrt(length(QClone_Output)))+1)+
            ggplot2::ggtitle("Criterion comparison")+
            ggplot2::xlim(c(0,1))+
            ggplot2::ylim(c(0,1))+
            ggplot2::xlab(Sample_names[1])+
            ggplot2::ylab(Sample_names[2])
    }
    
    else if(sum(grepl(pattern = "filtered.data",x = names(QClone_Output)))){
        message("Only one model identified")
        Cell <- QClone_Output$filtered.data
        M<-max(as.numeric(as.character(QClone_Output$cluster)))
        cluster<-as.factor(QClone_Output$cluster)
        if(is.null(Sample_names)){
            Sample_names<-unlist(lapply(X = QClone_Output$filtered.data,FUN = function(df){
                df[1,1]
            }))
            
        }
        ### Usual
        if(length(sample_selected)==2){
            message("Two samples identified...")
            result<-list()
            df<-data.frame(x=Cell[[sample_selected[1]]]$Cellularity,
                           y=Cell[[sample_selected[2]]]$Cellularity,
                           colour = cluster)
            if(!simulated){
                
                q<-ggplot2::ggplot(df,aes_string(x = "x",y="y",colour = "colour"))+
                    ggplot2::geom_point()+
                    ggplot2::xlab(paste('Cellular prevalence',Sample_names[sample_selected[1]]))+
                    ggplot2::ylab(paste('Cellular prevalence',Sample_names[sample_selected[2]]))+ 
                    ggplot2::scale_colour_discrete(name='Clone')+
                    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
                    ggplot2::theme_bw()
                
            }
            else{
                df$shape<-as.factor(Cell[[sample_selected[1]]]$Chr)
                
                q<-ggplot2::ggplot(df,aes_string(x = "x",y="y",colour = "colour",shape = "shape"))+
                    ggplot2::geom_point()+
                    ggplot2::xlab(paste('Cellular prevalence',Sample_names[sample_selected[1]]))+
                    ggplot2::ylab(paste('Cellular prevalence',Sample_names[sample_selected[2]]))+ 
                    ggplot2::scale_colour_discrete(name='Clone')+
                    ggplot2::scale_shape_discrete(name='Clone \n(simulated)')+
                    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
                    ggplot2::theme_bw()
            }
            return(q)
        }
        else if(length(sample_selected)==1){
            message("One sample identified...")
            df<-data.frame(x=Cell[[sample_selected[1]]]$Cellularity,
                           y=jitter(rep(0.5,times=length(Cell[[sample_selected[1]]]$Cellularity)),factor = 5),
                           colour = cluster)
            if(!simulated){
                
                result<-ggplot2::ggplot(df,aes_string(x = "x",y="y",colour = "colour") )+
                    ggplot2::geom_point()+
                    ggplot2::xlab(paste('cellularity',Sample_names[sample_selected[1]]))+
                    ggplot2::ylab('')+
                    ggplot2::scale_colour_discrete(name='Clone')+
                    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
                    ggplot2::theme_bw()+
                    ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                                   axis.ticks.y=ggplot2::element_blank(),
                                   panel.background  = ggplot2::element_blank(),
                                   axis.text.y = ggplot2::element_blank())
            }
            else{
                df$shape<-as.factor(Cell[[sample_selected[1]]]$Chr)
                
                result<-ggplot2::ggplot(df,aes_string(x = "x",y="y",colour = "colour",shape = "shape") )+
                    ggplot2::geom_point()+
                    xlab(paste('cellularity',Sample_names[sample_selected[1]]))+
                    ylab('')+
                    ggplot2::scale_colour_discrete(name='Clone')+
                    ggplot2::scale_shape_discrete("Clone \n (simulated)")+
                    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
                    ggplot2::theme_bw()+
                    ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                                   axis.ticks.y=ggplot2::element_blank(),
                                   panel.background  = ggplot2::element_blank(),
                                   axis.text.y = ggplot2::element_blank())
            }
            
        }
        else{
            # Get all possible combination of samples,without redundancy
            grid<-expand.grid(sample_selected,sample_selected)
            # facet on the grid row
            df<-NULL
            if(!simulated){
                for(row in 1:nrow(grid)){
                    df<-rbind(df,
                              data.frame(x=Cell[[sample_selected[grid[row,1]]]]$Cellularity,
                                         y=Cell[[sample_selected[grid[row,2]]]]$Cellularity,
                                         colour = cluster,
                                         facet_x = Sample_names[sample_selected[grid[row,1]]],
                                         facet_y = Sample_names[sample_selected[grid[row,2]]])
                    )
                }
                result<-ggplot2::ggplot(df,aes_string(x = "x",y="y",colour = "colour"))+
                    ggplot2::geom_point()+
                    ggplot2::xlab('Cellular prevalence')+
                    ggplot2::ylab('Cellular prevalence')+ 
                    ggplot2::scale_colour_discrete(name='Clone')+
                    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
                    ggplot2::facet_grid(facets = facet_x ~ facet_y)+
                    ggplot2::theme_bw()
            }
            else{
                shape<- as.factor(Cell[[sample_selected[1]]]$Chr)
                for(row in 1:nrow(grid)){
                    df<-rbind(df,
                              data.frame(x=Cell[[sample_selected[grid[row,1]]]]$Cellularity,
                                         y=Cell[[sample_selected[grid[row,2]]]]$Cellularity,
                                         shape = shape,
                                         colour = cluster,
                                         facet_x = Sample_names[sample_selected[grid[row,1]]],
                                         facet_y = Sample_names[sample_selected[grid[row,2]]])
                    )
                }
                result<-ggplot2::ggplot(df,aes_string(x = "x",y="y",colour = "colour"))+
                    ggplot2::geom_point()+
                    ggplot2::xlab('Cellular prevalence')+
                    ggplot2::ylab('Cellular prevalence')+ 
                    ggplot2::scale_colour_discrete(name='Clone')+
                    ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+
                    ggplot2::facet_grid(facets = facet_x ~ facet_y)+
                    ggplot2::theme_bw()
            }
            #stop("Number of samples can only be 1 or 2 for this function.Use sample_selected parameter.")
        }
    }
    else{
        stop("Incorrect input, please provide output from QuantumClone clustering")
    }
    return(result)
}

#' Evolution plot
#' 
#' Plots evolution in time of clones
#' @param QC_out : Output from One_step_clustering
#' @param Sample_names : character vector of the names of each sample (in the same order as the data)
#' @export evolution_plot
#' @examples 
#' require(ggplot2)
#' evolution_plot(QC_output)
evolution_plot<-function(QC_out,Sample_names=NULL){
    if(sum(grepl(x = names(QC_out$EM.output),pattern = "normalized.centers"))){
        L<-length(QC_out$EM.output$normalized.centers)
    }
    else{
        L<-length(QC_out$EM.output$centers)
    }
    if(is.null(Sample_names)){
        warning(paste("Samples_names is empty, will use 1 to",L))
        Sample_names<-1:L
        
    }
    x<-character()
    y<-numeric()
    col<-rep(1:length(unique(QC_out$cluster)),times = L)
    col<-as.factor(col)
    clone_width<-sapply(col,FUN = function(z){
        sum(as.factor(QC_out$cluster)==z)/length(QC_out$cluster)
    })
    
    for(i in 1:L){
        y<-c(y,QC_out$EM.output$centers[[i]])
        x<-c(x,rep(Sample_names[i],times = length(QC_out$EM.output$centers[[i]])))
    }
    
    df<-data.frame(row.names = 1:length(x))
    df$x<-x
    df$y<-y
    df$col<-col
    df$width<-clone_width
    
    #   q<-ggplot2::qplot(data = df,
    #                     x= x,
    #                     y= y,
    #                     colour = col,
    #                     xlab ="Sample",
    #                     ylab = "Cellularity",geom = "line")+ggplot2::theme_bw()+ggplot2::scale_colour_discrete("Clone") 
    
    q<-ggplot2::ggplot(df,ggplot2::aes_string(x ="x",y="y",
                                              group ="col",
                                              color = "col",
                                              size = "width"),
                       xlab = "Sample",
                       ylab = "Cellularity")+
        ggplot2::geom_line()+
        ggplot2::scale_color_discrete("Clone")+
        ggplot2::scale_size("Fraction of mutations",range = c(0.5,3))+
        ggplot2::xlab("Sample")+
        ggplot2::ylab("Cellularity")+
        ggplot2::theme_bw()
    
    return(q)
}

#' Plots multiple trees 
#'
#' Plots all trees created by the function Tree_generation. The red line means that mutations occured.
#' @param result_list List of lists (tree generated and the probability associated with each tree)
#' @param d Number of clusters found by QuantumClone
#' @param cex Coefficient of expansion for the texts in phylogenetic tree plots. Default is 0.8
#' @export
#' @keywords Clonal inference phylogeny
#' @examples multiplot_trees(QuantumClone::Tree, d= 4)

multiplot_trees<-function(result_list,d,cex=0.8){
    if(length(result_list)%%2==0){
        L<-length(result_list)%/%2
    }
    else{
        L<-length(result_list)%/%2+1
    }
    if(L>1){
        op<-par(mfrow = c(2,L),mar = rep(2, 4))
    }
    for(i in 1:length(result_list)){
        manual_plot_trees(result_list[[i]][[1]],d,cex,result_list[[i]][[2]])
    }
}

#' Plot tree 
#'
#' Creates a visual output for the phylogeny created by Tree_generation()
#' @param connexion_list Data frame of the concatenation of the interaction matrix and the cellularity of each clone at different time points.
#' @param d Number of clusters found by QuantumClone
#' @param cex Coefficient of expansion for the texts in phylogenetic tree plots. Default is 0.8
#' @param p Probability of a tree
#' @export
#' @examples # Extract one tree out of the 3 available trees:
#' Example_tree<-QuantumClone::Tree[[1]]
#' manual_plot_trees(Example_tree[[1]], d= 4,p = Example_tree[[2]])
#' @keywords Clonal inference phylogeny
manual_plot_trees<-function(connexion_list,d,cex=0.8,p){
    s<-dim(connexion_list[[1]][2])
    V<-numeric()
    X<-numeric()
    for(i in 1:(2*d-1)){
        V[i]<-longueur(connexion_list[1:(2*d-1),1:(2*d-1)],i)
        X[i]<-find_x_position(connexion_list[1:(2*d-1),1:(2*d-1)],i,d)
    }
    Y<-1-V/(max(V))
    plot(x=X,y=Y,xlim=c(-1,1),ylim=c(min(Y),1),cex=0, axes = F,xlab='',ylab='',main = paste('p = ',round(p,digits=5)))
    for(i in which(apply(X = connexion_list[1:(2*d-1),1:(2*d-1)],MARGIN = 1,FUN = sum)==2)){
        segments(x0=X[i],x1=X[i],y0=Y[i],y1=Y[i]-1/(max(V)))
        segments(x0=X[which(connexion_list[i,]==1)[1]],x1=X[i],y0=Y[i]-1/(max(V)),y1=Y[i]-1/(max(V)),col='red')
        segments(x0=X[i],x1=X[which(connexion_list[i,]==1)[2]],y0=Y[i]-1/(max(V)),y1=Y[i]-1/(max(V))) 
    }
    if(2*d<dim(connexion_list)[2]){
        LABELS<-apply(X = apply(X = connexion_list[1:(2*d-1),(2*d):(dim(connexion_list)[2])],2,FUN = round,digit=3),1,paste,collapse='\n')
        text(x=X,y=Y,labels = LABELS,pos = 3,cex = cex)
    }
    else{
        LABELS<-sapply(X = connexion_list[1:(2*d-1),(2*d)],FUN = round,digit=3)
        text(x=X,y=Y,labels = LABELS,pos = 3,cex = cex)
    }
}