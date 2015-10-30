#'Plot with margin densities
#'
#'Adapted from http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
#'Uses gridExtra package
#' @param QClone_Output Output from QuantumClone algorithm
#' @keywords Plot Densities
#' @export
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
#' @export
#' @examples
#' require(ggplot2)
#' message("Using preclustered data:")
#' QC_out<-QuantumClone::QC_output
#' plot_QC_out(QC_out)
#' 
plot_QC_out<-function(QClone_Output,Sample_names=NULL, simulated = FALSE,sample_selected = 1:2){
  Cell <- QClone_Output$filtered.data
  M<-max(as.numeric(as.character(QClone_Output$cluster)))
  cluster<-factor(QClone_Output$cluster)
  if(is.null(Sample_names)){
    Sample_names<-unlist(lapply(X = QClone_Output$filtered.data,FUN = function(df){
      df[1,1]
    }))
    
  }
  if(length(sample_selected)==2){
    result<-list()
    if(!simulated){
      q<-ggplot2::qplot(x=Cell[[sample_selected[1]]]$Cellularity,y=Cell[[sample_selected[2]]]$Cellularity, asp = 1,main=paste('Cellular prevalence',Sample_names[sample_selected[1]],Sample_names[sample_selected[2]]),
                        xlab=paste('Cellular prevalence',Sample_names[sample_selected[1]]),ylab=paste('Cellular prevalence',Sample_names[sample_selected[2]]), 
                        colour = cluster)+ggplot2::scale_colour_discrete(name='Clone')+ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+ggplot2::theme_bw()
      
    }
    else{
      q<-ggplot2::qplot(x=Cell[[sample_selected[1]]]$Cellularity,y=Cell[[sample_selected[2]]]$Cellularity, asp = 1,main=paste('Cellular prevalence plot',Sample_names[sample_selected[1]],Sample_names[sample_selected[2]]),
                        xlab=paste('Cellular prevalence',Sample_names[sample_selected[1]]),ylab=paste('Cellular prevalence',Sample_names[sample_selected[2]]),
                        colour = cluster,
                        shape=factor(Cell[[sample_selected[1]]]$Chr))+ggplot2::theme_bw()+ggplot2::scale_shape_discrete(factor(1:max(Cell[[sample_selected[1]]][,'Chr'])),
                                                                                                                        name='Clone \n(simulated)')+ggplot2::scale_colour_discrete(name='Cluster')+ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+ggplot2::theme_bw()
    }
    return(q)
  }
  else if(length(sample_selected)==1){
    if(!simulated){
      result<-ggplot2::qplot(x=Cell[[sample_selected[1]]]$Cellularity, y=jitter(rep(0.5,times=length(Cell[[sample_selected[1]]]$Cellularity)),factor = 5) , asp = 1,main=paste('Cellular prevalence',Sample_names[sample_selected[1]]),
                             xlab=paste('cellularity',Sample_names[sample_selected[1]]),ylab='', 
                             colour = cluster)+ggplot2::scale_colour_discrete(name='Clone')+ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+ggplot2::theme_bw()+ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                                                                                                                                                                                 axis.ticks.y=ggplot2::element_blank(),
                                                                                                                                                                                 panel.background  = ggplot2::element_blank(),
                                                                                                                                                                                 axis.text.y = ggplot2::element_blank())
    }
    else{
      result<-ggplot2::qplot(x=Cell[[sample_selected[1]]],y=jitter(rep(0.5,times=length(Cell[[sample_selected[1]]]$Cellularity)),factor = 5), asp = 1,main=paste('Cellular prevalence',Sample_names[sample_selected[1]]),
                             xlab=paste('Cellular prevalence',Sample_names[sample_selected[1]]),ylab='',
                             colour = cluster)+ggplot2::scale_colour_discrete(name='Cluster')+ggplot2::coord_cartesian(xlim=c(0,1),ylim=c(0,1))+ggplot2::scale_shape_discrete(factor(1:max(Cell[[1]][,'Chr'])),
                                                                                                                                                                              name='Clone \n(simulated)')+ggplot2::theme_bw()+ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                                                                                                                                                                                                                                             axis.ticks.y=ggplot2::element_blank(),
                                                                                                                                                                                                                                             panel.background  = ggplot2::element_blank(),
                                                                                                                                                                                                                                             axis.text.y = ggplot2::element_blank())
    }
  }
  else{
    stop("Number of samples can only be 1 or 2 for this function.Use sample_selected parameter.")
  }
  return(result)
}

#' Evolution plot
#' 
#' Plots evolution in time of clones
#' @param QC_out
#' @param Sample_names : character vector of the names of each sample (in the same order as the data)
#' @export
#' @examples 
#' require(ggplot2)
evolution_plot<-function(QC_out,Sample_names=NULL){
  if(is.null(Sample_names)){
    Sample_names<-unlist(lapply(X = QClone_Output$filtered.data,FUN = function(df){
      df[1,1]
    }))
    
  }
  y<-character()
  x<-numeric()
  col<-rep(1:length(unique(QC_out$cluster)),times = length(QC_out$EM.output$centers))
  col<-as.factor(col)
  for(i in 1:length(QC_out$EM.output$centers)){
    x<-c(x,QC_out$EM.output$centers[[i]])
    y<-c(y,rep(Sample_names[i],times = length(QC_out$EM.output$centers[[i]])))
  }
  q<-ggplot2::qplot(x = x,y =y,
                    xlab ="Cellularity",ylab = "Sample",
                    colour = col)+ggplot2::theme_bw()+ggplot2::scale_colour_discrete("Clone")
  return(q)
}