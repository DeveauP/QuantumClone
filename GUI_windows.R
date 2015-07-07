###########################################
#########GUI for QClone####################
###########################################
############Create file browser

create_browser_window<-function(button){
  dialog<-gtkFileChooserDialog(title = "Browse file", 
                               parent = NULL, action = "open",
                               "gtk-ok", GtkResponseType["ok"],
                               "gtk-cancel", GtkResponseType["cancel"],
                               show = T)
  gtkFileChooserSetSelectMultiple(object = dialog,select.multiple = TRUE)
  gSignalConnect(dialog, "response",f = function(dialog,response,data) {
    if(response == GtkResponseType["ok"]) {
      filename<-gtkFileChooserGetFilenames(object = dialog )
      if(button == bSNV){
        SNV.box$setText(paste(lapply(X = filename, FUN = function(z) gsub(pattern = "\\\\",replacement = "/",x = paste(z))),collapse = ","))
      }
      else if(button == bFREEC){
        FREEC.box$setText(gsub(pattern = "\\\\",replacement = "/",x = paste(filename)))
      }
    }
    dialog$destroy()
  })
}

####

##############Convert file for analysis

ConvertSNVFile<-function(path){
  paths<-strsplit(x = path,split=",")[[1]]
   result<-lapply(X=paths,FUN = function(z) read.table(file = z,header = TRUE,sep = "\t",as.is = TRUE))
#   for(i in 1:length(Names)){
#     result[[i]]<-file[grepl(pattern = Names[i],x = file[,1])]
#   }
  return(result)
}

##############Create window for outputs

create.output.window<-function(button,user.data){
  log = gtkWindow("modal")
  log$title<-"LOG"
  LOGbox<-gtkVBox(spacing = 0)
  log$add(LOGbox)
  box1<-gtkEntryNew()
  box2<-gtkEntryNew()
  
  LOGbox$add(box1)
  LOGbox$add(box2)
  box1$setText("Data pre-processing...")
  #################Recreate variables from input for analysis
  
  simulated<-Simulation$active
  save.data<-save.box$active
  if(simulated){
    box1$setText("Creating data...")
    SNV_list<-QuantumCat(number_of_clones = 4, number_of_mutations = 100 ,ploidy="AB",depth=100,
                         number_of_samples=2,Random_clones=F,contamination=c(0,0))
  }
  else{
    print(SNV.box$getText())
    SNV_list<-ConvertSNVFile(as.character(SNV.box$getText()))
  }
  if(sum(grepl(pattern = "Genotype",x = SNV_list))>0){
    FREEC_list<-NULL
  }
  else{
    FREEC_list<-ConvertSNVFile(as.character(SNV.box$getText()))
  }
  print(head(SNV_list[[1]]))
  conta<-eval(parse(text = paste("c(",Contamination$getText(),")",sep="")))
  nclone_range<-eval(parse(text = paste("c(",Clone_range$getText(),")",sep="")))
  clone_priors<-eval(parse(text = PriorClone$getText()))
  prior_weight<-eval(parse(text = PriorWeight$getText()))
  maxit<-eval(parse(text =Maxit$getText()))
  preclustering<-Preclust$active
  save_plot<-save.box$active
  ncores<-eval(parse(text= (Ncore$getText())))
  
  box1$setText("Data pre-processing... DONE")
  box2$setText("Clustering... MAY LAST SEVERAL MINUTES")
  result<-One_step_clustering(SNV_list = SNV_list,FREEC_list = FREEC_list,contamination = conta,nclone_range = nclone_range,
                              plot_3D = Plot3D$active,plot_3D_before_clustering = Plot3D$active,
                              clone_priors = clone_priors,prior_weight =prior_weight ,
                              maxit = maxit,preclustering = preclustering,
                              simulated = simulated,
                              save_plot = save_plot,ncores=ncores)
  
  box2$setText("Clustering... Finished")
  if(save.data){
    write.table(x = cbind(Reduce(f = cbind,x = result$filtered.data),result$cluster,result$EM.output$fik),file = paste(SNV_list[[1]][1,1],"/","clustering_",SNV_list[[1]][1,1],".tsv",sep=""),sep = "\t",row.names = FALSE,quote = FALSE)
  }
  Close<-gtkButtonNewFromStock("gtk-close")
  gSignalConnect(Close, "clicked", log$destroy)
  LOGbox$add(Close)
  
  
}
library(RGtk2)
library(QuantumClone)
library(ggplot2)
library(fpc)
library(doSNOW)
library(parallel)

###Graphical window to launch clonal heterogeneity research
## Inspired by http://tuxette.nathalievilla.org/?p=866&lang=en ; https://github.com/lawremi/RGtk2/blob/master/RGtk2/inst/doc/tutorial.sgml ; http://www.londonr.org/Presentations/Dec%202010%20-%20G%20Heywood%20Ruser%202010-12-06a.pdf
## http://www.r-bloggers.com/playing-with-guis-in-r-with-rgtk2/

window <- gtkWindow()
window["title"] <- "QuantumClone"
window$modifyBg(GtkStateType["normal"], "white")
frame <- gtkFrameNew("Files and options")
window$add(frame)
box0<-gtkVBox(spacing = 0)
frame$add(box0)
frame$modifyBg(GtkStateType["normal"], "orange")
### SNV file box
box0$add(File_dec_box<-gtkTextView())
File_desc <- File_dec_box$GetBuffer()
File_desc$SetText("Input: List of comma separated files
SNV_list is mandatory if \'Simulation\' is not ticked
FREEC_list is mandatory if genotype is not provided inside SNV_list dataframe")
gtkTextViewSetJustification(object = File_dec_box,justification = "center")
  
box1 <- gtkHBoxNew(spacing= 10)
box1$setBorderWidth(12)
box0$add(box1)   #add box1 to the frame

bSNV<-gtkButton("Browse")
gSignalConnect(bSNV, "clicked", f = create_browser_window)
SNV.box<-gtkEntryNew(show = TRUE)
SNV.box$setText("SNV file(s)...")
box1$packStart(SNV.box,fill=TRUE,expand=TRUE)
box1$packStart(bSNV,fill=FALSE,expand=FALSE)


### FREEC file box
box2<- gtkHBoxNew(spacing= 10)
box2$setBorderWidth(12)
box0$add(box2)

bFREEC<-gtkButton("Browse")
gSignalConnect(bFREEC, "clicked", f= create_browser_window)
FREEC.box<-gtkEntryNew()
FREEC.box$setText("FREEC file(s)...")
box2$packStart(FREEC.box,expand = TRUE)
box2$packStart(bFREEC,fill=FALSE,expand = FALSE)

#################################################Options
box0$add(Option_box_txt<-gtkTextView())
Option_txt<-Option_box_txt$GetBuffer()
Option_txt$SetText("*************************************
OPTIONS
*************************************")
gtkTextViewSetJustification(object = Option_box_txt,justification = "center")
# box0$add(File_dec_box<-gtkTextView())
# File_desc <- File_dec_box$GetBuffer()
# File_desc$SetText("Input: List of comma separated files
# SNV_list is mandatory if \'Simulation\' is not ticked
# FREEC_list is mandatory if genotype is not provided inside SNV_list dataframe")
# gtkTextViewSetJustification(object = File_dec_box,justification = "center")


box0$add(Conta_desc_box<-gtkTextView())
Conta_desc<-Conta_desc_box$GetBuffer()
Conta_desc$SetText("Contamination input: comma separated values of the fraction of normal cells in each samples,in the same order as the SNV_list input. (between 0 and 1)
Clone range: comma separated values of the suspected number of clones that will be tested by the algorithm
Save plot: Should the 2-dimensional plots be saved?
Genotype provided: is the genotype provided in the SNV list?")

#################################################Contamination and Clone range
box3<-gtkHBoxNew(spacing= 2)
box3$setBorderWidth(12)

box0$add(box3)   #add main box
Clone_range<-gtkEntryNew()
Clone_range$setText("2:5")
Contamination<-gtkEntryNew()
Contamination$setText("0,0")
Simulation<-gtkCheckButton(label = "Simulation",show = TRUE)
Align<-gtkAlignmentNew(xalign = 0.5)
#Simulation<-gcheckbox(text = "Simulation",checked = F,container = box0,toolkit = "RGtk2)

#setfont_hack(Simulation,list(scale=10))

####Separating the two frames in two sub-boxes

box4<-gtkHBoxNew()
box5<-gtkHBoxNew()
box3$add(box4)
box3$add(box5)
ContaFrame<-gtkFrame(label = "Contamination", show = TRUE)
box4$PackStart(ContaFrame,fill=TRUE,expand = TRUE)


#gtkTextViewSetJustification(object = Conta_desc,justification = "left")

ContaFrame$add(Contamination)
CloneFrame<-gtkFrame(label = "Clone range", show = TRUE)
box5$PackStart(CloneFrame,fill=TRUE,expand = TRUE)
CloneFrame$add(Clone_range)

####Add save plot on same line (box 3)
save.box<-gtkCheckButton(label = "Save plot",show=T)
save.data.box<-gtkCheckButton(label = "Save data in text",show=T)
save.data.box$active=T

save.box$active=T
save.data.box$active=T

box3$add(save.box)
box3$add(save.data.box)

############################## Priors and clustering options
AdvFrame<-gtkFrame(label = "Advanced options", show = TRUE)
AdvFrame$modifyBg(GtkStateType["normal"], "orange")
box0$PackStart(AdvFrame,fill=F)
OptBox<-gtkVBoxNew()
AdvFrame$add(OptBox)

Priors.box<-gtkHBoxNew()
OptBox$add(Priors.box)

PCbox<-gtkHBoxNew()
Priors.box$add(PCbox)
PCFrame<-gtkFrame(label = "Priors on clones", show = TRUE)
PCbox$add(PCFrame)
PWbox<-gtkHBoxNew()
PWFrame<-gtkFrame(label = "Priors on weight", show = TRUE)
PWbox$add(PWFrame)
PCbox$add(PWbox)

PriorClone<-gtkEntryNew()
PriorClone$setWidthChars(20)
PriorClone$setText("NULL")
PriorWeight<-gtkEntryNew()
PriorWeight$setWidthChars(20)
PriorWeight$setText("NULL")

PCFrame$add(PriorClone)
PWFrame$add(PriorWeight)

FinalOptions<-gtkHBoxNew()
OptBox$add(FinalOptions)
Preclust<-gtkCheckButton(label = "Preclustering",show = T)
Preclust$active<-TRUE
FinalOptions$add(Preclust)
Plot3D<-gtkCheckButton(label = "Plot 3D structures",show = T)
Plot3D$active<-F
FinalOptions$add(Plot3D)

MaxitBox<-gtkHBoxNew()
OptBox$add(MaxitBox)
MaxitFrame<-gtkFrame(label = "Number of initial conditions \n to be tested", show = TRUE)
MaxitBox$add(MaxitFrame)
Maxit<-gtkEntryNew()
Maxit$setWidthChars(2)
Maxit$setText("1")
MaxitFrame$add(Maxit)

NcoreBox<-gtkHBoxNew()
OptBox$add(NcoreBox)
NcoreFrame<-gtkFrame(label = "Number of cores", show = TRUE)
NcoreBox$add(NcoreFrame)
Ncore<-gtkEntryNew()
Ncore$setWidthChars(2)
Ncore$setText("1")
NcoreFrame$add(Ncore)

####################################################
####################################################
####################################################
##################Launch / Close
SNV.box$modifyBg(GtkStateType["normal"], "orange")
box0$add(Align)
Align$add(Simulation)
ButtonBox<-gtkHBoxNew()
box0$add(ButtonBox)
Start<-gtkButtonNewFromStock("gtk-ok")

Close<-gtkButtonNewFromStock("gtk-close")
gSignalConnect(Close, "clicked", window$destroy)
ButtonBox$add(Start)
ButtonBox$add(Close)


gSignalConnect(Start,signal = "clicked",create.output.window)