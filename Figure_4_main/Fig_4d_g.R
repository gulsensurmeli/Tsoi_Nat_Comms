#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Kernel plots of pooled cells for publication:
# Author: Sara Gomez-Arnaiz & Sau Yee Tsoi
# Last updated: 17/08/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##-- Load packages
library(tidyr)
library(ggplot2)
library(tidyverse)
library(MASS)
library(dplyr)
library(lattice)
library(gridExtra)
library(RColorBrewer)
library(png)
library(ggdensity)
library(grid)
library(magick)
library(imagefx)
library(openxlsx)


##-- Loading data:
borders_all_levels <-read.csv("cell_borders_for_publication.csv")
RSPORB <-read.csv("RSP_ORB-MO_cell_coordinates_01_07_24.csv")

##-- RSP data
# %>% is to pass the left hand side of the operator to the first argument 
# of the right hand side of the operator

RSP <- RSPORB %>%
  filter(inj_area == 'RSP') # filter rows that match a condition

##-- ORB and MOs data

ORB <- RSPORB %>%
   filter(inj_area == 'ORB')

MOs <- RSPORB %>%
   filter(inj_area == 'MOs')

ORB_MOs <- rbind(ORB,MOs)

ORB_MOs[ORB_MOs == "ORB"] <- "ORB|MOs"
ORB_MOs[ORB_MOs == "MOs"] <- "ORB|MOs"

ORB_RSP <- RSPORB %>%
  filter(inj_area == 'ORB'|inj_area == 'RSP')

coloc <- RSPORB %>%
  filter(inj_area == 'coloc')

RSPName <- 'RSP'
ORB_MOName <- 'ORB|MOs'
ORBName <- 'ORB'
MOName <- 'MOs'
ColocName <- 'Colocalised'

######################
## Preparing the data:
######################

colocplots <- list()  # new empty list to store the coloc scatter plots
colocplots_border <- list()
colorArrayColoc <- c("orange", "limegreen", "orangered3") #coloc, GFP and mCherry respectively
#colorArrayColoc <- c("orange", "orangered3", "limegreen") #coloc, GFP and mCherry respectively
colorArray <- c("limegreen", "orangered3") #coloc, GFP and mCherry respectively
#colorArray <- c("orangered3","limegreen") #coloc, GFP and mCherry respectively
colocplots_kernel <- list()

pax_levels <- c(4:4)
#pax_levels <- c(6:6) # if you want to plot only one level

for(layer in pax_levels){ 
  
  borders_layer <- borders_all_levels %>%
    filter(atlas_no == layer)
  
  RSP_layer <- RSP %>% 
    filter(atlas_no == layer)
  
  ORB_MOs_layer <- ORB_MOs %>% 
    filter(atlas_no == layer)
  
  ORB_RSP_layer <- ORB_RSP %>% 
    filter(atlas_no == layer)

  coloc_layer <- coloc %>%
    filter(atlas_no == layer)
  
  total_layer <- rbind(RSP_layer, ORB_MOs_layer)
  
  total_layer_with_coloc <- rbind(RSP_layer, ORB_MOs_layer,coloc_layer)

  # for kernel only:
  ORB_layer <- ORB %>% 
    filter(atlas_no == layer)
  
  MOs_layer <- MOs %>% 
    filter(atlas_no == layer)
  
  # Limits for kernel plots and scatter plot:
  minXkde <- 1500 
  maxXkde <- 4800 #4700 for 2? #4000
  minYkde <- 2000 #1750 #1750
  #maxYkde <- 5000 # level 1
  maxYkde <- 4800 #4500 #4800
  
  
  minXscat <- 1500 
  maxXscat <- 4800 #4700 #4800 #4000
  minYscat <- 1500
  maxYscat <- 4700
  
  #########################################
  ## Scatter plots of RSP and ORB|MO cells 
  # plotted at the end of the script
  # (but this scatter plots are not the publication ones):
  #########################################
  axisPlotsize = 14
  tickPlotsize = 12
  titlePlotsize = 16
  legendPlotsize = 12
  geom_point_size = 0.5 
  if (nrow(coloc_layer) > 1) {
    p_ggdensity4 <- ggplot(data = total_layer_with_coloc , mapping = aes ( x = X, y = Y,colour = factor(inj_area))) + 
      geom_point(size = geom_point_size) + scale_x_reverse() + scale_y_reverse() +
      scale_color_manual(values=colorArrayColoc) +
      xlim(maxXscat, minXscat) + ylim(maxYscat, minYscat) +
      #xlim(4000, 1500) + ylim(4500, 1750) +
      geom_path(data = borders_layer, color = "black") +
      theme(legend.justification=c(0,1), legend.position=c(0,1)) + theme(legend.title=element_blank()) +
      theme(axis.title.x=element_blank(), #remove x axis labels
            axis.title.y=element_blank()) +  #remove y axis labels+
      theme_minimal() +
      #ggtitle(paste0('Colocalised ',RSPName,' and ',ORB_MOName,' cells--\nLevel ', layer)) + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0('Level ', layer)) + 
      xlab("x- coordinate (px)") + ylab("y- coordinate (px)") +
      theme(axis.text.x = element_text(face="bold", color="black", size=tickPlotsize),
            axis.text.y = element_text(face="bold", color="black",size=tickPlotsize),
            plot.title = element_text(size=titlePlotsize, face= "bold", colour= "black" ),
            axis.title.x = element_text(size=axisPlotsize, face="bold", colour = "black"),    
            axis.title.y = element_text(size=axisPlotsize, face="bold", colour = "black")) +
      theme(legend.title = element_text(size=legendPlotsize, face="bold"),
            legend.text = element_text(size=legendPlotsize, face="bold"),
            legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
            legend.justification=c(0,1), legend.position=c(0,1),
            legend.margin = margin(c(0.5, 1.5, 0.5, 2.5),unit=c("mm","mm","mm","mm")),
            legend.spacing.x = unit(0.2, "mm"),
            legend.spacing.y = unit(0.1, "mm")) + # moving the legend inside the plot
      theme(legend.title = element_blank())
    
    #print(p_ggdensity4)
    
    colocplots[[layer]] <- p_ggdensity4  # add each plot into plot list
  } else {
    p_ggdensity4 <- ggplot(data = total_layer , mapping = aes ( x = X, y = Y,colour = factor(inj_area))) + 
      geom_point(size = geom_point_size) + scale_x_reverse() + scale_y_reverse() +
      scale_color_manual(values=colorArray) +
      #xlim(4000, 1500) + ylim(4500, 1750) +
      xlim(maxXscat, minXscat) + ylim(maxYscat, minYscat) +
      geom_path(data = borders_layer, color = "black") +
      theme(legend.justification=c(0,1), legend.position=c(0,1)) + theme(legend.title=element_blank()) +
      theme(axis.title.x=element_blank(), #remove x axis labels
            axis.title.y=element_blank()) +  #remove y axis labels+
      theme_minimal() +
      #ggtitle(paste0('Colocalised ',RSPName,' and ',ORB_MOName,' cells--\nLevel ', layer)) + theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0('Level ', layer)) + 
      xlab("x- coordinate (px)") + ylab("y- coordinate (px)") +
      theme(axis.text.x = element_text(face="bold", color="black", size=tickPlotsize),
            axis.text.y = element_text(face="bold", color="black",size=tickPlotsize),
            plot.title = element_text(size=titlePlotsize, face= "bold", colour= "black" ),
            axis.title.x = element_text(size=axisPlotsize, face="bold", colour = "black"),    
            axis.title.y = element_text(size=axisPlotsize, face="bold", colour = "black")) +
      theme(legend.title = element_text(size=legendPlotsize, face="bold"),
            legend.text = element_text(size=legendPlotsize, face="bold"),
            legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
            legend.justification=c(0,1), legend.position=c(0,1),
            legend.margin = margin(c(0.5, 1.5, 0.5, 2.5),unit=c("mm","mm","mm","mm")),
            legend.spacing.x = unit(0.2, "mm"),
            legend.spacing.y = unit(0.1, "mm")) + # moving the legend inside the plot
      theme(legend.title = element_blank())
    #print(p_ggdensity4)
    colocplots[[layer]] <- p_ggdensity4  # add each plot into plot list
  }
  
  # if (layer > 4) {
  # coloc_layer
  # }
  p_border_coloc <- ggplot (data = coloc_layer , mapping = aes ( x = X, y = Y)) + geom_point(size = geom_point_size) + scale_x_reverse() + scale_y_reverse() +
    #xlim(maxXs, minXs) + ylim(maxYs, minYs) + 
    #xlim(4000, 1500) + ylim(4500, 1750) + 
    xlim(maxXscat, minXscat) + ylim(maxYscat, minYscat) +
    geom_path(data = borders_layer, color = "black") +
    theme(legend.justification=c(0,1), legend.position=c(0,1)) + theme(legend.title=element_blank()) +
    theme(axis.title.x=element_blank(), #remove x axis labels
          axis.title.y=element_blank()) +  #remove y axis labels+
    #ggtitle(paste0('Colocalised ',RSPName,' and ',ORB_MOName,' cells-\nLevel ', layer)) + 
    ggtitle(paste0('Level ', layer)) + 
    theme(plot.title = element_text(hjust = 0.5))+
    xlab("x- coordinate (px)") + ylab("y- coordinate (px)") +
    theme_minimal() +
    theme(axis.text.x = element_text(face="bold", color="black", size=tickPlotsize),
          axis.text.y = element_text(face="bold", color="black",size=tickPlotsize),
          plot.title = element_text(size=titlePlotsize, face= "bold", colour= "black" ),
          axis.title.x = element_text(size=axisPlotsize, face="bold", colour = "black"),    
          axis.title.y = element_text(size=axisPlotsize, face="bold", colour = "black")) +
    theme(legend.title = element_text(size=legendPlotsize, face="bold"),
          legend.text = element_text(size=legendPlotsize, face="bold"),
          legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
          legend.justification=c(0,1), legend.position=c(0.7,1),
          legend.margin = margin(c(0.5, 1.5, 0.5, 2.5),unit=c("mm","mm","mm","mm")),
          legend.spacing.x = unit(0.2, "mm"),
          legend.spacing.y = unit(0.1, "mm"))  # moving the legend inside the plot
  
  #print(p_border_coloc)
  colocplots_border[[layer]] <- p_border_coloc  # add each plot into plot list
  
  ###############################
  ## Kernel density calculations:
  ###############################
 
  ## Kernel density of ORB region
  f1 <- kde2d(ORB_layer$X, ORB_layer$Y, n = 200, lims = c(minXkde,maxXkde,minYkde,maxYkde))
  df <- t(f1$z)
  df <- df[seq(from=nrow(df),to=1,by=-1),seq(from=ncol(df),to=1,by=-1)] #reversing the rows , SARA
  df1 <- t(df)
  
  dfORB_norm = (df1-min(df1))/(max(df1)-min(df1))
  
  ## Kernel density of MOs region
  f1 <- kde2d(MOs_layer$X, MOs_layer$Y, n = 200, lims = c(minXkde,maxXkde,minYkde,maxYkde))
  df <- t(f1$z)
  df <- df[seq(from=nrow(df),to=1,by=-1),seq(from=ncol(df),to=1,by=-1)] #reversing the rows , SARA
  df1 <- t(df)
  
  dfMOs_norm = (df1-min(df1))/(max(df1)-min(df1))

  ## Kernel density of joint/merged ORB and MOs regions
  f1 <- kde2d(ORB_MOs_layer$X, ORB_MOs_layer$Y, n = 200, lims = c(minXkde,maxXkde,minYkde,maxYkde))
  df <- t(f1$z)
  df <- df[seq(from=nrow(df),to=1,by=-1),seq(from=ncol(df),to=1,by=-1)] #reversing the rows , SARA
  df1 <- t(df)
  
  dfORB_MOs_norm = (df1-min(df1))/(max(df1)-min(df1))
  
  ## Kernel density of RSP
  f2 <- kde2d(RSP_layer$X, RSP_layer$Y, n = 200, lims = c(minXkde,maxXkde,minYkde,maxYkde))
  df <- t(f2$z)
  df <- df[seq(from=nrow(df),to=1,by=-1),seq(from=ncol(df),to=1,by=-1)] #reversing the rows , SARA
  df2 <- t(df)
  
  dfRSPnorm = (df2-min(df2))/(max(df2)-min(df2))
  
  ## Kernel density estimate for colocalised cells
  if (nrow(coloc_layer) > 1) { # calculate only if more than one colocalise cell 
  f2 <- kde2d(coloc_layer$X, coloc_layer$Y, n = 200, lims = c(minXkde,maxXkde,minYkde,maxYkde))
  df <- t(f2$z)
  df <- df[seq(from=nrow(df),to=1,by=-1),seq(from=ncol(df),to=1,by=-1)] #reversing the rows , SARA
  df2 <- t(df)
  
  dfColocnorm = (df2-min(df2))/(max(df2)-min(df2))
  }
  
  ## Kernel density SUBSTRACTION
  dfS_RSP_ORB     <- dfORB_norm-dfRSPnorm
  dfS_RSP_MOs     <- dfMOs_norm-dfRSPnorm
  dfS_RSP_ORB_MOs <- dfORB_MOs_norm-dfRSPnorm

  ################
  ## Kernel plots:
  ################
  pORBn <-levelplot(data.matrix(dfORB_norm),aspect="fill",xlab=NULL,ylab=NULL, colorkey=list(labels=list(cex=1.1),height=0.81, width=1.4,space="right"),scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ORBName,'- ','level ', layer),at=seq(0,1,length.out=100))
  leg.list <- pORBn$legend$right$args$key
  leg <-  draw.colorkey(leg.list)
  pORBn_noLeg <-levelplot(data.matrix(dfORB_norm),aspect="fill",xlab=NULL,ylab=NULL,  colorkey=FALSE,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ORBName,'- ','level ', layer),at=seq(0,1,length.out=100)) 
  
  pMOsn <-levelplot(data.matrix(dfMOs_norm),aspect="fill",xlab=NULL,ylab=NULL, colorkey=list(labels=list(cex=1.1),height=0.81, width=1.4,space="right"),scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(MOName,'- ','level ', layer),at=seq(0,1,length.out=100)) 
  leg.list <- pMOsn$legend$right$args$key
  leg <-  draw.colorkey(leg.list)
  pMOsn_noLeg <-levelplot(data.matrix(dfMOs_norm),aspect="fill",xlab=NULL,ylab=NULL,  colorkey=FALSE,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(MOName,'- ','level ', layer),at=seq(0,1,length.out=100)) 

  pORB_MOs_n <-levelplot(data.matrix(dfORB_MOs_norm),aspect="fill",xlab=NULL,ylab=NULL, colorkey=list(labels=list(cex=1.1),height=0.81, width=1.4,space="right"),scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ORB_MOName,'- ','level ', layer),at=seq(0,1,length.out=100)) #alternating=2 flips x axis labels up to top
  leg.list <- pORB_MOs_n$legend$right$args$key
  leg <-  draw.colorkey(leg.list)
  pORB_MOs_n_noLeg <-levelplot(data.matrix(dfORB_MOs_norm),aspect="fill",xlab=NULL,ylab=NULL,  colorkey=FALSE,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ORB_MOName,'- ','level ', layer),at=seq(0,1,length.out=100)) 
   
  pRSPn_noLeg <-levelplot(data.matrix(dfRSPnorm),aspect="fill",xlab=NULL,ylab=NULL,colorkey=FALSE,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(RSPName,'- ','level ', layer),at=seq(0,1,length.out=100)) 

  pColocn_noLeg <-levelplot(data.matrix(dfColocnorm),aspect="fill",xlab=NULL,ylab=NULL,colorkey=FALSE,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ColocName,'- ','level ', layer),at=seq(0,1,length.out=100)) 
  colocplots_kernel[[layer]] <- pColocn_noLeg
  
  brk <- seq(-1,1,length.out=100)
  cols <- colorRampPalette(brewer.pal(11, "PuOr"))
  pSusN_ORBonly <-levelplot(data.matrix(dfS_RSP_ORB),aspect="fill",col.regions = cols,colorkey = list(labels=list(cex=1.1),height=0.81, width=1.4,col = cols,at = brk,space="right"), xlab=NULL,ylab=NULL,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ORBName,'n - ',RSPName,'n'),at=brk) 
  leg.list <- pSusN_ORBonly$legend$right$args$key
  legS <-  draw.colorkey(leg.list)
  pSusN_ORBonly_noLeg <-levelplot(data.matrix(dfS_RSP_ORB),aspect="fill",col.regions = cols,colorkey = FALSE, xlab=NULL,ylab=NULL,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ORBName,'n - ',RSPName,'n'),at=brk) 
  
  pSusN_MOsonly <-levelplot(data.matrix(dfS_RSP_MOs),aspect="fill",col.regions = cols,colorkey = list(labels=list(cex=1.1),height=0.81, width=1.4,col = cols,at = brk,space="right"), xlab=NULL,ylab=NULL,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(MOName,'n - ',RSPName,'n'),at=brk) 
  leg.list <- pSusN_MOsonly$legend$right$args$key
  legS <-  draw.colorkey(leg.list)
  pSusN_MOsonly_noLeg <-levelplot(data.matrix(dfS_RSP_MOs),aspect="fill",col.regions = cols,colorkey = FALSE, xlab=NULL,ylab=NULL,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(MOName,'n - ',RSPName,'n'),at=brk) 

  pSusN <-levelplot(data.matrix(dfS_RSP_ORB_MOs),aspect="fill",col.regions = cols,colorkey = list(labels=list(cex=1.1),height=0.81, width=1.4,col = cols,at = brk,space="right"), xlab=NULL,ylab=NULL,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ORB_MOName,'n - ',RSPName,'n'),at=brk) 
  leg.list <- pSusN$legend$right$args$key
  legS <-  draw.colorkey(leg.list)
  pSusN_noLeg <-levelplot(data.matrix(dfS_RSP_ORB_MOs),aspect="fill",col.regions = cols,colorkey = FALSE, xlab=NULL,ylab=NULL,scales = list(tck = c(0,0),y=list(draw=FALSE),x=list(draw=FALSE)),main=paste0(ORB_MOName,'n - ',RSPName,'n'),at=brk) 
  
  ##########################
  ## Plotting with ggdensity
  ##########################
  probability = .80 # .99, .95, .80, .50
  p <- ggplot(ORB_RSP_layer, aes(X, Y, fill = inj_area))  +
    scale_fill_manual(values=colorArray) +
    #scale_color_manual(values=colorArray) +
    scale_x_reverse() + scale_y_reverse() +
    xlim(maxXkde, minXkde) + ylim(maxYkde, minYkde) + 
    #theme_classic() +
    theme_void() +
    theme(axis.title.x=element_blank(), #remove x axis labels
          axis.title.y=element_blank(),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(), 
          axis.text.y=element_blank(), 
          axis.ticks.y=element_blank()) +
    labs(fill = "Inj. area:") +
    #theme(legend.position="top")+
    theme(legend.justification=c(0,1), legend.position=c(0,1),
          legend.direction = "horizontal") +
    theme(legend.title = element_text(size=legendPlotsize, face="bold"),
          legend.text = element_text(size=legendPlotsize, face="bold")) +
  ggtitle(paste0('Level ', layer,' prob: ',probability)) 
  #ggtitle(paste0('Level ', layer,' prob ')) 
  
  p_dens_all <- p + geom_hdr() 
  # p_density

  p_dens_sub <- p + geom_hdr(probs = probability, alpha = .8) # best

  # p2 <- p + geom_hdr_lines()
  #p + geom_hdr(probs = .99, alpha = .5)
  
  #############################################
  # Kernel density figure
  #############################################  

  margW <-1/42
  legW <-1/15 #1/9
  figW <-1/3
  blank<-rectGrob(gp=gpar(col="white"))

  figure_g_density <- grid.arrange(blank,
                                  pORBn_noLeg, 
                                  pRSPn_noLeg,leg,
                                  pSusN_ORBonly_noLeg, legS, 
                                  blank, p_dens_sub, blank,
                                  nrow=1, widths = c(margW, figW,  figW, legW, figW, legW, margW,figW,margW),
                                  heights=unit(10, "cm"))
  
}

figure_g_densityColoc <- grid.arrange(colocplots_kernel[[1]])
figure_g_densityColoc <- grid.arrange(colocplots_kernel[[4]])
print(figure_g_densityColoc)

png(file=paste0('mice_pooled_L',4,'_coloc_cells_1200dpi_from01_07_24.png'),
    width=4000, height=4000,res=1200)
plot(figure_g_densityColoc)
dev.off()
figure_g_densityColoc <- grid.arrange(colocplots_kernel[[4]])

#############################################
# Scatter plot for all and colocalised cells
#############################################

## Colocalisation figures - Colocalised cells only
grid.arrange(colocplots_border[[1]])
grid.arrange(colocplots_border[[4]])
grid.arrange(colocplots_border[[6]])

## Colocalisation figures - For all cells, included colocalised cells
grid.arrange(colocplots[[1]])
grid.arrange(colocplots[[4]])
grid.arrange(colocplots[[6]])
