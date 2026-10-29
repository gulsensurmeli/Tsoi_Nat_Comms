#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Protocol to Assess Random Distribution of Cells (adapted from Surmeli et al., 
# 2011 Figure S4 
# https://www.sciencedirect.com/science/article/pii/S0092867411012104?via%3Dihub#figs2) 
# Author: Sara Gomez-Arnaiz; comments from Gulsen Surmeli protocol; edited by Sau Tsoi
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
library(ggpubr)
library(reshape)

#####clearing the environment#####
rm(list = ls())

#####prepping the datasets#####
RSPORB<-read.csv("all_combined_23_10_2025.csv")

colnames(RSPORB)<-c("X.1", "Label", "Mean", "X", "Y", "Ch", "inj_area", "animal", "section", "atlas_no")

#filling in the inj areas for each animal
inj_area<-data.frame(animal = c(659, 659, 1200, 1200, 1201, 1201, 1207, 1207, 1736895, 1736895),
                     Ch = c(1, 2, 1, 2, 1, 2, 1,2,1,2),
                     inj_area = c("RSP", "ORB", "ORB", "RSP", "RSP", "ORB", "ORB", "RSP", "RSP", "ORB"))

RSPORB<-left_join(RSPORB, inj_area, by = c("animal", "Ch"))
RSPORB<-RSPORB %>% 
  dplyr::select(-inj_area.x) %>% 
  dplyr::rename(Brain_Area = inj_area.y)

#adding coloc
RSPORB$Brain_Area[RSPORB$Ch==3]<-"coloc"

#for shuffling we are only using 659, 1200, 1201, 1207, as the numbers in 1736895 are too low
RSPORB<-RSPORB %>% 
  filter(animal==659|animal==1200|animal==1201|animal==1207)

######creating separate csv file for each animal######
# Get all unique combinations of animal and atlas_no
combinations <- RSPORB %>% distinct(animal, atlas_no)
# Loop through each combination
for (i in seq_len(nrow(combinations))) {
  this_animal <- combinations$animal[i]
  this_atlas <- combinations$atlas_no[i]
  
  # Filter the data for this combination
  subset_data <- RSPORB %>%
    filter(animal == this_animal, atlas_no == this_atlas)
  
  # Build a descriptive file name
  file_name <- sprintf("Datasets/shuffling/Results_%s_cell_counts_L%s.csv",
                       this_animal, this_atlas)
  
  # Save the filtered data
  write.csv(subset_data, file_name, row.names = FALSE)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Defining colocCellsXY and RefCellsXY data frames with coordinates of 
# colocalised and Ref cells for all animals and levels from excels
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#animal_number <- c(659,660, 1642, 1643) 
#animal_number <- c(659,660, 1642, 1643, 1736892, 1736895) 
animal_number <- c(659,1200, 1201,1207) 
#animal_number <- c(1642, 1643) 

area_ref <- "all" # "ORB_MOs", "RSP" & "all"

colocCellsXY <-data.frame()
RefCellsXY <-data.frame()
percentage_pAn_pL <-data.frame()
factorAnimal <-data.frame() # HERE OR OUTSIDE THE LOOP?
factorLayer <- data.frame() # HERE OR OUTSIDE THE LOOP?
zScoresAll <- data.frame()
z_scores_shufflingAll <- data.frame()
z_scores_pL <-data.frame()
actual_dist<-data.frame()
dist_all<-data.frame()
df <- c() # remove?
group <-c() #remove?
layer = 4 #rmb to change atlas no here
for(mouseN in animal_number){ 
  
  allCells <-read.csv(paste0("Datasets/shuffling/Results_",mouseN,"_cell_counts_L",layer,".csv"))
  
  colocCells <- allCells %>%
    filter(Brain_Area== "coloc")
  
  if(nrow(colocCells)>0){
    colocCellsXYaux <- colocCells[c("X","Y","Brain_Area")]
    RefCellsXYaux <- allCells[c("X","Y","Brain_Area")]
    
    if (area_ref == "ORB_MOs") {
      RefCellsXYaux <- RefCellsXYaux %>%
        filter(Brain_Area == 'ORB'|Brain_Area == 'MOs')
      
    } else if (area_ref == "RSP") {      
      RefCellsXYaux <- RefCellsXYaux %>%
        filter(Brain_Area== "RSP")
      
    } else if (area_ref == "all") {
      RefCellsXYaux <- RefCellsXYaux %>%
        filter(Brain_Area == 'ORB'|Brain_Area == 'MOs'|Brain_Area== "RSP")
      
    } 
    
    ###
    factorAnimal <- data.frame(animal = mouseN)
    factorLayer <- data.frame(level = layer)
    
    colocCellsXYaux <- cbind(factorAnimal,factorLayer,colocCellsXYaux)
    colocCellsXY <- rbind (colocCellsXY,colocCellsXYaux)
    
    RefCellsXYaux <- cbind(factorAnimal,factorLayer,RefCellsXYaux)
    RefCellsXY <- rbind (RefCellsXY,RefCellsXYaux)
  }
  
  
}

# Bar chart to observe the number of colocalised (pooled) and Ref cells
########################################################################
colocRef_long_pAn <-data.frame()
for(mouseN in animal_number){ 
  # start of for animal 
  
  index=1
  #coloc_number <- vector( "numeric" , 3 )
  #Ref_number <- vector( "numeric" , 3 )
  coloc_number <- vector( "numeric" , 1 )
  Ref_number <- vector( "numeric" , 1 )
  
  
  # colocCellsXY_level <- colocCellsXY %>%
  #   filter(level == layer)
  colocCellsXY_level <- colocCellsXY %>%
    filter(level == layer & animal == mouseN) # changed
  
  coloc_number[index] <- nrow(colocCellsXY_level)
  
  # RefCellsXY_level <- RefCellsXY %>%
  #   filter(level == layer)
  
  RefCellsXY_level <- RefCellsXY %>%
    filter(level == layer & animal == mouseN) # changed)
  Ref_number[index] <- nrow(RefCellsXY_level)-nrow(colocCellsXY_level) #this is the total number of cells in a section (double-counted cells accounted for)
  
  index=index+1
  
  
  ids <- unique(factor(colocCellsXY$level))
  values_bars <- data.frame(
    level = sort(unique(ids)),
    coloc_num = coloc_number,
    Ref_num = Ref_number
  )
  
  # Create data.frame with percentages, animals and levels
  
  factorAnimal <- data.frame(animal = mouseN)
  factorLayer <- data.frame(layer)
  percentage_vals <- coloc_number/Ref_number*100
  percentage_pAn_pLAux <- cbind(factorAnimal,factorLayer,percentage_vals)
  percentage_pAn_pL <- rbind (percentage_pAn_pL,percentage_pAn_pLAux)
  
  # convert values_bars to long format before plotting
  colocRef_long_pAnAux <- melt(values_bars, id = "level", variable_name = "cell_type")
  colocRef_long_pAnAux <- cbind(factorAnimal,colocRef_long_pAnAux)
  colocRef_long_pAn <- rbind (colocRef_long_pAn,colocRef_long_pAnAux)
  
  # plotting the bar chart of pooled colocalised versus Ref (total) cells per level
  p_bar <-ggplot(data=colocRef_long_pAnAux, aes( x = level, y = value, fill = cell_type)) +
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) + 
    geom_text(aes(label=value), vjust=-0.25, fontface = "bold", color="black", size=6,
              position=position_dodge(.9)) +
    ggtitle(paste0('Animal: ',mouseN, ' Colocalised and ',area_ref,' (reference) cell counts ')) +
    theme(plot.title = element_text(hjust = 0.85)) + # centering the plot title
    xlab("Atlas level (Dorsal --> Ventral)") + ylab("Number of neurons (pooled)") +
    theme_minimal() + 
    ylim(0,500) +
    #ylim(0,4500) +
    theme(axis.text.x = element_text(face="bold", color="black", 
                                     size=18),
          axis.text.y = element_text(face="bold", color="black", 
                                     size=18),
          plot.title = element_text(size=20, face= "bold", colour= "black" ),
          axis.title.x = element_text(size=20, face="bold", colour = "black"),    
          axis.title.y = element_text(size=20, face="bold", colour = "black")) +
    theme(legend.title = element_text(size=14, face="bold"),
          legend.text = element_text(size=12, face="bold"),
          legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
          legend.justification=c(0,1), legend.position=c(0,1))  # moving the legend inside the plot
  print(p_bar)
  # png(file=paste0('Plots/',mouseN,'bars_coloc_Ref',area_ref,'.png'),
  #     width=800, height=581)
  # plot(p_bar)
  # dev.off()
  
}
## end of for animal

#################################################
# Protocol to Assess Random Distribution of Cells
#################################################

# setting a random seed to obtain the same random split each time we run 
# the script for reproducibility 
histplotsLevels <- list()
zhistplotsLevels <- list()

pax_levels <- c(layer:layer)

shuffling_times <- 1000 # 200-fold iterative shuffling initially suggested
# start of animal 

for(mouseN in animal_number){ 
  ## Long format to wide format and allow only those sections where # of coloc cells is > 10%,
  ## add an if here when finished
  
  colocRef_long_pAnAux <- colocRef_long_pAn %>%
    filter(animal == mouseN) # changed
  
  colocRef_long_pAnAux <- subset(colocRef_long_pAnAux,select=-animal)
  
  colocRef_wide <- reshape(colocRef_long_pAnAux, idvar = "level", timevar = "cell_type", direction = "wide")
  
  percentages <- colocRef_wide$value.coloc_num/colocRef_wide$value.Ref_num*100
  #percentages[order(percentages$level,decreasing=FALSE),]
  set.seed(42)
  index_layer = 1
  
  if(percentages>=4 & !is.na(percentages)) {
    ## Sum Observed Pairwise Distances: 
    # Calculate the summed pairwise distances from the actual,
    # non-shuffled neuron locations (per level).
    
    colocCellsXY_level <- colocCellsXY %>%
      filter(level == layer & animal == mouseN) # changed
    eu_dist_act <- dist(colocCellsXY_level[,c("X","Y")]) 
    # summed pairwise distance of actual colocalised cells in pixels:
    sum_dist_act <- sum(eu_dist_act) 
    
    ## Shuffling Coordinates: 
    # Randomly shuffle the positions of the same number of neurons as in your 
    # double-labeled samples (colocalised cells), within the confines of the 
    # neurotrace cells in EC L5a (Ref cells). 
    # Repeat shuffling 200 times **(or as required)** to generate a robust distribution of random positions.
    # Calculate Pairwise Distances: For each shuffle, compute the pairwise Euclidean distances between all neurons.
    # Sum the Distances: For each shuffled dataset, calculate the summed pairwise distances.
    
    numColocL <- colocRef_long_pAn %>% # getting number of colocalised cells in that level
      filter(cell_type== "coloc_num")%>%
      filter(level== layer & animal == mouseN)
    
    RefCellsXY_level <- RefCellsXY %>% # Ref coordinated in that level
      filter(level == layer & animal == mouseN) # changed
    
    sum_dist <- vector( "numeric" , shuffling_times ) # initialise vector to store the summed distances
    for (ii in 1:shuffling_times){ # iterative shuffling
      
      # slice_sample() allows you to random select with or without replacement
      # https://web.ma.utexas.edu/users/parker/sampling/repl.htm
      # https://stats.stackexchange.com/questions/49480/sampling-with-or-without-replacement
      RefCellsXY_subset <- RefCellsXY_level %>% # sampling/selecting n= number of colocalised cells from Ref coordinates
        slice_sample(n = numColocL$value) 
      # calculating the Euclidean distance between rows of subset matrix   
      # https://www.geeksforgeeks.org/pairwise-distance-matrix-in-r/
      eu_dist <- dist(RefCellsXY_subset[,c("X","Y")]) 
      sum_dist[ii] <- sum(eu_dist) # in pixels
    }
    
    # Compute Z Scores: 
    # Determine how many standard deviations the observed summed distance is from 
    # the mean of the random distribution. 
    z_scores <- (sum_dist-mean(sum_dist))/sd(sum_dist)
    z_scores_act <- (sum_dist_act-mean(sum_dist))/sd(sum_dist)
    
    ## Plot histogram of pixel distances with mean and standard deviation:
    size_font <- 6
    histdata <-hist(sum_dist, breaks="FD")
    height_label <- max(histdata$counts) -max(histdata$counts)*5/100
    shuffle_index <- 1:shuffling_times
    sum_distDT <- data.frame(shuffle_index,sum_dist)
    p_simple <- ggplot(sum_distDT, aes(x=sum_dist)) + 
      #geom_histogram(color="black", fill="gray80",bins=round(shuffling_times/10)+10) + # add binwidth=1 argument
      geom_histogram(color="black", fill="gray80",binwidth=histdata$breaks[2]-histdata$breaks[1], boundary = histdata$breaks[2]-histdata$breaks[1]) + # add binwidth=1 argument
      geom_vline(aes(xintercept=mean(sum_dist)), color="red", linetype="solid", size=1) +
      geom_vline(aes(xintercept=mean(sum_dist)+sd(sum_dist)),color="black", linetype="dashed", size=1) +
      geom_vline(aes(xintercept=mean(sum_dist)-sd(sum_dist)),color="black", linetype="dashed", size=1) +
      geom_vline(aes(xintercept=sum_dist_act),color="blue", linetype="solid", size=1) +
      annotate("text", x = mean(sum_dist)+abs(mean(sum_dist))*1/100,y = height_label, label = "mean",colour = "red", fontface =2, size = size_font, angle = -90, hjust = 0.1)+
      annotate("text", x = mean(sum_dist)+sd(sum_dist)+abs(mean(sum_dist))*1/100,y = height_label, label = "SD",colour = "black", fontface =2, size = size_font, angle = -90, hjust = 0.1)+
      annotate("text", x = mean(sum_dist)-sd(sum_dist)+abs(mean(sum_dist))*1/100,y = height_label, label = "SD",colour = "black", fontface =2, size = size_font, angle = -90, hjust = 0.1)+
      annotate("text", x = sum_dist_act+abs(mean(sum_dist))*1/100,y = height_label, label = paste0("obs summed dist= ",round(mean(sum_dist),2)),colour = "blue", fontface =2, size = size_font, angle = -90, hjust = 0.1)+
      ggtitle(paste0("Level ", layer,": \ndistribution of randomly shuffled coordinates \n(n=",shuffling_times," shuffles; n=",numColocL$value," sampled cells, pooled; bin width = ",histdata$breaks[2]-histdata$breaks[1],")")) +
      xlab("summed pairwise distances\n(pixels)") + ylab("Frequency") +
      theme_minimal() + 
      ylim(0,max(histdata$counts)) +
      theme(axis.text.x = element_text(face="bold", color="black", 
                                       size=18),
            axis.text.y = element_text(face="bold", color="black", 
                                       size=18),
            plot.title = element_text(size=20, face= "bold", colour= "black" ),
            axis.title.x = element_text(size=20, face="bold", colour = "black"),    
            axis.title.y = element_text(size=20, face="bold", colour = "black")) 
    
    print(p_simple)
    
    
    # Plot with z-scores, mean and standard deviation:
    zhistdata <-hist(z_scores, breaks="FD")
    size_font <- 6
    height_label <- max(zhistdata$counts) -max(zhistdata$counts)*5/100
    # Sample plot but with ggplot
    shuffle_index <- 1:shuffling_times
    sum_distDTz <- data.frame(shuffle_index,z_scores)
    m <- cbind(factorAnimal,sum_distDTz)
    #z_scores_shufflingAll <- cbind(z_scores_shufflingAll,z_scores)
    factorAnimal <- data.frame(animal = mouseN)
    z_scores_Aux <- cbind(factorAnimal,z_scores)
    z_scores_pL <- rbind (z_scores_pL,z_scores_Aux)
    
    
    p_simpleZ <- ggplot(sum_distDTz, aes(x=z_scores)) + 
      #geom_histogram(color="black", fill="gray80",bins=round(shuffling_times/10)+10) + # add binwidth=1 argument
      geom_histogram(color="black", fill="gray80",binwidth=zhistdata$breaks[2]-zhistdata$breaks[1], boundary = zhistdata$breaks[2]-zhistdata$breaks[1]) + # add binwidth=1 argument
      geom_vline(aes(xintercept=mean(z_scores)),color="red", linetype="solid", size=1) +
      geom_vline(aes(xintercept=mean(z_scores)+sd(z_scores)),color="black", linetype="dashed", size=1) +
      geom_vline(aes(xintercept=mean(z_scores)-sd(z_scores)),color="black", linetype="dashed", size=1) +
      geom_vline(aes(xintercept=z_scores_act),color="blue", linetype="solid", size=1) +
      annotate("text", x = mean(z_scores)+0.15,y = height_label, label = "mean",colour = "red", fontface =2, size = size_font, angle = -90, hjust = 0.1)+
      annotate("text", x = mean(z_scores)+sd(z_scores)+0.15,y = height_label, label = "SD",colour = "black", fontface =2, size = size_font, angle = -90, hjust = 0.1)+
      annotate("text", x = mean(z_scores)-sd(z_scores)+0.15,y = height_label, label = "SD",colour = "black", fontface =2, size = size_font, angle = -90, hjust = 0.1)+
      annotate("text", x = z_scores_act+0.15,y = height_label, label = paste0("obs summed dist = ",round(z_scores_act,2)),colour = "blue", fontface =2, size = size_font, angle = -90, hjust = 0.1)+
      ggtitle(paste0("Level ", layer," animal ",mouseN,": \ndistribution of randomly shuffled coordinates \n(n=",shuffling_times," shuffles; n=",numColocL$value," sampled cells, pooled; bin width = ",zhistdata$breaks[2]-zhistdata$breaks[1],")")) +
      xlab("summed pairwise distances\n(z-score)") + ylab("Frequency") +
      theme_minimal() +  
      scale_x_continuous(breaks = seq(-6, 6, by = 1),limits = c(-6, 6)) +
      ylim(0,max(zhistdata$counts)) +
      theme(axis.text.x = element_text(face="bold", color="black", 
                                       size=18),
            axis.text.y = element_text(face="bold", color="black", 
                                       size=18),
            plot.title = element_text(size=20, face= "bold", colour= "black" ),
            axis.title.x = element_text(size=20, face="bold", colour = "black"),    
            axis.title.y = element_text(size=20, face="bold", colour = "black")) 
    
    print(p_simpleZ)
    # png(file=paste0('Plots/',mouseN,'_L',layer,'_bins_coloc.png'),
    #     width=900, height=581)
    # plot(p_simpleZ)
    # dev.off()
    print(paste0("Animal:",mouseN,"; perc=",percentages[index_layer],"--> Level: ",layer," --> sum_dist_act= ",round(sum_dist_act,2)," px & z_scores_act= ", round(z_scores_act,2))) 
    
    p_simpleZ <- ggplot(sum_distDTz, aes(x=z_scores)) + 
      #geom_histogram(color="black", fill="gray80",bins=round(shuffling_times/10)+10) + # add binwidth=1 argument
      geom_histogram(color="black", fill="gray80",binwidth=zhistdata$breaks[2]-zhistdata$breaks[1], boundary = zhistdata$breaks[2]-zhistdata$breaks[1]) + # add binwidth=1 argument
      geom_vline(aes(xintercept=mean(z_scores)),color="red", linetype="solid", size=1)
    
    
  }else{
    print(paste0("Animal:",mouseN,"; perc=",percentages[index_layer],"--> Level: ",layer,"; didn't make it  ")) 
  }
  
  # accumulate only if you have obtained a z-score
  if(percentages>=4 & !is.na(percentages)) {
    zScoresAll <- rbind(zScoresAll,c(mouseN, z_scores_act))
    actual_dist<- rbind(actual_dist, c(mouseN, sum_dist_act)) #creating data frame of the actual summed dist for coloc
    dist_all<- rbind(dist_all, data.frame(shuffled_distances = sum_dist, animal = rep(mouseN, length(sum_dist)))) #summed distances for each shuffle
  }
} # end of animal  

colnames(zScoresAll) <- c("animal", "z_scores")

#test <- zScoresAll[1:nrow(zScoresAll)-1,]
zScoresAll <- rbind(zScoresAll,c(00,  mean(zScoresAll$z_scores)))
zScoresAll$animal<-as.character(zScoresAll$animal)

zScoresAll[zScoresAll == "65912011207"] <- "Pooled"
#zScoresAll[zScoresAll == "16421643"] <- "Pooled"
zScoresAll[zScoresAll == "0"] <- "Mean"
#zScoresAll = zScoresAll[!zScoresAll == "660",] # REMOVING LINE WITH 660 BECAUSE IT IS JUST A COPY OF 659 AS IT DID NOT PASS THE 5% THRESHOLD

# Plotting z-scores all at the same time:
axisPlotsize = 16
tickPlotsize = 14
titlePlotsize = 18
legendPlotsize = 18
# Draw Plot with the nrow function
p_sumZscore <- ggplot(zScoresAll , aes (y = z_scores,x = 1)) + 
  geom_dotplot(binaxis='y', stackdir='center',binwidth = 1/20) 
print(p_sumZscore)
p_sumZscore <- ggplot(zScoresAll , aes (y = z_scores,x = 1,fill = factor(animal))) + 
  geom_dotplot(binaxis='y', stackdir='centerwhole',binwidth = 1/8,method="dotdensity",
               stackgroups = T,
               binpositions="all") + # binwidth modulates size of point
  #geom_point(aes(fill=factor(animal), colour=factor(animal)), pch= 16, size = 3) +
  ylim(-5,0) + 
  xlim(0.8,1.2) +
  theme(axis.title.x=element_blank(), #remove x axis labels
        axis.title.y=element_blank()) +  #remove y axis labels+
  theme_classic() +
  ggtitle(paste0('Level ', layer)) + ylab("z- scores") +
  theme(axis.text.x = element_text(face="bold", color="black", size=tickPlotsize),
        axis.text.y = element_text(face="bold", color="black",size=tickPlotsize),
        plot.title = element_text(size=titlePlotsize, face= "bold", colour= "black" ),
        axis.title.x = element_text(size=axisPlotsize, face="bold", colour = "black"),
        axis.title.y = element_text(size=axisPlotsize, face="bold", colour = "black")) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  theme(legend.title = element_text(size=legendPlotsize, face="bold"),
        legend.text = element_text(size=legendPlotsize, face="bold")) +
  theme(legend.title = element_text(size=legendPlotsize, face="bold"),
        legend.text = element_text(size=legendPlotsize, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.justification=c(1,1), legend.position=c(1,1),
        legend.margin = margin(c(0.5, 1.5, 0.5, 2.5),unit=c("mm","mm","mm","mm")),
        legend.spacing.x = unit(0.2, "mm"),
        legend.spacing.y = unit(0.1, "mm")) + # moving the legend inside the plot
  theme(legend.title = element_blank())
print(p_sumZscore)
# png(file=paste0('/Users/sara/Documents_2/R_scripts/level_4_analysis/figures/figures_shuffling/L',layer,'_z-score_coloc.png'),
#     width=581/2, height=900/2)
# plot(p_sumZscore)
# dev.off()


## HERE new plot with smoothed lines and summary
# Change density plot line colors by groups

# remove pooled data from show (65912001201)
z_scores_pL2 <- z_scores_pL %>%
  filter(animal != 65912011207)
# z_scores_pL2 <- z_scores_pL %>%
#   filter(animal != 16421643)


zScoresAll2 <- zScoresAll[c(1:3,5),]

p_all_pos <- ggplot(z_scores_pL2, aes(z_scores, colour = factor(animal) )) +
  geom_density(size=1.5) +
  geom_vline(data=zScoresAll, aes(xintercept=z_scores, color=animal),size=1.5,
             linetype="dashed",key_glyph = "rect") +
  theme_classic() +
  ggtitle(paste0('Shuffling protocol, atlas level ',layer))  +
  xlab("z-scores") + ylab("density") +
  scale_x_continuous(breaks = c(-4.:4),limits =c(-4,4)) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold", color="black", size=tickPlotsize),
        axis.text.y = element_text(face="bold", color="black",size=tickPlotsize),
        plot.title = element_text(size=titlePlotsize, face= "bold", colour= "black" ),
        axis.title.x = element_text(size=axisPlotsize, face="bold", colour = "black"),    
        axis.title.y = element_text(size=axisPlotsize, face="bold", colour = "black")) +
  theme(legend.title = element_text(size=legendPlotsize, face="bold"),
        legend.text = element_text(size=legendPlotsize, face="bold"),
        legend.background = element_rect(size=0.5, linetype="solid",colour ="black"),
        legend.justification=c(0.95,1), legend.position=c(0.95,1),
        legend.margin = margin(c(1.5, 1.5, 1.5, 2.5),unit=c("mm","mm","mm","mm")),
        legend.spacing.x = unit(0.5, "mm"),
        legend.spacing.y = unit(0.5, "mm")) + # moving the legend inside the plot
  theme(legend.title = element_blank()) #+
  # theme(axis.text.x=element_blank(), 
  #       axis.title.x =element_blank(), 
  #       #axis.ticks.x=element_blank(), 
  #       axis.text.y=element_blank(), 
  #       axis.title.y=element_blank(), 
  #       #axis.ticks.y=element_blank()
  # ) 
print(p_all_pos)
png(file=paste0('Plots/L',layer,'_shuffling_zscores.png'),
     width=863, height=542)
plot(p_all_pos)
dev.off()















