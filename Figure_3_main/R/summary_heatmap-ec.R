################################
## Heat plot to summarize HBMAP results
################################

library(ComplexHeatmap)
library(circlize)

## load results
load("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/BrainConnectivity/Results/ECdata_results_with_analysis.RData")

## Replace with results of the data!!!
mcmc_unique = mcmc_unique_EC8 #MCMC results for motif parameters given optimal clustering
mcmc_omega = omega_JM_mcmc #MCMC results for weights given optimal clustering

# Matrix of expected projection strength across clusters (regions x clusters)
df = data.frame(matrix(mcmc_unique$estimated.projection.df$projection.mean, R, mcmc_unique$J),row.names = mcmc_unique$estimated.projection.df$region[1:R])
names(df) = seq(1:mcmc_unique$J)
# Matrix of expected weight across clusters for the global and local weights for each mouse  (clusters x M+1)
weightmat = matrix(cbind(mcmc_omega$mean_omega_J, mcmc_omega$mean_omega_J_M), nrow = mcmc_unique$J, ncol = M +1)
# Matrix of expected projection strength for each mouse (regions x M)
psmat = matrix(cbind(ps_m1$eps,ps_m2$eps,ps_m3$eps,ps_m4$eps, ps_m5$eps, ps_m6$eps), nrow = R, ncol = M)
# vector of probability that the global weight is greater than a threshold (length = num of clusters)
pg = prob_greater_global
## variable motifs relative to the null
#locvar.prob = local_weights_analysis_ec$probability
#names(df) = sapply(1:mcmc_unique$J, function(j){
#  if(locvar.prob[j]>0.95){ cname = paste0(j,'*')}
#  else{cname = paste(j)}
#  cname
#    })

# Colorscale for heatmap
col_fun = colorRamp2(c(0, 1), c( "white", "blue"))

# Top annotation: expected global and local weights for each mouse
# Note: if more mice, add another color. Red indicates global weights
colnames(weightmat) = c("global", paste("mouse", c(1:M)))
mouse_colors = topo.colors(M)
column_ha  = HeatmapAnnotation(w = anno_barplot(weightmat, 
                                        beside = TRUE, attach = TRUE, 
                                        barwidth = 0.99,
                                        gp = gpar(fill = c("red",mouse_colors)),
                                        height = unit(6, "cm"))
)
# If you only want the global weights in the top annotation
#column_ha = HeatmapAnnotation(global = anno_barplot(mcmc_omega$mean_omega_J))

# Left annotation: expected projection strength for each mouse 
row_ha  = rowAnnotation(p = anno_barplot(psmat, 
                                         beside = TRUE, attach = TRUE, 
                                         barwidth = 0.99,
                                         gp = gpar(fill = mouse_colors),
                                         border = FALSE,
                                         width = unit(2, "cm"))
)


# Heat map of projection strength for each cluster with annotation
Heatmap(df,
        cluster_rows = FALSE, 
        row_names_side = "left",
        cluster_columns = FALSE, 
        show_column_names = TRUE,
        heatmap_legend_param = list(title = ""),
        col = col_fun,
        column_order = order(pg, decreasing = TRUE),
        column_split = as.factor(2*(pg<0.95&weightmat[,1]<0.01) + (pg<0.95&weightmat[,1]>0.01)),
        column_title_gp = gpar(fontsize=0),
        border = TRUE,
        top_annotation = column_ha,
        left_annotation = row_ha
)
lgd = Legend(labels = c("Global",paste("Mouse", c(1:M))), title = "", 
         legend_gp = gpar(fill = c("red",mouse_colors)),
         border = "black")
draw(lgd,x = unit(0.9, "npc"), y = unit(0.95, "npc"),just = c("right", "top"))

