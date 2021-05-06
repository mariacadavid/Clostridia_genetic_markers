#Compare gene distances (dnak1 vs. gyrb vs. 16s)
#Table of pairwise distances per pair of genes needed 
#By: Maria Cadavid Feb2021


#libraries
library(tidyr)#drop_na


#16s vs dnak1
distances_16s_dnka1 <- read.csv(file='/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_16s_dnak1.csv', sep = "\t", header = TRUE)
common_pairwise_16s_dnak1 <- drop_na(distances_16s_dnka1)
common_pairwise_16s_dnak1$simil_16s <- ((1 - common_pairwise_16s_dnak1$dist_16s)*100 )
common_pairwise_16s_dnak1$simil_dnak1 <- ((1 - common_pairwise_16s_dnak1$dist_dnak1)*100)
write.table(common_pairwise_16s_dnak1, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

#16s vs gyrb
distances_16s_gyrb <- read.csv(file='/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_16s_gyrb.csv', sep = "\t", header = TRUE)
common_pairwise_16s_gyrb <- drop_na(distances_16s_gyrb)
common_pairwise_16s_gyrb$simil_16s <- ((1 - common_pairwise_16s_gyrb$dist_16s)*100 )
common_pairwise_16s_gyrb$simil_gyrb <- ((1 - common_pairwise_16s_gyrb$dist_gyrb)*100)
write.table(common_pairwise_16s_gyrb, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)

#dnak1 vs gyrb
distances_dnak1_gyrb <- read.csv(file='/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/distances_dnak1_gyrb.csv', sep = "\t", header = TRUE)
common_pairwise_dnak1_gyrb <- drop_na(distances_dnak1_gyrb)
common_pairwise_dnak1_gyrb$simil_dnak1 <- ((1 - common_pairwise_dnak1_gyrb$dist_dnak1)*100 )
common_pairwise_dnak1_gyrb$simil_gyrb <- ((1 - common_pairwise_dnak1_gyrb$dist_gyrb)*100)
write.table(common_pairwise_dnak1_gyrb, file = "/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common.tsv", sep = "\t",
            row.names = FALSE, col.names = TRUE)



###Plot color by density: https://slowkow.com/notes/ggplot2-color-by-density/
library(MASS)
library(ggplot2)
library(viridis)
theme_set(theme_bw(base_size = 16))

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#Add r^2 to graph
library(ggpubr)
library(ggpmisc)



#16s vs dnak1
common_pairwise_16s_dnak1$density <- get_density(common_pairwise_16s_dnak1$dist_16s, common_pairwise_16s_dnak1$dist_dnak1, n = 100)
#distance plot 16s vs dnak1
plot_16s_dnak1_dist <- ggplot(common_pairwise_16s_dnak1, aes(x=dist_16s, y=dist_dnak1)) + 
  geom_point(aes(dist_16s, dist_dnak1, color = density), size= 0.5) + 
  scale_color_viridis() +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "8")) +
  theme(legend.position = "right", legend.text = element_text(size = 8, colour = "black")) +
  labs(title= "Pairwise distances 16s vs. dnak1") +
  xlab ("16s distance") + ylab ("dnak1 distance") +
  xlim(0,1) + ylim (0,1) +
  geom_smooth(method=lm, formula= y~x, color="red", size=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="black", size= 0.5)
plot_16s_dnak1_dist
#similarity plot 16s vs dnak1
plot_16s_dnak1_similiarity <- ggplot(common_pairwise_16s_dnak1, aes(x=simil_16s, y=simil_dnak1)) + 
  geom_point(aes(simil_16s, simil_dnak1, color = density), size= 0.5) + 
  scale_color_viridis() +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity 16s vs. dnak1") +
  xlab ("16s identity (%)") + ylab ("dnak1 identity (%)") +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity


#16s vs gyrb
common_pairwise_16s_gyrb$density <- get_density(common_pairwise_16s_gyrb$dist_16s, common_pairwise_16s_gyrb$dist_gyrb, n = 100)
#distance plot 16s vs gyrb
plot_16s_gyrb_dist <- ggplot(common_pairwise_16s_gyrb, aes(x=dist_16s, y=dist_gyrb)) + 
  geom_point(aes(dist_16s, dist_gyrb, color = density), size= 0.5) + 
  scale_color_viridis() +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise distances 16s vs. gyrb") +
  xlab ("16s distance") + ylab ("gyrb distance") +
  xlim(0,1) + ylim (0,1) +
  geom_smooth(method=lm, formula= y~x, color="red", size=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="black", size= 0.5)
plot_16s_gyrb_dist
#similarity plot 16s vs gyrB
plot_16s_gyrb_similiarity <- ggplot(common_pairwise_16s_gyrb, aes(x=simil_16s, y=simil_gyrb)) + 
  geom_point(aes(simil_16s, simil_gyrb, color = density), size= 0.5) + 
  scale_color_viridis() +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity 16s vs. gyrB") +
  xlab ("16s identity (%)") + ylab ("gyrB identity (%)") +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity

#dnak1 vs gyrb
common_pairwise_dnak1_gyrb$density <- get_density(common_pairwise_dnak1_gyrb$dist_dnak1, common_pairwise_dnak1_gyrb$dist_gyrb, n = 100)
#distance plot dnak1 vs gyrb
plot_dnak1_gyrb_dist <- ggplot(common_pairwise_dnak1_gyrb, aes(x=dist_dnak1, y=dist_gyrb)) + 
  geom_point(aes(dist_dnak1, dist_gyrb, color = density), size= 0.5) + 
  scale_color_viridis() +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise distances dnak1 vs. gyrb") +
  xlab ("dnak1 distance") + ylab ("gyrb distance") +
  xlim(0,1) + ylim (0,1) +
  geom_smooth(method=lm, formula= y~x, color="red", size=0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="black", size= 0.5)
plot_dnak1_gyrb_dist
#similarity plot dnak1 vs gyrB
plot_dnak1_gyrb_similiarity <- ggplot(common_pairwise_dnak1_gyrb, aes(x=simil_dnak1, y=simil_gyrb)) + 
  geom_point(aes(simil_dnak1, simil_gyrb, color = density), size= 0.5) + 
  scale_color_viridis() +
  labs(colour = "Number of pairwise \ncomparisons") +
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity dnaK1 vs. gyrB") +
  xlab ("dnaK1 identity (%)") + ylab ("gyrB identity (%)") +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity


#Same plots colored by taxonomy for level analysis

#Open tables with taxonomy analysis 
tax_analysis_16s_dnak1 <- read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_dnak1_common_taxonomy_analysis.tsv", sep = "\t", header = TRUE)
tax_analysis_16s_gyrb <- read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_16s_gyrb_common_taxonomy_analysis.tsv", sep = "\t", header = TRUE)
tax_analysis_dnak1_gyrb <- read.csv(file="/Users/mariacadavid/Desktop/Tesis_Clostridiales/Data/Compared_distances/identity_dnak1_gyrb_common_taxonomy_analysis.tsv", sep = "\t", header = TRUE)


#similarity plot 16s vs dnak1 color taxonomy analysis
plot_16s_dnak1_similiarity_tax <- ggplot(tax_analysis_16s_dnak1, aes(x=identity_16s, y=identity_dnak1) + 
  geom_point(aes(identity_16s, identity_dnak1, color = same_specie), size= 0.5) +  #Here change tax level as desired
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity 16s vs. dnaK1") +
  xlab ("16s identity (%)") + ylab ("dnaK1 identity (%)") +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_dnak1_similiarity_tax

#similarity plot 16s vs gyrB color taxonomy analysis
plot_16s_gyrb_similiarity_tax <- ggplot(tax_analysis_16s_gyrb, aes(x=identity_16s, y=identity_gyrb)) + 
  geom_point(aes(identity_16s, identity_gyrb, color = same_specie), size= 0.5) +  #Here change tax level as desired
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity 16s vs. gyrB") +
  xlab ("16s identity (%)") + ylab ("gyrB identity (%)") +
  xlim(65,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_16s_gyrb_similiarity_tax


#similarity plot dnak1 vs gyrB color taxonomy analysis
plot_dnak1_gyrb_similiarity_tax <- ggplot(tax_analysis_dnak1_gyrb, aes(x=identity_dnak1, y=identity_gyrb)) + 
  geom_point(aes(identity_dnak1, identity_gyrb, color = same_specie), size= 0.5) + #Here change tax level as desired
  theme(legend.title = element_text(size = "10")) +
  theme(legend.position = "right", legend.text = element_text(size = 10, colour = "black")) +
  labs(title= "Pairwise identity dnak1 vs. gyrB") +
  xlab ("dnak1 identity (%)") + ylab ("gyrB identity (%)") +
  xlim(0,100) + ylim (0,100) +
  geom_smooth(method=lm, formula= y~x, color="black", size=0.8) +
  stat_regline_equation(label.y = 10, aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 5, aes(label = ..rr.label..)) +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="red", size= 0.5)
plot_dnak1_gyrb_similiarity_tax



