library(nonlinearTseries)
library(TDAstats)
library(tidyverse)
library(TDA)
library(seewave)

set.seed(1)

##########################################
# music data
##########################################

# flute
flute_pd <- list()
Flute_full_note <- read.csv("Flute A full time series.csv")
for (i in 1:60) {
  # sample time series
  ai <- sample(10000:37999,1,replace = FALSE, prob = NULL)
  bi <- ai + 1000
  sample_sequence_ai <- Flute_full_note[ai:bi,]
  
  # Takens embedding
  flute_A4_matrix <- data.matrix(sample_sequence_ai)
  tak <- buildTakens(flute_A4_matrix,2,3)
  
  # compute PH
  flute_pd[[i]] <- calculate_homology(tak,return_df = TRUE) 
}

# clarinet
clarinet_pd <- list()
Clarinet_full_note <- read.csv("Clarinet A full time series.csv")
for (i in 1:60) {
  ai2 <- sample(10000:95000,1,replace = FALSE, prob = NULL)
  bi2 <- ai2 + 1000
  sample_sequence_ai2 <- Clarinet_full_note[ai2:bi2,]
  
  # Takens embedding
  Clarinet_A4_matrix <- data.matrix(sample_sequence_ai2)
  tak2 <- buildTakens(Clarinet_A4_matrix,2,3)
  
  clarinet_pd[[i]] <- calculate_homology(tak2,return_df = TRUE) 
}

#extracting one-dimensional features
flute_onedim_list = lapply(flute_pd,function(df){
  df[df$dimension ==1, c('birth','death')]
})
clarinet_onedim_list = lapply(clarinet_pd,function(df){
  df[df$dimension ==1, c('birth','death')]
})

for (i in seq_along(flute_onedim_list)) {
  write.csv(flute_onedim_list[[i]], paste0("flute_diagrams/flute_onedim_", i, ".csv"), row.names = FALSE)
}
for (i in seq_along(clarinet_onedim_list)) {
  write.csv(clarinet_onedim_list[[i]], paste0("clarinet_diagrams/clarinet_onedim_", i, ".csv"), row.names = FALSE)
}

##########################################
# plots
##########################################
library(ggplot2)
library(patchwork)
set.seed(12)

# # sample time series
# ai <- sample(10000:37999,1,replace = FALSE, prob = NULL)
# bi <- ai + 1000
# sample_sequence_ai <- Flute_full_note[ai:bi,]
# sample_sequence_ai2 <- Clarinet_full_note[ai:bi,]
# 
# 
# npoint=1001
# plotdata=data.frame(times = 50*rep(ai:(ai+npoint-1)/1000000,2),
#                     y=c(sample_sequence_ai2[1:npoint],sample_sequence_ai[1:npoint]),
#                     Instrument=c(rep("Clarinet",npoint),rep("Flute",npoint)))
# figure1 = ggplot(plotdata, aes(x=times,y=y,color=Instrument,linetype=Instrument)) +
#   geom_line()+
#   labs(x="Second",y="")+
#   scale_color_manual(values = c("blue","red"), 
#                      labels = c( "Clarinet", "Flute"))+
#   scale_linetype_manual(values = c(1,5), 
#                         labels = c("Clarinet", "Flute"))+
#   theme_minimal()+
#   theme(panel.grid.minor = element_blank(),
#                          panel.background = element_blank(),
#                          legend.position = "top", 
#                          legend.background=element_rect(colour = "white"),
#                          legend.title = element_text(size=15),
#                          legend.text = element_text(size=15))
# ggsave("../figures/sound_waves.png", plot = figure1, dpi = 300,bg = "white" )
# 
# # flute
# x<- data.matrix(sample_sequence_ai)
# a<- buildTakens(x,2,3)
# figure2=ggplot(data.frame(a,col="a"))+
#   geom_point(aes(x=X1,y=X2,color=col))+
#   labs(x=expression(Z[t]),y=expression(Z[t+3]))+
#   scale_color_manual(values = c("red"))+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   ylim(-0.4,0.4)+
#   xlim(-0.4,0.4)
# ggsave("../figures/flute_cloud.png", plot = figure2, dpi = 300,bg = "white" )
# # clarinet
# x<- data.matrix(sample_sequence_ai2)
# a<- buildTakens(x,2,3)
# figure3=ggplot(data.frame(a,col="a"))+
#   geom_point(aes(x=X1,y=X2,color=col))+
#   labs(x=expression(Z[t]),y=expression(Z[t+3]))+
#   scale_color_manual(values = c("blue"))+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   ylim(-0.4,0.4)+
#   xlim(-0.4,0.4)
# 
# ggsave("../figures/clarinet_cloud.png", plot = figure3, dpi = 300,bg = "white" )
# 


set.seed(12)

# sample time series
ai <- sample(10000:37999,1,replace = FALSE, prob = NULL)
bi <- ai + 1000
sample_sequence_ai <- Flute_full_note[ai:bi,]
sample_sequence_ai2 <- Clarinet_full_note[ai:bi,]

npoint=1001
plotdata=data.frame(
  times = 50*rep(ai:(ai+npoint-1)/1000000,2),
  y=c(sample_sequence_ai2[1:npoint], sample_sequence_ai[1:npoint]),
  Instrument=c(rep("Clarinet",npoint), rep("Flute",npoint))
)

### ---- figure1 ----
figure1 = ggplot(plotdata, aes(x=times,y=y,color=Instrument,linetype=Instrument)) +
  geom_line()+
  labs(x="Second",y="")+
  scale_color_manual(values = c("#41b6c4","#253494"), 
                     labels = c("Clarinet", "Flute"))+
  scale_linetype_manual(values = c(1,5), 
                        labels = c("Clarinet", "Flute"))+
  scale_linewidth_manual(values = c(1.4, 1.2))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "top", 
        legend.background=element_rect(colour = "white"),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))+
 theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank()
)

### ---- figure2 ----
x<- data.matrix(sample_sequence_ai)
a<- buildTakens(x,2,3)
figure2=ggplot(data.frame(a,col="a"))+
  geom_point(aes(x=X1,y=X2,color=col))+
  labs(x=expression(Z[t]),y=expression(Z[t+3]))+
  scale_color_manual(values = c("#253494"))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  ylim(-0.4,0.4)+
  xlim(-0.4,0.4)+
  theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank()
)

### ---- figure3 ----
x<- data.matrix(sample_sequence_ai2)
a<- buildTakens(x,2,3)
figure3=ggplot(data.frame(a,col="a"))+
  geom_point(aes(x=X1,y=X2,color=col))+
  labs(x=expression(Z[t]),y=expression(Z[t+3]))+
  scale_color_manual(values = c("#41b6c4"))+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        legend.position = "none")+
  ylim(-0.4,0.4)+
  xlim(-0.4,0.4)+
 theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank()
)


tag_theme <- theme(
  plot.tag = element_text(
    size = 20,
    face = "bold",
    colour = "black"
  ),
  plot.tag.position = c(0.01, 0.99),
  plot.margin = margin(t = 12, r = 8, b = 8, l = 12)
)

figure1 <- figure1 + labs(tag = "A") + tag_theme
figure2 <- figure2 + labs(tag = "B") + tag_theme
figure3 <- figure3 + labs(tag = "C") + tag_theme

combined_plot <- figure1 / (figure2 | figure3) +
  plot_layout(heights = c(1, 1))

### ---- SAVE ----
ggsave("../results/figures/figure_file/sound_wave_clouds.pdf",
       plot = combined_plot,
       dpi = 300,
       width = 8,
       height = 8,
       bg = "white")
