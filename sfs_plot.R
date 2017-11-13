require(ggplot2)

setwd("~/Documents/PhD/dutch_carriage/")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

sfs_plot_input <- read.delim("~/Documents/PhD/dutch_carriage/sfs_plot_input.txt", 
                             header=FALSE)
colnames(sfs_plot_input) = c("MAF", "Consequence", "Niche")
sfs_plot_input$Consequence <- factor(sfs_plot_input$Consequence,
                                     levels=c("LoF", "Intergenic", "Missense", "Synonymous"))

ggplot(sfs_plot_input) + geom_histogram(aes(x=MAF, (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..],
                                            fill=Consequence), binwidth = 0.05,position='dodge') +
  facet_grid(. ~ Niche, drop=T) + 
  theme_bw(base_size = 14) +
  ylab("Normalised frequency") + 
  scale_fill_manual(values=cbPalette)