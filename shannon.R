require(vegan)
require(ggplot2)

setwd("~/Documents/PhD/dutch_carriage/ivr")
num_samples = 1729

means <- readRDS("means.Rdata")

div_mat <- matrix(0, nrow = num_samples, ncol = 6)
for (i in 1:6)
{
  pattern = paste0("pi\\[[0-9]+,",i,"\\]")
  div_mat[,i] = means[grep(pattern, names(means))]
}
div_mat = as.data.frame(div_mat)

shannon = diversity(div_mat)
shannon_div = data.frame(Shannon=shannon, Niche=c(rep("Invasive", 1052), 
                                                  rep("Carriage", 677)))

ggplot(shannon_div) + geom_histogram(aes(x=Shannon, fill=factor(Niche)), 
                                     position = "identity", alpha=0.4, binwidth = 0.15) + 
  xlab("Shannon diversity") + 
  theme_bw(base_size = 16, base_family = "Helvetica") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_manual(values=cbbPalette,name="Niche")

median(shannon_div[shannon_div$Niche=="Carriage",1])
median(shannon_div[shannon_div$Niche=="Invasive",1])

wilcox.test(shannon_div[shannon_div$Niche=="Carriage",1], shannon_div[shannon_div$Niche=="Invasive",1])
