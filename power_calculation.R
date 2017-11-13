# Should be able to get a ef <-> OR formula
# Should be able to do 3x2 properly

require(pwr)
require(ggplot2)

get_w = function(maf1, maf2, ef)
{
  p0 <- matrix(c(maf1*maf2, maf1*(1-maf2), maf2*(1-maf1), (1-maf1)*(1-maf2)), nrow = 2)
  p1 <- matrix(c(maf1*maf2*ef, maf1*(1-maf2*ef), maf*(1-maf), (1-maf)^2), nrow = 2)
  return(ES.w1(p0, p1))
}

power_test = function(N, w)
{
  return(pwr.chisq.test(w=w, N = N, sig.level = 1e-11, df = 2)$power)
}

# Test
maf1 = 0.25
maf2 = 0.25
p0 <- matrix(c(maf1*maf2, maf1*(1-maf2), maf2*(1-maf1), (1-maf1)*(1-maf2)), nrow = 2)

ef = 2.27
p1 <- matrix(c(maf1*maf2*ef, maf1*(1-maf2*ef), maf*(1-maf), (1-maf)^2), nrow = 2)
or = (p1[1,1] / p1[1,2]) / (p1[2,1]/p1[2,2])
#or

w = ES.w1(p0, p1)
pwr.chisq.test(w=w, N = 460, sig.level = 1e-11, df = 2)

pwr_3 = unlist(lapply(seq(10,2000,10), power_test, w=get_w(0.25,0.25,2)))
pwr_4 = unlist(lapply(seq(10,2000,10), power_test, w=get_w(0.25,0.25,2.28)))
pwr_5 = unlist(lapply(seq(10,2000,10), power_test, w=get_w(0.25,0.25,2.5)))

power_plot = as.data.frame(cbind(rep(seq(10,2000,10), 3), c(pwr_3, pwr_4, pwr_5), c(rep(3,200), rep(4,200), rep(5,200))))
colnames(power_plot) = c("Samples", "Power", "OR")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(power_plot) + geom_line(aes(x=Samples, y=Power, colour=factor(OR)), size=1.5) +
  geom_segment(aes(x = 460, xend = 460, y = 0, yend = 1), linetype=2, size=0.1, show.legend = F) +
  scale_colour_manual(values=cbPalette, name="OR") +
  theme_bw(base_size = 14) +
  xlim(0,1100)
