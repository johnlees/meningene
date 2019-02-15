pnas_1613937114_sd01_1_ <- read_excel("Downloads/pnas.1613937114.sd01 (1).xlsx",skip = 1)
View(pnas_1613937114_sd01_1_)
core_omega = as.numeric(pnas_1613937114_sd01_1_[pnas_1613937114_sd01_1_$`Sequences in filtered alignment`==616,13]$ω)
core_omega = core_omega[!is.na(core_omega) & core_omega < 900]
hist(core_omega, breaks = c(seq(0,1,0.1),10))
hits_omega = c(0.2770, 0.1485, 0.2788, 0.1780, 0.2788, 0.6085)
wilcox.test(core_omega, hits_omega)
median(core_omega)
median(hits_omega)
