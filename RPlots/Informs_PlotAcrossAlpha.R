#Plots across Alpha for a single K Run

library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
loadfonts()

###The old data
dat = read_csv("../Results/InformsPlots/singleKAcrossAlpha_10000_.95.csv")
# dat = read_csv("../Results/InformsPlots/singleKAcrossAlphaBad_999_.95.csv")


SAA0_val <- filter(dat, alpha < 1e-6) %>% select(SAAEstim)
dat <- dat %>% mutate(SAA0 = SAA0_val[[1]], 
               SAASubOpt = SAAEstim - SAA0, 
               Stability = LOOEstim_Sc - SAAEstim)

g<- ggplot(dat, aes(alpha, OR),) + 
  geom_point(color="#00AFBD") + geom_line(color="#00AFBD") + 
  theme_bw(base_size=18) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .9)) + 
  ylab("Cost") + xlab("")

# add a line for oracle value
alphaOR = max(filter(dat, OR == min(OR)) %>% select(alpha))

g<- g + geom_vline(xintercept = alphaOR, linetype="dotted")

#ggsave("./Figures/informs_Newsvendor_AcrossAlpha_10000_.95.pdf", height=5.34, width=9, units="in")
#ggsave("./Figures/informs_Newsvendor_AcrossAlphaBad_999_.95.pdf", height=5.34, width=9, units="in")

## Now make a plot with the LOO estimate on it too.  
alphaLOO = max(filter(dat, LOOEstim_Sc == min(LOOEstim_Sc)) %>% select(alpha))
g2 <- g + geom_line(aes(alpha, LOOEstim_Sc), color="#FBBEBA")
  
#ggsave("./Figures/informs_Newsvendor_AcrossAlphaBoth_10000_.95.pdf",height=5.34, width=9, units="in")

########
#  Plot the breakdown by stability

g3 <- dat %>% ggplot(aes(alpha, SAASubOpt + Stability )) + 
  geom_line(color="#00AFBD") + 
  geom_line(aes(alpha, SAASubOpt), linetype="dotted") + 
  geom_line(aes(alpha, Stability), linetype="dotted") + 
  theme_bw(base_size=18) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .9)) + 
  ylab("") + xlab("")

#ggsave("./Figures/informs_Newsvendor_Stability_10000_.95.pdf", height=5.34, width=9, units="in")
ggsave("./Figures/informs_Newsvendor_StabilityBad_999_.95.pdf", height=5.34, width=9, units="in")
