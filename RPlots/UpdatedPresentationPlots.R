### Updated Presentation Plots Sept 2019
#Plots across Alpha for a single K Run
library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
loadfonts()

dat = read_csv("../Results/single_path_curves_presentation.csv")
dat = read_csv("../Results/single_path_curves_bad_presentation.csv")
dat.melt <- gather(dat, key=Type, value=value, -alpha)

#teal #00AFBD
#pink #FBBEBA

(
  g<- ggplot(dat, aes(alpha, OR),) + 
  geom_point(color="#00AFBD") + geom_line(color="#00AFBD") + 
  theme_bw(base_size=18) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .9)) + 
  ylab("Cost") + xlab("") + 
  ylim(31.5, 33.1))
)

#save down a "clean one.
ggsave("../RPlots/Figures/updated_Good_OR_curve.pdf", 
       g, height=5.34, width=9, units="in")


# add a line for oracle value
alphaOR = max(filter(dat, OR == min(OR)) %>% select(alpha))
g<- g + geom_vline(xintercept = alphaOR, linetype="dotted")
ggsave("../RPlots/Figures/updated_Good_OR_curve2.pdf", 
       g, height=5.34, width=9, units="in")

#ggsave("./Figures/informs_Newsvendor_AcrossAlpha_10000_.95.pdf", height=5.34, width=9, units="in")
#ggsave("./Figures/informs_Newsvendor_AcrossAlphaBad_999_.95.pdf", height=5.34, width=9, units="in")

## Now make a plot with the LOO estimate on it too.  
alphaLOO = max(filter(dat, LOO == min(LOO)) %>% select(alpha))
g2 <- g + geom_line(aes(alpha, LOO), color="#FBBEBA") + 
          geom_point(aes(alpha, LOO), color="#FBBEBA", shape="square")

#no line for the LOO minimum because same as OR
#Save it down
ggsave("../RPlots/Figures/updated_Good_ORandLOO_curve.pdf",
       g2, height=5.34, width=9, units="in")


########
#  Plot the breakdown by stability
###
#First one that's just SAA Subopt
g3 <- ggplot(dat, aes(alpha, SAASubopt)) + 
  geom_line(linetype="dashed") + 
  theme_bw(base_size=18) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .9)) + 
  ylab("Cost") + xlab("") + ylim(0, 8)
g3

ggsave("../RPlots/Figures/updated_Good_SAASubopt.pdf",
       g3, height=5.34, width=9, units="in")  

  
g3 <- g3 +   geom_line(aes(alpha, Instab), linetype="dashed") 
ggsave("../RPlots/Figures/updated_Good_SAASubopt_and_Instab.pdf",
       g3, height=5.34, width=9, units="in")

g3 <- g3 +
  geom_line(aes(alpha, SAASubopt + Instab), color="#FBBEBA") + 
  geom_point(aes(alpha, SAASubopt + Instab), color="#FBBEBA", shape="square")

ggsave("../RPlots/Figures/updated_Good_SAASubopt_and_Instab_and_LOO.pdf",
       g3, height=5.34, width=9, units="in")


##### Instability Plots for bad example
########
#  Plot the breakdown by stability
###
#First one that's just SAA Subopt
g3 <- ggplot(dat, aes(alpha, SAASubopt)) + 
  geom_line(linetype="dashed") + 
  theme_bw(base_size=18) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .9)) + 
  ylab("Cost") + xlab("") + ylim(0, .42)
g3

ggsave("../RPlots/Figures/updated_Bad_SAASubopt.pdf",
       g3, height=5.34, width=9, units="in")  

g3 <- g3 +   geom_line(aes(alpha, Instab), linetype="dashed") 
ggsave("../RPlots/Figures/updated_Bad_SAASubopt_and_Instab.pdf",
       g3, height=5.34, width=9, units="in")

g3 <- g3 +
  geom_line(aes(alpha, SAASubopt + Instab), color="#FBBEBA") + 
  geom_point(aes(alpha, SAASubopt + Instab), color="#FBBEBA", shape="square")

ggsave("../RPlots/Figures/updated_Bad_SAASubopt_and_Instab_and_LOO.pdf",
       g3, height=5.34, width=9, units="in")

