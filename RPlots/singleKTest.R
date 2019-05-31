### Creates a series of plots for INFORMS 2018 Presentation

library(tidyverse)
library(ggplot2)
library(forcats)
library(extrafont)
loadfonts()

#pink, blue, green
cbPallette = c("#FBBEBA", "#B5CEFF", "#9ADDA5")

#first is teal, then pink, 
cbPallette2 = c("#00AFBD", "#FBBEBA")


###
#First create a plot of 1 newsvendor SAA
#####
dat = read_csv("../Results/singleKUnifNewsvendor_1_1000.csv")

##add pretty labels
dat <- dat %>% mutate(Method = factor(Method), 
              Label = fct_recode(Method, Pooled="LOO_avg"))

dat.sum<- dat %>% group_by(Label) %>%
  summarise(avg = mean(TruePerf), 
            std = sd(TruePerf), 
            avgAlpha0 = mean(alpha)) 

g <- dat %>% filter(Method %in% c("SAA")) %>%
ggplot(aes(x=TruePerf, group=Label, fill=Label)) + 
  geom_histogram(binwidth=.25, size=0) 
  
#add lines for means
g <- g + 
  geom_vline(data = filter(dat.sum, Label %in% c("FullInfo")),  
             aes(xintercept=avg), linetype="dashed" )
g <- g + theme_bw() + 
  theme(legend.title=element_blank(), 
        legend.position=c(.8, .9)) + 
#  ylim(0, 90) + 
  ylab("") + xlab("Cost") + 
  scale_fill_manual(values = cbPallette2)

ggsave("./Figures/informs_singleNewsvendor_1.pdf", height=3.8, width=5.3, units="in")


g <- g + 
  geom_vline(data = filter(dat.sum, Label %in% c("SAA")),  
             aes(xintercept=avg), linetype="dotted" )

ggsave("./Figures/informs_singleNewsvendor_1_withSAALine.pdf", height=3.8, width=5.3, units="in")


####
# Similar graph but for many K's
####
#dat = read_csv("../Results/singleKUnifNewsvendor_1000_1000.csv")
dat = read_csv("../Results/InformsPlots/singleKUnifNewsvendor_1000_1000.csv")

##add pretty labels
dat <- dat %>% mutate(Method = factor(Method), 
                      Label = fct_recode(Method, Pooled="LOO_avg"), 
                      Label = fct_relevel(Label, "FullInfo", "SAA", "Oracle", "Pooled", "LOO_unif"))

dat.sum<- dat %>% group_by(Label) %>%
  summarise(avg = mean(TruePerf), 
            std = sd(TruePerf), 
            avgAlpha0 = mean(alpha))

###For just SAA
g <- dat %>% filter(Method %in% c("SAA")) %>%
  ggplot(aes(x=TruePerf, group=Label, fill=Label)) + 
  geom_histogram(size=0, binwidth = .01) 

#add lines for means
g<- g + 
  geom_vline(data = filter(dat.sum, Label %in% c("FullInfo")),  
             aes(xintercept=avg, linetype="dashed")) + 
  geom_vline(data = filter(dat.sum, Label %in% c("SAA")),  
            aes(xintercept=avg, linetype="dotted"))
                          
g <-  g + theme_bw() + 
  theme(legend.title=element_blank(), 
        legend.position=c(.5, .85)) + 
  xlim(4.3, 4.87) + ylim(0, 475) +
ylab("") + xlab("Cost") + 
  scale_fill_manual(values = cbPallette2) + 
  scale_linetype_manual(values=c("dashed", "dotted"), guide="none")

ggsave("./Figures/informs_singleNewsvendor_1000.pdf", height=3.8, width=5.3, units="in")


###For both on one graph
g <- dat %>% filter(Method %in% c("SAA", "LOO_avg")) %>%
  ggplot(aes(x=TruePerf, group=Label, fill=Label)) + 
  geom_histogram(size=0, binwidth = .01) 

#add lines for means
g<- g + 
  geom_vline(data = filter(dat.sum, Label %in% c("FullInfo", "SAA", "Pooled")),  
             aes(xintercept=avg, linetype=Label)) 
g <- 
  g + theme_bw() + 
  theme(legend.title=element_blank(), 
        legend.position=c(.5, .85)) + 
  xlim(4.3, 4.87) + ylim(0, 475) +
  ylab("") + xlab("Cost") + 
  scale_fill_manual(values = cbPallette2) + 
  scale_linetype_manual(values=c("dashed", "dotted", "dotted"), guide="none")

ggsave("./Figures/informs_singleNewsvendor_1000_both.pdf", height=3.8, width=5.3, units="in")



