##Updated plots for illustrating data-pooling

### Creates Plot from Introduction
library(tidyverse)
library(ggplot2)
library(forcats)
library(latex2exp)
library(stringr)
library(extrafont)
loadfonts()

library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

#first is teal, then pink, 
cbPallette2 = c("#00AFBD", "#FBBEBA")

#UPDATED Sept 2019
dat = read_csv("../Results/UPDATED_singleKUnifNewsvendor_10000_1000.csv")

dat <- dat %>% mutate(Method = factor(Method), 
                      Label = fct_recode(Method, 'Shrunken-SAA'="LOO_avg"), 
                      Label = fct_relevel(Label, "FullInfo", "SAA", "Oracle", 'Shrunken-SAA'))

dat.sum<- dat %>% group_by(Label) %>%
  summarise(avg = mean(TruePerf), 
            std = sd(TruePerf), 
            avgAlpha0 = mean(alpha))

g <- dat %>% filter(Method %in% c("SAA", "LOO_avg")) %>%
  ggplot(aes(x=TruePerf, group=Label, fill=Label)) + 
  geom_histogram(size=0, binwidth = .01) 

#Just the SAA
# g <- dat %>% filter(Method %in% c("SAA")) %>%
#   ggplot(aes(x=TruePerf, group=Label, fill=Label)) + 
#   geom_histogram(size=0, binwidth = .01) 


#add lines for means
g<- g +
  geom_vline(data = filter(dat.sum, Label %in% c("FullInfo")),
             aes(xintercept=avg), linetype="dotted")

g <- 
  g + theme_minimal(base_family = "Times New Roman", base_size=10) + 
  theme(legend.title=element_blank(), 
        legend.position=c(.25, .85)) + 
  ylab("") + xlab("Cost") + 
  scale_fill_manual(values = cbPallette2) + 
  xlim(28.5, 34) + ylim(0, 90)


ggsave("../RPlots/Figures/IllustratingDataPooling_SAA_presentation.pdf",
       g, height=5/7 * 4, width=4, units="in")



ggsave("../RPlots/Figures/IllustratingDataPooling_SAA_only_presentation.pdf",
       g, height=5/7 * 4, width=4, units="in")





