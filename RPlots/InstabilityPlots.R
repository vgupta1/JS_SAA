#### 
#Generates the SubOptimality-Instability Tradeoff Curves 
# for the Simple Newsvendor examples in Fig. 5 and EC.1
####
library(tidyverse)
library(ggplot2)
library(tikzDevice)

# dat = read_csv("../Results/TradeoffBadtrue.csv")
# graph_path = "InstabilityBad.tex"

# dat = read_csv("../Results/TradeoffGoodP0true.csv")
# graph_path = "InstabilityGoodP0.tex"

dat = read_csv("../Results/TradeoffGoodStrue.csv")
graph_path = "InstabilityGoodS.tex"


#melt it down for simplicity
dat.melt <- dat %>% select(-LOO,-OR) %>%
  gather(Type, Value,-Alpha)

g1 <- dat.melt %>% 
  ggplot(aes(Alpha, Value, group=Type, color=Type)) + 
  geom_line(aes(linetype=Type)) + 
  theme_minimal(base_size=8) + 
  xlab("$\\alpha$") + ylab("") + 
  theme(legend.title=element_blank(), 
        legend.position=c(.4, .8) )

if (graph_path == "InstabilityGoodS.tex"){
  g1 <- g1 + theme(legend.position = c(.6, .5))
}

tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/", graph_path), 
     width = 2, height = 2)
g1
dev.off()


### The LOO version
#melt it down again
N = 10 #(for scaling purposes)
dat.melt <- dat %>% 
#  mutate(LOO = LOO  / N) %>%
  select(Alpha, LOO,OR) %>%
  gather(Type, Value,-Alpha)

g2 <-   dat.melt %>% 
  ggplot(aes(Alpha, Value, group=Type, color=Type)) + 
  geom_line(aes(linetype=Type)) + 
  theme_minimal(base_size=8) + 
  xlab("$\\alpha$") + ylab("") + 
  theme(legend.title=element_blank(), 
        legend.position=c(.4, .8) )

tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/LOO", graph_path), 
     width = 2, height = 2)
g2
dev.off()


