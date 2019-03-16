#### 
#Generates the SubOptimality-Instability Tradeoff Curves 
# for Simple Newsvendor example
####
library(tidyverse)
library(ggplot2)
library(latex2exp)
library(stringr)
library(extrafont)

library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()


# dat = read_csv("../Results/PaperPlots/TradeoffBad.csv")
# graph_path = "InstabilityBad.pdf"

dat = read_csv("../Results/PaperPlots/TradeoffGoodP0.csv")
graph_path = "InstabilityGoodP0.pdf"

# dat = read_csv("../Results/PaperPlots/TradeoffGoodS.csv")
# graph_path = "InstabilityGoodS.pdf"


#melt it down for simplicity
dat.melt <- dat %>% select(-LOO,-OR) %>%
  gather(Type, Value,-Alpha)

g1 <- dat.melt %>% 
  ggplot(aes(Alpha, Value, group=Type, color=Type)) + 
  geom_line(aes(linetype=Type)) + 
  theme_minimal(base_size=10, base_family = "Times New Roman") + 
  xlab(TeX("$\\alpha$")) + ylab("") + 
  theme(legend.title=element_blank(), 
        legend.position=c(.4, .8) )

if (graph_path == "InstabilityGoodS.pdf"){
  g1 <- g1 + theme(legend.position = c(.6, .5))
}
ggsave(str_c("../../DataPoolingTex/Paper_V1/Figures/", graph_path), 
       g1, height=2, width=2, units="in")


### The LOO version
#melt it down agai
N = 10 #(for scaling purposes)
dat.melt <- dat %>% 
  mutate(LOO = LOO  / N) %>%
  select(Alpha, LOO,OR) %>%
  gather(Type, Value,-Alpha)

g2 <-   dat.melt %>% 
  ggplot(aes(Alpha, Value, group=Type, color=Type)) + 
  geom_line(aes(linetype=Type)) + 
  theme_minimal(base_size=10, base_family = "Times New Roman") + 
  xlab(TeX("$\\alpha$")) + ylab("") + 
  theme(legend.title=element_blank(), 
        legend.position=c(.4, .8) )

ggsave(str_c("../../DataPoolingTex/Paper_V1/Figures/LOO", graph_path), 
       g2, height=2, width=2, units="in")

