##Plotting the synthetic Rossman Experiment 
#
library(tidyverse)
library(ggplot2)

dat <- read_csv("../Results/res_RossN_1115_50_0.95_true_100.csv")
dat <- read_csv("../Results/res_RossN_1115_50_0.95_false_100.csv")

full_info <- dat %>% 
  filter(Method=="FullInfo", Run==1) %>% 
  select(K, d, N, TruePerf)

dat <- left_join(dat, full_info, by=c("K", "d", "N")) %>%
    rename(TruePerf = TruePerf.x,
           FullInfo = TruePerf.y) %>%
  mutate(RelLoss = TruePerf / FullInfo - 1)


#Summarize
dat.sum <- dat %>% group_by(K, d, N, Method) %>%
    summarise(avgPerf = mean(TruePerf), 
              avgTime = mean(time), 
              avgAlpha = mean(alpha), 
              avgRelLoss = mean(RelLoss), 
              sdRelLoss = sd(RelLoss), 
              upRelLoss = quantile(RelLoss, .9), 
              downRelLoss = quantile(RelLoss, .1)
    )

dat.sum %>% filter(Method != "FullInfo", K==1115) %>%
  ggplot(aes(N, avgRelLoss, group=Method, color=Method)) + 
  geom_point() + geom_line() + 
  theme_bw() + 
  xlim(10, 50)

##By K 
dat.sum %>% filter(Method != "FullInfo", N==10) %>%
  ggplot(aes(K, avgRelLoss, group=Method, color=Method)) + 
  geom_point() + geom_line() + 
  theme_bw()

dat.sum %>% filter(Method != "FullInfo", K==1115) %>%
  arrange(N, Method) %>% 
  select(Method, avgRelLoss) %>% 
  View()

dat.sum %>% ungroup() %>%
  filter(Method != "FullInfo") %>%
  select(K,d, avgRelLoss, Method) %>% 
  spread(Method, avgRelLoss)

dat.sum %>% ungroup() %>%
  filter(Method != "FullInfo") %>%
  select(K,d, avgAlpha, Method) %>% 
  spread(Method, avgAlpha)

#######
dat.sum %>% filter(K > 600) %>%
  arrange(K, d)
