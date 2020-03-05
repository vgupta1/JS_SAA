### Plots to provide extra comparisons of Nhat Poisson approx

library(tidyverse)
library(ggplot2)
library(forcats)
library(stringr)
library(tikzDevice)
library(scales) #To properly format percent axis


dat <- read_csv("../src/tempNVarying_syntheticRossman_0.95_true_8.csv")


##Break people up
dat <- dat %>% mutate(Method = as.factor(Method),
               Method = fct_relevel(Method, 
                                    "BetaOptOR", "BetaOptLOO", 
                                    "OraclePhat", "LOO_avg", "MSE_GM", "CV5_avg", 
                                    "Oracle", "LOO_unif", "MSE", "CV5_unif", 
                                    "SAA", "KS", "FullInfo"), 
               Labels = fct_recode(Method, 
                                   `JS-Fixed`= "MSE", 
                                   `JS-GM` = "MSE_GM",
                                   `S-SAA-Fixed`="LOO_unif", 
                                   `S-SAA-GM`="LOO_avg", 
                                   `Oracle-Fixed`="Oracle", 
                                   `Oracle-GM`="OraclePhat", 
                                   `Oracle-Beta`="BetaOptOR", 
                                   `S-SAA-Beta`="BetaOptLOO"))
#Throw away SAA and full-info bc you won't need them
dat <- filter(dat, ! Method %in% c("FullInfo", "SAA"))
dat <- dat %>% mutate(Anchor = fct_collapse(Method, 
                                     Fixed=c("Oracle", "LOO_unif"), 
                                     GM = c("OraclePhat", "LOO_avg"), 
                                     Beta = c("BetaOptOR", "BetaOptLOO")
                                     ), 
               Type = fct_collapse(Method, 
                                   Oracle=c("Oracle", "OraclePhat", "BetaOptOR"), 
                                   LOO=c("LOO_unif", "LOO_avg", "BetaOptLOO"))
               )

#Match to Appropriate Oracle
dat.oracle <- filter(dat, Type == "Oracle") %>% 
  select(Run, K, d, N, Anchor, TruePerf)
dat <- 
  left_join(dat, dat.oracle, by=c("Run", "K", "d", "N", "Anchor")) %>%
  rename(TruePerf = TruePerf.x,
         Oracle = TruePerf.y) %>%
  mutate(RelGap =  1 - Oracle/TruePerf)

##summarize
dat.sum <- dat %>% group_by(K, d, N, Method, Anchor, Type) %>%
  summarise(avgPerf = mean(TruePerf), 
            sdPerf = sd(TruePerf),
            avgTime = mean(time), 
            avgAlpha = mean(alpha), 
            avgRelGap = mean(RelGap), 
            sdRelGap = sd(RelGap),
            upRelGap = quantile(RelGap, .9), 
            downRelGap = quantile(RelGap, .1), 
            numReps = n(), 
            stdErr = sdRelGap/sqrt(numReps)
  )

dat.sum <- mutate(dat.sum, 
                  MethodN = str_c(Method, N))


dat.sum %>% 
  filter(Type=="LOO") %>%
  ggplot(aes(K, avgRelGap, group=MethodN, color=Anchor)) + 
  geom_point(aes(shape=Anchor)) + 
  geom_line(aes(linetype=as.factor(N))) + 
  theme_minimal(base_size=8) +
  scale_x_log10() + 
  scale_y_continuous(labels=scales::percent_format(suffix="/%"))+ 
  ylab("Relative Gap to Oracle (/%)") + 
  theme(legend.title = element_blank(), 
        legend.position="top")

