library(tidyverse)
library(ggplot2)
library(forcats)
library(stringr)
library(tikzDevice)
library(scales) #To properly format percent axis

usePoisson = TRUE

dat = read_csv(str_c("../Results/paperv2__cv_synthetic_0.95_", usePoisson, "_200.csv"))

#split and relabel everybody for ease
dat <- dat %>% mutate(Method = as.factor(Method), 
                              Method2 = fct_recode(Method, 
                                                   CV10_GM = "CV10_avg", 
                                                   CV10_Fixed="CV10_unif", 
                                                   CV2_GM = "CV2_avg", 
                                                   CV2_Fixed="CV2_unif", 
                                                   CV5_GM = "CV5_avg", 
                                                   CV5_Fixed="CV5_unif", 
                                                   LOO_GM="LOO_avg", 
                                                   LOO_Fixed ="LOO_unif", 
                                                   OR_Fixed="Oracle", 
                                                   OR_GM="OraclePhat",
                                                   SAA_Decoupled="SAA",
                                                   FullInfo_Decoupled="FullInfo"
                              )                                      
) 

getCVType <- function(s){ str_split(s, "_", simplify=TRUE)[, 1]}
getAnchorType <- function(s){ str_split(s, "_", simplify=TRUE)[, 2]}

dat <- dat %>% mutate(Folds=getCVType(Method2), 
                              Anchor=getAnchorType(Method2))


###rescale everyoen relative to their own oracle
dat_fixed <- filter(dat, Anchor == "Fixed")
dat_GM <- filter(dat, Anchor == "GM")

#rescale relative to oracle
rescale_OR <- function(dat){
  or_perf <- dat %>%
    filter(Folds =="OR") %>%
    select(Run, K, d, N, TruePerf)
  
  dat <- left_join(dat, or_perf, by=c("Run", "K", "d", "N")) %>%
    rename(TruePerf = TruePerf.x,
           OR = TruePerf.y) %>%
    mutate(RelSubOpt =  1 - TruePerf / OR)
  return(dat)
}
dat_fixed <- rescale_OR(dat_fixed)
dat_GM <- rescale_OR(dat_GM)

my_summarize <- function(dat){
  dat.sum <- dat %>% group_by(K, d, N, Method2, Folds, Anchor) %>%
    summarise(avgPerf = mean(TruePerf), 
              sdPerf = sd(TruePerf),
              avgTime = mean(time), 
              avgAlpha = mean(alpha), 
              avgRelSubOpt = mean(RelSubOpt), 
              sdRelSubOpt = sd(RelSubOpt), 
              upRelSubOpt = quantile(RelSubOpt, .9), 
              downRelSubOpt = quantile(RelSubOpt, .1), 
              stdErrRelSubOpt = sdRelSubOpt/ sqrt(n())
    )
  return(dat.sum)
}
dat_fixed <- my_summarize(dat_fixed)
dat_GM <- my_summarize(dat_GM)


genPlot<- function(dat)
{
  g<- dat %>% filter (Folds != "OR") %>%
    ggplot(aes(K, avgRelSubOpt, group=Folds, color=Folds, linetype=Folds)) + 
    geom_line() + 
    geom_point(aes(shape=Folds)) + 
    theme_minimal(base_size = 8) + 
    scale_x_log10() + 
    scale_y_continuous(labels=function(x){percent(x, suffix="\\%", accuracy = .1)}) + 
    theme(legend.text=element_text(size=6),
          legend.title=element_blank(), 
          legend.position=c(.8, .2)) +
    ylab("Relative SubOpt (\\%)") 
  return(g)
}

##Save down the two plots
tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/compareCV_Fixed_Poisson", usePoisson, ".tex"),
     width = 3, height = 3)
print(genPlot(dat_fixed))
dev.off()

tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/compareCV_GM_Poisson", usePoisson, ".tex"),
     width = 3, height = 3)
print(genPlot(dat_GM))
dev.off()





