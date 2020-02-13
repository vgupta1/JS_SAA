### Used to analyze results from teh back-testing experiments
library(tidyverse)
library(ggplot2)
library(stringr)
library(forcats)
library(scales) #To properlyf ormat percent axis
library(tikzDevice)

## Load up 3 files of various d's
dat1 = read_csv("../Results/PaperPlots/paperv1_rolling_RossRolling_0.95__20_AdjSales_NoWeekends_ShuffleWithinCol.csv")
dat2 = read_csv("../Results/PaperPlots/paperv1_rolling_RossRolling_0.95__50_AdjSales_NoWeekends_ShuffleWithinCol.csv")
dat3 = read_csv("../Results/PaperPlots/paperv1_rolling_RossRolling_0.95__1000_AdjSales_NoWeekends_ShuffleWithinCol.csv")

dat1 = read_csv("../Results/paperv2_rolling_RossRolling_0.95__20_AdjSales_NoWeekends_ShuffleWithinCol.csv")

dat = rbind(dat1, dat2, dat3)
rm(dat1, dat2, dat3)

###############
#Just look at it in terms of improvement too SAA too.  
dat.saa <- dat %>% 
  filter(Method=="SAA") %>% 
  select(StartDate, K, d, N, TruePerf)

dat <- 
left_join(dat, dat.saa, by=c("StartDate", "K", "d", "N")) %>%
  rename(TruePerf = TruePerf.x,
         SAA = TruePerf.y) %>%
  mutate(RelBenefit =  1 - TruePerf/SAA)

#Summarize
dat.sum <- dat %>% group_by(K, d, N, Method) %>%
  summarise(avgPerf = mean(TruePerf), 
            sdPerf = sd(TruePerf),
            avgTime = mean(time), 
            avgAlpha = mean(alpha), 
            avgRelBenefit = mean(RelBenefit), 
            sdRelBenefit = sd(RelBenefit),
            upRelBenefit = quantile(RelBenefit, .9), 
            downRelBenefit = quantile(RelBenefit, .1), 
            numReps = n(), 
            stdErr = sdRelBenefit/sqrt(numReps)
  )

####
#Relabel the methods
dat.sum <- dat.sum %>% mutate(Method = as.factor(Method),
                              Method = fct_relevel(Method, "OraclePhat", "LOO_avg", "MSE_GM", "Oracle", "LOO_unif", "MSE", "SAA", "FullInfo"), 
                              Labels = fct_recode(Method, `JS-Fixed`= "MSE", 
                                                  `JS-GM` = "MSE_GM",
                                                  `S-SAA-Fixed`="LOO_unif", 
                                                  `S-SAA-GM`="LOO_avg", 
                                                  `Oracle-Fixed`="Oracle", 
                                                  `Oracle-GM`="OraclePhat") )


#############################



#############################
#Plot Relative to SAA for base case different values of $d$
saveHistoricalPlot <- function(d_val){
  pd <- position_dodge(width=.2)
  g<- dat.sum %>% filter(! Method %in%c("SAA","FullInfo"), 
                         N == 10, d ==d_val) %>%
    ggplot(aes(K, avgRelBenefit, group=Labels, color=Labels, linetype=Labels, shape=Labels)) + 
    geom_point(position=pd) + geom_line(position=pd) +
    geom_errorbar(aes(ymax = avgRelBenefit + stdErr, 
                      ymin=avgRelBenefit-stdErr), position=pd) +
    theme_minimal(base_size = 8) + 
    scale_y_continuous(labels=function(x){percent(x, suffix="\\%", accuracy = .1)}) + 
    scale_x_log10()+
    ylab("Benefit over SAA (\\%)") +
    theme(legend.title=element_blank(), 
          legend.position = c(.5, .9), 
          legend.direction= "horizontal") 

  # tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/backtest_", d_val, ".tex"), 
  #      width = 3, height = 3)
  print(g)
  # dev.off()
}

saveHistoricalPlot(20)
saveHistoricalPlot(50)
saveHistoricalPlot(1000)

### Similar graph but varying d and only the ShrunkenSAA
#First filter down to the Shrunken SAA and add a column renaming
pd = position_dodge(.1)
g <- dat.sum %>% 
  filter(Method %in% c("LOO_avg", "LOO_unif"), 
         N == 10) %>%
  arrange(K, d) %>%
  mutate(Labels2 = str_c(Labels, ", $d=", if_else(d==1000, "\\infty", as.character(d)), "$"), 
         ) %>%
  ggplot(aes(K, avgRelBenefit, group=Labels2, color=Labels2, linetype=Labels2, shape=Labels2)) + 
  geom_point(position=pd) + geom_line(position=pd) +
  theme_minimal(base_size = 8) + 
  scale_y_continuous(labels=function(x){percent(x, suffix="\\%", accuracy = .1)}) + 
  scale_x_log10()+
  ylab("Benefit over SAA (\\%)") +
  theme(legend.title=element_blank(), 
        legend.position = c(.7, .3), 
        legend.direction= "horizontal") + 
  guides(color=guide_legend(ncol=1,byrow=TRUE))

# ggsave("../../DataPoolingTex/Paper_V1/Figures/comparing_d.pdf", 
#        g, units="in", width=3.25, height=3.25)

tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/comparing_d.tex", 
     width = 3, height = 3)
g
dev.off()


##########################################
### Comparing what happens as N increases, d = infinity
###############
genHistPlot_N<- function( Nval ){
  dat.sum %>% filter(d == 1000, N == Nval, !Method %in% c("SAA", "FullInfo")) %>%
    ggplot(aes(K, avgRelBenefit, group=Labels, color=Labels, linetype=Labels, shape=Labels)) + 
    geom_point(position=pd) + geom_line(position=pd) +
    geom_errorbar(aes(ymax = avgRelBenefit + stdErr, 
                      ymin=avgRelBenefit-stdErr), position=pd) +
    theme_minimal(base_size = 8) + 
    scale_y_continuous(labels=function(x){percent(x, suffix="\\%", accuracy = .1)}) + 
    scale_x_log10()+
    ylab("Benefit over SAA (\\%)") +
    theme(legend.title=element_blank(), 
          legend.position = c(.5, .9), 
          legend.direction= "horizontal")  
}

g <- genHistPlot_N(20)
# ggsave("../../DataPoolingTex/Paper_V1/Figures/HistPlotN20.pdf", 
#        g, units="in", width=3.25, height=3.25)

tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/HistPlotN20.tex", 
     width = 3, height = 3)
g
dev.off()


g <- genHistPlot_N(40)
# ggsave("../../DataPoolingTex/Paper_V1/Figures/HistPlotN40.pdf", 
#        g, units="in", width=3.25, height=3.25)

tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/HistPlotN40.tex", 
     width = 3, height = 3)
g
dev.off()
