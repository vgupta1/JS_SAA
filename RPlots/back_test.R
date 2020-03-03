### Used to analyze results from teh back-testing experiments
library(tidyverse)
library(ggplot2)
library(stringr)
library(forcats)
library(scales) #To properly format percent axis
library(tikzDevice)
library(xtable)

#Used to make the tables a bit sparser
K_grid = c(10, 32,  64, 128,  256, 362,  431,  512,  609,  724,  861, 1024, 1115)



## Load up 3 files of various d's
dat1 = read_csv("../Results/PaperPlots/paperv1_rolling_RossRolling_0.95__20_AdjSales_NoWeekends_ShuffleWithinCol.csv")
dat2 = read_csv("../Results/PaperPlots/paperv1_rolling_RossRolling_0.95__50_AdjSales_NoWeekends_ShuffleWithinCol.csv")
dat3 = read_csv("../Results/PaperPlots/paperv1_rolling_RossRolling_0.95__1000_AdjSales_NoWeekends_ShuffleWithinCol.csv")


dat1 = read_csv("../Results/paperv2_Roll2_Ross_0.95__20_AdjSales_NoWeekends_ShuffleWithinCol.csv")
dat2 = read_csv("../Results/paperv2_Roll2_Ross_0.95__50_AdjSales_NoWeekends_ShuffleWithinCol.csv")
dat3 = read_csv("../Results/paperv2_Roll2_Ross_0.95__1000_AdjSales_NoWeekends_ShuffleWithinCol.csv")

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
                                                  `S-SAA-Beta`="BetaOptLOO") )


#############################



#############################
#Plot Relative to SAA for base case different values of $d$
saveHistoricalPlot <- function(d_val){
  pd <- position_dodge(width=.2)
  g<- dat.sum %>% filter(! Method %in%c("SAA","FullInfo", "KS", "CV5_avg", "CV5_unif", "BetaOptLOO", "BetaOptOR"), 
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

  tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/backtest_", d_val, ".tex"),
       width = 3, height = 3)
  print(g)
  dev.off()
}

saveHistoricalPlot(20)
saveHistoricalPlot(50)
saveHistoricalPlot(1000)

#######
#Dump the info for a full table in excel that uses all the policies.
#######

saveHistoricalTable<- function(d_val)
{
  dat.sum<- ungroup(dat.sum)
  tdat <- dat.sum %>% filter(d == d_val, N==10, !Method %in% c("FullInfo", "CV5_avg", "CV5_unif")) %>% 
    filter(K %in% K_grid) %>%
    select(K, Labels, avgRelBenefit) %>%
    mutate(avgRelBenefit = avgRelBenefit * 100) %>%
    spread(Labels, avgRelBenefit) 
  xdat <- xtable(tdat, auto=TRUE, digits=2)
  
  print(xdat, str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/HistDataTrueAllPolicies", d_val, ".tex"), 
        type="latex",
        booktabs=TRUE, 
        include.rownames=FALSE, 
        floating = FALSE)  
}
saveHistoricalTable(20)
saveHistoricalTable(50)
saveHistoricalTable(1000)


### Similar graph but varying d and only the ShrunkenSAA
#First filter down to the Shrunken SAA and add a column renaming
pd = position_dodge(.1)
g <- dat.sum %>% 
  filter(Method %in% c("LOO_avg", "LOO_unif", "BetaOptLOO"), 
         N == 10) %>%
  arrange(K, d) %>%
  mutate(AnchorType = str_remove(Labels, "S-SAA-"), 
         dType = str_c("$d=", if_else(d==1000, "\\infty", as.character(d)), "$"), 
           Labels2= str_c(AnchorType, ", ", dType), 
         ) %>%
  mutate( dType = as.factor(dType), 
          dType = fct_relevel(dType, "$d=20$", "$d=50$") 
          ) %>%
  ggplot(aes(K, avgRelBenefit, group=Labels2, color=AnchorType, 
             linetype=dType, shape=AnchorType)) + 
  geom_point(position=pd, size=1) + geom_line(position=pd) +
  theme_minimal(base_size = 8) + 
  scale_y_continuous(labels=function(x){percent(x, suffix="\\%", accuracy = .1)}) + 
  scale_x_log10()+
  ylab("Benefit over SAA (\\%)") +
  theme(legend.title=element_blank(), 
        legend.position = c(.7, .3), 
        legend.direction= "horizontal") + 
  guides(color="none")


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
tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/HistPlotN20.tex", 
     width = 3, height = 3)
print(g)
dev.off()

g <- genHistPlot_N(40)
tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/HistPlotN40.tex", 
     width = 3, height = 3)
g
dev.off()


genHistTable_N <- function(Nval) 
{
  dat.sum<- ungroup(dat.sum)
  tdat <- dat.sum %>% filter(d == 1000, N==Nval,
                             !Method %in% c("FullInfo", "CV5_avg", "CV5_unif")) %>% 
    filter(K %in% K_grid) %>%
    select(K, Labels, avgRelBenefit) %>%
    mutate(avgRelBenefit = avgRelBenefit * 100) %>%
    spread(Labels, avgRelBenefit) 
  xdat <- xtable(tdat, auto=TRUE, digits=2)
  
  print(xdat, str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/HistDataLargeN_", Nval, ".tex"), 
        type="latex",
        booktabs=TRUE, 
        include.rownames=FALSE, 
        floating = FALSE)  
}
genHistTable_N(20)
genHistTable_N(40)
