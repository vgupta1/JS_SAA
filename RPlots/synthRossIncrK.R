##Plotting results out of the syntheticData Experiment
#
library(tidyverse)
library(ggplot2)
library(forcats)
library(stringr)
library(tikzDevice)
library(scales) #To properly format percent axis
library(xtable)
K_grid = c(10, 32,  64, 128,  256, 362,  431,  512,  609,  724,  861, 1024, 1115)


usePoisson = FALSE

dat <- read_csv(str_c("../Results/paperv2_syntheticRossman_0.95_", usePoisson, "_200.csv"))


#express everything as benefit over SAA
saa_perf <- dat %>%
  filter(Method =="SAA") %>%
  select(Run, K, d, N, TruePerf)

dat <- left_join(dat, saa_perf, by=c("Run", "K", "d", "N")) %>%
    rename(TruePerf = TruePerf.x,
           SAA = TruePerf.y) %>%
  mutate(RelBenefit =  1 - TruePerf / SAA)

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
              stdErrRelBenefit = sdRelBenefit/ sqrt(n())
    )


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
                                                      `CV-5-Fixed`="CV5_unif", 
                                                      `CV-5-GM`="CV5_avg", 
                                                      `Oracle-Beta` = "BetaOptOR", 
                                                      `S-SAA-Beta` = "BetaOptLOO") )

##By K 
pd = position_dodge(.2)
g <- dat.sum %>% filter(!Method %in% c("SAA", "FullInfo", "KS", "CV5_unif", "CV5_avg", "BetaOptOR", "BetaOptLOO"), N==10) %>%
  ggplot(aes(K, avgRelBenefit, group=Labels, color=Labels, shape=Labels, linetype=Labels)) + 
  geom_point(position=pd) + geom_line(position=pd) + 
  geom_errorbar(aes(ymin=avgRelBenefit-stdErrRelBenefit, ymax=avgRelBenefit + stdErrRelBenefit), position=pd) + 
  theme_minimal(base_size = 8) + 
  scale_x_log10() + 
    scale_y_continuous(limits=c(NA, .175), 
                       labels=function(x){percent(x, suffix="\\%", accuracy = .1)}) + 
    theme(legend.title=element_blank(), 
          legend.position = c(.55, .85)) +
    ylab("Benefit over SAA (\\%)") + 
  guides(color=guide_legend(nrow=3,byrow=TRUE))


tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/synthData", usePoisson, ".tex"), 
     width = 3.2, height = 3.2)
g
dev.off()


##Standard Deviation Plot
g <- dat.sum %>% filter(!Method %in% c("FullInfo", "KS", "CV5_unif", "CV5_avg"), N==10) %>%
  ggplot(aes(K, sdPerf, group=Labels, color=Labels, shape=Labels, linetype=Labels)) + 
  geom_point(position=pd) + geom_line(position=pd) +
  theme_minimal(base_size = 8) + 
  scale_x_log10() + scale_y_log10(limits=c(NA, 3000)) + 
  theme(legend.text=element_text(size=6),
    legend.title=element_blank(), 
        legend.position=c(.7, .8)) +
  ylab("Std. of Performance") + 
  guides(color=guide_legend(ncol=2,byrow=TRUE))
g

tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/synDataStdDev_", usePoisson, ".tex"), 
     width = 3.2, height = 3.2)
g
dev.off()


##### Make a big table for the appendix of all the methods
dat.sum<- ungroup(dat.sum)
tdat <- dat.sum %>% filter(!Method %in% c("FullInfo", "CV5_avg", "CV5_unif")) %>% 
  filter(K %in% K_grid) %>%
  select(K, Labels, avgRelBenefit) %>%
  mutate(avgRelBenefit = avgRelBenefit * 100) %>%
  spread(Labels, avgRelBenefit) 
xdat <- xtable(tdat, auto=TRUE, digits=2)

print(xdat, str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/SynthDataPoisson", usePoisson, "_20_10.tex"), 
      type="latex",
      booktabs=TRUE, 
      include.rownames=FALSE, 
      floating = FALSE)  

#####
#Plots of alpha Convergence
####
g <- dat.sum %>% filter(Method %in% c("MSE", "MSE_GM", "LOO_unif", "Oracle")) %>%
    ggplot(aes(K, avgAlpha, group=Labels, color=Labels, shape=Labels, linetype=Labels)) + 
  geom_point(position=pd) + geom_line(position=pd) + 
  theme_minimal(base_size = 8) + 
  scale_x_log10() + scale_y_continuous() + 
  theme(legend.title=element_blank(), legend.position=c(.6, .8)) +
  ylab("$\\alpha$") + 
  guides(color=guide_legend(nrow=3,byrow=TRUE))

if(!usePoisson)
{
  g <- g + theme(legend.position=c(.6, .4))
}
g

tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/synDataAlphaNonGM_", usePoisson, ".tex"), 
     width = 3.2, height = 3.2)
g
dev.off()


g <- dat.sum %>% filter(Method %in% c("OraclePhat", "LOO_avg", "BetaOptOR", "BetaOptLOO"), N==10) %>%
  ggplot(aes(K, avgAlpha, group=Labels, color=Labels, shape=Labels, linetype=Labels)) + 
  geom_point(position=pd) + geom_line(position=pd) + 
  theme_minimal(base_size = 8) + 
  scale_x_log10() + scale_y_continuous() + 
  theme(legend.title=element_blank(), legend.position=c(.6, .8)) +
  ylab("$\\alpha$") + 
  guides(color=guide_legend(nrow=4,byrow=TRUE))
g

if(!usePoisson)
{
  g <- g + theme(legend.position=c(.8, .2))
}
g

tikz(file = str_c("../../DataPoolingTex/MS_Submission_R2/Paper/Figures/synDataAlphaGM_", usePoisson, ".tex"), 
     width = 3.2, height = 3.2)
g
dev.off()

