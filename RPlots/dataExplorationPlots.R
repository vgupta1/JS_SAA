####
# Some Data Exploration Plots
####
library(tidyverse)
library(ggplot2)
library(tikzDevice)

dat = read_csv("../RossmanKaggleData/CleanedData/AdjSales_NoWeekends.csv")

#####
##How much missing data?
#####
dat <- dat %>% gather(Store, AdjSales, -Date)
dat %>% group_by(Store) %>%
  summarise(numBlank = sum(is.na(AdjSales)), 
            num = n(), 
            frac = numBlank/num) %>%
  summarise( sum(frac > 0) )

#####
##Look at Avg Demand
#####
g <- dat %>% group_by(Store) %>%
  summarise(avg = mean(AdjSales, na.rm=TRUE), 
            std = sd(AdjSales, na.rm=TRUE)) %>%
  ggplot(aes(x=avg)) + geom_histogram(bins=75) + 
  theme_minimal(base_size=8) + 
  scale_x_continuous(labels=scales::comma) + 
  xlab("") + 
  ylab("")

tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/AvgDailyDemand.tex", 
     width = 3, height = 3)
g
dev.off()


dat %>% group_by(Store) %>%
  summarise(avg = mean(AdjSales, na.rm=TRUE), 
            std = sd(AdjSales, na.rm=TRUE)) %>% 
  summarise(max(avg), 
            min(avg))

#####
##Plot densities for a few choice stores
#####
g2<- dat %>% filter(!is.na(AdjSales), 
                Store %in% c(208, 817, 261, 403, 700)) %>%
  ggplot(aes(x=AdjSales, group=Store, color=as.factor(Store))) + 
    geom_density() + 
  theme_minimal(base_size = 8) + 
  ylab("") + xlab("") + 
  theme(legend.position="none", 
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) + 
  scale_x_continuous(labels=scales::comma)

tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/DensitiesByStore.tex", 
     width = 3, height = 3)
g2
dev.off()

########
#  Try to look at scale adjusted histograms
#########
dat <- read_csv("../RossmanKaggleData/CleanedData/ps_full20.csv", 
                col_names =FALSE)
dat$supp = seq(1, 20)
dat.melt <- dat %>% gather(Store, prob, -supp)
dat.melt <- dat.melt %>% mutate(Store = as.factor(str_replace(Store, "X", "")))
  
g <- dat.melt %>% 
  filter(Store %in% c("208", "817", "261", "403", "700") ) %>%
  ggplot(aes(supp, prob, group=Store, shape=Store, color=Store, linetype=Store))  + 
  geom_point() + geom_line()

g<- g + theme_minimal(base_size=8, base_family="Times New Roman") + 
  theme(legend.position="none") + 
  xlab("Bin Number ($i$)") + ylab("Probability ($p_{ki})$")
  
# ggsave("../../DataPoolingTex/Paper_V1/Figures/SamplePks.pdf", 
#        g, height=3, width=3, units="in")

tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/SamplePks.tex", 
     width = 2.9, height = 2.9)
g
dev.off()


# Compute the .95 quantile for each of them
# Plot across stores
wQuant <- function(probs){
    which(cumsum(probs) >= .95)[1]
}
    
g <- dat.melt %>% group_by(Store) %>%
  summarise(nv_indx = wQuant(prob)) %>%
  ggplot(aes(nv_indx), data=.) + geom_histogram(bins=50) + 
  theme_minimal(base_size=8) + 
  xlab("Bin Number (i)") + ylab("Number of Stores")

# ggsave("../../DataPoolingTex/Paper_V1/Figures/quantileDistribution.pdf", 
#        g, height=3, width=3, units="in")

tikz(file = "../../DataPoolingTex/MS_Submission_R2/Paper/Figures/quantileDistribution.tex", 
     width = 2.9, height = 2.9)
g
dev.off()

dat <- dat.melt %>% spread(supp, prob)
dat <- dat %>% mutate(Store = as.double(str_replace(Store, "X", ""))) %>%
    arrange(Store)

dat %>% 
  filter(Store %in% c(208, 817, 261, 403, 700) ) %>%
  ggplot(aes(supp, prob, group=Store, fill=Store), data=.) + 
  geom_bar()
