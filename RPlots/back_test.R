### Used to analyze results from teh back-testing experiments
library(tidyverse)
library(ggplot2)
library(forcats)

## the base case: Dataset 1 
dat = read_csv("../Results/permuted_10_20_RossBacktest_1115_20_0.95_10.csv")

## Med N:  Dataset 2
dat = read_csv("../Results/permuted_20_20_RossBacktest_1115_20_0.95_10.csv")

## Large N:  Dataset 3
dat = read_csv("../Results/permuted_40_20_RossBacktest_1115_20_0.95_40.csv")


## Medium d: Dataset 4
dat = read_csv("../Results/permuted_10_50_RossBacktest_1115_50_0.95_10.csv")

## Large d:  Dataset 5 (only partial so far!)
dat = read_csv("../Results/permuted_10_1000_RossBacktest_1115_1000_0.95_10.csv")

#historical:  Dataset 6
dat = read_csv("../Results/historical_10_20_RossBacktest_1115_20_0.95_10.csv")
dat = read_csv("../Results/hist_10_50_RossBacktest_1115_50_0.95_10.csv")

#New:  HIstorical with d = 1000
dat = read_csv("../Results/historical_10_1000_RossBacktest_1115_1000_0.95_10.csv")




###############


#Express things as percentage of full Info
full_info <- dat %>% 
  filter(Method=="FullInfo") %>% 
  select(StartDate, K, d, N, TruePerf)

#Just look at it in terms of improvement oto SAA too.  
dat.saa <- dat %>% 
  filter(Method=="SAA") %>% 
  select(StartDate, K, d, N, TruePerf)


dat <- left_join(dat, full_info, by=c("StartDate", "K", "d", "N")) %>%
  rename(TruePerf = TruePerf.x,
         FullInfo = TruePerf.y) %>%
  mutate(RelLoss = TruePerf / FullInfo - 1)

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
            avgRelLoss = mean(RelLoss), 
            sdRelLoss = sd(RelLoss), 
            avgRelBenefit = mean(RelBenefit), 
            sdRelBenefit = sd(RelBenefit),
            upRelLoss = quantile(RelLoss, .9), 
            downRelLoss = quantile(RelLoss, .1), 
            numReps = n()
  )

####
#Relabel the methods
dat.sum<- dat.sum %>% mutate(Labels = fct_recode(Method, 
                                       `LOO`="LOO_unif", 
                                       `Oracle GM`="OraclePhat", 
                                       `LOO GM` = "LOO_avg"), 
                             Labels= fct_relevel(Labels, 
                                    c("FullInfo",  "SAA", "MSE", "LOO", "Oracle", "LOO GM", "Oracle GM")
                                    )
                   ) 
#############################



#############################
pd <- position_dodge(width=40)
dat.sum %>% filter(Method !="FullInfo", 
                   K >= 100, N==10, d==20) %>%
  ggplot(aes(K, avgRelLoss, group=Labels, color=Labels)) + 
  geom_point(position=pd) + geom_line(position=pd) + 
  geom_errorbar(aes(ymin=avgRelLoss + sdRelLoss/sqrt(numReps),
                    ymax=avgRelLoss- sdRelLoss/sqrt(numReps)), 
                position=pd) +
  theme_bw()

#Plot Relative to SAA
pd <- position_dodge(width=20)
g<- dat.sum %>% filter(Method != "SAA", Method != "FullInfo",
                   !K %in% c(1, 20, 30, 40, 50, 60, 70, 80, 90)) %>%
  ggplot(aes(K, avgRelBenefit, group=Labels, color=Labels)) + 
  geom_point(position=pd) + geom_line(position=pd) + 
  geom_errorbar(aes(ymin=avgRelBenefit + sdRelBenefit/sqrt(numReps),
                    ymax=avgRelBenefit- sdRelBenefit/sqrt(numReps)), 
                position=pd)
g <- g +
  theme_bw(base_size = 16) + scale_y_continuous(labels = scales::percent) + 
  ylab("Benefit over SAA") +
  theme(legend.title=element_blank(), 
        legend.position = c(.5, .9), 
        legend.direction= "horizontal") 

# ggsave("Figures/columbiaIEOR_permute_10_10.pdf", g, 
#        width=9, height=5.34, units = "in")

# ggsave("Figures/columbiaIEOR_permute_20_20.pdf", g,
#        width=9, height=5.34, units = "in")

# ggsave("Figures/columbiaIEOR_permute_40_20.pdf", g,
#        width=9, height=5.34, units = "in")


ggsave("Figures/columbiaIEOR_permute_10_1000.pdf", g,
       width=9, height=5.34, units = "in")




###############################################
###############################################
###
# Things below this line mostly for exploration. 
###############################################
###############################################
###########################################
### Things to study:
  #Is there at time-series effect in the true performance?
  #Possible sources of correlations:  across days, across stores
  #Possible sources of systematic bias:   across days, across stores
  #Try permutation tests?

dat %>% filter(K==1115, N==10, d==20, Method=="SAA") %>%
  ggplot(aes(StartDate, TruePerf, group=Method, color=Method)) + 
  geom_point() + geom_line()

#Possibility 1: It's seasonal.  Possibility 2:  It's promotions
#No obvious pattern in dates of spikes
dat %>% filter(K==1115, N==10, d==20, Method=="SAA", 
               TruePerf > 6000) %>%
  select(StartDate, TruePerf) %>% View()

#Bring in the promotions data and store data
dat.stores <- read_csv("../RossmanKaggleData/store.csv")
dat.sales = read_csv("../RossmanKaggleData/train.csv", 
               col_types = cols(
                 Store = col_integer(),
                 DayOfWeek = col_integer(),
                 Date = col_date(format = ""),
                 Sales = col_integer(),
                 Customers = col_integer(),
                 Open = col_integer(),
                 Promo = col_integer(),
                 StateHoliday = col_character(),
                 SchoolHoliday = col_integer()
               ))
#drop december and de-trend the way they did before.  
dat.sales <- filter(dat.sales, Open == 1)
dat.sales <- dat.sales %>% mutate(Month = lubridate::month(Date), 
                      Year = lubridate::year(Date) )
t <- lm(Sales~Year + as.factor(Month), data=dat.sales)
beta <- t$coefficients["Year"]
dat.sales <- mutate(dat.sales, 
                    AdjSales = Sales - beta * (Year - 2013) )

#################
##Now look at firm level sumamries across the day
## VG:  This is a bad way to do it, becasue you'd like mean's across stores AND across the window
dat.daily <- dat.sales %>% group_by(Date) %>%
              summarise_all(mean)

#Join the SAA Performance to this bad boy...
#:( Essentially throwing away info during the window, only keeps starts. 
dat.daily <- dat %>% filter(K == 1115, d == 20, N==10, Method=="SAA") %>%
  select(StartDate, TruePerf) %>%
  inner_join(dat.daily, ., by=c("Date"="StartDate")) %>%
  arrange(Date) 

#can we find a structural pattern in perf of SAA now?
#Not an obvious one...  blergh.
dat.daily %>%
  ggplot(aes(SchoolHoliday, TruePerf)) + 
  geom_point() + geom_line() +
  theme_bw()


##############################