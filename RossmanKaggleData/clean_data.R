#Playing with the Rossman Data
library(tidyverse)
library(ggplot2)
library(broom)

dat = read_csv("../RossmanKaggleData/train.csv", 
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

#Throw away days they are closed because zero sales
dat <- filter(dat, Open == 1)

#How to treat december
dat <- dat %>% mutate(Month = lubridate::month(Date), 
                      Year = lubridate::year(Date) )

#Should we throw away december?
#Seemingly in some stores there is a strong effect, others not
dat %>% filter(Store == 307,
  ggplot(aes(Date, Sales)) + geom_point()

dat %>% filter(Store < 10) %>% 
  ggplot(aes(Date, Sales, group=Store, color=as.factor(Store))) + 
  geom_smooth(method="lm") + 
  geom_smooth(linetype="dotted")

salesByTime.lm <- lm(Sales~Date, data=dat)

#double check with a linear gression
salesByMonth.lm <- lm(Sales ~ Year + as.factor(Month), data=dat)



#There are time trends, but most significant seems to be in December
#Also a 2-3% growh per year, maybe due to GDP?

#Correct the linear trend, and drop december
dat <- filter(dat, Month != 12)
beta <- salesByMonth.lm$coefficients["Year"]
dat <- dat %>% mutate(AdjSales = Sales - beta * (Year - 2013) )

#Notes:  After filtering
# Store with most/least data: 85 / 644
# Store with most/least avg sales 817 / 307
# average store has 712 data points, 6,832 Sales Volume
# Num Stores: 1115
dat %>% group_by(Store) %>% 
  summarise(num=n(), sales = mean(Sales)) %>%
  arrange((sales))


#Steps:  
# Fix a discretization grid for each store 
# Compute the true p for each store over entire time period and dump to file
# Do this for several different discretizations to assess (via simulation) the effect of d


#how many bins would a histogram a representation use?
dat %>% group_by(Store) %>% 
  summarise(numBins = nclass.scott(Sales)) %>%
  with(., unique(numBins))



#To do:
# Need to encode actual counts somehow still
computeSupport <- function(df, d){
  #Take range, blow up by .1% divy into d pieces and round everyone to the left. 
  rng = range(df$AdjSales)
  bin_width <- (1.001*rng[2] - .999*rng[1]) / d
  supp <- seq(.999 * rng[1], 1.001 * rng[2], bin_width)[1:d]
  return(data.frame(t(supp)))
}

## VG: need to revisit 
#This does not handle missing values quite correctly.  
computeCounts <- function(df, d){
  supp <- computeSupport(df, d)
  counts <- findInterval(df$AdjSales, supp, rightmost.closed=TRUE) %>% 
    table() %>%
    as.numeric() 
  return(data.frame(t(counts)))
}



###############################
#  Compute relevant support and counts for d = 10, 20, 50, 100
#not the most beautiful way to do this, but gets it done
d = 20
store.supp<- dat %>% group_by(Store) %>%
            do(computeSupport(., d))
store.counts <- dat %>% group_by(Store) %>%
            do(computeCounts(., d))

#Just shuffle them once for safety....
# N <- dim(store.supp)[1]
# store_order <- sample(N)
# 
# store.supp[store_order, ] %>%
#   write_csv(str_c("../RossmanKaggleData/Results/support", d, ".csv"))
# store.counts[store_order, ] %>%
#   write_csv(str_c("../RossmanKaggleData/Results/counts", d, ".csv"))

#################################

#####
#What consitutes a reasonble anchor?
######
dat %>% ggplot(aes(AdjSales)) + geom_density() + 
  stat_function(fun= dexp, geom="line", color="red", args = list(rate = 1/mean(dat$AdjSales)))


#Attempt 1:  Fit a density to the overall distribution
rng = range(dat$AdjSales)
mu = mean(dat$AdjSales)
std = sd(dat$AdjSales)

dat %>% filter(AdjSales > 0) %>%
  mutate(logAdjSales = log(AdjSales)) %>%
  summarise(mean(logAdjSales), sd(log(AdjSales))
  )

#a lognormal fits the entire dataset very well....
dat %>% ggplot(aes(AdjSales)) + geom_density() + 
  stat_function(fun= dlnorm, geom="line", color="red", 
                args = list(meanlog = 8.71, sdlog = .433))

#compare to the distributions by stores
#Looks awful.. :)
iStore = 878
dat %>% filter(Store == iStore) %>%
  ggplot(aes(AdjSales)) + geom_density() + 
  stat_function(fun= dlnorm, geom="line", color="red", 
                args = list(meanlog = 8.71, sdlog = .433))

#########
# Preprocess some data for the backtesting experiments
#########
#add a column which indicates which suppport point demand belongs to. 
#not the most indiomatic or efficient way to do this, but whatever....
dat$xi_indx = 0
for( ix in 1:nrow(dat) ){
  ix_store = which(store.supp$Store == dat$Store[ix] )
  dat[ix, "xi_indx"] = findInterval(dat[ix, "AdjSales"], 
                                     store.supp[ix_store, 2:(d-1)], 
                                    rightmost.closed = TRUE)
}

#Look at subset of dates where everyone has data...
dat %>% group_by(Date) %>% 
  summarise(nStoresData = n()) %>%
  filter(nStoresData == 1115) %>%
  summarise(n())




t(store.supp[1, ])



findInterval(dat1$AdjSales, store.supp[dat1$Store])



counts <- findInterval(df$AdjSales, supp, rightmost.closed=TRUE) %>% 


#In Julia, will subset into groups of size N = 10, 20, 30 ..., 100 data points...
#and compute relevant mhat for





#Save another version of datasets where you shuffle columns (break dependence between stores)

