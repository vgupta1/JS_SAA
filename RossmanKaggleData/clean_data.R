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
               Date <= "2015-01-01", 
               Date >= "2014-01-01") %>%
  ggplot(aes(Date, Sales)) + geom_point()

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
computeCounts(dat1, 10)
computeSupport(dat1, 10)

#not the most beautiful way to do this, but gets it done
d = 50
store.supp<- dat %>% group_by(Store) %>%
            do(computeSupport(., d))
store.counts <- dat %>% group_by(Store) %>%
            do(computeCounts(., d))

write_csv(store.supp, str_c("../RossmanKaggleData/Results/support", d, ".csv"))
write_csv(store.counts, str_c("../RossmanKaggleData/Results/counts", d, ".csv"))


