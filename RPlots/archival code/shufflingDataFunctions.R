#####
#  This file has been deprecated.  All Shuffling done the cleanData.R file now.  
#####

library(tidyverse)
library(ggplot2)
library(lubridate)
library(forcats)

###  Shuffling the raw data to test correlation structure
set.seed(8675309)

#load in the original counts
dat <- read_csv("../RossmanKaggleData/Results/AdjSalesByStore_Binned20.csv")

#shuffle across entire time horizon for each store
full_shuffle <- dat
full_shuffle[, 2:ncol(dat)] <- lapply(full_shuffle[, 2:ncol(dat)], sample)
write_csv(full_shuffle, 
          "../RossmanKaggleData/Results/AdjSalesByStore_Binned20_FullShuffle.csv")

#shuffle within seasons.  (3 months)  (use grouping and do)
shuff_all_cols <- function(df){
  df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df)], sample)
  return(df)
}

season_shuffle <- dat %>% mutate(Season = quarter(Date, with_year = TRUE, fiscal_start=12))
season_shuffle <- season_shuffle %>% group_by(Season) %>%
  do(shuff_all_cols(.)) %>%
  ungroup() %>%
  select(-Season)
write_csv(season_shuffle, 
          "../RossmanKaggleData/Results/AdjSalesByStore_Binned20_SeasonShuffle.csv")

#shuffle within 1 month 
month_shuffle <- dat %>% mutate(Month = month(Date), Year=year(Date))
month_shuffle <- month_shuffle %>% group_by(Year, Month) %>%
  do(shuff_all_cols(.)) %>%
  ungroup() %>%
  select(-Month, -Year)
write_csv(month_shuffle, 
          "../RossmanKaggleData/Results/AdjSalesByStore_Binned20_MonthShuffle.csv")

#shuffle within week
week_shuffle <- dat %>% mutate(Month = month(Date), Year=year(Date), Week = week(Date))
week_shuffle <- week_shuffle %>% group_by(Year, Month, Week) %>%
  do(shuff_all_cols(.)) %>%
  ungroup() %>%
  select(-Month, -Year, -Week)
write_csv(week_shuffle, 
          "../RossmanKaggleData/Results/AdjSalesByStore_Binned20_WeekShuffle.csv")

#shuffles all columns in a consistent way so that correlations between columns untouched
#but breaks correlations between rows
shuff_days <- function(df){ 
  df[, 2:ncol(df)] <- df[sample(nrow(df)), 2:ncol(df)]  
  return(df)
  }

#shuffle across entire time horizon for each store
full_shuffle <- dat
full_shuffle <- shuff_days(full_shuffle)
write_csv(full_shuffle, 
          "../RossmanKaggleData/Results/AdjSalesByStore_Binned20_FullShuffle_KeepStore.csv")

#shuffle within seasons.  (3 months)  (use grouping and do)
season_shuffle <- dat %>% mutate(Season = quarter(Date, with_year = TRUE, fiscal_start=12))
season_shuffle <- season_shuffle %>% group_by(Season) %>%
  do(shuff_days(.)) %>%
  ungroup() %>%
  select(-Season)
write_csv(season_shuffle, 
          "../RossmanKaggleData/Results/AdjSalesByStore_Binned20_SeasonShuffle_KeepStore.csv")

#shuffle within 1 month 
month_shuffle <- dat %>% mutate(Month = month(Date), Year=year(Date))
month_shuffle <- month_shuffle %>% group_by(Year, Month) %>%
  do(shuff_days(.)) %>%
  ungroup() %>%
  select(-Month, -Year)
write_csv(month_shuffle, 
          "../RossmanKaggleData/Results/AdjSalesByStore_Binned20_MonthShuffle_KeepStore.csv")


#shuffle within week
week_shuffle <- dat %>% mutate(Month = month(Date), Year=year(Date), Week = week(Date))
week_shuffle <- week_shuffle %>% group_by(Year, Month, Week) %>%
  do(shuff_days(.)) %>%
  ungroup() %>%
  select(-Month, -Year, -Week)
write_csv(week_shuffle, 
          "../RossmanKaggleData/Results/AdjSalesByStore_Binned20_WeekShuffle_KeepStore.csv")




##################
#  Analyze rough results
##################
dat.orig <- read_csv("../Results/backtest_RossBacktest_1115_20_0.95_160.csv") %>%
  filter(K==1115, N==20, d==20)
dat.full_shuffle <- read_csv("../Results/tempFullShuffle_RossBacktest_1115_20_0.95_20.csv")
dat.season_shuffle <- read_csv("../Results/tempSeasonShuffle_RossBacktest_1115_20_0.95_20.csv")
dat.month_shuffle <- read_csv("../Results/tempMonthShuffle_RossBacktest_1115_20_0.95_20.csv")
dat.week_shuffle <- read_csv("../Results/tempWeekShuffle_RossBacktest_1115_20_0.95_20.csv")

#bind them together to make a facet graph
dat.orig$Type = "Original"
dat.full_shuffle$Type = "FullShuffle"
dat.season_shuffle$Type = "SeasonShuffle"
dat.month_shuffle$Type = "MonthShuffle"
dat.week_shuffle$Type = "WeekShuffle"

dat = rbind(dat.orig, dat.full_shuffle, dat.season_shuffle, dat.month_shuffle, dat.week_shuffle)
dat <- dat %>% 
  mutate(Type = as.factor(Type)) %>%
  mutate(Type = fct_relevel(Type, c("Original", "WeekShuffle", "MonthShuffle", "SeasonShuffle", "FullShuffle")
  ) )

dat %>% 
  filter(Method=="SAA", K==1115, N==20, d==20) %>%
  ggplot(aes(StartDate, TruePerf, group=Type, color=Type)) + 
  geom_point() + 
  theme_bw() + scale_x_date(date_breaks="6 month") +
  facet_grid(Type~.) + 
  theme(legend.position="None") + 
  ylim(4200, 5800)

########  
#### analyze the results when keeping stores fixed
dat.orig <- read_csv("../Results/backtest_RossBacktest_1115_20_0.95_160.csv") %>%
  filter(K==1115, N==20, d==20)
dat.full_shuffle <- read_csv("../Results/tempFullShuffleStore_RossBacktest_1115_20_0.95_20.csv")
dat.season_shuffle <- read_csv("../Results/tempSeasonShuffleStore_RossBacktest_1115_20_0.95_20.csv")
dat.month_shuffle <- read_csv("../Results/tempMonthShuffleStore_RossBacktest_1115_20_0.95_20.csv")
dat.week_shuffle <- read_csv("../Results/tempWeekShuffleStore_RossBacktest_1115_20_0.95_20.csv")

#bind them together to make a facet graph
dat.orig$Type = "Original"
dat.full_shuffle$Type = "FullShuffle"
dat.season_shuffle$Type = "SeasonShuffle"
dat.month_shuffle$Type = "MonthShuffle"
dat.week_shuffle$Type = "WeekShuffle"

dat = rbind(dat.orig, dat.full_shuffle, dat.season_shuffle, dat.month_shuffle, dat.week_shuffle)
dat <- dat %>% 
  mutate(Type = as.factor(Type)) %>%
  mutate(Type = fct_relevel(Type, c("Original", "WeekShuffle", "MonthShuffle", "SeasonShuffle", "FullShuffle")
                            ) )
dat %>% 
  filter(Method=="SAA", K==1115, N==20, d==20) %>%
  ggplot(aes(StartDate, TruePerf, group=Type, color=Type)) + 
  geom_point() + 
  theme_bw() + scale_x_date(date_breaks="6 month") +
  facet_grid(Type~.) + 
  theme(legend.position="None") + 
  ylim(4200, 5800)
