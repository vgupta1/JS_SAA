#exploring Store Features

library(tidyverse)
library(ggplot2)

dat <- read_csv("../RossmanKaggleData/store.csv")
#Recode some things as factors for safety
dat <- dat %>% mutate(StoreType= as.factor(StoreType), 
              Assortment = as.factor(Assortment), 
              hasPromo2 = !is.na(PromoInterval), 
              hasCompStart = !is.na(CompetitionOpenSinceMonth)
              )

dat %>% select(-PromoInterval, -CompetitionOpenSinceMonth, -CompetitionOpenSinceYear) %>% 
  summarise_all(funs(sum(is.na(.))))

clusters <- dat %>% select(-PromoInterval) %>% 
            kmeans(centers = 10, nstart=10)

kmeans(dat, )
