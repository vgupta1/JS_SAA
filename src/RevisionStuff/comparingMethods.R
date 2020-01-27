library(stringr)
library(scales)
library(tidyverse)
library(ggplot2)
library(forcats)
library(latex2exp)
library(purrr)
library(broom)
library(extrafont)
loadfonts()
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

#high idiosyncratic attempt
#dat = read_tsv("../src/RevisionStuff/Results/highIdioResults_10_5_1.0_100_8675309.tsv")
#dat = read_tsv("../src/RevisionStuff/Results/highIdioResults_20_5_1.0_100_8675309.tsv")
dat = read_tsv("../src/RevisionStuff/Results/highIdioResults_10_3_10.0_100_8675309.tsv")

#Low idiosyncratic effect attempts
#dat = read_tsv("../src/RevisionStuff/temp3Results_10_0.0_50_8675309.tsv")
dat = read_tsv("../src/RevisionStuff/Results/lowIdioResults_10_0.0_10.0_100_8675309.tsv")

#medium effects
dat = read_tsv("../src/RevisionStuff/temp3Results_10_3_10.0_50_8675309.tsv")
dat = read_tsv("../src/RevisionStuff/medIdioResults_10_2_5.0_100_8675309.tsv")

dat <- dat %>% mutate(Method = fct_relevel(Method, c("SAA", "S-SAA", "NW", "HS-NW", "NWbar", "S-NW")), 
                      Label = fct_recode(Method, `NW-f`="NW", 
                                                  `NW-fbar`="NWbar"))

dat.sum <- dat %>% group_by(Label) %>%
  summarise(avg = mean(TruePerf), 
            std = sd(TruePerf))

zstar = filter(dat.sum, Label=="FullInfo") %>% select(avg)
dat$FullInfo <- zstar[[1]]
dat <- mutate(dat, RelPerf = TruePerf/FullInfo - 1)
dat.sum <- mutate(dat.sum, rel_avg = avg/zstar[[1]] -1 )

# dat %>% filter(Label != "FullInfo") %>%
#   ggplot(aes(RelPerf, group=Label, color=Label, fill=Label)) + 
#   geom_density(position="identity", alpha=.2) + 
#   geom_vline(aes(xintercept=rel_avg, color=Label), linetype="dashed", 
#              data = filter(dat.sum, Label!="FullInfo")) + 
#   theme_minimal(base_size = 10, base_family = "Times New Roman") + 
#   ylab("") + xlab("% Loss") + 
#   theme(legend.title=element_blank()) + 
#   scale_x_continuous(labels= scales::percent)

#build custom labels based on the means
dat.sum <- dat.sum %>% mutate(pLabel = str_c(Label, " - ", percent(rel_avg)))
v <- as.character(dat.sum$Label)
names(v) <- as.character(dat.sum$pLabel)
dat.sum <- dat.sum %>% mutate(pLabel = fct_recode(Label, !!!v))
dat <- dat %>% mutate(pLabel = fct_recode(Label, !!!v))


# dat %>% filter(!Method %in% c("FullInfo", "SAA", "HS-NW")) %>%
#   ggplot(aes(Method, RelPerf)) + 
#   geom_boxplot()

### Publication quality pLot
##Convert the densities into line plots first
# dat.density <- dat %>% filter(Method != "FullInfo") %>%
#     nest(-Label) %>%
#   mutate( smooth_vals = map(data, ~ density(.$RelPerf)), 
#           smooth_table = map(smooth_vals, tidy) ) %>%
#   unnest(smooth_table, .drop=TRUE)
# dat.density<- dat.density %>% mutate(pLabel = fct_recode(Label, !!!v))

### for the Low-
(
g<- 
    dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "S-NW")) %>%
    ggplot(aes(Label, RelPerf, color=pLabel)) +
    geom_boxplot() + 
    theme_minimal(base_size = 8, base_family = "Times New Roman") +
    xlab("") + ylab("(%) Loss") +
    theme(legend.title=element_blank(),
          legend.position=c(.75, .6),
          plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    scale_y_continuous(labels= scales::percent)
)  
#ggsave("~/Dropbox/NSF Career 2019/V1/Figures/LowIdiosyncracy_full.pdf", g, height=2.16, width=2.16, units="in")


### for the High
(
  g<- dat.density %>% filter(! Label %in% c("HS-NW")) %>%
    ggplot(aes(x,y, group=pLabel, color=pLabel, linetype=pLabel)) +
    geom_line() +
    theme_minimal(base_size = 8, base_family = "Times New Roman") +
    ylab("") + xlab("(%) Loss") +
    theme(legend.title=element_blank(), 
          legend.position=c(.75, .6), 
          plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    scale_x_continuous(labels= scales::percent, limits = c(0, .5))
)  

ggsave("~/Dropbox/NSF Career 2019/V1/Figures/HighIdiosyncracy.pdf", g, height=2.16, width=2.16, units="in")


### for the Medium-
(
  g<- dat.density %>% filter(! Label %in% c("HS-NW")) %>%
    ggplot(aes(x,y, group=pLabel, color=pLabel, linetype=pLabel)) +
    geom_line() +
    theme_minimal(base_size = 8, base_family = "Times New Roman") +
    ylab("") + xlab("(%) Loss") +
    theme(legend.title=element_blank(), 
          legend.position=c(.75, .6), 
          plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    scale_x_continuous(labels= scales::percent)
)  
ggsave("~/Dropbox/NSF Career 2019/V1/Figures/MedIdiosyncracy.pdf", g, height=2.16, width=2.16, units="in")




