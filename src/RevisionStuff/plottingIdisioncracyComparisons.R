library(stringr)
library(scales)
library(tidyverse)
library(ggplot2)
library(forcats)
library(purrr)
library(broom)
#loadfonts()
library(showtext)
font_add("Times New Roman", "Times New Roman.ttf")
showtext_auto()

library(tikzdevice)

#high idiosyncratic attempt  - VG Complete!
dat = read_tsv("../src/RevisionStuff/Results/highIdioResults_20_5_1.0_100_8675309.tsv")

#Low idiosyncratic effect attempts
dat = read_tsv("../src/RevisionStuff/Results/low_50.0_0.0_10.0_100_8675309.tsv")

#medium effects
dat = read_tsv("../src/RevisionStuff/Results/med_5.0_1.0_10.0_100_8675309.tsv")


### Do some analysis 
{
  dat <- dat %>% mutate(Method = fct_relevel(Method, c("SAA", "S-SAA", "NW", "HS-NW", "NWbar", "S-NW")), 
                        Label = fct_recode(Method, `NW-f`="NW", 
                                           `NW-fbar`="NWbar", 
                                           `Proposal`="S-NW"))
  
  #First compute relative values
  dat.sum <- dat %>% group_by(Label) %>%
    summarise(avg = mean(TruePerf), 
              std = sd(TruePerf))
  
  zstar = filter(dat.sum, Label=="FullInfo") %>% select(avg)
  dat$FullInfo <- zstar[[1]]
  dat <- mutate(dat, RelPerf = TruePerf/FullInfo - 1)
  dat.sum <- mutate(dat.sum, rel_avg = avg/zstar[[1]] -1 )
  
  #build custom labels based on the means
  dat.sum <- dat.sum %>% mutate(pLabel = str_c(Label, " - ", percent(rel_avg), "-BOO"),
                                pRel_avg = str_c(percent(rel_avg)))

  v <- as.character(dat.sum$Label)
  names(v) <- as.character(dat.sum$pLabel)
  dat.sum <- dat.sum %>% mutate(pLabel = fct_recode(Label, !!!v))
  dat <- dat %>% mutate(pLabel = fct_recode(Label, !!!v))
}




latex_labels <- list("S-SAA", "$NW-\alpha$", "$NW-\\mathbf f$", "Proposal")
(
  
  g<- dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal")) %>%
    ggplot(aes(Label, RelPerf, color=Label)) +
    geom_boxplot()+ 
    geom_text(aes(Label, label=pRel_avg), y=.0, 
              size=2.5, color="black", 
              data=filter(dat.sum, Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal"))) +
    theme_minimal(base_size = 8) +
    xlab("") + ylab("(%) Loss") + scale_y_continuous(labels=scales::percent) + 
    scale_x_discrete(labels=latex_labels) + 
    ggtitle("Informative Features /\nLow-Idiosyncratic Effect") + 
    theme(legend.title=element_blank(),
          legend.position="none",
          plot.margin=grid::unit(c(.5,0,0,0), "mm"), 
          plot.title = element_text(size = 8, hjust=.5))     +
    coord_cartesian(ylim=c(-.005, .20))
)

tikzDevice::tikz(file = "~/Dropbox/NSF Career 2019/V2/LowIdiosyncracy_play.tex", height=2.16, width=2.16)
g
dev.off()

ggsave("~/Dropbox/NSF Career 2019/V2/Figures/LowIdiosyncracy_play.pdf", g, height=2.16, width=2.16, units="in")






#Low idiosyncracy/ Informative FEatures
# latex_labels <- list("S-SAA", "$NW-\\mathbf f$", "$NW-\\mathbf f$", "Proposal")
# (
#   
#   g<- dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal")) %>%
#     ggplot(aes(Label, RelPerf, color=Label)) +
#     geom_boxplot()+ 
#     geom_text(aes(Label, label=pRel_avg), y=.0, 
#               family="Times New Roman", size=2.5, color="black", 
#               data=filter(dat.sum, Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal"))) +
#     theme_minimal(base_size = 8, base_family = "Times New Roman") +
#     xlab("") + ylab("(%) Loss") + scale_y_continuous(labels=scales::percent) + 
#     scale_x_discrete(labels=latex_labels) + 
#     ggtitle("Informative Features /\nLow-Idiosyncratic Effect") + 
#     theme(legend.title=element_blank(),
#           legend.position="none",
#           plot.margin=grid::unit(c(.5,0,0,0), "mm"), 
#           plot.title = element_text(size = 8, family="Times New Roman", hjust=.5))     +
#     coord_cartesian(ylim=c(-.005, .20))
# )
# 
# tikzDevice::tikz(file = "~/Dropbox/NSF Career 2019/V2/LowIdiosyncracy_play.tex", height=2.16, width=2.16)
# g
# dev.off()
# 
# ggsave("~/Dropbox/NSF Career 2019/V2/Figures/LowIdiosyncracy_play.pdf", g, height=2.16, width=2.16, units="in")
# 

### for the High
(
  g<- 
    dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal")) %>%
    ggplot(aes(Label, RelPerf, color=pLabel)) +
    geom_boxplot()+ 
    geom_text(aes(Label, label=pRel_avg), y=.05, 
              family="Times New Roman", size=2.5, color="black", 
              data=filter(dat.sum, Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal"))) +
    theme_minimal(base_size = 8, base_family = "Times New Roman") +
    xlab("") + ylab("(%) Loss") + scale_y_continuous(labels=scales::percent) + 
    ggtitle("Uninformative Features /\nHigh-Idiosyncratic Effect") + 
    theme(legend.title=element_blank(),
          legend.position="none",
          plot.margin=grid::unit(c(.5,0,0,0), "mm"), 
          plot.title = element_text(size = 8, family="Times New Roman", hjust=.5))  +
    coord_cartesian(ylim=c(.045, .85))
  )  
ggsave("~/Dropbox/NSF Career 2019/V2/Figures/HighIdiosyncracy.pdf", g, height=2.16, width=2.16, units="in")


### for the Medium-
(
  g<- 
    dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal")) %>%
    ggplot(aes(Label, RelPerf, color=pLabel)) +
    geom_boxplot()+ 
    geom_text(aes(Label, label=pRel_avg), y=.04, 
              family="Times New Roman", size=2.5, color="black", 
              data=filter(dat.sum, Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal"))) +
    theme_minimal(base_size = 8, base_family = "Times New Roman") +
    xlab("") + ylab("(%) Loss") + scale_y_continuous(labels=scales::percent) + 
    ggtitle("Realistic-Setting") + 
    theme(legend.title=element_blank(),
          legend.position="none",
          plot.margin=grid::unit(c(.5,0,0,0), "mm"), 
          plot.title = element_text(size = 8, family="Times New Roman", hjust=.5))  +
     coord_cartesian(ylim=c(.03, .35))
)  
ggsave("~/Dropbox/NSF Career 2019/V2/Figures/MedIdiosyncracy.pdf", g, height=2.16, width=2.16, units="in")




