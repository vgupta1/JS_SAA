library(stringr)
library(tidyverse)
library(ggplot2)
library(forcats)

library(tikzDevice)

latex_labels <- list("S-SAA", "\\sf NW-$\\mathbf f$", "\\sf NW-$\\mathbf{f}$-ID", "Proposal")

#high idiosyncratic attempt  - VG Complete!
dat = read_tsv("../src/RevisionStuff/Results/highIdioResults_20_5_1.0_100_8675309.tsv")
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
  dat.sum <- dat.sum %>% mutate(pRel_avg = scales::percent_format(suffix="\\%")(rel_avg),
                                pLabel = str_c(Label, " - ", pRel_avg)
                                )
  
  v <- as.character(dat.sum$Label)
  names(v) <- as.character(dat.sum$pLabel)
  dat.sum <- dat.sum %>% mutate(pLabel = fct_recode(Label, !!!v))
  dat <- dat %>% mutate(pLabel = fct_recode(Label, !!!v))
}

##Build graph
{
  g<- 
    dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal")) %>%
    ggplot(aes(Label, RelPerf, color=pLabel)) +
    geom_boxplot()+ 
    geom_text(aes(Label, label=pRel_avg), y=.05, 
              size=2.5, color="black", 
              data=filter(dat.sum, Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal"))) +
    theme_minimal(base_size = 7) +
    xlab("") + ylab("(\\%) Loss") + scale_y_continuous(labels=scales::percent_format(suffix="\\%")) +
    scale_x_discrete(labels=latex_labels) + 
    ggtitle("Uninformative Features /\nHigh-Idiosyncratic Effect") + 
    theme(legend.title=element_blank(),
          legend.position="none",
          plot.margin=grid::unit(c(.5,0,0,.5), "mm"), 
          plot.title = element_text(size = 7, hjust=.5))  +
    coord_cartesian(ylim=c(.045, .85))
  
}
#Send it to tikz
tikzDevice::tikz(file = "~/Dropbox/NSF Career 2019/V2_newnotation/highIdiosyncracy.tex", height=2.16, width=2.16)
g
dev.off()


#Low idiosyncratic effect attempts
dat = read_tsv("../src/RevisionStuff/Results/low_50.0_0.0_10.0_100_8675309.tsv")
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
  dat.sum <- dat.sum %>% mutate(pRel_avg = scales::percent_format(suffix="\\%")(rel_avg),
                                pLabel = str_c(Label, " - ", pRel_avg)
  )
  
  v <- as.character(dat.sum$Label)
  names(v) <- as.character(dat.sum$pLabel)
  dat.sum <- dat.sum %>% mutate(pLabel = fct_recode(Label, !!!v))
  dat <- dat %>% mutate(pLabel = fct_recode(Label, !!!v))
}
##Build graph
{
  g<- dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal")) %>%
    ggplot(aes(Label, RelPerf, color=Label)) +
    geom_boxplot()+ 
    geom_text(aes(Label, label=pRel_avg), y=.0, 
              size=2.5, color="black", 
              data=filter(dat.sum, Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal"))) +
    theme_minimal(base_size = 7) +
    xlab("") + ylab("(\\%) Loss") + scale_y_continuous(labels= scales::percent_format(suffix="\\%")) +
    scale_x_discrete(labels=latex_labels) + 
    ggtitle("Informative Features /\nLow-Idiosyncratic Effect") + 
    theme(legend.title=element_blank(),
          legend.position="none",
          plot.margin=grid::unit(c(.5,0,0,.5), "mm"), 
          plot.title = element_text(size = 7, hjust=.5))     +
    coord_cartesian(ylim=c(-.005, .20))
}
#Send it to tikz
tikzDevice::tikz(file = "~/Dropbox/NSF Career 2019/V2_newnotation/LowIdiosyncracy.tex", height=2.16, width=2.16)
g
dev.off()







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
  dat.sum <- dat.sum %>% mutate(pRel_avg = scales::percent_format(suffix="\\%")(rel_avg),
                                pLabel = str_c(Label, " - ", pRel_avg)
  )
  
  v <- as.character(dat.sum$Label)
  names(v) <- as.character(dat.sum$pLabel)
  dat.sum <- dat.sum %>% mutate(pLabel = fct_recode(Label, !!!v))
  dat <- dat %>% mutate(pLabel = fct_recode(Label, !!!v))
}
##Build graph
{
  g<- 
    dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal")) %>%
    ggplot(aes(Label, RelPerf, color=pLabel)) +
    geom_boxplot()+ 
    geom_text(aes(Label, label=pRel_avg), y=.04, 
              size=2.5, color="black", 
              data=filter(dat.sum, Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal"))) +
    theme_minimal(base_size = 7) +
    xlab("") + ylab("(\\%) Loss") + scale_y_continuous(labels=scales::percent_format(suffix="\\%")) + 
    scale_x_discrete(labels=latex_labels) + 
    ggtitle("Realistic-Setting") + 
    theme(legend.title=element_blank(),
          legend.position="none",
          plot.margin=grid::unit(c(.5,0,0,.5), "mm"), 
          plot.title = element_text(size = 7, hjust=.5))  +
    coord_cartesian(ylim=c(.03, .35))
  
}
#Send it to tikz
tikzDevice::tikz(file = "~/Dropbox/NSF Career 2019/V2_newnotation/MedIdiosyncracy.tex", height=2.16, width=2.16)
g
dev.off()










