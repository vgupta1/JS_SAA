library(stringr)
library(scales)
library(tidyverse)
library(ggplot2)
library(forcats)

#Tikz does all the font stuff for you.  so dn't worry about it.  
#just worry about size and makign sure special characters are escaped correctly
# This is a much better way to do this going forwards!!!

##Some "Gotchyas"  Be careful with labels and unintentional special characters, e.g., Loss (%)  needs to be Loss (\%)
## Building on above, if you have a percentage scale, can't simply use scales::percent.  Instead pass formating function directly.  labels= function(x) str_c(x*100, "\\%")
## Can't use macros you defined elsewhere

library(tikzDevice)

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
                                pRel_avg = percent_format(suffix="\\%")(rel_avg))
  
  v <- as.character(dat.sum$Label)
  names(v) <- as.character(dat.sum$pLabel)
  dat.sum <- dat.sum %>% mutate(pLabel = fct_recode(Label, !!!v))
  dat <- dat %>% mutate(pLabel = fct_recode(Label, !!!v))
}

latex_labels <- list("S-SAA", "$\\sf NW-\\mathbf f$", "$\\sf NW-\\mathbf{\\overline{f}}$", "Proposal")

#Creata a really simple plot first to mess with
g<- dat %>% filter(Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal")) %>%
  ggplot(aes(Label, RelPerf, color=Label)) +
  geom_boxplot()+ 
  theme_minimal(base_size = 8) +
  xlab("") + ylab("(\\%) Loss") + 
  scale_y_continuous(labels= percent_format(suffix="\\%")) + 
  scale_x_discrete(labels=latex_labels) +
  ggtitle("Informative Features /\nLow-Idiosyncratic Effect") + 
  geom_text(aes(Label, label=pRel_avg), y=.0, 
              size=2.5, color="black", 
              data=filter(dat.sum, Label %in% c("S-SAA", "NW-f", "NW-fbar", "Proposal"))) +
    theme(legend.title=element_blank(),
        legend.position="none",
        plot.margin=grid::unit(c(.5,0,0,.5), "mm"), 
        plot.title = element_text(size = 8, hjust=.5))     +
  coord_cartesian(ylim=c(-.005, .20))

#Send it to tikz
tikzDevice::tikz(file = "~/Dropbox/NSF Career 2019/V2_newnotation/tikz_play.tex", height=2.16, width=2.16)
g
dev.off()

