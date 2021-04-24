#     Script visualizes results from experiment on the impact of 
#     accounting for births/deaths, outputs figure.
#     Copyright (C) 2021  Kathyrn R Fair
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

library(ggplot2); library(patchwork); library(RColorBrewer); library(gmodels); library(dplyr);

### Read in files
tags<-c("NOnpi", "npi", "NObirthdeath", "birthdeath")

for (i in 1:2) {
  for (j in 3:4) {
    tmp<-read.csv(sprintf("COVIDsim_birthdeathexp_%s_%s.csv", tags[i], tags[j]))
    tmp<-tmp[,colnames(tmp) %in% c("ts","infected", "initpop")]
 
    tmp$tag1<-tags[i]
    tmp$tag2<-tags[j]
    if (i==1 && j==3) {
      df<-tmp
    } else {
      df<-rbind(df, tmp)
    }
  }
}

df.sum <- df %>% group_by(initpop,ts,tag1,tag2) %>% summarise(mean = mean(infected), min=min(infected), max=max(infected))

labs<-c("30,000", "50,000", "100,000", "1,000,000")

### Viz

p.no <- ggplot(data=df.sum[(df.sum$tag1=="NOnpi"  & df.sum$ts<121),], aes(x=ts, group=interaction(tag1,tag2,initpop), colour=as.factor(initpop), fill=as.factor(initpop), linetype=tag2)) +
  geom_line(aes(y=log(mean))) +
  geom_ribbon( aes(ymin=log(min), ymax=log(max)), alpha=0.1, colour=NA) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,120,30)) +
  scale_colour_brewer(expression(P[0]), labels=labs, palette="RdBu", direction=-1, guide = guide_legend(title.position = "top", nrow = 1)) +
  scale_fill_brewer(expression(P[0]), labels=labs, palette="RdBu", direction=-1, guide = guide_legend(title.position = "top", nrow = 1)) +
  scale_linetype_manual("Accounting for births, non-COVID-19 deaths?", labels=c("Yes", "No"), values=c("solid", "longdash"), guide = guide_legend(title.position = "top", nrow = 1)) + 
  xlab("Day") +
  ylab("Log # infected") +
  ggtitle("No NPIs") +
  theme_light()

p.yes <- ggplot(data=df.sum[(df.sum$tag1=="npi" & df.sum$ts<10*365),], aes(x=ts/365, group=interaction(tag1,tag2,initpop), colour=as.factor(initpop), fill=as.factor(initpop), linetype=tag2)) +
  geom_line(aes(y=log(mean))) +
  geom_ribbon( aes(ymin=log(min), ymax=log(max)), alpha=0.1, colour=NA) +
  scale_x_continuous(expand=c(0,0),breaks=seq(0,10,1)) +
  scale_colour_brewer(expression(P[0]), labels=labs, palette="RdBu", direction=-1, guide = guide_legend(title.position = "top", nrow = 1)) +
  scale_fill_brewer(expression(P[0]), labels=labs,palette="RdBu", direction=-1, guide = guide_legend(title.position = "top", nrow = 1)) +
  scale_linetype_manual("Accounting for births, non-COVID-19 deaths?", labels=c("Yes", "No"), values=c("solid", "longdash"), guide = guide_legend(title.position = "top", nrow = 1)) + 
  xlab("Year") +
  ylab("Log # infected") +
  ggtitle("Baseline NPIs") +
  theme_light()

p.comb<-wrap_plots(p.no, p.yes, nrow = 2) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a') & theme(legend.direction = "vertical", 
                                            legend.position = "bottom",
                                            legend.box = "horizontal", 
                                            legend.background = element_blank(),
                                            legend.box.background = element_rect(colour =NA, fill = NA, 
                                                                                 linetype='solid'),
                                            plot.margin = margin(0, 7, 2, 0, "pt"))



png("COVIDsim_birthdeath.png", width=20, height=20, units="cm", res=500)
print(p.comb)
dev.off()