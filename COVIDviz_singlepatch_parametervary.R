#     Script visualizes results from experiment varying a single parameter 
#     at a time, outputs figure.
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

library(ggplot2); library(RColorBrewer); library(patchwork); library(scales); library(ggtext);

parm.list<-c("beta", "epsilon", "gamma", "omega", "tau_I", "b", "xi", "d", "delta_C")

for (i in 1:length(parm.list)) {
  
  ### Select parameter of interest, some other variables for data viz
  parm<-parm.list[i]
  spacer<-0.1
  ymax<-0.9
  xmin<-0.75e6
  xmax<-2.25e6
  
  ### Read in files
  df.sum.comb<-read.csv(sprintf("COVIDsim_CCS_fadeoutdat_%s.csv", parm),
                        header=TRUE,stringsAsFactors = FALSE);
  df.ccs.comb<-read.csv(sprintf("COVIDsim_CCS_ccsdat_%s.csv", parm),
                        header=TRUE,stringsAsFactors = FALSE);
  df.ccs.comb$labs<-format(df.ccs.comb$ccs.smart, big.mark=",",trim = 
                             "TRUE")
  
  
  if (parm=="beta") {
    p.group<-"beta"
    p.colour<-"as.factor(2*beta)"
    legend.title<-expression("Transmission probability (symptomatic), "*beta[0]^{I}*"  ")
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$beta))
    labs<-c("0.65", "0.75", "0.85", "0.95")
    custom.pal<-brewer.pal(length(labs), "RdBu")
  }
  
  if (parm=="b") {
    p.group<-"b"
    p.colour<-"as.factor(b)"
    legend.title<-"Birth probability, b  "
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$b))
    labs<-c("0", expression("2 x"*10^{-5}), expression("5 x"*10^{-5}), expression("2 x"*10^{-4}))
    custom.pal<-brewer.pal(length(labs), "RdBu")
  }
  
  if (parm=="omega") {
    p.group<-"omega"
    p.colour<-"as.factor(omega)"
    legend.title<-expression("Risk perception coefficient, "*omega*"  ")
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$omega))
    labs<-labs<-c(expression("4.4 x"*10^{2}), expression("4.4 x"*10^{3}), expression("4.4 x"*10^{4}), expression("4.4 x"*10^{5}))
    custom.pal<-brewer.pal(length(labs), "RdBu")
  }
  
  if (parm=="tau_I") {
    p.group<-"tau_I"
    p.colour<-"as.factor(tau_I)"
    legend.title<-expression("Symptomatic testing probability, "*tau[I]*"  ")
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$tau_I))
    labs<-c(0.1, 0.3, 0.5, 0.7)
    custom.pal<-brewer.pal(length(labs), "RdBu")
  }
  
  if (parm=="epsilon") {
    p.group<-"epsilon"
    p.colour<-"as.factor(epsilon)"
    legend.title<-expression("NPI efficacy, "*epsilon*"  ")
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$epsilon))
    labs<-c(0.5, 0.65, 0.8, 0.95)
    custom.pal<-brewer.pal(length(labs), "RdBu")
  }
  
  if (parm=="d") {
    p.group<-"d"
    p.colour<-"as.factor(d)"
    legend.title<-"Death probability, d  "
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$d))
    labs<-c("0", expression("3 x"*10^{-6}), expression("2 x"*10^{-5}), expression("2 x"*10^{-4}))
    custom.pal<-brewer.pal(length(labs), "RdBu")
    
  }
  
  if (parm=="xi") {
    p.group<-"xi"
    p.colour<-"as.factor(xi)"
    legend.title<-expression("Propotion of mass-action incidence, "*xi*"  ")
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$xi))
    labs<-c("0", expression("1 x"*10^{-5}), expression("1 x"*10^{-3}), expression("1 x"*10^{-1}))
    custom.pal<-brewer.pal(length(labs), "RdBu")
  }
  
  if (parm=="gamma") {
    p.group<-"gamma"
    p.colour<-"as.factor(gamma)"
    legend.title<-expression("Trigger prevalence, "*gamma*"  ")
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$gamma))
    labs<-c(expression("1.5 x"*10^{-5}), expression("1.5 x"*10^{-4}), expression("1.5 x"*10^{-3}), expression("1.5 x"*10^{-2}))
    custom.pal<-brewer.pal(length(labs), "RdBu")
  }
  
  if (parm=="delta_C") {
    p.group<-"delta_C"
    p.colour<-"as.factor(delta_C)"
    legend.title<-expression("Closure duration (days), "*delta[C]*"  ")
    ann.min<-ymax-spacer*length(unique(df.ccs.comb$delta_C))
    labs<-c(15, 30, 60, 120)
    custom.pal<-brewer.pal(length(labs), "RdBu")
  }
  
  ### Viz 
  
  theme.size = 10
  geom.text.size = (1/3) * theme.size
  
  breaks <- 10^(4:6)
  minor_breaks <- rep(1:9, 21)*(10^rep(4:6, each=9))
  
  p1<-ggplot(data=df.sum.comb[df.sum.comb$ts==1825,], aes_string(x="initpop", group=p.group, fill=p.colour, colour=p.colour)) +
    geom_line(aes(y=zeroevent/N), linetype="solid", key_glyph = "rect") +
    ylab(bquote(atop(bold('5 years'),"Proportion of fade-outs"))) + 
    xlab(expression("Initial population size, "*P[0])) +
    coord_cartesian(ylim=c(0,1), clip = 'off') +
    scale_colour_manual(legend.title, values=custom.pal, labels=labs, guide = guide_legend(nrow = 1)) +
    scale_fill_manual(legend.title, values=custom.pal, labels=labs, guide = guide_legend(nrow = 1)) +
    geom_text(data=df.ccs.comb[df.ccs.comb$ts==1825,], aes(x=xmin+10000, y = ymax+0.025, label = "CCS:"), colour="black", hjust = 0,  size=geom.text.size) +
    geom_text(data=df.ccs.comb[df.ccs.comb$ts==1825,], aes_string(x="xmin+10000", y = "rev(seq(ann.min,ymax-spacer,spacer))", label = "labs", colour=p.colour), hjust = 0, size=geom.text.size) +
    scale_x_log10(labels=comma, expand=c(0,0), breaks=breaks, minor_breaks=minor_breaks) +
    scale_y_continuous(expand=c(0.01,0)) +
    theme_light() + 
    guides(fill=FALSE, colour=FALSE) +
    theme(legend.position = 'right', plot.margin = margin(2, 3, 2, 2, "pt"))
  
  
  p2<-ggplot(data=df.sum.comb[df.sum.comb$ts==365,], aes_string(x="initpop", group=p.group, fill=p.colour, colour=p.colour)) +
    geom_line(aes(y=zeroevent/N), linetype="solid", key_glyph = "rect") +
    ylab(bquote(atop(bold('1 year'),"Proportion of fade-outs"))) + 
    xlab(expression("Initial population size, "*P[0])) +
    coord_cartesian(ylim=c(0,1), clip = 'off') +
    scale_colour_manual(legend.title, values=custom.pal, labels=labs, guide = guide_legend(nrow = 1)) +
    scale_fill_manual(legend.title, values=custom.pal, labels=labs, guide = guide_legend(nrow = 1)) +
    geom_text(data=df.ccs.comb[df.ccs.comb$ts==365,], aes(x=xmin+10000, y = ymax+0.025, label = "CCS:"), colour="black", hjust = 0,  size=geom.text.size) +
    geom_text(data=df.ccs.comb[df.ccs.comb$ts==365,], aes_string(x="xmin+10000", y = "rev(seq(ann.min,ymax-spacer,spacer))", label = "labs", colour=p.colour), hjust = 0, size=geom.text.size) +
    scale_x_log10(labels=comma, expand=c(0,0), breaks=breaks, minor_breaks=minor_breaks) +
    scale_y_continuous(expand=c(0.01,0)) +
    theme_light() + 
    guides(fill=FALSE, colour=FALSE) +
    theme(legend.position = 'right', plot.margin = margin(2, 3, 2, 2, "pt"))
  
  
  p.inf.1<-ggplot(df.sum.comb[df.sum.comb$ts==1825,], aes_string(x="initpop", y="mean.propinf.boot", ymin="min.propinf", ymax="max.propinf", colour=p.colour, fill=p.colour)) +
    scale_colour_manual(legend.title, values=custom.pal, labels=labs, guide = guide_legend(nrow = 1)) +
    scale_fill_manual(legend.title, values=custom.pal, labels=labs, guide = guide_legend(nrow = 1)) +
    geom_ribbon(alpha=0.1, size=0.25, key_glyph = "rect") +
    geom_line(size=0.25, linetype="dashed", key_glyph = "rect") +
    scale_x_log10(labels=comma, expand=c(0,0), breaks=breaks, minor_breaks=minor_breaks) +
    scale_y_continuous(expand=c(0.01,0), limits=c(0,1)) +
    xlab(expression("Initial population size, "*P[0])) +
    ylab("Proportion infected") +
    coord_cartesian(ylim=c(0,1), clip = 'off') +
    theme_light() +
    theme(plot.title=element_text(hjust=0.5))
  
  p.inf.2<-ggplot(df.sum.comb[df.sum.comb$ts==365,], aes_string(x="initpop", y="mean.propinf.boot", ymin="min.propinf", ymax="max.propinf", colour=p.colour, fill=p.colour)) +
    scale_colour_manual(legend.title, values=custom.pal, labels=labs, guide = guide_legend(nrow = 1)) +
    scale_fill_manual(legend.title, values=custom.pal, labels=labs, guide = guide_legend(nrow = 1)) +
    scale_x_log10(labels=comma, expand=c(0,0), breaks=breaks, minor_breaks=minor_breaks) +
    scale_y_continuous(expand=c(0.01,0), limits=c(0,1)) +
    xlab(expression("Initial population size, "*P[0])) +
    ylab("Proportion infected") +
    coord_cartesian(ylim=c(0,1), clip = 'off') +
    geom_ribbon(alpha=0.1, size=0.25, key_glyph = "rect") +
    geom_line(size=0.25, linetype="dashed", key_glyph = "rect") +
    theme_light()+
    theme(plot.title=element_text(hjust=0.5))
  
  
  
  p.comb.total <- wrap_plots(p2,  p.inf.2, p1, p.inf.1, nrow = 2) +
    plot_layout(guides = 'collect', heights = c(1, 1)) +
    plot_annotation(tag_levels = 'a') & theme(legend.position = "bottom",
                                              legend.background = element_blank(),
                                              legend.margin = margin(0, 0, 0, 0, "pt"),
                                              plot.margin = margin(1, 1, 1, 1, "pt"))
  
  
  png(sprintf("COVIDsim_parametervary_%s.png", parm), width=20, height=10, units="cm", res=500)
  print(p.comb.total)
  dev.off()
  
}
