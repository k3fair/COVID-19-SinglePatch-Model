#     Script runs simulations for experiment varying a single parameter at 
#     a time, outputs results of these simulations.
#     Copyright (C) 2021  Kathyrn R Fair, Vadim A Karatayev
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

library(abind); library(fields); library(Hmisc); library(dplyr);
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

tstart<-Sys.time()

SpatiotempPlot=function (spacetime, Xax = 1:dim(spacetime)[2], Yax = 1:dim(spacetime)[1], 
                         XaxN = "X", YaxN = "Y", figtitle = "Title", Zlim = c(0, max(spacetime, 
                                                                                     na.rm = TRUE)), cont = NULL, cexAx = 1, contPlot = spacetime, 
                         cexCont = 1.5 * cexAx, Conts = NULL, contSpan = 1, palette = 1) 
{
  require(fields)
  spacetime[is.na(spacetime)] = Zlim[1] - 1
  COL = rev(rainbow(1000, start = 0, end = 0.7))
  if (palette > 1) {
    require(pals)
    COL = list(parula(1000), head(tail(parula(1000), -50, -50)))[[palette - 1]]
  }
  if (length(Zlim) == 1) {
    Zlim = quantile(spacetime, c(Zlim, 1 - Zlim), na.rm = T)
    Rng = range(spacetime)
    if (Zlim[1] < Rng[1]) 
      Zlim[1] = Rng[1]
    if (Zlim[2] > Rng[2]) 
      Zlim[2] = Rng[2]
  }
  spacetime[which(is.na(spacetime), arr.ind = TRUE)] = max(Zlim) + 1
  image.plot(x = Xax, y = Yax, z = t(spacetime), zlim = Zlim, 
             xlab = XaxN, ylab = YaxN, cex.axis = cexAx, cex.lab = cexAx, 
             legend.cex = cexAx, main = figtitle, col = COL)
  box()
  if (!is.null(cont)) {
    if (abs(log(max(contPlot, na.rm = T), 10)) > 2) 
      options(scipen = -10)
    if (contSpan > 1) {
      smoo1 = t(apply(contPlot, 1, function(x) supsmu(1:ncol(contPlot), 
                                                      x, span = contSpan/ncol(contPlot))$y))
      smoo2 = t(apply(smoo1, 1, function(x) supsmu(1:ncol(smoo1), 
                                                   x, span = contSpan/ncol(smoo1))$y))
      contPlot = smoo2
    }
    if (is.null(Conts)) 
      contour(x = Xax, y = Yax, z = t(contPlot), add = T, 
              col = cont, lwd = 1.5 * cexCont, labcex = cexCont)
    if (!is.null(Conts)) 
      contour(x = Xax, y = Yax, z = t(contPlot), levels = Conts, 
              add = T, col = cont, lwd = 1.5 * cexCont, labcex = cexCont)
    if (abs(log(max(contPlot, na.rm = T), 10)) > 2) 
      options(scipen = 0)
  }
}

#Notes:
# one parameter renamed here, we have xi=dlt
#pi here (not in parms matrix, fixed at 0.2) is fraction of people who never develop symptoms. So asymptomatic stage is a bit longer and a fraction pi go A->R
#Mfact scales how much travel happens compared to reality. Mfact=1 is just movement from commuting data but really other people travel too, so baseline travel is Mfact=1.5 (irrelevant here because we don't include travel)
#epsP allows voluntary distancing efficacy to differ from that of closures (epsP=eps in paper)
#Tg and Tl are the closure thresholds, n0 is fraction initially infected
#Msave is the travel matrix (currently set so there is no connectivity between patches)

#Take cmd line input to determine which parameter is being varied (or, to run directly in Rstudio/R comment out the args and manually select a variable)
var.list<-c("beta", "epsilon", "gamma", "omega", "tau_I", "b", "xi", "d", "delta_C")
args <- commandArgs(trailingOnly = TRUE)
#Select with parameter you want to vary (xi, beta, epsilon, gamma, omega, tau_I, b, d, delta_C)
select.var<-var.list[as.numeric(args[1])]
print(select.var)

### GENERIC SETUP
years<- 5 
samples<-500
popsample<-c(1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,1e5,1.5e5,2e5,2.5e5,3e5,3.5e5,4e5,4.5e5,5e5,5.5e5,6e5,6.5e5,7e5,7.5e5,8e5,8.5e5,9e5,9.5e5,1e6,1.1e6,1.2e6,1.3e6,1.4e6,1.5e6,1.6e6,1.7e6,1.8e6,1.9e6,2e6,2.1e6,2.2e6,2.3e6,2.4e6,2.5e6)

if (select.var=="beta") var.sample<- c(0.325,0.375,0.425,0.475)
if (select.var=="epsilon") var.sample<- c(0.5,0.65,0.8,0.95);
if (select.var=="gamma") var.sample<- c(1.5e-5, 1.5e-4, 1.5e-3,1.5e-2)
if (select.var=="omega") var.sample<-c(4.4e2, 4.4e3, 4.4e4, 4.4e5)
if (select.var=="tau_I") var.sample<-seq(0.1,0.7,0.2)
if (select.var=="b") var.sample<-c(0, 2e-5, 5e-5, 2e-4)
if (select.var=="xi") var.sample<-c(0, 1e-5, 1e-3, 1e-1)
if (select.var=="d") var.sample<-c(0, 3e-6, 2e-5, 2e-4)
if (select.var=="delta_C") var.sample<-c(15,30,60,120)

counter<-0
for (z1 in 1:length(var.sample)) 
{
  
  for (z2 in 1:length(popsample))
  {
    
    pop<-popsample[z2]

    parms=cbind(N=rep(pop, samples),births=0.000051084, deaths=0.000020737, coviddeaths=0.0066, Mfact=1.5,s=0.2,Tg=1,Tl=1.5e-4,beta=0.43,tauI=0.3,tauEAf=0,epsP=0.85,eps=0.85,w=0.45,omg=44356,r=0.19,eta=0.8,dlt=0,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,M=diag(pop, samples, samples));

    inClen=30*3 #Duration of initial closure
    Resurge=30 #Duration of additional closures
    
    if (select.var=="beta") parms[,"beta"]<- var.sample[z1]
    if (select.var=="epsilon") {parms[,"eps"]<- var.sample[z1]; parms[,"epsP"]<-var.sample[z1];}
    if (select.var=="gamma") parms[,"Tl"] <- var.sample[z1];
    if (select.var=="omega") parms[,"omg"]<- var.sample[z1]
    if (select.var=="tau_I") parms[,"tauI"] <- var.sample[z1]
    if (select.var=="b") parms[,"births"]<- var.sample[z1]
    if (select.var=="xi") parms[,"dlt"] <- var.sample[z1]
    if (select.var=="d") parms[,"deaths"] <- var.sample[z1]
    if (select.var=="delta_C") Resurge<-var.sample[z1];
    
    inCstart <- parms[1,"Tl"]; 

    NP=ncol(parms)-nrow(parms)
    
    #Defining model state space. Tn, Tk, Da, and Di are all untested, tested, asymptomatics, and infecteds.
    #S prefix indicates superspreaders, k suffix indicates an individual who has tested positive
    #Rd_covid and Rd are removed individuals who died (due to COVID-19 and non-COVID-19 related reasons, respectively)
    #Nt tracks cumulative # positive cases (including those recovered)
    #C tracks # days since last closure, but in output msim converts all positive C entries into eps, Ci is for the intial closure, Ct combines the two
    #CN is newly exposed individuals
    Tn=c("A","SA","I","SI"); Tk=c("Ak","SAk","Ik","SIk"); Da=c("A","Ak","SA","SAk"); Di=c("I","Ik","SI","SIk"); Snm=c("S","E",Da,Di,"R", "Rd_covid", "Rd", "Nt", "C", "Ci", "Ct", "B", "D", "CN");
    
    
    #these functions handle all the state transitions. Input x is a matrix where 1st half of columns are source states and 2nd half are destination states
    gpTransB=function(x,Prs,seed,nc=ncol(x)){ 
      xvp=cbind(as.vector(x[,1:(nc/2)]),as.vector(Prs))
      if(max(xvp[,1])==0) return(x);  #if no untested, nobody to move over to tested, quit
      nz=(xvp[,1]*xvp[,2])>0; xvp[!nz,1]=0; #set any sample sizes where the prob is zero to zero so we don't need to run a binom for them
      set.seed(seed); xvp[nz,1]=apply(matrix(xvp[nz,],ncol=2),1,function(y) rbinom(1,y[1],y[2]));#run binom to figure out transfers
      return(x+matrix(c(-xvp[,1],xvp[,1]),ncol=nc))
    }
    #gpTrans is a simplified version where one transition probability for all states. 
    #If recovery=TRUE have individuals from 1st columns of x all transitioning into the last column
    gpTrans=function(x,Pr,seed,Recovery=FALSE, nc=ncol(x)){
      xv=as.vector(x[,1:c(nc/2,nc-1)[1+Recovery]])
      if(max(xv)==0) return(x); #if nobody to potentially transfer, quit
      set.seed(seed); xv[xv>0]=sapply(xv[xv>0], function(y) rbinom(1,y,Pr)); #decides how many from each compartment get shifted within each patch
      if(Recovery){ Tr=matrix(xv,ncol=nc-1); return(x+cbind(-Tr,rowSums(Tr))); };  #If for recovery, all get shifted to recovered
      return(x+matrix(c(-xv,xv),ncol=nc)); #Not for recovery so they get shifted out of first 1/2 of matrix into 2nd
    }
    
    #Transition probabilities for travel. Multinomial faster than many binomials. rs=TRUE returns just total # people going to each province
    reshfl2=function(x,M,seed,rs=TRUE,L=ncol(M),xnm=diag(x)){
      #if there's no one to transfer just leave as is, otherwise calc transport
      set.seed(seed); if(max(x)>0) xnm[,x>0]=apply(matrix(rbind(x,M)[,x>0],nrow=L+1), 2, function(y) rmultinom(1,y[1],y[-1])); #use columns where top row is number in comparment, rows below are probs of staying in/leaving home patch
      if(rs) return(rowSums(xnm)); return(xnm); 
    }
    
    #birth function, input total population and output additions to succeptibles
    gpBirth=function(x,Pr,seed,nc=ncol(x)){
      xv=as.vector(x)
      if(max(xv)==0) return(x); #if pop is zero, just quit bc no births
      set.seed(seed);
      #decides how many births to add to succeptible compartment for each patch
      xv[xv>0]=sapply(xv[xv>0], function(y) rbinom(1,y,Pr));
      #Add births to compartment in each patch
      Tr=matrix(xv); return(Tr); 
    }
    
    #death function, input information on a compartment, outputs deaths for that compartment across patches
    gpDeath=function(x,Pr,seed,nc=ncol(x)){
      xv=as.vector(x)
      if(max(xv)==0) return(x); #if no one in compartment, quit
      set.seed(seed);
      #decides how many deaths to remove from compartment for each patch
      xv[xv>0]=sapply(xv[xv>0], function(y) rbinom(1,y,Pr));
      #Remove deaths from the compartment in each patch
      Tr=matrix(xv); return(x-Tr); 
    }
    
    
    #Modifier of travel matrix as people sick and/or tested positive less likely to travel by a proportion pstay
    Mstay=function(M,pstay,Mod=M*(-(diag(ncol(M))-1)*(1-pstay))) Mod+diag(1-colSums(Mod))
    
    meansd=function(x,wtR=1+x*0,wts=wtR/sum(wtR)){ xm=weighted.mean(x,wts); return(c(xm,sqrt(sum(wts*(x-xm)^2)))); }
    
    #Main function handling all within-day transitions over state space S
    m3iter=function(S,parms,seed,Ns=parms[,"N"]){
      
      #Introduce births and deaths
      #births
      Bk=S[,c("S")];
      S[,c("S")]=S[,c("S")] + gpBirth(rowSums(S[,c("S","E",Da,Di,"R")]),parms[1,"births"],seed+27);
      S[,c("B")] = S[,c("S")]-Bk;
      # non-covid deaths
      Nold=rowSums(S[,c("R","E",Da,Di, "S")]);
      Rold=S[,"R"]
      S[,c("R","E",Da,Di, "S")]=apply(rbind(seed+(16:26), S[,c("R","E",Da,Di, "S")]), 2, function(x) gpDeath(x[-(1)],parms[1,"deaths"],x[1]))
      S[,"Rd"]=Rold-S[,"R"]; #Calculate how many individuals who recovered from COVID-19 died of unrelated causes at this timestep
      S[,c("D")] = Nold - rowSums(S[,c("R","E",Da,Di, "S")]); #Calculate total (non-covid) deaths at this timestep
      #Count up number of people in each patch after biths and deaths, to use for thresholds, etc
      Ns=rowSums(S[,c("S","E",Da,Di,"R")])
      
      
      #Record all previously known tests AFTER any birth/death - this avoids underestimating number of new cases, as if we record the previously known before death we're getting a "net" new including deaths, not actual new
      S0k=S[,Tk]; #These are positive tests we already knew about
      #Then testing results come back
      S[,c(Tn,Tk)]=gpTransB(S[,c(Tn,Tk)],as.vector(parms[,"tauI"]%*%t(c(parms[1,c("tauEAf","tauEAf")],1,1))),seed+3); 
      #update cumulative total, and local active case prevlanence based on our new testing
      S[,"Nt"]=S[,"Nt"]+rowSums(S[,Tk]-S0k); pos=rowSums(S[,Tk])/Ns; 
      
      #Disease progression: pi is the fraction of people who never show symptoms (modeled implicitly)
      pi=0.2; 
      ## Removal without death from I
      S[,c(Di,"R")]=gpTrans(S[,c(Di,"R")],parms[1,"rho"]*(1-parms[1,"coviddeaths"]),seed,TRUE); 
      ## Removal via death from I
      S[,c(Di,"Rd_covid")]=gpTrans(S[,c(Di,"Rd_covid")],parms[1,"rho"]*parms[1,"coviddeaths"],seed,TRUE); 
      ## All other transitions
      S[,c(Da,"R")]=gpTrans(S[,c(Da,"R")],pi*prod(parms[1,c("sig","rho")]),seed,TRUE); 
      S[,c(Da,Di)]=gpTrans(S[,c(Da,Di)],parms[1,"sig"]*(1-pi),seed+1); 
      S[,c("E","A")]=gpTrans(S[,c("E","A")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+2); S[,c("E","SA")]=gpTrans(S[,c("E","SA")],parms[1,"alpha"]*parms[1,"s"],seed+2.5);
      
      
      #Closure decisions: C tracks # days dC since last closure (ie units 1/eps rather than integers).
      ### Replace dC, dCt with -1 to get rid of closures
      dC=round(S[,"C"]/parms[1,"eps"]);
      dCt=round(S[,"Ct"]/parms[1,"eps"]);
      past_init<-(dCt>=inClen); in_init <-((S[,"Nt"]/Ns >= inCstart & dCt==0) | (dCt<inClen & dCt>0));
      decide= past_init & (dC==0 | dC>Resurge-1); #determines if we're past initial and any subsequent closures
      closed=parms[1,"eps"]*((pos > parms[,"Tl"]*decide & parms[,"Tl"]<1 & past_init==TRUE) | in_init | (pos==0 & decide==FALSE & parms[,"Tl"]<1 & past_init==TRUE));
      S[,"C"]=S[,"C"]*(!decide)+(past_init*closed); 
      S[,"Ci"]=in_init*closed; 
      S[,"Ct"]=S[,"Ct"] + closed;
      
      ###Note: section below does calculations for travel, as this was included in a previous version of the model, however travel matrix ("M" in the "parms" data frame) has been modified to eliminate travel
      #Make modified travel matrices for people feeling sick and/or tested positive, then implement travel
      M=Mc=parms[,-(1:NP)]; Mc=M[1:nrow(M),]*(1-closed); diag(Mc)=diag(Mc)+1-colSums(Mc); McA=abind(Mc,Mstay(Mc,parms[1,"eta"]),Mstay(Mc,parms[1,"r"]),Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])),along=3);
      Ss=reshfl2(S[,"S"],Mc,seed+4,FALSE); Rearr=apply(rbind(seed+(5:14), c(1,1,1:2,1:2,3:4,3:4), S[,c("R","E",Da,Di)]), 2, function(x) reshfl2(x[-(1:2)],McA[,,x[2]],x[1]))
      Pstar=rowSums(cbind(Ss,Rearr)); 
      
      #Calculate impact of populatin size on transmission probability
      dlt=parms[1,"dlt"]; 
      Pmod= dlt + (1-dlt)/Pstar;

      #Baseline infection probability. The last term with rowSums is how many people end up in each county for the day.
      Infect = parms[1,"beta"]*(parms[1,"w"]*(1-closed) + (1-parms[1,"w"])*(1-parms[1,"epsP"]*(1-exp(-parms[1,"omg"]*pos)))) * Pmod
      #Class-specific infection modifiers. Net infection Pr is then 1-Pr(avoid infection by anyone).
      modK=1-parms[1,"eta"]; modS=1/parms[1,"s"]-1; modA=cbind(Infect,modK*Infect,modS*Infect,modS*modK*Infect);
      Infects=1 - apply((1-cbind(modA,2*modA))^Rearr[,-(1:2)], 1, prod)
      #Implement infection and move susceptibles and newly exposeds back to home county
      
      #Record all previously exposed
      C0=S[,c("E")];
      S[,c("S","E")]=cbind(0,S[,"E"]) + colSums(gpTransB(cbind(Ss,0*Ss),round(Infects,10),seed+15))
      #Update with newly exposed
      S[,"CN"]=S[,c("E")]-C0;
      
          return(S)
    }; FUN=m3iter
    
    closein=function(threshs,tsNt,tm){ 
      return(threshs) #Initial closure over
    }
    #Implements decline in testing time from initial 10 days, and allows parameters (if specified in retest) to change after re-opening.
    closeinappl=function(parms,TS,tm=dim(TS)[3],delayInit=10){
      Nta=colSums(t(t(TS[,"Nt",]))); parms[,c("Tl","Tg")]=closein(parms[,c("Tl","Tg")],Nta,tm);
      return(parms);
    }
    Resagg=function(TS,parms,plotgive=TRUE){  tmp=t(apply(TS,2:3,sum))/sum(parms[,"N"]); tmp2=cbind(tmp[,"R"],rowSums(tmp[,2:10])*1e2,rowSums(tmp[,Tk])*1000,tmp[,"C"]); if(plotgive) matplot(tmp2,type="l"); return(tmp2); }
    #Function to implement simulations. InitInf sets which stages the initially sick people are in, InitInf=2 default is all sick initially exposed
    msim3=function(FUN,parms,Trun=365,seed0=11,plotgive=FALSE,InitInf=c(2)){
      L=nrow(parms); nI=length(InitInf); Ns=parms[,"N"]; 
      #assign initial infections
      set.seed(seed0); Infs=NI0=apply(parms[,c("N","n0")], 1, function(x) rbinom(n=nI,x[1],prob=x[2]/nI)); if(nI>1){ Infs=t(Infs); NI0=rowSums(Infs); };
      
      TS=array(0,dim=c(L,length(Snm),1)); TS[,c(1,InitInf),1]=cbind(Ns-NI0,Infs); TS[,13,1]=0; colnames(TS)=Snm; set.seed(seed0+2); Seed=rnorm(Trun,1e6,1e6);
      Mn=round(parms[1,"Mfact"]*parms[,-(1:NP)]*(1-diag(L))); diag(Mn)=Ns-colSums(Mn); parms[,-(1:NP)]=Mn%*%diag(1/colSums(Mn));
      for(i in 2:Trun) TS=abind(TS,FUN(TS[,,i-1],closeinappl(parms,TS),seed=Seed[i]),along=3); TS[,"C",]=parms[1,"eps"]*(TS[,"C",]>0);
      #Different levels of aggregation in model output. plotgive=3 or 3.5 give shortest output form (tracking only infections and costs) in integer format to reduce output size
      if(plotgive=="TS") return(TS); if(plotgive==TRUE) return(Resagg(TS,parms,plotgive==1));
      if(plotgive%in%c(3,3.5)){ TS[,"C",]=TS[,"C",]*matrix(Ns,nrow(parms),Trun); TSn=apply(TS,c(2,3),sum); out=matrix(as.integer(TSn),nrow(TSn),ncol(TSn)); if(plotgive==3.5) return(out[c(1,13),]); return(rbind(out,colMeans(TS[,"C",]>0))); }
      #In fitting also tracked mean and variance in proportion distancing (omgs) and proportion of infections in Toronto (fracTor)
      if(plotgive=="fit"){ omgs=t(apply(apply(TS[,Tk,],c(1,3),sum)*matrix(1/Ns,L,Trun), 2, function(x) meansd(1-exp(-parms[1,"omg"]*x),Ns)));
      propCits=t(TS[,"Nt",]); states=t(apply(TS,c(2,3),sum)); return(cbind(states,omgs,propCits)); }
    }

    x=(msim3(m3iter,parms,(years*365)+1,plotgive="TS",seed0=1,InitInf=c(2)));

    for (i in 1:length(x[,1,1])) {
      
      tmp<-as.data.frame(t(x[i,,]))
      tmp$ts <- as.numeric(row.names(tmp))-1
      tmp$run<-i
      
      if (i==1) {
        df<-tmp
      } else {
        df<-rbind(df,tmp)
      }
    }
    
    df$sick<-rowSums(df[,colnames(df) %in% c("E",Di,Da)])
    df$know.current<-rowSums(df[,colnames(df) %in% Tk])
    
    df<- df %>% group_by(run) %>% mutate(new.K = Nt - lag(Nt), cum.D = cumsum(D), cum.Rd=cumsum(Rd))
    
    
    if (select.var=="beta") df$beta<-parms[1,"beta"]
    if (select.var=="epsilon") df$epsilon<-parms[1,"eps"]
    if (select.var=="gamma") df$gamma<-inCstart
    if (select.var=="omega") df$omega<-parms[1,"omg"]
    if (select.var=="tau_I") df$tau_I<-parms[1,"tauI"]
    if (select.var=="b") df$b<-parms[1,"births"]
    if (select.var=="xi")   df$xi<-parms[1,"dlt"]
    if (select.var=="d") df$d<-parms[1,"deaths"]
    if (select.var=="delta_C") df$delta_C<-Resurge
    
    df$initpop<-pop
    
    df.init<-df[df$ts==min(df$ts),]
    outbreaks<-df.init$run[df.init$sick>0]
    ### Only keep runs that had an outbreak (matters for very small inital pops)
    df.outbreaks<-df[df$run %in% outbreaks,]
    
    df.end<-df.outbreaks[df.outbreaks$ts==max(df.outbreaks$ts),]
    df.end$inf<-rowSums(df.end[,c("E",Da,Di,"R","Rd_covid", "cum.Rd")]) #Everyone who currently has covid + everyone who has previously had covid
    df.end$currentpop<-rowSums(df.end[,c("S","E",Da,Di,"R")]) #current population size
    df.end$alltimepop<-rowSums(df.end[,c("S","E",Da,Di,"R","Rd_covid", "cum.D")]) #All individuals who ever lived, don't need to include cum.Rd here because it's already counted in cum.D
    
    df.sum.end <- df.end %>% group_by_at(vars(one_of(c("initpop", select.var, "ts")))) %>% summarise(
      mean.propinf.boot= smean.cl.boot(inf/alltimepop, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[1],
      lowCI95.propinf.boot= smean.cl.boot(inf/alltimepop, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[2],
      highCI95.propinf.boot= smean.cl.boot(inf/alltimepop, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[3],
      min.propinf=min(inf/alltimepop), max.propinf=max(inf/alltimepop),
      mean.inf=mean(inf), min.inf=min(inf), max.inf=max(inf),
      mean.curpop=mean(currentpop), min.curpop=min(currentpop), max.curpop=max(currentpop),
      mean.alltimepop=mean(alltimepop), min.alltimepop=min(alltimepop), max.alltimepop=max(alltimepop),
      zeroevent = sum(sick==0),
      N=length(unique(df.end$run)))
    
    df.365<-df.outbreaks[df.outbreaks$ts==365,]
    df.365$inf<-rowSums(df.365[,c("E",Da,Di,"R","Rd_covid", "cum.Rd")])
    df.365$currentpop<-rowSums(df.365[,c("S","E",Da,Di,"R")])
    df.365$alltimepop<-rowSums(df.365[,c("S","E",Da,Di,"R","Rd_covid","cum.D")])
    
    df.sum.365 <- df.365 %>% group_by_at(vars(one_of(c("initpop", select.var, "ts")))) %>% summarise(
      mean.propinf.boot= smean.cl.boot(inf/alltimepop, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[1],
      lowCI95.propinf.boot= smean.cl.boot(inf/alltimepop, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[2],
      highCI95.propinf.boot= smean.cl.boot(inf/alltimepop, conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)[3],
      min.propinf=min(inf/alltimepop), max.propinf=max(inf/alltimepop),
      mean.inf=mean(inf), min.inf=min(inf), max.inf=max(inf),
      mean.curpop=mean(currentpop), min.curpop=min(currentpop), max.curpop=max(currentpop),
      mean.alltimepop=mean(alltimepop), min.alltimepop=min(alltimepop), max.alltimepop=max(alltimepop),
      zeroevent = sum(sick==0),
      N=length(unique(df.end$run)))
    
    if (z1==1 && z2==1) {
      df.sum.end.comb<-df.sum.end
      df.sum.365.comb<-df.sum.365
    } else {
      df.sum.end.comb<-rbind(df.sum.end.comb, df.sum.end)
      df.sum.365.comb<-rbind(df.sum.365.comb, df.sum.365)
    }

    counter<-counter+1
    #Displays progress
    print(counter/(length(var.sample)*length(popsample)))
    
  } 
  }


df.sum.comb<-rbind(df.sum.end.comb, df.sum.365.comb)
write.csv(df.sum.comb, sprintf("COVIDsim_CCS_fadeoutdat_%s.csv", select.var), row.names = FALSE)


df.ccs.end<-df.sum.end.comb %>% group_by_at(vars(one_of(c(select.var)))) %>% summarize(ccs=min(initpop[zeroevent == 0]), #ccs finds the smallest population size with no "zero-events"
                                                             ccs.smart=initpop[length(zeroevent)-Position(function(x) x > 0, rev(zeroevent))+2]) #ccs.smart finds the smallest population size with no "zero-events" for which all larger population sizes also experience no "zero-events"
df.ccs.365<-df.sum.365.comb %>% group_by_at(vars(one_of(c(select.var)))) %>% summarize(ccs=min(initpop[zeroevent == 0]),
                                                             ccs.smart=initpop[length(zeroevent)-Position(function(x) x > 0, rev(zeroevent))+2])

df.ccs.end[is.infinite(df.ccs.end$ccs)==TRUE,c(2:5)]<-NA
df.ccs.end$ts<-max(df$ts)
df.ccs.365[is.infinite(df.ccs.365$ccs)==TRUE,c(2:5)]<-NA
df.ccs.365$ts<-365

df.ccs<-rbind(df.ccs.end,df.ccs.365)
write.csv(df.ccs, sprintf("COVIDsim_CCS_ccsdat_%s.csv", select.var), row.names = FALSE)

tfin<-Sys.time()

print(tfin-tstart)