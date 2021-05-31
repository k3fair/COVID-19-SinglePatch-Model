#     covidHier0.35Lite_singlepatch_scenarios.R runs simulations for NPI scenarios, visualizes results and outputs figure.
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

library(abind); library(fields); library(ggplot2); library(RColorBrewer); library(patchwork);library(dplyr);

#Select the scenario you want to consider; moderate level of NPIs ("moderate"), high trigger prevalence ("hightrig"), low percieved risk ("lowomega"), high percieved risk ("highomega"), no NPIs ("donothing"), long closures ("longclose")
scenario.select<-"moderate"

#Notes:
# one parameter renamed here, we have xi=dlt
#pi here (not in parms matrix, fixed at 0.2) is fraction of people who never develop symptoms. So asymptomatic stage is a bit longer and a fraction pi go A->R
#Mfact scales how much travel happens compared to reality. Mfact=1 is just movement from commuting data but really other people travel too, so baseline travel is Mfact=1.5 (irrelevant here because we don't include travel)
#epsP allows voluntary distancing efficacy to differ from that of closures (epsP=eps in paper)
#Tg and Tl are the closure thresholds, n0 is fraction initially infected
#Msave is the travel matrix (currently set so there is no connectivity between patches)

### GENERIC SETUP
years<- 1
samples<-25
pop<-50000
parms=cbind(N=rep(pop, samples),births=0.000051084, deaths=0.000020737, coviddeaths=0.0066, Mfact=1.5,s=0.2,Tg=1,Tl=1.5e-4,beta=0.43,tauI=0.3,tauEAf=0,epsP=0.85,eps=0.85,w=0.45,omg=44356,r=0.19,eta=0.8,dlt=0,alpha=0.4,sig=0.4,rho=0.67,n0=0.0001,M=diag(pop, samples, samples));

if (scenario.select=="hightrig") {parms[,"Tl"]<-2e-3;}
if (scenario.select=="lowomega") {parms[,"omg"]<-4435.6;}
if (scenario.select=="highomega") {parms[,"omg"]<-443560;}
if (scenario.select=="donothing") {parms[, c("tauI")]<-0;}


inCstart= parms[1,"Tl"]; #set initial trigger prevalence equal to dynamic trigger prevalence
inClen=30*3 #Duration of initial province closure
Resurge=30 #Duration of additional closures
NP=ncol(parms)-nrow(parms)

if (scenario.select=="longclose") {Resurge=120;}

#Defining model state space. Tn, Tk, Da, and Di are all untested, tested, asymptomatics, and infecteds.
#Nt tracks cumulative # positive cases (including those recovered)
#C tracks # days since last closure, but in output msim converts all positive C entries into eps
Tn=c("A","SA","I","SI"); Tk=c("Ak","SAk","Ik","SIk"); Da=c("A","Ak","SA","SAk"); Di=c("I","Ik","SI","SIk"); Snm=c("S","E",Da,Di,"R", "Rd", "Nt","C", "Ci", "Ct", "B", "D", "CN");


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
  Ncurr=rowSums(S[,c("R","E",Da,Di, "S")]);
  S[,c("R","E",Da,Di, "S")]=apply(rbind(seed+(16:26), S[,c("R","E",Da,Di, "S")]), 2, function(x) gpDeath(x[-(1)],parms[1,"deaths"],x[1]))
  S[,c("D")] = Ncurr - rowSums(S[,c("R","E",Da,Di, "S")]);
  #Count up number of people in each patch after biths and deaths, to use for thresholds, etc
  Ns=rowSums(S[,c("S","E",Da,Di,"R")])

  S0k=S[,Tk]; #These are positive tests we already knew about (for individuals who are still alive)
	#Then testing results come back
	S[,c(Tn,Tk)]=gpTransB(S[,c(Tn,Tk)],as.vector(parms[,"tauI"]%*%t(c(parms[1,c("tauEAf","tauEAf")],1,1))),seed+3);
	#update cumulative total, and local active case prevalence based on our new testing
	S[,"Nt"]=S[,"Nt"]+rowSums(S[,Tk]-S0k); pos=rowSums(S[,Tk])/Ns;

	#Disease progression: pi is the fraction of people who never show symptoms (modeled implicitly)
	pi=0.2;
	## Removal without death from I
	S[,c(Di,"R")]=gpTrans(S[,c(Di,"R")],parms[1,"rho"]*(1-parms[1,"coviddeaths"]),seed,TRUE);
	## Removal via death from I
	S[,c(Di,"Rd")]=gpTrans(S[,c(Di,"Rd")],parms[1,"rho"]*parms[1,"coviddeaths"],seed,TRUE);
	## All other transitions
	S[,c(Da,"R")]=gpTrans(S[,c(Da,"R")],pi*prod(parms[1,c("sig","rho")]),seed,TRUE);
	S[,c(Da,Di)]=gpTrans(S[,c(Da,Di)],parms[1,"sig"]*(1-pi),seed+1);
	S[,c("E","A")]=gpTrans(S[,c("E","A")],parms[1,"alpha"]*(1-parms[1,"s"]),seed+2); S[,c("E","SA")]=gpTrans(S[,c("E","SA")],parms[1,"alpha"]*parms[1,"s"],seed+2.5);

	#Closure decisions: C tracks # days dC since last closure (ie units 1/eps rather than integers).
  ### Replace dC, dCt with -1 to get rid of closures
  dC=round(S[,"C"]/parms[1,"eps"]);
  dCt=round(S[,"Ct"]/parms[1,"eps"]);
  past_init<-(dCt>=inClen); in_init <-((S[,"Nt"]/Ns >= inCstart & dCt==0) | (dCt<inClen & dCt>0));
	decide= past_init & (dC==0 | dC>Resurge-1); #determines if we're past initial and any subesequent closures
  closed=parms[1,"eps"]*((pos > parms[,"Tl"]*decide & parms[,"Tl"]<1 & past_init==TRUE) | in_init | (pos==0 & decide==FALSE & parms[,"Tl"]<1 & past_init==TRUE));
	S[,"C"]=S[,"C"]*(!decide)+(past_init*closed);
	S[,"Ci"]=in_init*closed;
	S[,"Ct"]=S[,"Ct"] + closed;

	### NOTE: right now the below travel stuff is OK because we have no travel and the initial matrix (which has static population size) which is being used for caculating travel proportions is using the initial static population size to do so --- we will need to fix this to incorporate dynamic pop size if we add travel back in
	#Make modified travel matrices for people feeling sick and/or tested positive, then implement travel
	M=Mc=parms[,-(1:NP)]; Mc=M[1:nrow(M),]*(1-closed); diag(Mc)=diag(Mc)+1-colSums(Mc); McA=abind(Mc,Mstay(Mc,parms[1,"eta"]),Mstay(Mc,parms[1,"r"]),Mstay(Mc,1-(1-parms[1,"eta"])*(1-parms[1,"r"])),along=3);
	Ss=reshfl2(S[,"S"],Mc,seed+4,FALSE); Rearr=apply(rbind(seed+(5:14), c(1,1,1:2,1:2,3:4,3:4), S[,c("R","E",Da,Di)]), 2, function(x) reshfl2(x[-(1:2)],McA[,,x[2]],x[1]))
  Pstar=rowSums(cbind(Ss,Rearr));
	dlt=parms[1,"dlt"];
	Pmod= dlt + (1-dlt)/Pstar;

	#Baseline infection probability. The last term with rowSums is how many people end up in each county for the day.
	Infect = parms[1,"beta"]*(parms[1,"w"]*(1-closed) + (1-parms[1,"w"])*(1-parms[1,"epsP"]*(1-exp(-parms[1,"omg"]*pos)))) * Pmod
	#Class-specific infection modifiers. Net infection Pr is then 1-Pr(avoid infection by anyone).
	modK=1-parms[1,"eta"]; modS=1/parms[1,"s"]-1; modA=cbind(Infect,modK*Infect,modS*Infect,modS*modK*Infect);
	Infects=1 - apply((1-cbind(modA,2*modA))^Rearr[,-(1:2)], 1, prod)
	if (max(Infects)>1) print(max(Infects));

	#Implement infection and move susceptibles and newly exposeds back to home county
	#Record all previously exposed
	C0=S[,c("E")];
	S[,c("S","E")]=cbind(0,S[,"E"]) + colSums(gpTransB(cbind(Ss,0*Ss),round(Infects,10),seed+15))
	#Update with newly exposed
	S[,"CN"]=S[,c("E")]-C0;

	return(S)
}; FUN=m3iter

closein=function(threshs,tsNt,tm){
#   #uncomment lines below to revert to initial closure in the style of Karateyev et al. (2020), right now this function doesn't do anything
# 	if(max(tsNt)<inCstart) return(1+0*threshs) #Closures not yet started
# 	if((tm-which(tsNt>=inCstart)[1])<inClen) return(0*threshs) #Initial closure started
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

df<- df %>% group_by(run) %>% mutate(new.K = Nt - lag(Nt), cum.D = cumsum(D))

df.fin<-df[df$ts==max(df$ts),]
df.fin$propinf<-rowSums(df.fin[,c("E",Da,Di,"R","Rd")])/rowSums(df.fin[,c("S","E",Da,Di,"R","Rd","cum.D")])

df$sick[df$sick==0]<-NA

transp<-0.25

p1<-ggplot(data=df, aes(x=ts, y=new.K)) +
  geom_line(aes(group=run), colour="dodgerblue3", alpha=transp) +
  stat_summary(geom="line", fun=mean, colour="black", linetype="solid") +
  xlab("Day") +
  ylab("# New confirmed cases") +
  scale_x_continuous(expand = c(0,0), breaks=seq(0,390,30)) +
  scale_y_continuous(expand=c(0.015,0)) +
  theme_bw()

c.breaks <- seq(0, 390, 30)
c.labels <- as.character(c.breaks)
c.labels[!(c.breaks %% 60 == 0)] <- ''

p2<-ggplot(data=df, aes(x=ts, y=run)) +
  geom_tile(aes(group=run, fill=as.factor((C+Ci)/parms[1,"eps"]))) +
  xlab("Day") +
  ylab("Realization") +
  scale_x_continuous(expand = c(0,0), breaks=c.breaks, labels=c.labels) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual("Closed?", values=c("dodgerblue3", "firebrick2"), labels=c("No", "Yes"), drop=FALSE) +
  theme_bw()

if (scenario.select=="donothing") {
  p1<-ggplot(data=df, aes(x=ts, y=CN)) +
    geom_line(aes(group=run), colour="dodgerblue3", alpha=transp) +
    stat_summary(geom="line", fun=mean, colour="black", linetype="solid") +
    xlab("Day") +
    ylab("# New infections") +
    scale_x_continuous(expand = c(0,0), breaks=seq(0,390,30)) +
    theme_bw()

  p2<-ggplot(data=df, aes(x=ts, y=run)) +
    geom_tile(aes(group=run, fill=factor((C+Ci)/parms[1,"eps"], levels=c("0", "1")))) +
    xlab("Day") +
    ylab("Realization") +
    scale_x_continuous(expand = c(0,0), breaks=c.breaks, labels=c.labels) +
    scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual("Closed?", values=c("dodgerblue3", "firebrick2"), labels=c("No", "Yes"), drop=FALSE) +
    theme_bw()
}


p3<-ggplot(data=df, aes(x=ts, y=run)) +
  geom_tile(aes(group=run, fill=sick)) +
  xlab("Day") +
  ylab("Realization") +
  scale_fill_distiller("# Infected", palette="YlOrRd", direction = 1, na.value = "black") +
  scale_x_continuous(expand = c(0,0), breaks=c.breaks, labels=c.labels) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()


p.comb<-wrap_plots( wrap_plots(p3, p2, nrow=1), p1 + guides(colour=FALSE), nrow = 2) +
  plot_annotation(tag_levels = 'a') & theme(legend.position = 'right', legend.margin = margin(0, 0, 0, -5, "pt"), plot.margin = margin(0, 5, 0, 0, "pt"))

print(p.comb)

png(sprintf("COVIDsim_scenario_%s.png",scenario.select), width=20, height=10, units="cm", res=500)
print(p.comb)
dev.off()
