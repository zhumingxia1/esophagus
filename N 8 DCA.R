library(survival)
library(ggDCA)

#for crt
rt=read.table("crt.txt", header=T, sep="\t", check.names=F, row.names=1)
rt[,"Age"]=ifelse(rt[,"Age"]>68, 1, 0)
rt$futime=rt$futime/365
predictTime=2    #predicted time
Risk<-coxph(Surv(futime,fustat)~riskScore,rt)
BT<-coxph(Surv(futime,fustat)~LMR + SII + Prealbumin,rt)
Clinical<-coxph(Surv(futime,fustat)~Age + Gender + Stage,rt)
RnB<-coxph(Surv(futime,fustat)~riskScore + LMR + SII + Prealbumin,rt)
Nomogram <-coxph(Surv(futime,fustat)~ riskScore + LMR + SII + Prealbumin + Age + Gender + Stage,rt)

pdf(file="DCA crt 2year.pdf", width=6.5, height=5.2)
d_train=dca(Risk,BT,Clinical,RnB,Nomogram, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()




#for ici
rt=read.table("ici.txt", header=T, sep="\t", check.names=F, row.names=1)
rt[,"Age"]=ifelse(rt[,"Age"]>68, 1, 0)
rt$futime=rt$futime/365


predictTime=2    #predicted time
Risk<-coxph(Surv(futime,fustat)~riskScore,rt)
BT<-coxph(Surv(futime,fustat)~LMR + CRP + Albumin + PNI + Prealbumin,rt)
Clinical<-coxph(Surv(futime,fustat)~Age + Gender + Stage,rt)
RnB<-coxph(Surv(futime,fustat)~riskScore + LMR + CRP + Albumin + PNI + Prealbumin,rt)
Nomogram <-coxph(Surv(futime,fustat)~ riskScore + LMR + CRP + Albumin + PNI + Prealbumin + Age + Gender + Stage,rt)

pdf(file="DCA ici 2year.pdf", width=6.5, height=5.2)
d_train=dca(Risk,BT,Clinical,RnB,Nomogram, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()


