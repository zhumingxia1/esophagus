library(survival)
library(ggDCA)

################################### CRT cohort ##################################
# Read CRT dataset
rt=read.table("crt.txt", header=T, sep="\t", check.names=F, row.names=1)
# Dichotomize Age at 68 years
rt[,"Age"]=ifelse(rt[,"Age"]>68, 1, 0)
# Convert follow-up time from days to years
rt$futime=rt$futime/365
predictTime=2    # Predicted time for DCA
# Define Cox models used in DCA
Risk<-coxph(Surv(futime,fustat)~riskScore,rt)
BT<-coxph(Surv(futime,fustat)~LMR + SII + Prealbumin,rt)
Clinical<-coxph(Surv(futime,fustat)~Age + Gender + Stage,rt)
RnB<-coxph(Surv(futime,fustat)~riskScore + LMR + SII + Prealbumin,rt)
Nomogram <-coxph(Surv(futime,fustat)~ riskScore + LMR + SII + Prealbumin + Age + Gender + Stage,rt)
# Draw DCA curves and save to PDF
pdf(file="DCA crt 2year.pdf", width=6.5, height=5.2)
d_train = dca(Risk, BT, Clinical, RnB, Nomogram, times = predictTime)  # compute net benefit
ggplot(d_train, linetype = 1)  # plot decision curves (linetype fixed as 1)
dev.off()




################################### ICI cohort ##################################
# Read ICI dataset
rt=read.table("ici.txt", header=T, sep="\t", check.names=F, row.names=1)
rt[,"Age"]=ifelse(rt[,"Age"]>68, 1, 0)
rt$futime=rt$futime/365


predictTime=2    # Predicted time for DCA
# Define Cox models used in DCA
Risk<-coxph(Surv(futime,fustat)~riskScore,rt)
BT<-coxph(Surv(futime,fustat)~LMR + CRP + Albumin + PNI + Prealbumin,rt)
Clinical<-coxph(Surv(futime,fustat)~Age + Gender + Stage,rt)
RnB<-coxph(Surv(futime,fustat)~riskScore + LMR + CRP + Albumin + PNI + Prealbumin,rt)
Nomogram <-coxph(Surv(futime,fustat)~ riskScore + LMR + CRP + Albumin + PNI + Prealbumin + Age + Gender + Stage,rt)
# Draw DCA curves and save to PDF
pdf(file="DCA ici 2year.pdf", width=6.5, height=5.2)
d_train=dca(Risk,BT,Clinical,RnB,Nomogram, times=predictTime)
ggplot(d_train, linetype=1)
dev.off()


