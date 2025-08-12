
library(survival)
library(pec)


dat <- read.table("crt score.txt", header=T, sep="\t", check.names=F,row.names=1)

dat <- as.data.frame(na.omit(dat)) 
dat$OS.time <- dat$OS.time/365 
head(dat)

#for crt group
cox1 <- coxph(Surv(OS.time,OS)~riskScore,data = dat,x=TRUE,y=TRUE) 
cox2 <- coxph(Surv(OS.time,OS)~LMR + SII + Prealbumin,data = dat,x=TRUE,y=TRUE) 
cox3 <- coxph(Surv(OS.time,OS)~Age + Gender + Stage,data = dat,x=TRUE,y=TRUE) 
cox4 <- coxph(Surv(OS.time,OS)~riskScore + LMR + SII + Prealbumin,data = dat,x=TRUE,y=TRUE) # 
cox5 <- coxph(Surv(OS.time,OS)~riskScore + LMR + SII + Prealbumin + Age + Gender + Stage,data = dat,x=TRUE,y=TRUE)


#for ici group
cox1 <- coxph(Surv(OS.time,OS)~riskScore,data = dat,x=TRUE,y=TRUE) 
cox2 <- coxph(Surv(OS.time,OS)~LMR + CRP + Albumin + PNI + Prealbumin,data = dat,x=TRUE,y=TRUE) 
cox3 <- coxph(Surv(OS.time,OS)~Age + Gender + Stage,data = dat,x=TRUE,y=TRUE) 
cox4 <- coxph(Surv(OS.time,OS)~riskScore + LMR + CRP + Albumin + PNI + Prealbumin,data = dat,x=TRUE,y=TRUE) 
cox5 <- coxph(Surv(OS.time,OS)~riskScore + LMR + CRP + Albumin + PNI + Prealbumin + Age + Gender + Stage,data = dat,x=TRUE,y=TRUE)


set.seed(123456) 
eval.time <- seq(1,floor(max(dat$OS.time)),0.5) 

#for crt group
obj <- list("Risk"=cox1,
            "BT"=cox2,
            "Clinical"=cox3,
            "RnB"=cox4,
            "All"=cox5)


#for ici group
obj <- list("Risk"=cox1,
            "BT"=cox2,
            "Clinical"=cox3,
            "RnB"=cox4,
            "All"=cox5)
			
			

timeC <- pec::cindex(object = obj,
                      formula=Surv(OS.time,OS)~.,
                      data=dat,
                      eval.times=eval.time, 
                      splitMethod = "BootCv") # bootstrap cross validation


timeC.mat <- do.call(cbind,timeC$BootCvCindex) 
ymin <- min(timeC.mat) 

#color
mycol <- RColorBrewer::brewer.pal(n = ncol(timeC.mat), name = 'Set1')
pdf("crt time-dependent Cindex.pdf",width = 6,height = 5.5)
par(bty="l", 
    mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)

#plot
for (i in 1:ncol(timeC.mat)) { 
  if(i == 1){ 
    plot(eval.time,timeC.mat[,i],
         type="l",
         col = mycol[i],
         lwd = 2,
         ylim = c(ymin,1),xlim = range(dat$OS.time),
         xaxt = "n",
         xlab="Time (Years)",ylab = "Concordance index")
    axis(side = 1,
         at = seq(0,max(eval.time),1),
         labels = seq(0,max(eval.time),1))
  } else { 
    lines(eval.time,timeC.mat[,i],
          col = mycol[i],
          lwd = 2)
  }
}

 if(ymin < 0.5) {abline(h = 0.5,lty = 4,col = "grey50",lwd = 2)} 


legend("topright",  
       legend = colnames(timeC.mat),
       col = mycol,
       lty = 1,
       lwd = 2,
       y.intersp = 1, x.intersp = 0.5, 
       bty = "o") 
invisible(dev.off()) 