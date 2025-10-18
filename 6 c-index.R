
library(survival)
library(pec)

############################# Read data ##########################
dat <- read.table("crt score.txt", header=T, sep="\t", check.names=F,row.names=1)
# Remove rows containing NA values
dat <- as.data.frame(na.omit(dat))
# Convert survival time from days to years
dat$OS.time <- dat$OS.time/365 
head(dat)

############################# Cox model definitions ############################
# ---------- For CRT group ----------
cox1 <- coxph(Surv(OS.time,OS)~riskScore,data = dat,x=TRUE,y=TRUE) 
cox2 <- coxph(Surv(OS.time,OS)~LMR + SII + Prealbumin,data = dat,x=TRUE,y=TRUE) 
cox3 <- coxph(Surv(OS.time,OS)~Age + Gender + Stage,data = dat,x=TRUE,y=TRUE) 
cox4 <- coxph(Surv(OS.time,OS)~riskScore + LMR + SII + Prealbumin,data = dat,x=TRUE,y=TRUE) # 
cox5 <- coxph(Surv(OS.time,OS)~riskScore + LMR + SII + Prealbumin + Age + Gender + Stage,data = dat,x=TRUE,y=TRUE)


# ---------- For ICI group ----------
cox1 <- coxph(Surv(OS.time,OS)~riskScore,data = dat,x=TRUE,y=TRUE) 
cox2 <- coxph(Surv(OS.time,OS)~LMR + CRP + Albumin + PNI + Prealbumin,data = dat,x=TRUE,y=TRUE) 
cox3 <- coxph(Surv(OS.time,OS)~Age + Gender + Stage,data = dat,x=TRUE,y=TRUE) 
cox4 <- coxph(Surv(OS.time,OS)~riskScore + LMR + CRP + Albumin + PNI + Prealbumin,data = dat,x=TRUE,y=TRUE) 
cox5 <- coxph(Surv(OS.time,OS)~riskScore + LMR + CRP + Albumin + PNI + Prealbumin + Age + Gender + Stage,data = dat,x=TRUE,y=TRUE)

############################# Evaluation time setup #############################
# Set seed for reproducibility
set.seed(123456) 
# Define evaluation times (from 1 year to max OS.time, step = 0.5 years)
eval.time <- seq(1,floor(max(dat$OS.time)),0.5) 

############################# Define model list #################################
# For CRT group
obj <- list("Risk"=cox1,
            "BT"=cox2,
            "Clinical"=cox3,
            "RnB"=cox4,
            "All"=cox5)


# For ICI group
obj <- list("Risk"=cox1,
            "BT"=cox2,
            "Clinical"=cox3,
            "RnB"=cox4,
            "All"=cox5)
			
			
############################# Time-dependent C-index ############################
timeC <- pec::cindex(object = obj,
                      formula=Surv(OS.time,OS)~.,
                      data=dat,
                      eval.times=eval.time, 
                      splitMethod = "BootCv") # bootstrap cross validation

# Extract C-index matrix
timeC.mat <- do.call(cbind,timeC$BootCvCindex) 
# Find minimal C-index value for setting y-axis limits
ymin <- min(timeC.mat) 

############################# Plot settings #####################################
# Define color palette for each model
mycol <- RColorBrewer::brewer.pal(n = ncol(timeC.mat), name = 'Set1')
pdf("crt time-dependent Cindex.pdf",width = 6,height = 5.5)
par(bty="l", 
    mgp = c(1.9,.33,0), mar=c(4.1,4.1,2.1,2.1)+.1, las=1, tcl=-.25)

############################# Plot time-dependent C-index #######################
# Loop over models and plot time-dependent C-index curves
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
# Reference line at C-index = 0.5 (random performance)
 if(ymin < 0.5) {abline(h = 0.5,lty = 4,col = "grey50",lwd = 2)} 

# Add legend with model names and colors
legend("topright",  
       legend = colnames(timeC.mat),
       col = mycol,
       lty = 1,
       lwd = 2,
       y.intersp = 1, x.intersp = 0.5, 
       bty = "o") 
invisible(dev.off()) 