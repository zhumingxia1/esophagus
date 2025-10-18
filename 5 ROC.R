library(survival)
library(survminer)
library(timeROC)

############################# Read data ##########################
risk=read.table("score.txt", header=T, sep="\t", check.names=F,row.names=1)
# Convert follow-up time from days to years
risk$futime=risk$futime/365

# Generate distinct colors for each variable column
bioCol=rainbow(ncol(risk)-1, s=0.9, v=0.9)

############################# ROC for 1, 2, 3 years #############################
# Time-dependent ROC analysis for 1, 2, and 3 years using the risk score
ROC_rt = timeROC(
  T       = risk$futime,     # survival time (in years)
  delta   = risk$fustat,     # event status (1 = event, 0 = censored)
  marker  = risk$riskScore,  # risk score variable
  cause   = 1,               # define event of interest
  weighting = 'aalen',       # weighting method (Aalen estimator)
  times     = c(1, 2, 3),    # time points (1, 2, and 3 years)
  ROC       = TRUE           # compute ROC curves
)
pdf(file="ROC ici score.pdf", width=5, height=5)
# Plot ROC curve for 1 year
plot(ROC_rt, time = 1, col = bioCol[3], title = FALSE, lwd = 2)
# Add 2-year ROC curve
plot(ROC_rt, time = 2, col = bioCol[4], add = TRUE, title = FALSE, lwd = 2)
# Add 3-year ROC curve
plot(ROC_rt, time = 3, col = bioCol[5], add = TRUE, title = FALSE, lwd = 2)
# Add legend with AUC values for each time point
legend('bottomright',
	   c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	     paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	     paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	   col=bioCol[3:5], lwd=2, bty = 'n')
dev.off()


############################# Clinical ROC (2-year) #############################
predictTime=2     #Define the predicted years
aucText=c()
pdf(file="cliROC ici 2year.pdf", width=6, height=6)

# Plot ROC for risk score first
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
# Add diagonal reference line
abline(0,1)

# Loop over other clinical variables for ROC comparison
for(i in 4:ncol(risk)){
	ROC_rt=timeROC(T=risk$futime,
				   delta=risk$fustat,
				   marker=risk[,i], cause=1,
				   weighting='aalen',
				   times=c(predictTime),ROC=TRUE)
	plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
	aucText=c(aucText, paste0(colnames(risk)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(risk)-1)])
dev.off()


