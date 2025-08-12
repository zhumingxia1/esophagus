library(survival)
library(survminer)
library(timeROC)

risk=read.table("score.txt", header=T, sep="\t", check.names=F,row.names=1)
risk$futime=risk$futime/365

#color
bioCol=rainbow(ncol(risk)-1, s=0.9, v=0.9)

######AUC at 1,2,3 years######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
	           marker=risk$riskScore,cause=1,
	           weighting='aalen',
	           times=c(1,2,3),ROC=TRUE)
pdf(file="ROC ici score.pdf", width=5, height=5)
plot(ROC_rt,time=1,col=bioCol[3],title=FALSE,lwd=2)
plot(ROC_rt,time=2,col=bioCol[4],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[5],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
	   c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	     paste0('AUC at 2 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	     paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	   col=bioCol[3:5], lwd=2, bty = 'n')
dev.off()


######cli ROC######
predictTime=2     #Define the predicted years
aucText=c()
pdf(file="cliROC ici 2year.pdf", width=6, height=6)

i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)

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


