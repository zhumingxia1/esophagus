library(limma)
library(survival)

## Load data
rt=read.table("score.txt", header=T, sep="\t", check.names=F,row.names=1)
## Convert follow-up time from days to years
rt$futime=rt$futime/365

##########################  Univariate Cox analysis  ##########################
# Create empty objects to store results
outTab=data.frame()
sigGenes=c()
# Loop through each feature column
for(i in colnames(rt[,3:ncol(rt)])){
	 # Build a univariate Cox proportional hazards model for each feature
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	# Summarize Cox regression results
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<1){
		sigGenes=c(sigGenes,i)
		 # Combine the results: HR, 95% CI, and p-value
		outTab=rbind(outTab,
				         cbind(id=i,
				         HR=coxSummary$conf.int[,"exp(coef)"],
				         HR.95L=coxSummary$conf.int[,"lower .95"],
				         HR.95H=coxSummary$conf.int[,"upper .95"],
				         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
				        )
	}
}
###########################  Save Cox results  #################################
write.table(outTab,file="uniCox score cli.txt",sep="\t",row.names=F,quote=F)
sigGeneExp=data[sigGenes,]
sigGeneExp=rbind(id=colnames(sigGeneExp), sigGeneExp)
write.table(sigGeneExp, file="uni scorecli.txt", sep="\t", quote=F, col.names=F)

