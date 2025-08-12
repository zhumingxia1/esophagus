library(limma)
library(survival)

rt=read.table("score.txt", header=T, sep="\t", check.names=F,row.names=1)
rt$futime=rt$futime/365


outTab=data.frame()
sigGenes=c()
for(i in colnames(rt[,3:ncol(rt)])){
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<1){
		sigGenes=c(sigGenes,i)
		outTab=rbind(outTab,
				         cbind(id=i,
				         HR=coxSummary$conf.int[,"exp(coef)"],
				         HR.95L=coxSummary$conf.int[,"lower .95"],
				         HR.95H=coxSummary$conf.int[,"upper .95"],
				         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
				        )
	}
}

write.table(outTab,file="uniCox score cli.txt",sep="\t",row.names=F,quote=F)

sigGeneExp=data[sigGenes,]
sigGeneExp=rbind(id=colnames(sigGeneExp), sigGeneExp)
write.table(sigGeneExp, file="uni scorecli.txt", sep="\t", quote=F, col.names=F)

