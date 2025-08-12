library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

##################################GSEA##########################################
#combination
mrna=read.table("mRNA selection.txt",sep="\t",header=T,check.names=F, row.names = 1)
risk=read.table("riskscore selection.txt",sep="\t",header=T,check.names=F, row.names = 1)
risk<-t(risk)
sameSample=intersect(colnames(mrna), colnames(risk))
mrna=mrna[,sameSample]
risk=risk[3:4,sameSample]
rt=rbind(mrna, risk)
rt=rt[-19499,]
gene="riskScore"                               
gmtFile="c5.go.v2024.1.Hs.symbols.gmt"          

#read gmt file
gmt=read.gmt(gmtFile)


	
#take the average of duplicate genes
rt=as.matrix(rt)
rownames(rt)=rownames(rt)
exp=rt[,1:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#define groups
dataL=data[,(data[gene,]<=median(data[gene,]))]   
dataH=data[,(data[gene,]>median(data[gene,]))]     
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH/meanL)
logFC=sort(logFC,decreasing=T)

#GSEA
kk=GSEA(logFC,TERM2GENE=gmt, nPerm=100,pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
write.table(kkTab,file=paste0("GO",".txt"),sep="\t",quote=F,row.names = F)
	
termNum=10
if(nrow(kkTab)>=termNum){
		gseaplot=gseaplot2(kk, row.names(kkTab)[1:termNum],base_size =12)
		pdf(file=paste0("GO",".pdf"),width=12,height=11)
		print(gseaplot)
		dev.off()
	}



####################################ssGSEA################################
library(GSVA)
library(limma)
library(GSEABase)


#define function
immuneScore=function(expFile=null, gmtFile=null, socreFile=null){
	rt=read.table(expFile, header=T, sep="\t", check.names=F)
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	mat=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	mat=avereps(mat)
	mat=mat[rowMeans(mat)>0,]
	
	geneSet=getGmt(gmtFile, geneIdType=SymbolIdentifier())
	
	#ssgsea analyse
	ssgseaScore=gsva(mat, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
	normalize=function(x){
	  return((x-min(x))/(max(x)-min(x)))}
	ssgseaOut=normalize(ssgseaScore)
	ssgseaOut=rbind(id=colnames(ssgseaOut),ssgseaOut)
	write.table(ssgseaOut, file=socreFile, sep="\t", quote=F, col.names=F)
}

immuneScore(expFile="mRNA.txt", gmtFile="Hallmark.gmt", socreFile="Hallmark score.txt")	


#######comaprison between different groups
rt=read.table("Hallmark score group.txt",sep="\t",header=T,check.names=F,row.name=1)
my_comparisons <- list( c("low", "high") )

for(j in colnames(rt)[2:51]) {
  data = rt[c("group", j)]
  colnames(data) = c("group", "value")
  data = subset(data, value != "unknow")
  if(nrow(data) == 0) next
  stat_res <- compare_means(
    formula     = value ~ group,
    data        = data,
    method      = "wilcox.test",
    comparisons = my_comparisons
  )
  
  if(any(stat_res$p < 0.05)) {
    boxplot = ggboxplot(
      data,
      x = "group", y = "value", 
      legend.title = "group",
      color = "group",
      palette = c("green", "red", "yellow", "blue", "purple"),
      add = "jitter"
    ) +
      stat_compare_means(comparisons = my_comparisons)  
    
    pdf(file = paste0(j, ".pdf"), width = 5.5, height = 5)
    print(boxplot)
    dev.off()
  }
}


################################cibersort############################################
library("limma")                                                      
             
rt=read.table("symbol.txt",sep="\t",header=T,check.names=F)          
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=data[,group==0]

data=data[rowMeans(data)>0,]


v <-voom(data, plot = F, save.plot = F)
out=v$E
out=rbind(ID=colnames(out),out)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)      

#run CIBERSORT
source("ssGSEA18.CIBERSORT.R")
results=CIBERSORT("imunecells.txt", "uniq.symbol.txt", perm=100, QN=TRUE)
	