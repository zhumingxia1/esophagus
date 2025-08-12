
library(survival)
library(regplot)


rt=read.table("crt score cli for nomo.txt", header=T, sep="\t", check.names=F, row.names=1)
rt$futime=rt$futime/365

#plot
res.cox=coxph(Surv(futime, fustat) ~ . , data = rt)
nom1<-regplot(res.cox,
              plots = c("density", "boxes"),
              clickable=F,
              title="",
              points=TRUE,
              droplines=TRUE,
              observation=rt[1,],
              rank="sd",
              failtime = c(1,2,3),
              prfail = T)

