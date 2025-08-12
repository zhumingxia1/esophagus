library(survival)
library(survminer)

# cohort
rt   <- read.table("expTime.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
test <- read.table("expTimetest.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

rt$futime   <- rt$futime/365
test$futime <- test$futime/365

riskScore <- colnames(rt)[3]

#best cutoff from training cohort
res.cut <- surv_cutpoint(rt, time = "futime", event = "fustat", variables = riskScore)
res.cat <- surv_categorize(res.cut)  

#grouping by cutoff
cut.df <- res.cut$cutpoint
cutoff <- as.numeric(cut.df$cutpoint)

res.cat$group <- res.cat[[riskScore]]


## ---------- training KM ----------
diff.tr <- survdiff(Surv(futime, fustat) ~ group, data = res.cat)
p.tr <- 1 - pchisq(diff.tr$chisq, df = 1)
lab.tr <- if (p.tr < 0.001) "p<0.001" else paste0("p=", sprintf("%.03f", p.tr))

fit.tr <- survfit(Surv(futime, fustat) ~ group, data = res.cat)
plt.tr <- ggsurvplot(
  fit.tr, data = res.cat,
  title = paste0(title=paste0("Cohort: Training")),
  pval = lab.tr, pval.size = 6,
  legend.labs = c("low","high"),
  legend.title = paste0(riskScore, " levels"),
  font.legend = 12,
  xlab = "Time (years)", ylab = "Overall survival",
  break.time.by = 1, palette = c("red","blue"),
  conf.int = TRUE, fontsize = 4, surv.median.line = "hv",
  risk.table = TRUE, risk.table.title = "", risk.table.height = .25
)
  print(plt.tr)
  dev.off()


## ----------testing KM（cutoff from training） ----------
test$group <- factor(ifelse(test[[riskScore]] <= cutoff, "low", "high"), levels = c("low","high"))
diff.te <- survdiff(Surv(futime, fustat) ~ group, data = test)
p.te <- 1 - pchisq(diff.te$chisq, df = 1)
lab.te <- if (p.te < 0.001) "p<0.001" else paste0("p=", sprintf("%.03f", p.te))

fit.te <- survfit(Surv(futime, fustat) ~ group, data = test)
plt.te <- ggsurvplot(
  fit.te, data = test,
  title = paste0(title=paste0("Cohort: Testing")),
  pval = lab.te, pval.size = 6,
  legend.labs = c("low","high"),
  legend.title = paste0(riskScore, " levels"),
  font.legend = 12,
  xlab = "Time (years)", ylab = "Overall survival",
  break.time.by = 1, palette = c("red","blue"),
  conf.int = TRUE, fontsize = 4, surv.median.line = "hv",
  risk.table = TRUE, risk.table.title = "", risk.table.height = .25
)
  print(plt.te)
  dev.off()


