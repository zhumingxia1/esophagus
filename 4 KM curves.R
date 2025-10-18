library(survival)
library(survminer)

################################# Read cohorts #################################
rt   <- read.table("expTime.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
test <- read.table("expTimetest.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Convert follow-up time from days to years
rt$futime   <- rt$futime/365
test$futime <- test$futime/365

riskScore <- colnames(rt)[3]

########################## Determine best cutoff (train) ########################
res.cut <- surv_cutpoint(rt, time = "futime", event = "fustat", variables = riskScore)
res.cat <- surv_categorize(res.cut)  

# Extract numeric cutoff value from surv_cutpoint() result
cut.df <- res.cut$cutpoint
cutoff <- as.numeric(cut.df$cutpoint)

res.cat$group <- res.cat[[riskScore]]


############################ KM curve (training set) ############################
# Log-rank test on training cohort
diff.tr <- survdiff(Surv(futime, fustat) ~ group, data = res.cat)
p.tr <- 1 - pchisq(diff.tr$chisq, df = 1)
lab.tr <- if (p.tr < 0.001) "p<0.001" else paste0("p=", sprintf("%.03f", p.tr))
# Kaplan–Meier fit on training cohort
fit.tr <- survfit(Surv(futime, fustat) ~ group, data = res.cat)

# KM plot for training cohort
plt.tr <- ggsurvplot(
  fit.tr, data = res.cat,
  title = paste0(title = paste0("Cohort: Training")),  # plot title
  pval = lab.tr, pval.size = 6,                        # display p-value
  legend.labs = c("low","high"),                       # legend labels (groups)
  legend.title = paste0(riskScore, " levels"),         # legend title
  font.legend = 12,
  xlab = "Time (years)", ylab = "Overall survival",    # axis labels
  break.time.by = 1, palette = c("red","blue"),        # time breaks & colors
  conf.int = TRUE, fontsize = 4, surv.median.line = "hv", # CI & median lines
  risk.table = TRUE, risk.table.title = "", risk.table.height = .25 # risk table
)
print(plt.tr)  # print training KM plot to current device
dev.off()      # close current graphic device if previously opened


######################## KM curve (testing set; fixed cutoff) ###################
# Apply the training-derived cutoff to the testing cohort
test$group <- factor(ifelse(test[[riskScore]] <= cutoff, "low", "high"), levels = c("low","high"))
diff.te <- survdiff(Surv(futime, fustat) ~ group, data = test)
p.te <- 1 - pchisq(diff.te$chisq, df = 1)
lab.te <- if (p.te < 0.001) "p<0.001" else paste0("p=", sprintf("%.03f", p.te))

fit.te <- survfit(Surv(futime, fustat) ~ group, data = test)
# KM plot for testing cohort (using the same visual settings)
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


