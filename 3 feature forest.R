# ============================ User parameters =================================
inputFile = "feature forest.txt"   # input table containing HR, 95% CI, p-value
outFile   = "feature forest.pdf"   # output PDF for the forest plot
# Working directory
setwd("F:\\radiomics")

##########load data################
rt = read.table(inputFile, header=T, sep="\t", row.names=1, check.names=F)
gene = rownames(rt)

#Convert character types to numeric types
hrNum     = as.numeric(rt$"HR")
hrLowNum  = as.numeric(rt$"HR.95L")
hrHighNum = as.numeric(rt$"HR.95H")

# Format numbers for display in the textual columns
hr        = sprintf("%.3f", hrNum)
hrLow     = sprintf("%.3f", hrLowNum)
hrHigh    = sprintf("%.3f", hrHighNum)

# Construct "Hazard ratio" text column
Hazard.ratio = paste0(hr,"(", hrLow, "-", hrHigh, ")")
# p-value display
pVal         = ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))

# Open PDF device; width/height control the overall figure size
pdf(file=outFile, width = 18, height = 7)
n = nrow(rt)   # number of features
nRow = n + 1    # one extra row for column headers in the left panel
ylim = c(1, nRow)   # y-axis limits
layout(matrix(c(1,2), nc=2), width=c(3,2))   # panel layout

# ============================ Left panel (labels) =============================
xlim = c(0,3)   # X-axis range for the left panel where text columns are drawn
par(mar=c(4,2,1.5,1.5)) # Margins
plot(1, xlim=xlim, ylim=ylim, type="n", axes=F, xlab="", ylab="")
text.cex = 0.8
text(0, n:1, gene, adj=0, cex=text.cex)
text(1.7, n:1, pVal, adj=1, cex=text.cex); text(1.7, n+1, 'pvalue', cex=text.cex, font=2, adj=1)     #1.7就是 pVal 那一列所在的 x 坐标。如果你想让这一列向右移动，可以直接增大这个 x 值即可
text(3, n:1, Hazard.ratio, adj=1, cex=text.cex); text(3, n+1, 'Hazard ratio', cex=text.cex, font=2, adj=1)

# ============================ Right panel (forest) ============================
par(mar=c(4,1,1.5,1), mgp=c(2,0.5,0))

##range of the X-axis
xLowLim  = min(hrLowNum, na.rm=TRUE)
xHighLim = max(hrHighNum, na.rm=TRUE)
# Bound the limits to avoid extreme/invalid values on a log scale
if(xLowLim < 1e-6)  xLowLim  = 1e-6
if(xHighLim < 1e-6) xHighLim = 1e-6 
if(xHighLim > 1e6) xHighLim = 1e6
# Initialize the forest plot area (log-scaled x-axis for HR)
plot(
    1,
    xlim = c(xLowLim, xHighLim),
    ylim = ylim,
    type = "n",
    axes = FALSE,
    ylab = "",
    xlab = "Hazard ratio",
    log  = "x"       # log-scale is standard for HR visualization   
)

# Draw confidence interval bars as horizontal arrows from lower to upper CI
arrows(
    hrLowNum, n:1,
    hrHighNum, n:1,
    angle = 90,
    code  = 3,
    length= 0.03,
    col   = "darkblue",
    lwd   = 2.5
)

# Reference line at HR = 1
abline(v=1, col="black", lty=2, lwd=2)


boxcolor = ifelse(hrNum > 1, "red", "blue")
points(hrNum, n:1, pch=15, col=boxcolor, cex=1.3)

axis(1)
# Close the PDF device
dev.off()