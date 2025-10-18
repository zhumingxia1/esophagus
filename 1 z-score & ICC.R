#############################Z score normalization####################

# Read train/test files
train <- read.csv("trainorig.csv", header = TRUE, check.names = FALSE)
test <- read.csv("testorig.csv", header = TRUE, check.names = FALSE)

# Column index of the ID
idcol <- 1

# Extract numeric feature matrices (exclude ID column)
Xtr  <- as.matrix(train[ , -idcol, drop = FALSE])
Xtes <- as.matrix(test[ , -idcol, drop = FALSE])

# Compute training mean and sd
mu    <- colMeans(Xtr, na.rm = TRUE)
sigma <- apply(Xtr,  2, sd, na.rm = TRUE)

# Scale
Ztr  <- scale(Xtr,  center = mu, scale = sigma)
Ztes <- scale(Xtes, center = mu, scale = sigma)

# Reattach ID column; preserve original feature names
train_z   <- data.frame(train[ , idcol, drop = FALSE], as.data.frame(Ztr),  check.names = FALSE)
test_z<- data.frame(test[ , idcol, drop = FALSE], as.data.frame(Ztes), check.names = FALSE)

# Save normalized datasets
write.csv(train_z,"train.csv",   row.names = FALSE)
write.csv(test_z,"test.csv",row.names = FALSE)





###########################ICC computation##########################################
# Read files
T1 <- read.csv("tumor1.csv",header = F)
T2 <- read.csv("tumor2.csv",header = F)
#dim(T1);dim(T2)

# Combine side-by-side: 
x <- dim(T1)[1]   
y <- dim(T1)[2]   
T12 <- cbind(T1,T2) 

#dim(T12) #View(T12)

library(psych)
t= 2      
icc <- c(1:y)
for(i in 1:y) {icc[i] <- ICC(T12[,c(i,i+y)])$results$ICC[t]}   
#icc

#mean(icc)
#median(icc)

#l <- length(which(icc >= 0.75))  
#s <- length(which(icc <= 0.75))

#save feature ICCs
write.csv(icc, file="icc_inter.csv")


