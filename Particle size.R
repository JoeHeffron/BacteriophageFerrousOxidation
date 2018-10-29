library(dplyr)
require(akima)
library(rgl)


data <- read.csv("Particle size.csv")


data <- mutate(data, time = time *60)


allmodel <- lm(Z.Average ~ dosesqrt + timesqrt, data = data)

allmat <- data.frame(matrix(nrow = 100, ncol = 3))

colnames(allmat) <- c("Dose", "time", "Z.pred")

allrange <- seq(0.1,12, length.out = 10)

allmat$Dose <- rep(allrange, each = 10)

allmat$time <- seq(1, 60, length.out = 10)

allmat$timesqrt <- sqrt(allmat$time)

allmat$dosesqrt <- sqrt(allmat$Dose)

allmat$Z.pred <- predict(allmodel, allmat)


### using a cutoff to model low and high doses

dosecutoff <- 3

lowdose <- data[data$Dose < dosecutoff, ]

highdose <- data[data$Dose > dosecutoff, ]



lowmodel <- lm(Z.Average ~ Dose + timesqrt, data = lowdose)


highmodel <- lm(Z.Average ~ Dose + timesqrt, data = highdose)



# construct matrix of predicted values based on low and high models


lowmat <- data.frame(matrix(nrow = 100, ncol = 3))

colnames(lowmat) <- c("Dose", "time", "Z.pred")

lowrange <- seq(0.1,dosecutoff -0.5, length.out = 10)

lowmat$Dose <- rep(lowrange, each = 10)

lowmat$time <- seq(1, 60, length.out = 10)

lowmat$timesqrt <- sqrt(lowmat$time)

lowmat$Z.pred <- predict(lowmodel, lowmat)

highrange <- seq(dosecutoff+0.5, 12, length.out = 10)

highmat <- data.frame(matrix(nrow = 100, ncol = 3))

colnames(highmat) <- c("Dose", "time", "Z.pred")

highmat$Dose <- rep(highrange, each = 10)

highmat$time <- seq(1, 60, length.out = 10)

highmat$timesqrt <- sqrt(highmat$time)

highmat$Z.pred <- predict(highmodel, highmat)

allmat <- rbind(lowmat, highmat)



# interpolate a regular matrix from the predicted data frame

s <- interp(allmat$time, allmat$Dose, allmat$Z.pred, nx =50, ny = 50)


open3d()
par3d(cex = 1.1, font =2)
persp3d(s, col = "blue", xlab = " ", ylab = " ", zlab = " ", alpha = 0.3, specular = "white", box = FALSE)
points3d(x = data$time, y = data$Dose, z = data$Z.Average)