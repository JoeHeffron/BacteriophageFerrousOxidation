library(dplyr)
library(ggplot2)

data = read.csv("FeCl2 DO.csv")

## remove controls

data <- data[data$Sample != "control", ]

## cut off low counts with high variance (95% CI should not include 0)
data <- mutate(data, meanC= (C1+C2+C3+C4+C5+C6+C7+C8+C9+C10)/10)
data <- mutate(data, varC = 1/9*((C1 - meanC)^2 + (C2 - meanC)^2 + (C3 - meanC)^2 + (C4 - meanC)^2 + (C5 - meanC)^2 + (C6 - meanC)^2 + (C7 - meanC)^2 + (C8 - meanC)^2 + (C9 - meanC)^2 + (C10 - meanC)^2))
data <- mutate(data, SEC = sqrt(varC/10))
data <- mutate(data, lowerC = meanC - 1.96*SEC)

data <- data[data$lowerC > 0, ]

## transform covariates

data <- mutate(data, dOH = 10^(-pHi) - 10^(pHf))


## transform independent variable (dissolved oxygen, mg/L)

data <- mutate(data, dDO = DOf - DOi)

data <- mutate(data, invDO = DOi^-1)
data <- mutate(data, invDO2 = DOi^-2)



MS2 <- data[data$Phage == "MS2", ]
P22 <- data[data$Phage == "P22", ]

plot <- ggplot(data, aes(DOi, Log.R))+ 
  geom_point(size = 2, stroke = 1.5, aes(shape = Phage), na.rm = TRUE) +
  scale_shape_manual(values = c(1, 17)) +
  theme_classic(base_size = 18)


## Develop models 

# MS2

Minvlm2 <- lm(Log.R ~ invDO + invDO2, data = MS2)

# remove leverage points
MS2lev <- MS2[-c(7,8), ]

Minvlm2 <- lm(Log.R ~ invDO + invDO2, data = MS2lev)

# inverse square not significant

Minvlm <- lm(Log.R ~ invDO, data = MS2lev)


# evaluate model including changes in hydroxide concentration and DO during test

Minvlmplus <- lm(Log.R ~ invDO + dOH + dDO, data = MS2lev)



# P22


Pinvlm2 <- lm(Log.R ~ invDO + invDO2, data = P22)

#remove leverage point
P22lev <- P22[-c(25), ]
Pinvlm2 <- lm(Log.R ~ invDO + invDO2, data = P22lev)

#inverse square term not significant

Pinvlm <- lm(Log.R ~ invDO, data = P22lev)

Pinvlmplus <- lm(Log.R ~ invDO + dOH + dDO, data = P22lev)


### add linear regressions to plot

predicted_MS2 <- data.frame(MS2_pred = predict(Minvlm, MS2lev), DO = MS2lev$DOi)

predicted_P22 <- data.frame(P22_pred = predict(Pinvlm, P22), DO = P22$DOi)

plot + geom_line(linetype = 2, data = predicted_MS2, aes(x=DO, y=MS2_pred)) +
  geom_line(linetype = 2, data = predicted_P22, aes(x=DO, y=P22_pred)) +
  labs(x = "Dissolved Oxygen (mg/L)", y = "Log Inactivation") +
  expand_limits(y = 0) 