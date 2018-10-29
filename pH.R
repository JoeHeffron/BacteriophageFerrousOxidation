library(dplyr)
library(ggplot2)

data = read.csv("FeCl2 pH.csv")

## remove controls

data <- data[data$Sample != "control", ]


## cut off low counts with high variance (95% CI should not include 0)
data <- mutate(data, meanC= (C1+C2+C3+C4+C5+C6+C7+C8+C9+C10)/10)
data <- mutate(data, varC = 1/9*((C1 - meanC)^2 + (C2 - meanC)^2 + (C3 - meanC)^2 + (C4 - meanC)^2 + (C5 - meanC)^2 + (C6 - meanC)^2 + (C7 - meanC)^2 + (C8 - meanC)^2 + (C9 - meanC)^2 + (C10 - meanC)^2))
data <- mutate(data, SEC = sqrt(varC/10))
data <- mutate(data, lowerC = meanC - 1.96*SEC)

data <- data[data$lowerC > 0, ]



## transform independent variable ([OH], mol/L)

data <- mutate(data, pOH = 14-pHi)
data <- mutate(data, OHi = 10^-(pOH))

data <- mutate(data, invOH = OHi^-1)
data <- mutate(data, invOH2 = OHi^-2)
data <- mutate(data, invOH3 = OHi^-3)



MS2 <- data[data$Phage == "MS2", ]
P22 <- data[data$Phage == "P22", ]

plot <- ggplot(data, aes(pHi, Log.R))+ 
  geom_point(size = 3, stroke = 1.5, aes(shape = Phage), na.rm = TRUE) +
  scale_shape_manual(values = c(1, 17)) +
  theme_classic(base_size = 18)


## Develop models 

# MS2

Minvlm3 <- lm(Log.R ~ invOH + invOH2 + invOH3, data = MS2)


# inverse cube not significant
Minvlm2 <- lm(Log.R ~ invOH + invOH2, data = MS2)


# remove leverage points
MS2lev <- MS2[-c(7), ]
MS2lev <- MS2lev[-c(6), ]


Minvlm2 <- lm(Log.R ~ invOH + invOH2, data = MS2lev)

#inverse square not significant

Minvlm <- lm(Log.R ~ invOH, data = MS2lev)


# P22

Pinvlm3 <- lm(Log.R ~ invOH + invOH2 + invOH3, data = P22)


Pinvlm2 <- lm(Log.R ~ invOH + invOH2, data = P22)


# remove leverage points
P22lev <- P22[-c(17), ]

Pinvlm3 <- lm(Log.R ~ invOH + invOH2 + invOH3, data = P22lev) 
#inverse cube not significant

# remove leverage points
P22lev <- P22lev[-c(16), ]
P22lev <- P22lev[-c(7), ]
P22lev <- P22lev[-c(14), ]

P22lev <- P22lev[-c(5,13), ]

Pinvlm2 <- lm(Log.R ~ invOH + invOH2, data = P22lev)

# inverse square not significant
Pinvlm <- lm(Log.R ~ invOH, data = P22lev)


### add linear regressions to plot

predicted_MS2 <- data.frame(MS2_pred = predict(Minvlm, MS2lev), pH = MS2lev$pHi)

predicted_P22 <- data.frame(P22_pred = predict(Pinvlm, P22lev), pH = P22lev$pHi)

plot + geom_line(linetype = 2, data = predicted_MS2, aes(x=pH, y=MS2_pred)) +
  geom_line(linetype = 4, data = predicted_P22, aes(x=pH, y=P22_pred)) +
  labs(x = "pH", y = "Log Inactivation") +
  expand_limits(y = 0) 
