library(dplyr)
library(ggplot2)

data = read.csv("FeCl2 dose.csv")

## remove controls

data <- data[data$Sample != "control", ]


## cut off low counts with high variance (95% CI should not include 0)
data <- mutate(data, meanC= (C1+C2+C3+C4+C5+C6+C7+C8+C9+C10)/10)
data <- mutate(data, varC = 1/9*((C1 - meanC)^2 + (C2 - meanC)^2 + (C3 - meanC)^2 + (C4 - meanC)^2 + (C5 - meanC)^2 + (C6 - meanC)^2 + (C7 - meanC)^2 + (C8 - meanC)^2 + (C9 - meanC)^2 + (C10 - meanC)^2))
data <- mutate(data, SEC = sqrt(varC/10))
data <- mutate(data, lowerC = meanC - 1.96*SEC)

data <- data[data$lowerC > 0, ]


## transform independent variable ([Fe], umol/L)

# transform from OD510 to mg/L Fe

data <- mutate(data, Total = Total * 5)
data <- mutate(data, Ferrous = Ferrous * 5)

data <- mutate(data, dDose = (Total - Ferrous)/55.85*1000)

data <- mutate(data, Total2 = Total^2)


MS2 <- data[data$Phage == "MS2", ]
P22 <- data[data$Phage == "P22", ]


plot <- ggplot(data, aes(Total, Log.R))+ 
  geom_point(size = 3, stroke = 1.5, aes(shape = Phage), na.rm = TRUE) +
  scale_shape_manual(values = c(1, 17)) +
  theme_classic(base_size = 18)


## Develop models 

# MS2

Mlm2 <- lm(Log.R ~ Total + Total2 - 1, data = MS2)

Mlm <- lm(Log.R ~ Total -1, data = MS2)

# remove leverage points

MS2lev <- MS2[-c(8,9), ]

Mlm <- lm(Log.R ~ Total, data = MS2lev)



# P22

# consider doses less than 3 mg/L (same as removal of leverage points for MS2)
P22 <- P22[P22$Total < 3, ]

Plm2 <- lm(Log.R ~ Total + Total2 - 1, data = P22)

#remove leverage points

P22lev <- P22[-c(2), ]

Plm2 <- lm(Log.R ~ Total + Total2 - 1, data = P22lev)

# square not significant
Plm <- lm(Log.R ~ Total -1, data = P22lev)


### add linear regressions to plot

predicted_MS2 <- data.frame(MS2_pred = predict(Mlm, MS2lev), Total = MS2lev$Total)

predicted_P22 <- data.frame(P22_pred = predict(Plm, P22lev), Total = P22lev$Total)

plot + geom_line(linetype = 2, data = predicted_MS2, aes(x=Total, y=MS2_pred)) +
  geom_line(linetype = 4, data = predicted_P22, aes(x=Total, y=P22_pred)) +
  labs(x = "Ferrous dose (mg/L Fe)", y = "Log Inactivation") +
  expand_limits(y = 0) +
  geom_vline(xintercept = 3, color = "red", linetype = "twodash") +
  scale_x_continuous(breaks = c(0,3,6,9))
