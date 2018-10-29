library(dplyr)
library(ggplot2)

data = read.csv("FeCl2 time.csv")

## cut off low counts with high variance (95% CI should not include 0)
data <- mutate(data, meanC= (C1+C2+C3+C4+C5+C6+C7+C8+C9+C10)/10)
data <- mutate(data, varC = 1/9*((C1 - meanC)^2 + (C2 - meanC)^2 + (C3 - meanC)^2 + (C4 - meanC)^2 + (C5 - meanC)^2 + (C6 - meanC)^2 + (C7 - meanC)^2 + (C8 - meanC)^2 + (C9 - meanC)^2 + (C10 - meanC)^2))
data <- mutate(data, SEC = sqrt(varC/10))
data <- mutate(data, lowerC = meanC - 1.96*SEC)

data <- data[data$lowerC > 0, ]


## transform independent variable (Percent oxidation)

data <- mutate(data, Percent.Oxidation = Percent.Oxidation * 100)


MS2 <- data[data$Phage == "MS2", ]
P22 <- data[data$Phage == "P22", ]

plot <- ggplot(data, aes(Percent.Oxidation, Log.R))+ 
  geom_point(size = 3, stroke = 1.5, aes(shape = Phage), na.rm = TRUE) +
  scale_shape_manual(values = c(1, 17)) +
  labs(x = "Percent Iron Oxidation", y = "Log Inactivation") +
  theme_classic(base_size = 18)



## Develop models 

# MS2

Mlm <- lm(Log.R ~ Percent.Oxidation, data = MS2)



# P22

Plm <- lm(Log.R ~ Percent.Oxidation, data = P22)


### add linear regressions to plot

predicted_MS2 <- data.frame(MS2_pred = predict(Mlm, MS2), Percent.Oxidation = MS2$Percent.Oxidation)

predicted_P22 <- data.frame(P22_pred = predict(Plm, P22), Percent.Oxidation = P22$Percent.Oxidation)

plot + geom_line(linetype = 2, data = predicted_MS2, aes(x=Percent.Oxidation, y=MS2_pred)) +
  geom_line(linetype = 4, data = predicted_P22, aes(x=Percent.Oxidation, y=P22_pred)) + 
  labs(x = "Percent Ferrous Oxidation (time)", y = "Log Inactivation") + 
  expand_limits(y=0)