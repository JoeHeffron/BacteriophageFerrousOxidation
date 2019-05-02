library(dplyr)
library(ggplot2)

data = read.csv("Fe oxidation.csv")

pHdata <- data[data$Test == "pH", ]
DOdata <- data[data$Test == "DO", ]


#Transform independent variables
pHdata <- mutate(pHdata, OH2 = OH^2)

DOdata <- mutate(DOdata, DO2 = DO^2)



plot <- ggplot(pHdata, aes(OH, k))+
  geom_point(size = 3, stroke = 1.5, shape = 1, na.rm = TRUE) +
  theme_classic(base_size = 18)

plot2 <- ggplot(DOdata, aes(DO, k))+
  geom_point(size = 3, stroke = 1.5, shape = 1, na.rm = TRUE) +
  theme_classic(base_size = 18)



## Develop models 

# pH

pHlm2 <- lm(k ~ OH + OH2, data = pHdata)

pHlm <- lm(k ~ OH, data = pHdata)

# DO

DOlm2 <- lm(k ~ DO + DO2, data = DOdata)
# second-order variable not significant: 

DOlm <- lm(k ~ DO, data = DOdata)


### add linear regressions to plot

predicted_pH <- data.frame(pH_pred = predict(pHlm2, pHdata), OH = pHdata$OH)

predicted_DO <- data.frame(DO_pred = predict(DOlm, DOdata), DO = DOdata$DO)

plot + geom_line(linetype = 2, data = predicted_pH, aes(x=OH, y=pH_pred)) +
  labs(x = bquote("[" *OH^'-' *"] (M)"), y = bquote("k ("*min^'-1' *")")) +
#  xlim(c(1e-9, 1.2e-6))
#  expand_limits(y = 0)+
  scale_x_continuous(breaks = c(1e-9,2e-7,4e-7,6e-7, 8e-7, 1e-6))

plot2 + geom_line(linetype = 2, data = predicted_DO, aes(x=DO, y=DO_pred)) +
  labs(x = "Dissolved Oxygen (mg/L)", y = bquote("k ("*min^'-1' *")")) #+
#  expand_limits(y = 0)+
#  scale_x_continuous(breaks = c(0,3,6,9))



#
##
###
#### Virus inactivation as a function of pH
###
##
#


pHdata = read.csv("FeCl2 pH.csv")

## remove controls

pHdata <- pHdata[pHdata$Sample != "control", ]


## cut off low counts with high variance (95% CI should not include 0)
pHdata <- mutate(pHdata, meanC= (C1+C2+C3+C4+C5+C6+C7+C8+C9+C10)/10)
pHdata <- mutate(pHdata, varC = 1/9*((C1 - meanC)^2 + (C2 - meanC)^2 + (C3 - meanC)^2 + (C4 - meanC)^2 + (C5 - meanC)^2 + (C6 - meanC)^2 + (C7 - meanC)^2 + (C8 - meanC)^2 + (C9 - meanC)^2 + (C10 - meanC)^2))
pHdata <- mutate(pHdata, SEC = sqrt(varC/10))
pHdata <- mutate(pHdata, lowerC = meanC - 1.96*SEC)

pHdata <- pHdata[pHdata$lowerC > 0, ]



## transform independent variable ([OH], mol/L)

pHdata <- mutate(pHdata, pOH = 14-pHi)
pHdata <- mutate(pHdata, OHi = 10^-(pOH))

pHdata <- mutate(pHdata, invOH = OHi^-1)
pHdata <- mutate(pHdata, invOH2 = OHi^-2)
pHdata <- mutate(pHdata, invOH3 = OHi^-3)

pHdata <- mutate(pHdata, OH = OHi)
pHdata <- mutate(pHdata, OH2 = OH^2)

colnames(pHdata)[3] <- "pH"


pHdata <- mutate(pHdata, k_pred = predict(pHlm2, pHdata))
pHdata <- mutate(pHdata, k_inv = k_pred^-1)


pHMS2 <- pHdata[pHdata$Phage == "MS2", ]
pHP22 <- pHdata[pHdata$Phage == "P22", ]

plot <- ggplot(pHdata, aes(pHi, Log.R))+ 
  geom_point(size = 3, stroke = 1.5, aes(shape = Phage), na.rm = TRUE) +
  scale_shape_manual(values = c(1, 17)) +
  theme_classic(base_size = 18)





## Develop models 

# MS2

Mkinvlm <- lm(Log.R ~ k_inv, data = pHMS2)


# P22

Pkinvlm <- lm(Log.R ~ k_inv, data = pHP22)

### add linear regressions to plot

predicted_MS2 <- data.frame(MS2_pred = predict(Mkinvlm, pHMS2), pH = pHMS2$pHi)

predicted_P22 <- data.frame(P22_pred = predict(Pkinvlm, pHP22), pH = pHP22$pHi)

plot + geom_line(linetype = 2, data = predicted_MS2, aes(x=pH, y=MS2_pred)) +
  geom_line(linetype = 4, data = predicted_P22, aes(x=pH, y=P22_pred)) +
  labs(x = "pH", y = "Log Inactivation") +
  expand_limits(y = 0) 

#
##
###
#### Virus inactivation as a function of DO
###
##
#

DOdata = read.csv("FeCl2 DO.csv")

## remove controls

DOdata <- DOdata[DOdata$Sample != "control", ]

## cut off low counts with high variance (95% CI should not include 0)
DOdata <- mutate(DOdata, meanC= (C1+C2+C3+C4+C5+C6+C7+C8+C9+C10)/10)
DOdata <- mutate(DOdata, varC = 1/9*((C1 - meanC)^2 + (C2 - meanC)^2 + (C3 - meanC)^2 + (C4 - meanC)^2 + (C5 - meanC)^2 + (C6 - meanC)^2 + (C7 - meanC)^2 + (C8 - meanC)^2 + (C9 - meanC)^2 + (C10 - meanC)^2))
DOdata <- mutate(DOdata, SEC = sqrt(varC/10))
DOdata <- mutate(DOdata, lowerC = meanC - 1.96*SEC)

DOdata <- DOdata[DOdata$lowerC > 0, ]





## transform independent variable (dissolved oxygen, mg/L)

DOdata <- mutate(DOdata, dDO = DOf - DOi)

DOdata <- mutate(DOdata, DO = DOi)

DOdata <- mutate(DOdata, invDO = DOi^-1)
DOdata <- mutate(DOdata, invDO2 = DOi^-2)

## predict ferrous oxidation rate under given conditions

DOdata <- mutate(DOdata, OH = 10^(pHi-14))
DOdata <- mutate(DOdata, OH2 = OH^2)

DOdata <- mutate(DOdata, k_pred = predict(DOlm, DOdata))
DOdata <- mutate(DOdata, k_inv = k_pred^-1)


DOMS2 <- DOdata[DOdata$Phage == "MS2", ]
DOP22 <- DOdata[DOdata$Phage == "P22", ]

DOplot <- ggplot(DOdata, aes(DOi, Log.R))+ 
  geom_point(size = 2, stroke = 1.5, aes(shape = Phage), na.rm = TRUE) +
  scale_shape_manual(values = c(1, 17)) +
  theme_classic(base_size = 18)


## Develop models 

# MS2

DOMkinvlm <- lm(Log.R ~ k_inv, data = DOMS2)


# remove leverage points
DOMS2lev <- DOMS2[-c(7,8), ]

DOMkinvlm <- lm(Log.R ~ k_inv, data = DOMS2lev)



# P22

DOPkinvlm <- lm(Log.R ~ k_inv, data = DOP22)

#remove leverage point
DOP22lev <- DOP22[-c(25), ]

DOPkinvlm <- lm(Log.R ~ k_inv, data = DOP22lev)


### add linear regressions to plot

predicted_MS2_DO <- data.frame(MS2_pred = predict(DOMkinvlm, DOMS2lev), DO = DOMS2lev$DOi)

predicted_P22_DO <- data.frame(P22_pred = predict(DOPkinvlm, DOP22lev), DO = DOP22lev$DOi)

DOplot + geom_line(linetype = 2, data = predicted_MS2_DO, aes(x=DO, y=MS2_pred)) +
  geom_line(linetype = 2, data = predicted_P22_DO, aes(x=DO, y=P22_pred)) +
  labs(x = "Dissolved Oxygen (mg/L)", y = "Log Inactivation") +
  expand_limits(y = 0) 