

# Stage 2 analysis of the Multi100 project Barreca2016

# The analysis is based on the following instructions receivedf from the Multi100 team:
# "Your analysis should produce a single, main result in terms of statistical families 
# of z-, t-, F-, or χ² tests (or their alternative or non-parametric versions).
# You should use 90 °F as a threshold for temperature extremes. You should compute 
# the decline between the 1931–59 and 1960–2004 periods. You should conduct your 
# analysis without geographical, or precipitation frequency specification of the data. 
# You should use a 2 months temperature exposure window in your analysis. You should 
# disregard socio-economic factors and technological and public health advancements 
# over time in your analysis. You should not specify the causes of death in your analysis."


##########################################
#             Load packages              #
##########################################

library(haven)
library(tidyverse)
library(zoo) # for as.yearmon
library(splines) # for flexible splines
library(lme4)
library(cdlTools) # for fips
library(cAIC4)

library(lmerTest)

##########################################
#             Load data files            #
##########################################

data_folder <- "C:\\Users\\User\\Documents\\Temp\\Multi100\\2012608data\\"
#data_folder <- "C:\\Users\\kekec\\Documents\\Temp\\Multi100\\"

data_1900_2004 <- read_dta(paste0(data_folder, "DATA_1900_2004.dta"))


##########################################
#           Data processing              #
##########################################


data_1900_2004$time = as.yearmon(paste(data_1900_2004$year, data_1900_2004$month), "%Y %m")
data_1900_2004$deathper100k = data_1900_2004$dths_tot/(data_1900_2004$totalpop/100000)
data_1900_2004$state = fips(data_1900_2004$stfips, to = "Name")

data_1900_2004 = data_1900_2004 %>% 
  mutate(after1960 = case_when(year < 1960 ~ 0,
                               year >= 1960 ~ 1))

data_1900_2004_list = split(data_1900_2004, f = data_1900_2004$state)

### Computing the b10_10 for the current and the previous month to get two month temperature exposure window

for(i in 1:length(data_1900_2004_list)){
  data_1900_2004_list[[i]]$b10_10_lastmonth = c(0, data_1900_2004_list[[i]]$b10_10[1:length(data_1900_2004_list[[i]]$b10_10)-1])
  data_1900_2004_list[[i]]$b10_10_thisandlastmonth = data_1900_2004_list[[i]]$b10_10 + data_1900_2004_list[[i]]$b10_10_lastmonth
}


data_1900_2004_withlaggedb10_10 = do.call("rbind", data_1900_2004_list)



data_1931_2004_withlaggedb10_10 = 
  data_1900_2004_withlaggedb10_10 %>% 
  filter(year > 1930)

##########################################
#           Data analysis                #
##########################################

mod2 = lmer(deathper100k ~ b10_10_thisandlastmonth + after1960 + sh_0000 + sh_4564 + sh_6599 + i1_statemove + ns(month, df = 4):factor(year) + ns(time, df = 10) + (1 | state), data = data_1931_2004_withlaggedb10_10)
mod2_int = lmer(deathper100k ~ b10_10_thisandlastmonth * after1960 + sh_0000 + sh_4564 + sh_6599 + i1_statemove + ns(month, df = 4):factor(year) + ns(time, df = 10) + (1 | state), data = data_1931_2004_withlaggedb10_10)


### Main analysis inference: whether the coefficient of b10_10_thisandlastmonth:after1960 is significant
### the coefficient should be negative

summary(mod2_int)

### Alternative tests of the same hypothesis: difference of 

AIC(mod2)
AIC(mod2_int)

anova(mod2, mod2_int)




