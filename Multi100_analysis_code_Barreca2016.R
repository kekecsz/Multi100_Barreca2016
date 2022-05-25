
##########################################
#             Load packages              #
##########################################

library(haven)
library(tidyverse)
library(zoo) # for as.yearmon
library(splines) # for flexible splines
library(lme4)
library(cdlTools) # for fips

##########################################
#          Custom functions              #
##########################################


model_function = function(main_predictor, data, yearmin, yearmax, decade_effect){
  
  outcome <- "deathper100k"
  if(decade_effect == T){
    variables <- c(main_predictor, "sh_0000", "sh_4564", "sh_6599", "devp25", "devp75", "lri", "i1_physician", "i1_urban", "i1_black", "i1_statemove", "ns(month,df=4):factor(year)", "ns(time,df=10)", "(1|state)")
  } else if(decade_effect == F){
    variables <- c(main_predictor, "sh_0000", "sh_4564", "sh_6599", "devp25", "devp75", "lri", "i1_physician", "i1_urban", "i1_black", "i1_statemove", "ns(month,df=4):factor(year)", "time", "(1|state)")
  } else {print("no decade effect")}
  
  f <- as.formula(
    paste(outcome, 
          paste(variables, collapse = " + "), 
          sep = " ~ "))
  
  data_for_model = data %>% 
    filter(year >  (yearmin-1)) %>% 
    filter(year < (yearmax+1))
  
  model <- lmer(f, data = data_for_model)
  
  output_row = which(row.names(summary(model)[[10]]) == main_predictor)
  
  rowname = row.names(summary(model)[[10]])[output_row]
  estimate = summary(model)[[10]][output_row,1]
  std_error = summary(model)[[10]][output_row,2]
  t_val = summary(model)[[10]][output_row,3]
  aic = AIC(model)
  
  return(c(rowname, yearmin, yearmax, decade_effect, estimate, std_error, t_val, aic))
}



model_function_interaction = function(main_predictor, data, yearmin, yearmax, decade_effect){
  outcome <- "deathper100k"
  main_predictor_with_time_int = paste0(main_predictor, " + ", main_predictor, ":time")
  if(decade_effect == T){
    variables <- c(main_predictor_with_time_int, "sh_0000", "sh_4564", "sh_6599", "devp25", "devp75", "lri", "i1_physician", "i1_urban", "i1_black", "i1_statemove", "ns(month,df=4):factor(year)", "ns(time,df=10)", "(1|state)")
  } else if(decade_effect == F){
    variables <- c(main_predictor_with_time_int, "sh_0000", "sh_4564", "sh_6599", "devp25", "devp75", "lri", "i1_physician", "i1_urban", "i1_black", "i1_statemove", "ns(month,df=4):factor(year)", "time", "(1|state)")
  } else {print("no decade effect")}
  
  f <- as.formula(
    paste(outcome, 
          paste(variables, collapse = " + "), 
          sep = " ~ "))
  
  data_for_model = data %>% 
    filter(year >  (yearmin-1)) %>% 
    filter(year < (yearmax+1))
  
  model <- lmer(f, data = data_for_model)
  
  output_row = which(row.names(summary(model)[[10]]) == paste0(main_predictor, ":time"))
  
  rowname = row.names(summary(model)[[10]])[output_row]
  estimate = summary(model)[[10]][output_row,1]
  std_error = summary(model)[[10]][output_row,2]
  t_val = summary(model)[[10]][output_row,3]
  aic = AIC(model)
  
  return(c(rowname, yearmin, yearmax, decade_effect, estimate, std_error, t_val, aic))
}

##########################################
#             Load data files            #
##########################################

#data_folder <- "C:\\Users\\User\\Documents\\Temp\\Multi100\\2012608data\\"
data_folder <- "C:\\Users\\kekec\\Documents\\Temp\\Multi100\\"

data_1900_2004 <- read_dta(paste0(data_folder, "DATA_1900_2004.dta"))


##########################################
#           Data processing              #
##########################################


data_1900_2004$time = as.yearmon(paste(data_1900_2004$year, data_1900_2004$month), "%Y %m")
data_1900_2004$deathper100k = data_1900_2004$dths_tot/(data_1900_2004$totalpop/100000)
data_1900_2004$state = fips(data_1900_2004$stfips, to = "Name")


##########################################
#           Data analysis                #
##########################################



main_predictor_variables = paste0("b10_", 1:10)


output_table_maineffect = data.frame(rowname = NA, yearmin = NA, yearmax = NA, decade_effect = NA, estimate = NA, std_error = NA, t_val = NA, aic = NA)

for(i in 1:length(main_predictor_variables)){
  print(paste0("currently analyzing ", main_predictor_variables[i]))
  output_table_maineffect[i,] = t(model_function(main_predictor = main_predictor_variables[i], data = data_1900_2004, yearmin = 1900, yearmax = 2004, decade_effect = T))
}

output_table_interactioneffect = data.frame(rowname = NA, yearmin = NA, yearmax = NA, decade_effect = NA, estimate = NA, std_error = NA, t_val = NA, aic = NA)

for(i in 1:length(main_predictor_variables)){
  print(paste0("currently analyzing ", main_predictor_variables[i]))
  output_table_interactioneffect[i,] = t(model_function_interaction(main_predictor = main_predictor_variables[i], data = data_1900_2004, yearmin = 1900, yearmax = 2004, decade_effect = T))
}

### Main hypothesis test
### comparing AIC of models with and without the interaction of days in the given temperature bin and time (e.g. b10_1:time)
### AIC difference of at least two is considered significant

output_table_interactioneffect$AIC_difference = as.numeric(output_table_maineffect$aic) - 
as.numeric(output_table_interactioneffect$aic)

##########################################
#           Visualization                #
##########################################


### Visualization of the change of the effect over time

main_predictor_variables = paste0("b10_", 1:10)
yearbin_minimums = seq(from = 1930, to = 1990, by = 10)

output_table_maineffect_bydecade = data.frame(rowname = NA, yearmin = NA, yearmax = NA, decade_effect = NA, estimate = NA, std_error = NA, t_val = NA, aic = NA)

k = 0

for(i in 1:length(main_predictor_variables)){
  for(j in 1:length(yearbin_minimums)){
    k = k + 1
    print(paste0("currently analyzing ", main_predictor_variables[i], " ", yearbin_minimums[j], "-", yearbin_minimums[j]+9))
    output_table_maineffect_bydecade[k,] = t(model_function(main_predictor = main_predictor_variables[i], data = data_1900_2004, yearmin = yearbin_minimums[j], yearmax = yearbin_minimums[j]+9, decade_effect = F))
  }
}



output_table_maineffect_bydecade$yearmin = as.numeric(output_table_maineffect_bydecade$yearmin)
output_table_maineffect_bydecade$estimate = as.numeric(output_table_maineffect_bydecade$estimate)

### change rowname in the code below to get estimates for a specific temperature bin


output_table_maineffect_bydecade %>% 
  filter(rowname == "b10_9") %>% 
  ggplot() +
  aes(x = yearmin, y = estimate, color = rowname) +
  geom_point() +
  geom_line()




### Visualization of the change of the effect over time for each climate region


main_predictor_variables = paste0("b10_", 1:10)
yearbin_minimums = seq(from = 1930, to = 1990, by = 10)

output_table_maineffect_uscr = data.frame(uscr = NA, rowname = NA, yearmin = NA, yearmax = NA, decade_effect = NA, estimate = NA, std_error = NA, t_val = NA, aic = NA)


k = 0

climate_regions_to_analyze = c(1, 2, 3, 4, 5, 7, 8, 9) ### region 6 did not have enough data to perfomr this analysis

for(h in climate_regions_to_analyze){
  data_uscr = data_1900_2004 %>% 
    filter(uscr == h)
  for(i in 1:length(main_predictor_variables)){
    for(j in 1:length(yearbin_minimums)){
      k = k + 1
      print(paste0("currently analyzing ", main_predictor_variables[i], " ", yearbin_minimums[j], "-", yearbin_minimums[j]+9))
      output_table_maineffect_uscr[k,2:9] = t(model_function(main_predictor = main_predictor_variables[i], data = data_uscr, yearmin = yearbin_minimums[j], yearmax = yearbin_minimums[j]+9, decade_effect = F))
      output_table_maineffect_uscr[k,"uscr"] = h
    }
  }
}




output_table_maineffect_uscr$yearmin = as.numeric(output_table_maineffect_uscr$yearmin)
output_table_maineffect_uscr$estimate = as.numeric(output_table_maineffect_uscr$estimate)
output_table_maineffect_uscr$uscr = as.numeric(output_table_maineffect_uscr$uscr)

### change rowname in the code below to get estimates for a specific temperature bin

output_table_maineffect_uscr %>% 
  filter(rowname == "b10_10") %>% 
  ggplot() +
  aes(x = yearmin, y = estimate, color = factor(uscr)) +
  geom_point() +
  geom_line()

table(data_1900_2004$uscr, data_1900_2004$state)



