# Load required libraries
library(deSolve)
library(pracma)
library(ggplot2)
library(dplyr)
library(tidyr)

lambda <- mean(rpois(1000, 2))

# Define the differential equations function for men and women
DiffEq <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    # For men
    DS_men <- -lambda*beta*S_men*I_women/(V_men+I_men+S_men) - S_men*alpha_men 
    DI_men <- lambda*beta*S_men*I_women/(V_men+I_men+S_men) - (I_men*0.9)/24 + beta*V_men*(1-0.998)
    DV_men <- S_men*alpha_men - beta*V_men*(1-0.998)
    
    # For women
    DS_women <- -beta*S_women*I_men/(V_women+I_women+S_women) - S_women*alpha_women
    DI_women <- beta*S_women*I_men/(V_women+I_women+S_women) - (I_women*0.9)/24 + beta*V_women*(1-0.998)
    DV_women <- S_women*alpha_women - beta*V_women*(1-0.998)
    
    
    # Introduce a pulse of susceptible individuals every 12 steps for both men and women
    if (round(t) %% 12 == 0) {
      DS_men <- DS_men + pulse_amount*0.41 - 1/4*DS_men
      DV_men <- DV_men - 1/4*DV_men + pulse_amount*0.59
      DI_men <- DI_men - 1/4*DI_men
      
      DS_women <- DS_women + pulse_amount*0.36 - 1/4*DS_women
      DV_women <- DV_women - 1/4*DV_women + pulse_amount*0.64
      DI_women <- DI_women - 1/4*DI_women 
    }
    
    list(c(DS_men, DI_men, DV_men, DS_women, DI_women, DV_women))
  })
}

# Define the differential equations function for declining teen vaccination rates
DiffEq1 <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    # For men
    DS_men <- -lambda*beta*S_men*I_women/(V_men+I_men+S_men) - S_men*alpha_men 
    DI_men <- lambda*beta*S_men*I_women/(V_men+I_men+S_men) - (I_men*0.9)/24 + beta*V_men*(1-0.998)
    DV_men <- S_men*alpha_men - beta*V_men*(1-0.998)
    
    # For women
    DS_women <- -beta*S_women*I_men/(V_women+I_women+S_women) - S_women*alpha_women
    DI_women <- beta*S_women*I_men/(V_women+I_women+S_women) - (I_women*0.9)/24 + beta*V_women*(1-0.998)
    DV_women <- S_women*alpha_women - beta*V_women*(1-0.998)
    
    # Introduce a pulse of susceptible individuals every 12 steps for both men and women
    if (round(t) %% 12 == 0) {
      # Current Adolescent Vaccination Rates 3.3 percentage pts every year
      vrate_men <- 0.59-0.033
      vrate_women <- 0.64-0.033
      DS_men <- DS_men + pulse_amount*(1-vrate_men) - 1/4*DS_men
      DV_men <- DV_men - 1/4*DV_men + pulse_amount*vrate_men
      DI_men <- DI_men - 1/4*DI_men
      
      DS_women <- DS_women + pulse_amount*(1-vrate_women) - 1/4*DS_women
      DV_women <- DV_women - 1/4*DV_women + pulse_amount*vrate_women
      DI_women <- DI_women - 1/4*DI_women
    }
    
    list(c(DS_men, DI_men, DV_men, DS_women, DI_women, DV_women))
  })
}

# Define the differential equations function for matching MMR adolescent vaccination rates
DiffEq2 <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    # For men
    DS_men <- -lambda*beta*S_men*I_women/(V_men+I_men+S_men) - S_men*alpha_men 
    DI_men <- lambda*beta*S_men*I_women/(V_men+I_men+S_men) - (I_men*0.9)/24 + beta*V_men*(1-0.998)
    DV_men <- S_men*alpha_men - beta*V_men*(1-0.998)
    
    # For women
    DS_women <- -beta*S_women*I_men/(V_women+I_women+S_women) - S_women*alpha_women
    DI_women <- beta*S_women*I_men/(V_women+I_women+S_women) - (I_women*0.9)/24 + beta*V_women*(1-0.998)
    DV_women <- S_women*alpha_women - beta*V_women*(1-0.998)
    
    
    # Introduce a pulse of susceptible individuals every 12 steps for both men and women
    if (round(t) %% 12 == 0) {
      DS_men <- DS_men + pulse_amount*0.087 - 1/4*DS_men
      DV_men <- DV_men - 1/4*DV_men + pulse_amount*0.913
      DI_men <- DI_men - 1/4*DI_men
      
      DS_women <- DS_women + pulse_amount*0.087 - 1/4*DS_women
      DV_women <- DV_women - 1/4*DV_women + pulse_amount*0.913
      DI_women <- DI_women - 1/4*DI_women 
    }
    
    list(c(DS_men, DI_men, DV_men, DS_women, DI_women, DV_women))
  })
}

# Set initial state for men and women given 60000 students with 54.6% women
# 59% of men and 64% of women up to date on vaccine
state1 <- c(S_men = 11168, I_men = 1000, V_men = 16072,
           S_women = 11794, I_women = 1000, V_women = 20966)
# 91% of men and women up to date
state2 <- c(S_men = 2370, I_men = 1000, V_men = 24870,
            S_women = 2850, I_women = 1000, V_women = 29910)
pulse_amount <- 9000
# Set time points for simulation
times <- seq(0, 4*12, by = 1)

### Case Scenario: Current rates of vaccination
parameters <- c(beta = 0.8795*0.479, alpha_men = 0.34, alpha_women = 0.552, lambda=lambda)
out <- ode(y = state1, times = times, func = DiffEq, parms = parameters)

### Scenario 1: Declining adolescent vacc
out1 <- ode(y = state1, times = times, func = DiffEq1, parms = parameters)

### Scenario 2: Match MMR adolescent vaccination
out2 <- ode(y = state2, times = times, func = DiffEq2, parms = parameters)


# Initialize vectors to store total infected at the end of each year (12, 24, 36, 48 months)
total_infected <- numeric(16)
total_infected_1 <- numeric(16)
total_infected_2 <- numeric(16)
total_men <- numeric(16)
total_men_1 <- numeric(16)
total_men_2 <- numeric(16)
total_women <- numeric(16)
total_women_1 <- numeric(16)
total_women_2 <- numeric(16)

# Count cumulative infected individuals (currently infected + recovered) for Case Scenario 1
for (i in 1:16) {
  time_index <- i * 3 # Time point at the end of each year
  total_infected[i] <- out[time_index, "I_men"] + out[time_index, "I_women"]
  total_men[i] <- out[time_index, "I_men"]
  total_women[i] <- out[time_index, "I_women"]
}

# Count cumulative infected individuals (currently infected + recovered) for Scenario 1
for (i in 1:16) {
  time_index <- i * 3 # Time point at the end of each year
  total_infected_1[i] <- out1[time_index, "I_men"] + out1[time_index, "I_women"]
  total_men_1[i] <- out1[time_index, "I_men"]
  total_women_1[i] <- out1[time_index, "I_women"]
}

# Count cumulative infected individuals (currently infected + recovered) for Scenario 2
for (i in 1:16) {
  time_index <- i * 3 # Time point at the end of each year
  total_infected_2[i] <- out2[time_index, "I_men"] + out2[time_index, "I_women"]
  total_men_2[i] <- out2[time_index, "I_men"]
  total_women_2[i] <- out2[time_index, "I_women"]
}


# Print the cumulative infected counts for Scenario 1
print("Cumulative Infected Counts current Vaccination rate (Every 3 Months):")
print(round(sum(total_infected),0))
print("Cumulative Infected Counts decreasing Vaccination rate (Every 3 Months):")
print(round(sum(total_infected_1),0))
print("Cumulative Infected Counts 91% Vaccination rate (Every 3 Months):")
print(round(sum(total_infected_2),0))
print("Cumulative Infected men Counts Case Scenario (Every 3 Months):")
print(round(sum(total_men),0))
print(round(sum(total_women), 0))
print("Cumulative Infected men Counts Scenario 1 (Every 3 Months):")
print(round(sum(total_men_1),0))
print(round(sum(total_women_1),0))
print("Cumulative Infected men Counts Scenario 2 (Every 3 Months):")
print(round(sum(total_men_2),0))
print(round(sum(total_women_2),0))


library(ggplot2)
library(tidyr)
# Prepare the data for combined Case Scenario and Scenario 1 (only Infected)
data_combined1 <- data.frame(
  time = out[, "time"],
  I_men_current = out[, "I_men"],         # Case Scenario
  I_women_current = out[, "I_women"],    
  I_men_dropping = out1[, "I_men"],       # Scenario 1
  I_women_dropping = out1[, "I_women"]
)

# Convert to long format
data_combined_long1 <- gather(data_combined1, key = "Population", value = "Count", -time)

# Define colors and line types
colors_combined1 <- c("I_men_current" = "blue", "I_women_current" = "purple", 
                      "I_men_dropping" = "red", 'I_women_dropping' = 'darksalmon') 

# Update the line types to valid ggplot2 types
line_types_combined1 <- c("I_men_current" = "solid", "I_women_current" = "solid", 
                          "I_men_dropping" = "dashed", 'I_women_dropping' = 'dashed')

# Plot for combined Scenario 1 and Scenario 2 (Infected)
ggplot(data_combined_long1, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined1) +
  scale_linetype_manual(values = line_types_combined1) +
  labs(x = "Time", y = "Infected Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")

#Prepare the data for combined Case Scenario, Scenario 1, and Scenario 2 (Infected Only)
data_combined2 <- data.frame(
  time = out[, "time"],
  I_men_current = out[, "I_men"],         # Case Scenario
  I_women_current = out[, "I_women"],    
  I_men_dropping = out1[, "I_men"],       # Scenario 1
  I_women_dropping = out1[, "I_women"],
  I_men_0.91 = out2[, "I_men"],        # Scenario 2
  I_women_0.91 = out2[, 'I_women']
)

# Convert to long format
data_combined_long2 <- gather(data_combined2, key = "Population", value = "Count", -time)

# Define colors and line types
colors_combined2 <- c("I_men_current" = "blue", "I_women_current" = "purple", 
                     "I_men_dropping" = "red", 'I_women_dropping' = 'darksalmon',
                     'I_men_0.91' = 'darkgreen', "I_women_0.91" = 'firebrick') 

# Update the line types to valid ggplot2 types
line_types_combined2 <- c("I_men_current" = "solid", "I_women_current" = "solid", 
                         "I_men_dropping" = "dashed", 'I_women_dropping' = 'dashed',
                         "I_men_0.91" = 'dotted', "I_women_0.91" = "dotted")

# Plot for combined Scenario 1 and Scenario 2 (Infected)
ggplot(data_combined_long2, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined2) +
  scale_linetype_manual(values = line_types_combined2) +
  labs(x = "Time", y = "Infected Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Print the cumulative infected counts for Scenario 1
print("Cumulative Infected Counts Scenario 1 (Every 12 Months):")
print(round(sum(total_infected),0))
print("Cumulative Infected Counts Scenario 2 (Every 12 Months):")
print(round(sum(total_infected_1),0))
print("Cumulative Infected Counts Scenario 3 (Every 12 Months):")
print(round(sum(total_infected_2),0))
 
