# Load required libraries
library(deSolve)
library(pracma)
library(ggplot2)
library(dplyr)
library(tidyr)

lambda <- mean(rpois(1000, 2))

# Define the differential equations function for men and women
# Define the differential equations function for men and women
DiffEq <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    # For men
    DS_men <- -lambda*beta*S_men*I_women/(V_men+I_men+S_men) - S_men*alpha_men 
    DI_men <- lambda*beta*S_men*I_women/(V_men+I_men+S_men) - (I_men*0.9)/24 + (beta*V_men*(1-0.998))/(V_men+I_men+S_men)
    DV_men <- S_men*alpha_men - beta*V_men*(1-0.998)
    
    # For women
    DS_women <- -beta*S_women*I_men/(V_women+I_women+S_women) - S_women*alpha_women
    DI_women <- beta*S_women*I_men/(V_women+I_women+S_women) - (I_women*0.9)/24 + (beta*V_women*(1-0.998))/(V_women+I_women+S_women)
    DV_women <- S_women*alpha_women - beta*V_women*(1-0.998)
    
    
    # Introduce a pulse of susceptible individuals every 12 steps for both men and women
    if (round(t) %% 12 == 0) {
      DS_men <- DS_men + 0.485*pulse_amount*0.41 - 1/4*DS_men   #<- you need to change this number (0.0.41)
      DV_men <- DV_men - 1/4*DV_men + 0.485*pulse_amount*0.59   #<- you need to change this number (0.59)
      DI_men <- DI_men - 1/4*DI_men                       # these need to be out of 1 tho so make sure that is the case
      
      DS_women <- DS_women + 0.515*pulse_amount*0.36 - 1/4*DS_women #<- you need to change this number (0.36)
      DV_women <- DV_women - 1/4*DV_women + 0.515*pulse_amount*0.64 #<- you need to change this number (0.64)
      DI_women <- DI_women - 1/4*DI_women 
    }
    
    list(c(DS_men, DI_men, DV_men, DS_women, DI_women, DV_women))
  })
}

# Set initial state for men and women
state <- c(S_men = 60000*0.485, I_men = 1000, V_men = 0,
           S_women = 60000*0.515, I_women = 1000, V_women = 0)
pulse_amount <- 13665
# Set time points for simulation
times <- seq(0, 4*12, by = 1)

### Scenario 1: Separate alpha values for men and women
parameters1 <- c(beta = 0.8795*0.479, alpha_men = 0.34, alpha_women = 0.552, lambda=lambda)
out1 <- ode(y = state, times = times, func = DiffEq, parms = parameters1)

### Scenario 2: Same alpha for men and women
parameters2 <- c(beta = 0.8795*0.479, alpha_men = 0.552, alpha_women = 0.552, lambda=lambda)
out2 <- ode(y = state, times = times, func = DiffEq, parms = parameters2)

### Scenario 3: both sexes have vaccination rates of 0.75
parameters3 <- c(beta = 0.8795*0.479, alpha_men = 0.75, alpha_women = 0.75, lambda=lambda)
out3 <- ode(y=state, times = times, func=DiffEq, parms = parameters3)


### Scenario 1: Separate alpha values for men and women
parameters4 <- c(beta = 0.8795*0.479, alpha_men = 0.2, alpha_women = 0.2, lambda=lambda)
out4 <- ode(y = state, times = times, func = DiffEq, parms = parameters4)


# Initialize vectors to store total infected at the end of each year (12, 24, 36, 48 months)
total_infected_1 <- numeric(16)
total_infected_2 <- numeric(16)
total_infected_3 <- numeric(16)
total_infected_4 <- numeric(16)
total_men_1 <- numeric(16)
total_men_2 <- numeric(16)
total_men_3 <- numeric(16)
total_men_4 <- numeric(16)
total_women_1 <- numeric(16)
total_women_2 <- numeric(16)
total_women_3 <- numeric(16)
total_women_4 <- numeric(16)


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

# Count cumulative infected individuals (currently infected + recovered) for Scenario 2
for (i in 1:16) {
  time_index <- i * 3  # Time point at the end of each year
  total_infected_3[i] <- out3[time_index, "I_men"] + out3[time_index, "I_women"]
  total_men_3[i] <- out3[time_index, "I_men"]
  total_women_3[i] <- out3[time_index, "I_women"]
}

# Count cumulative infected individuals (currently infected + recovered) for Scenario 2
for (i in 1:16) {
  time_index <- i * 3  # Time point at the end of each year
  total_infected_4[i] <- out4[time_index, "I_men"] + out4[time_index, "I_women"]
  total_men_4[i] <- out4[time_index, "I_men"]
  total_women_4[i] <- out4[time_index, "I_women"]
}


# Print the cumulative infected counts for Scenario 1
print("Cumulative Infected Counts current Vaccination rate (Every 3 Months):")
print(round(sum(total_infected_1),0))
print("Cumulative Infected Counts matched Vaccination rate (Every 3 Months):")
print(round(sum(total_infected_2),0))
print("Cumulative Infected Counts 75% Vaccination rate (Every 3 Months):")
print(round(sum(total_infected_3),0))
print("Cumulative Infected men Counts Scenario 1 (over 4 years):")
print(round(sum(total_men_1),0))
print(round(sum(total_women_1), 0))
print("Cumulative Infected men Counts Scenario 2 (Every 3 Months):")
print(round(sum(total_men_2),0))
print(round(sum(total_women_2),0))
print("Cumulative Infected men Counts Scenario 2 (Every 3 Months):")
print(round(sum(total_men_3),0))
print(round(sum(total_women_3),0))
print(round(sum(total_men_4),0))
print(round(sum(total_women_4),0))


library(ggplot2)
library(tidyr)
# Prepare the data for combined Scenario 1 and Scenario 2 (only Infected)
data_combined1 <- data.frame(
  time = out1[, "time"],
  I_men_true_rate = out1[, "I_men"],         # Scenario 1
  I_women_true_rate = out1[, "I_women"],      # Scenario 1
  I_men_equal_rate = out2[, "I_men"],
  I_women_equal_rate = out2[, 'I_women']
)

# Convert to long format
data_combined_long1 <- gather(data_combined1, key = "Population", value = "Count", -time)

# Define colors and line types
colors_combined1 <- c("I_men_true_rate" = "blue", "I_women_true_rate" = "purple", 
                     "I_men_equal_rate" = "red", 'I_women_equal_rate' = 'darksalmon')

# Update the line types to valid ggplot2 types
line_types_combined1<- c("I_men_true_rate" = "solid", "I_women_true_rate" = "dashed", 
                         "I_men_equal_rate" = "solid",'I_women_equal_rate' = "dashed")

# Plot for combined Scenario 1 and Scenario 2 (Infected)
ggplot(data_combined_long1, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined1) +
  scale_linetype_manual(values = line_types_combined1) +
  labs(x = "Time in Months", y = "Infected Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom",axis.title.x = element_text(size = 20),  # Change size of x-axis label
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 11),   # Change size of legend title
        legend.text = element_text(size = 12))

data_combined2 <- data.frame(
  time = out1[, "time"],
  I_men_true_rate = out1[, "I_men"],         # Scenario 1
  I_women_true_rate = out1[, "I_women"],      # Scenario 1
  I_men_equal_rate = out2[, "I_men"],
  I_women_equal_rate = out2[, 'I_women'],
  I_women_75 = out3[, 'I_women'],
  I_men_75 = out3[, 'I_men']
)

# Convert to long format
data_combined_long2 <- gather(data_combined2, key = "Population", value = "Count", -time)

# Define colors and line types
colors_combined2 <- c("I_men_true_rate" = "blue", "I_women_true_rate" = "purple", 
                     "I_men_equal_rate" = "red", 'I_women_equal_rate' = 'darksalmon', 'I_women_75' = 'firebrick',
                     "I_men_75" = 'darkgreen') 

# Update the line types to valid ggplot2 types
line_types_combined2 <- c("I_men_true_rate" = "solid", "I_women_true_rate" = "dashed", 
                         "I_men_equal_rate" = "solid", 'I_women_75' = 'solid', "I_men_75" = 'solid'
                         , "I_women_equal_rate" = "dashed")

# Plot for combined Scenario 1 and Scenario 2 (Infected)
ggplot(data_combined_long2, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined2) +
  scale_linetype_manual(values = line_types_combined2) +
  labs(x = "Time in Months", y = "Infected Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom",axis.title.x = element_text(size = 20),  # Change size of x-axis label
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 12),   # Change size of legend title
        legend.text = element_text(size = 13))

# Prepare the data for combined Scenario 1 and Scenario 2 (only Infected)
data_combined3 <- data.frame(
  time = out1[, "time"],
  I_men_true_rate = out1[, "I_men"],         # Scenario 1
  I_women_true_rate = out1[, "I_women"],      # Scenario 1
  I_men_equal_rate = out2[, "I_men"],
  I_women_equal_rate = out2[, "I_women"],
  I_men_20 = out4[, 'I_men'],
  I_women_20 = out4[,"I_women"]
)

  # Convert to long format
data_combined_long3 <- gather(data_combined3, key = "Population", value = "Count", -time)

# Define colors and line types
colors_combined3 <- c("I_men_true_rate" = "blue", "I_women_true_rate" = "purple", 
                     "I_men_equal_rate" = "red", 'I_women_equal_rate' = 'darksalmon',"I_women_20" = 'orange', "I_men_20" = 'green')

# Update the line types to valid ggplot2 types
line_types_combined3 <- c("I_men_true_rate" = "solid", "I_women_true_rate" = "dashed", 
                         "I_men_equal_rate" = "solid", "I_women_20" ='dashed',"I_men_20" = 'dashed'
                         ,"I_women_equal_rate" = "dashed")
# Plot for combined Scenario 1 and Scenario 2 (Infected)
ggplot(data_combined_long3, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined3) +
  scale_linetype_manual(values = line_types_combined3) +
  labs(x = "Time in Months", y = "Infected Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom",axis.title.x = element_text(size = 20),  # Change size of x-axis label
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 12),   # Change size of legend title
        legend.text = element_text(size = 13))

# Prepare the data for combined Scenario 1 and Scenario 2 (only Infected)
data_combined4 <- data.frame(
  time = out1[, "time"],
  S_men_true_rate = out1[, "S_men"],         # Scenario 1
  S_women_true_rate = out1[, "S_women"],     # Scenario 1
  S_men_equal_rate = out2[, "S_men"],
  S_women_equal_rate = out2[, "S_women"],
  S_men_20 = out4[, 'S_men'],
  S_men_75 = out3[, 'S_men'],
  S_women_75 = out3[,'S_women'],
  S_women_20 = out4[, "S_women"]
)

# Convert to long format
data_combined_long4 <- gather(data_combined4, key = "Population", value = "Count", -time)

# Define colors and line types with correct column names
colors_combined4 <- c("S_men_true_rate" = "blue", "S_women_true_rate" = "purple", 
                      "S_men_equal_rate" = "red", 'S_women_equal_rate' = 'darksalmon', "S_women_20" = 'orange', "S_men_20" = 'green',
                      'S_women_75' = 'firebrick', "S_men_75" = 'darkgreen')
line_types_combined4 <- c("S_men_true_rate" = "solid", "S_women_true_rate" = "dashed", 
                          "S_men_equal_rate" = "solid", "S_women_20" ='dashed', "S_men_20" = 'dashed',
                          'S_women_75' = 'solid', "S_men_75" = 'solid', "S_women_equal_rate" = "dashed")

# Plot for combined Scenario 1 and Scenario 2 (Infected)
ggplot(data_combined_long4, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined4) +
  scale_linetype_manual(values = line_types_combined4) +
  labs(x = "Time in Months", y = "Susceptible Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom",axis.title.x = element_text(size = 20),  # Change size of x-axis label
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 12),   # Change size of legend title
        legend.text = element_text(size = 13))

# Prepare the data for combined Scenario 1 and Scenario 2 (only Vaccinated)
data_combined_vaccinated <- data.frame(
  time = out1[, "time"],
  V_men_true_rate = out1[, "V_men"],         # Scenario 1
  V_women_true_rate = out1[, "V_women"],     # Scenario 1
  V_men_equal_rate = out2[, "V_men"],
  V_women_equal_rate = out2[,"V_women"],     # Scenario 2
  V_men_75 = out3[,'V_men'],
  V_women_75 = out3[,'V_women'],
  V_men_20 = out4[, 'V_men'],                # Another scenario (20)
  V_women_20 = out4[,"V_women"]              # Another scenario (20)
)

# Convert to long format
data_combined_long_vaccinated <- gather(data_combined_vaccinated, key = "Population", value = "Count", -time)

# Define colors and line types for vaccinated data
colors_combined_vaccinated <- c("V_men_true_rate" = "blue", "V_women_true_rate" = "purple", 
                                "V_men_equal_rate" = "red", "V_women_20" = 'orange', "V_men_20" = "green",
                                "V_women_75" = "firebrick", 'V_men_75' = 'darkgreen', "V_women_equal_rate" = 'darksalmon')
line_types_combined_vaccinated <- c("V_men_true_rate" = "solid", "V_women_true_rate" = "dashed", 
                                    "V_men_equal_rate" = "solid", "V_women_20" = 'dashed', "V_men_20" = 'dashed',
                                    "V_men_75" = 'solid', "V_women_75" = "solid", "V_women_equal_rate" = "dashed")

# Plot for combined Scenario 1 and Scenario 2 (Vaccinated)
ggplot(data_combined_long_vaccinated, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined_vaccinated) +
  scale_linetype_manual(values = line_types_combined_vaccinated) +
  labs(x = "Time", y = "Vaccinated Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.x = element_text(size = 20),  # Change size of x-axis label
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 12),   # Change size of legend title
        legend.text = element_text(size = 13))



# Print the cumulative infected counts for Scenario 1
print("Cumulative Infected Counts Scenario 1 (Every 12 Months):")
print(round(sum(total_infected_1),0))
print("Cumulative Infected Counts Scenario 2 (Every 12 Months):")
print(round(sum(total_infected_2),0))
print("Cumulative Infected Counts Scenario 3 (Every 12 Months):")
print(round(sum(total_infected_3),0))
print("Cumulative Infected Counts Scenario 4 (Every 12 Months):")
print(round(sum(total_infected_4),0))
 
