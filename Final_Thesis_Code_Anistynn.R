# Install all required libraries
install.packages('deSolve')
install.packages('pracma')
install.packages('ggplot2')
install.packages('dplyr')
install.packages('tidyr')
install.packages('patchwork')

# Load required libraries
library(deSolve)
library(pracma)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

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

# Set initial state for men and women
state <- c(S_men = 60000*0.485, I_men = 1000, V_men = 0,
           S_women = 60000*0.515, I_women = 1000, V_women = 0)
pulse_amount <- 9000
# Set time points for simulation
times <- seq(0, 4*12, by = 1)



### Scenario 1: Separate alpha values for men and women (47.9)
parameters1a <- c(beta = 0.8795*0.479, alpha_men = 0.34, alpha_women = 0.552, lambda=lambda)
out1a <- ode(y = state, times = times, func = DiffEq, parms = parameters1a)

### Scenario 1: Separate alpha values for men and women (40.3%)
parameters1b <- c(beta = 0.8795*0.403, alpha_men = 0.34, alpha_women = 0.552, lambda=lambda)
out1b <- ode(y = state, times = times, func = DiffEq, parms = parameters1b)

### Scenario 1: Separate alpha values for men and women (32.5%)
parameters1c <- c(beta = 0.8795*0.325, alpha_men = 0.34, alpha_women = 0.552, lambda=lambda)
out1c <- ode(y = state, times = times, func = DiffEq, parms = parameters1c)


### Scenario 3: both sexes have vaccination rates of 0.75 (47.9%)
parameters3a <- c(beta = 0.8795*0.479, alpha_men = 0.75, alpha_women = 0.75, lambda=lambda)
out3a <- ode(y=state, times = times, func=DiffEq, parms = parameters3a)

### Scenario 3: both sexes have vaccination rates of 0.75 (40.3%)
parameters3b <- c(beta = 0.8795*0.403, alpha_men = 0.75, alpha_women = 0.75, lambda=lambda)
out3b <- ode(y=state, times = times, func=DiffEq, parms = parameters3b)

### Scenario 3: both sexes have vaccination rates of 0.75 (32.5%)
parameters3c <- c(beta = 0.8795*0.325, alpha_men = 0.75, alpha_women = 0.75, lambda=lambda)
out3c <- ode(y=state, times = times, func=DiffEq, parms = parameters3c)



# Initialize vectors to store total infected at the end of each year (12, 24, 36, 48 months)
total_infected_1a <- numeric(16)
total_infected_1b <- numeric(16)
total_infected_1c <- numeric(16)

total_infected_3a <- numeric(16)
total_infected_3b <- numeric(16)
total_infected_3c <- numeric(16)

total_men_1a <- numeric(16)
total_men_1b <- numeric(16)
total_men_1c <- numeric(16)

total_men_3a <- numeric(16)
total_men_3b <- numeric(16)
total_men_3c <- numeric(16)

total_women_1a <- numeric(16)
total_women_1b <- numeric(16)
total_women_1c <- numeric(16)

total_women_3a <- numeric(16)
total_women_3b <- numeric(16)
total_women_3c <- numeric(16)


# Count cumulative infected individuals (currently infected + recovered) for Scenario 1 (47.9%)
for (i in 1:16) {
  time_index <- i * 3 # Time point at the end of each year
  total_infected_1a[i] <- out1a[time_index, "I_men"] + out1a[time_index, "I_women"]
  total_men_1a[i] <- out1a[time_index, "I_men"]
  total_women_1a[i] <- out1a[time_index, "I_women"]
}

# Count cumulative infected individuals (currently infected + recovered) for Scenario 1 (40.3%)
for (i in 1:16) {
  time_index <- i * 3 # Time point at the end of each year
  total_infected_1b[i] <- out1b[time_index, "I_men"] + out1b[time_index, "I_women"]
  total_men_1b[i] <- out1b[time_index, "I_men"]
  total_women_1b[i] <- out1b[time_index, "I_women"]
}

# Count cumulative infected individuals (currently infected + recovered) for Scenario 1 (32.5%)
for (i in 1:16) {
  time_index <- i * 3 # Time point at the end of each year
  total_infected_1c[i] <- out1c[time_index, "I_men"] + out1c[time_index, "I_women"]
  total_men_1c[i] <- out1c[time_index, "I_men"]
  total_women_1c[i] <- out1c[time_index, "I_women"]
}


# Count cumulative infected individuals (currently infected + recovered) for Scenario 3
for (i in 1:16) {
  time_index <- i * 3  # Time point at the end of each year
  total_infected_3a[i] <- out3a[time_index, "I_men"] + out3a[time_index, "I_women"]
  total_men_3a[i] <- out3a[time_index, "I_men"]
  total_women_3a[i] <- out3a[time_index, "I_women"]
}

# Count cumulative infected individuals (currently infected + recovered) for Scenario 3
for (i in 1:16) {
  time_index <- i * 3  # Time point at the end of each year
  total_infected_3b[i] <- out3b[time_index, "I_men"] + out3b[time_index, "I_women"]
  total_men_3b[i] <- out3b[time_index, "I_men"]
  total_women_3b[i] <- out3b[time_index, "I_women"]
}

# Count cumulative infected individuals (currently infected + recovered) for Scenario 3
for (i in 1:16) {
  time_index <- i * 3  # Time point at the end of each year
  total_infected_3c[i] <- out3c[time_index, "I_men"] + out3c[time_index, "I_women"]
  total_men_3c[i] <- out3c[time_index, "I_men"]
  total_women_3c[i] <- out3c[time_index, "I_women"]
}

## INFECTED
# Prepare the data for combined Scenario 1 and Scenario 2 (only Infected)
data_combined1 <- data.frame(
  time = out1a[, "time"],
  I_men_47.9 = out1a[, "I_men"],         # Scenario 1
  I_women_47.9 = out1a[, "I_women"],      # Scenario 1
  I_men_40.3 = out1b[, "I_men"],
  I_women_40.3 = out1b[, 'I_women'],
  I_men_32.5 = out1c[, "I_men"],
  I_women_32.5 = out1c[, 'I_women']
)

# Convert to long format
data_combined_long1 <- gather(data_combined1, key = "Population", value = "Count", -time)

# Define colors and line types
colors_combined <- c("I_men_47.9" = "darkblue", "I_women_47.9" = "purple", 
                      "I_men_40.3" = "blue", 'I_women_40.3' = 'magenta',
                      "I_men_32.5" = "lightblue", 'I_women_32.5' = 'pink')

# Update the line types to valid ggplot2 types
line_types_combined<- c("I_men_47.9" = "solid", "I_women_47.9" = "dashed", 
                         "I_men_40.3" = "solid", 'I_women_40.3' = 'dashed',
                         "I_men_32.5" = "solid", 'I_women_32.5' = 'dashed')

# Plot for combined Scenario 1 and Scenario 2 (Infected)
iplot1 <- ggplot(data_combined_long1, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined) +
  scale_linetype_manual(values = line_types_combined) +
  labs(title='Given Current Rate of Vaccination', x = "Time", y = "Infected Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("plot1_all_encouter_rates.png", width = 8, height = 6, dpi = 300)

# Prepare the data for combined Scenario 1 and Scenario 2 (only Infected)
data_combined3 <- data.frame(
  time = out3a[, "time"],
  I_men_47.9 = out3a[, "I_men"],         # Scenario 1
  I_women_47.9 = out3a[, "I_women"],      # Scenario 1
  I_men_40.3 = out3b[, "I_men"],
  I_women_40.3 = out3b[, 'I_women'],
  I_men_32.5 = out3c[, "I_men"],
  I_women_32.5 = out3c[, 'I_women']
)

# Convert to long format
data_combined_long3 <- gather(data_combined3, key = "Population", value = "Count", -time)

# Plot for combined Scenario 1 and Scenario 2 (Infected)
iplot2 <- ggplot(data_combined_long3, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined) +
  scale_linetype_manual(values = line_types_combined) +
  labs(title='Given Target Rate of Vaccination', x = "Time", y = "Infected Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("plot2_all_encouter_rates.png", width = 8, height = 6, dpi = 300)

## SUSCEPTIBLE
# Prepare the data for combined Scenario 1 and Scenario 2 (only Vaccinated)
data_combined_susceptible1 <- data.frame(
  time = out1a[, "time"],
  S_men_47.9 = out1a[, "S_men"],         # Scenario 1
  S_women_47.9 = out1a[, "S_women"],      # Scenario 1
  S_men_40.3 = out1b[, "S_men"],
  S_women_40.3 = out1b[, 'S_women'],
  S_men_32.5 = out1c[, "S_men"],
  S_women_32.5 = out1c[, 'S_women']
)

# Convert to long format
data_combined_long_susceptible1 <- gather(data_combined_susceptible1, key = "Population", value = "Count", -time)

# Define colors and line types
colors_combined_susceptible <- c("S_men_47.9" = "darkblue", "S_women_47.9" = "purple", 
                                "S_men_40.3" = "blue", 'S_women_40.3' = 'magenta',
                                "S_men_32.5" = "lightblue", 'S_women_32.5' = 'pink')

# Update the line types to valid ggplot2 types
line_types_combined_susceptible <- c("S_men_47.9" = "solid", "S_women_47.9" = "dashed", 
                                    "S_men_40.3" = "solid", 'S_women_40.3' = 'dashed',
                                    "S_men_32.5" = "solid", 'S_women_32.5' = 'dashed')

# Plot for combined Scenario 1 and Scenario 2 (Vaccinated)
splot1 <- ggplot(data_combined_long_susceptible1, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined_susceptible) +
  scale_linetype_manual(values = line_types_combined_susceptible) +
  labs(title='Given Current Rate of Vaccination', x = "Time", y = "Susceptible Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("plot1_susceptible_all_encouter_rates.png", width = 8, height = 6, dpi = 300)


# Prepare the data for combined Scenario 1 and Scenario 2 (only Vaccinated)
data_combined_susceptible3 <- data.frame(
  time = out1a[, "time"],
  S_men_47.9 = out3a[, "S_men"],         # Scenario 1
  S_women_47.9 = out3a[, "S_women"],      # Scenario 1
  S_men_40.3 = out3b[, "S_men"],
  S_women_40.3 = out3b[, 'S_women'],
  S_men_32.5 = out3c[, "S_men"],
  S_women_32.5 = out3c[, 'S_women']
)

# Convert to long format
data_combined_long_susceptible3 <- gather(data_combined_susceptible3, key = "Population", value = "Count", -time)

# Plot for combined Scenario 1 and Scenario 2 (Vaccinated)
splot2 <- ggplot(data_combined_long_susceptible3, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined_susceptible) +
  scale_linetype_manual(values = line_types_combined_susceptible) +
  labs(title='Given Target Rate of Vaccination', x = "Time", y = "Susceptible Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("plot2_susceptible_all_encouter_rates.png", width = 8, height = 6, dpi = 300)


## VACCINATED
# Prepare the data for combined Scenario 1 and Scenario 2 (only Vaccinated)
data_combined_vaccinated1 <- data.frame(
  time = out1a[, "time"],
  V_men_47.9 = out1a[, "V_men"],         # Scenario 1
  V_women_47.9 = out1a[, "V_women"],      # Scenario 1
  V_men_40.3 = out1b[, "V_men"],
  V_women_40.3 = out1b[, 'V_women'],
  V_men_32.5 = out1c[, "V_men"],
  V_women_32.5 = out1c[, 'V_women']
)

# Convert to long format
data_combined_long_vaccinated1 <- gather(data_combined_vaccinated1, key = "Population", value = "Count", -time)

# Define colors and line types
colors_combined_vaccinated <- c("V_men_47.9" = "darkblue", "V_women_47.9" = "purple", 
                     "V_men_40.3" = "blue", 'V_women_40.3' = 'magenta',
                     "V_men_32.5" = "lightblue", 'V_women_32.5' = 'pink')

# Update the line types to valid ggplot2 types
line_types_combined_vaccinated <- c("V_men_47.9" = "solid", "V_women_47.9" = "dashed", 
                        "V_men_40.3" = "solid", 'V_women_40.3' = 'dashed',
                        "V_men_32.5" = "solid", 'V_women_32.5' = 'dashed')

# Plot for combined Scenario 1 and Scenario 2 (Vaccinated)
vplot1 <- ggplot(data_combined_long_vaccinated1, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined_vaccinated) +
  scale_linetype_manual(values = line_types_combined_vaccinated) +
  labs(title='Given Current Rate of Vaccination', x = "Time", y = "Vaccinated Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("plot1_vaccinated_all_encouter_rates.png", width = 8, height = 6, dpi = 300)


# Prepare the data for combined Scenario 1 and Scenario 2 (only Vaccinated)
data_combined_vaccinated3 <- data.frame(
  time = out1a[, "time"],
  V_men_47.9 = out3a[, "V_men"],         # Scenario 1
  V_women_47.9 = out3a[, "V_women"],      # Scenario 1
  V_men_40.3 = out3b[, "V_men"],
  V_women_40.3 = out3b[, 'V_women'],
  V_men_32.5 = out3c[, "V_men"],
  V_women_32.5 = out3c[, 'V_women']
)

# Convert to long format
data_combined_long_vaccinated3 <- gather(data_combined_vaccinated3, key = "Population", value = "Count", -time)

# Plot for combined Scenario 1 and Scenario 2 (Vaccinated)
vplot2 <- ggplot(data_combined_long_vaccinated3, aes(x = time, y = Count, color = Population, linetype = Population)) +
  geom_line(size = 1) +
  scale_color_manual(values = colors_combined_vaccinated) +
  scale_linetype_manual(values = line_types_combined_vaccinated) +
  labs(title='Given Target Rate of Vaccination', x = "Time", y = "Vaccinated Population", color = "Population") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("plot2_vaccinated_all_encouter_rates.png", width = 8, height = 6, dpi = 300)


# Combine the plots horizontally
combined_plot1 <- iplot1 | splot1 | vplot1
# Save the combined plot
ggsave("combined_plot1_horizontal.png", plot = combined_plot1, width = 18, height = 6, dpi = 300)

# Combine the plots horizontally
combined_plot2 <- iplot2 | splot2 | vplot2
# Save the combined plot
ggsave("combined_plot2_horizontal.png", plot = combined_plot2, width = 18, height = 6, dpi = 300)

# Print the cumulative infected counts for Scenario 1
print("Cumulative Infected Counts at Current Rate of Vaccination Annually at 47.9% of Sexual Encouters being Unprotected:")
print(round(sum(total_infected_1a),0))
print("Cumulative Infected Counts at Current Rate of Vaccination Annually at 40.3% of Sexual Encouters being Unprotected:")
print(round(sum(total_infected_1b),0))
print("Cumulative Infected Counts at Current Rate of Vaccination Annually at 32.5% of Sexual Encouters being Unprotected:")
print(round(sum(total_infected_1c),0))
print("Cumulative Infected Counts at Target Rate of Vaccination Annually at 47.9% of Sexual Encouters being Unprotected:")
print(round(sum(total_infected_3a),0))
print("Cumulative Infected Counts at Target Rate of Vaccination Annually at 40.3% of Sexual Encouters being Unprotected:")
print(round(sum(total_infected_3b),0))
print("Cumulative Infected Counts at Target Rate of Vaccination Annually at 32.5% of Sexual Encouters being Unprotected:")
print(round(sum(total_infected_3c),0))

