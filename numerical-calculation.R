# setup ----

library(tidyverse)
library(NonCompart) # install.packages('NonCompart')

# https://pharmacy.ufl.edu/files/2013/01/5127-28-equations.pdf

conc_time_steady_state <- function(time, dose, bioavailability, ka, ke, V){
  (bioavailability * dose * ka) / (V * (ka-ke)) * {exp(-ke*time)/(1-exp(-ke*tau)) - exp(-ka*time)/(1-exp(-ka*tau))} 
}

conc_time <- function(time, dose, bioavailability){
  (bioavailability * dose * ka) / (V * (ka-ke)) * {exp(-ke*time) - exp(-ka*time)} 
}

dose_group <- c(1500, 3000, 4500, 6000)
bioavailability_group <- c(1, 1, 1, 1.4)
cmax <- conc_time_steady_state(tmax, dose_group, bioavailability_group)

tau <- 12

V <- 236
ka <- 3.12
CL <- 120 #CL <- V * ke

up <- ka * (1 - exp(-ke*tau))
down <- ke * (1 - exp(-ka*tau))

tmax <- log(up/down) / (ka-ke)
tmax
FullCov <- matrix(c(0.42^2,0,0,
                    0, 0.62^2,0,
                    0,0, 1.57^2), nrow = 3, byrow=TRUE)

FullCov <- matrix(c(1, 0, 0,
                    0, 1, 0,
                    0,0, 1), nrow = 3, byrow=TRUE)

sim_n <- 120
rpk <- MASS::mvrnorm(sim_n, rep(0, 3), FullCov)
rpk_value <- matrix(rep(c(CL, V, ka), sim_n), nrow=sim_n, byrow=TRUE) * exp(rpk)


ipk <- rpk_value %>% 
  as.data.frame() %>% 
  dplyr::select(CL = V1, V = V2, ka = V3) %>% 
  dplyr::mutate(ke = CL/V) %>% 
  dplyr::mutate(SUBJ = dplyr::row_number()) %>% 
  mutate(dose = rep(dose_group, 30)) %>% 
  as_tibble() %>% 
  print()

time_dose <- expand.grid(x= seq(0, 24, by = 0.1), dose = dose_group) %>% 
  mutate(bioavailability = ifelse(dose == 6000, 1.4, 1)) %>% 
  mutate(dose_label = sprintf('%s mg', dose) %>% as_factor()) %>% 
  as_tibble() %>% 
  print()

intermediate_df <- ipk %>% 
  left_join(time_dose) %>% 
  mutate(yss = conc_time_steady_state(x, dose, bioavailability, ka, ke, V)) %>%
  arrange(SUBJ) %>% 
  print()

final_df <- intermediate_df %>% 
  select(SUBJ, dose_label, x, yss) %>% 
  group_by(dose_label, x) %>% 
  summarise(yss = median(yss)) %>% 
  print()

plot_conc_time <- ggplot(final_df, aes(color = dose_label)) +
  geom_line(aes(x,yss)) +
  facet_grid(.~dose_label)

plot_conc_time

# df ----


df <- time_dose %>% 
  arrange(dose, x) %>% 
  mutate(y = conc_time(x, dose, bioavailability)) %>% 
  mutate(yss = conc_time_steady_state(x, dose, bioavailability)) %>% 
  print()

plot_conc_time_combined <- ggplot(df, aes(group = dose_label, color = dose_label)) +
  geom_line(aes(x,yss)) 

plot_conc_time_combined


# test ----




# NCA ----

tblNCA(as.data.frame(df), 'dose', 'x', 'yss')

exp(rpk)


