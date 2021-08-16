# Code Description ############################################################
# Date: 09_08_2021
# Written by: Paul McMenemy, Matthew Baister, Ewan McTaggart
# Script to model COVID dynamics using meta-populations of 
# VULNERABLE population segment - CARE HOME RESIDENTS (C)
# SHIELDING population segment - CARE HOME STAFF (S) and
# REST OF population segment - REST OF POPULATION (R)
# 
# This script calls the SEIRD functions from script 1
# and generates data 
#
# 1) Preamble - Source any script with SEIRD functions
# 2) Define parameter values/vectors, and generate matrices/lists to use in SEIRD functions
# 3) Generate data - results of solving ODE system


# 1) Preamble ---------------------------------------------------

# wd_path <-
# setwd(wd_path)

# import the SEIRD functions from the R script:
source('1_SEIRD_FUNCTIONS.R')

# 2) Parameters -------------------------------------------------

# Reproduction rate values for R_matrix(ces) in form:

#        C           S          R
#   - - - - - - - - - - - - - - - - -
#  |           |          |          |
#C |           |          |          |
#  |           |          |          |
#   - - - - - - - - - - - - - - - - -
#  |           |          |          |
#S |           |          |          |
#  |           |          |          |
#   - - - - - - - - - - - - - - - - -
#  |           |          |          |
#R |           |          |          |
#  |           |          |          |
#   - - - - - - - - - - - - - - - - -
seeded_4home_E_best_fit$All_params

baseRHigh_pars <- list(
  
  beta_c = 4.7 , # internal Care home reproduction number
  beta_cc = 0,  # between Care home reproduction number
  
  beta_cs = 4.7, # Care home to Shielder/Staff reproduction number
  beta_sc = 4.7, # Staff to Care home reproduction number (assume equal to beta_cs)
  
  
  # Remaining Shielder/Staff transmission rates 
  
  beta_s = (4.7+4.1)/2,   # internal Staff reproduction number
  beta_ss = (4.7 +4.1)/2,  # between Staff reproduction number (could be higher than beta_rr due to mixing in workplace?)
  
  # Remaining Rest of Population transmission rates 
  beta_rr = 4.1 # rest of population to rest of population reproduction number    #4.56 is latency = 5.8 instead of 8
)

baseRLow_pars <- list(
  
  beta_c  = 0.6, # internal Care home reproduction number
  beta_cc = 0,  # between Care home reproduction number
  
  beta_cs = 0.6,# Care home to Shielder/Staff treproduction number
  beta_sc = 0.6,# Staff to Care home reproduction number (assume equal to beta_cs)
  
  # Remaining Shielder/Staff reproduction number
  
  beta_s = 0.6, # internal Staff reproduction number
  beta_ss = 0.6,# between Staff reproduction number (could be higher than beta_rr due to mixing in workplace?)
  
  # Remaining Rest of Population reproduction numbers
  beta_rr = 0.6 # rest of population to rest of population reproduction number    #0.7 if latency = 5.8 instead of 8
)

baseStart_pars <- list(
  
  c = -82,    # internal Care home logi function start time
  cc = -82,   # between Care home logi function start time
  cs = -82,   # Care home to Staff logi function start time
  sc = -82,   # Staff to Care home logi function start time
  
  # Remaining Shielder/Staff logi function start times
  
  s = -82,    # internal Staff logi function start time
  ss = -82,   # between Staff logi function start time
  
  # Note - Staff to rest of population logi function start time = rr
  rr = -82    # rest of population to rest of population logi function start time
)

baseEnd_pars <- list(
  
  c = 42,   # internal Care home logi function end time
  cc = 42,  # between Care home logi function end time
  cs = 42,  # Care home to Shielder/Staff logi function end time
  sc = 42,  # Staff to Care home logi function end time
  
  # Remaining Shielder/Staff logi function end times
  
  s = 42,   # internal Staff logi function end time
  ss = 42,  # between Staff logi function end time
  # # Note - Staff to rest of population logi function end time = rr
  rr = 22  # rest of population to rest of population logi function end time
)


baseRate_pars <- list(
  
  c = 0.5,  # internal Care home logi function rate value
  cc =0.5, # between Care home logi function rate value
  cs =0.5, # Care home to Shielder/Staff logi function rate value
  sc =0.5, # Staff to Care home logi function rate value
  
  # Remaining Shielder/Staff logi function rate value
  
  s = 0.5,  # internal Staff logi function rate value
  ss =0.5, # between Staff logi function rate value
  rr = 0.5  # rest of population to rest of population logi function rate value
)

## SECOND WAVE parameters ######################################## 

pulseRHigh_pars <- list(
  
  # NOTE: the pulseRHigh value is summed with the respective baseRLow values,
  # so take this into account when parameterising these values
  
  beta_c = 0, # internal Care home reproduction number
  beta_cc = 0,# between Care home reproduction number
  beta_cs = 0,# Care home to Shielder/Staff reproduction number
  beta_sc = 0,# Staff to Care home reproduction number (assume equal to beta_cs)
   
  # Remaining Shielder/Staff transmission rates 
  
  beta_s = 0, # internal staff reproduction number
  beta_ss = 0,# between Staff reproduction number (could be higher than beta_rr due to mixing in workplace?)
  
  # Remaining Rest of Population transmission rates 
  beta_rr = 0 # rest of population to rest of population reproduction number
)

pulseRLow_pars <- list(
  
  # note that, as we sum the base and pulse values, then we need only zero values for the pulse,
  # otherwise we would be summing the (baseRLow + pulseRLow) values outside the main waves of infection
  
  beta_c = 0, # internal Care home reproduction number
  beta_cc = 0,# between Care home reproduction number
  beta_cs = 0,# Care home to Shielder/Staff treproduction number
  beta_sc = 0,# Staff to Care home reproduction number (assume equal to beta_cs)
  
  # Remaining Shielder/Staff reproduction number
  
  beta_s = 0, # internal staff reproduction number
  beta_ss = 0,# between Staff reproduction number (could be higher than beta_rr due to mixing in workplace?)

  # Remaining Rest of Population reproduction numbers
  beta_rr = 0 # rest of population to rest of population reproduction number
)

pulseStart_pars <- list(
  
  c = 54,    # internal Care home logi function start time
  cc = 54,   # between Care home logi function start time
  cs = 54,   # Care home to Staff logi function start time
  sc = 54, # Staff to Care home logi function start time
  
  # Remaining Shielder/Staff logi function start times
  
  s = 54,  # internal Staff logi function start time
  ss = 54, # between Staff logi function start time
  
  # Note - Staff to rest of population logi function start time = rr
  rr = 54  # rest of population to rest of population logi function start time
  
)

pulseEnd_pars <- list(
  
  c = 206,   # internal Care home logi function end time
  cc = 206,  # between Care home logi function end time
  cs = 206,  # Care home to Shielder/Staff logi function end time
  sc = 206,  # Staff to Care home logi function end time
  
  # Remaining Shielder/Staff logi function end times
  
  s = 206,   # internal Staff logi function end time
  ss = 206,  # between Staff logi function end time
  
  # Note - Staff to rest of population logi function end time = rr
  rr = 206   # rest of population to rest of population logi function end time
)

pulseRate_pars <- list(
  
  c = 0.1,  # internal Care home logi function rate value
  cc = 0.1, # between Care home logi function rate value
  cs = 0.1, # Care home to Shielder/Staff logi function rate value
  sc = 0.1, # Staff to Care home logi function rate value
  
  # Remaining Shielder/Staff logi function rate value
  
  s = 0.1,  # internal Staff logi function rate value 
  ss = 0.1, # between Staff logi function rate value
  rr = 0.2  # rest of population to rest of population logi function rate value
)
 
## BOTH WAVE parameters ######################

# Parameters which influence reproductive rates between 
# populations over both waves

between_pars <- list( 

  alpha = 0         # Staff shielding from Rest of population. Value between 0 and 1 where 0 is no difference and 1 is full shielding from Rest of population
)


## OTHER parameters ######################################## 
 
 
simTime = seq(from = 0, by = 1, length.out = 102)   # EWISOTT
 
# number of sub-populations per segment
seg_pars <- list(   # EWISOTT
  
  C_pops = 109, 
  S_pops = 109,       # = C_pops, there is one shielder/staff sub-population per care home.
  R_pops = 1
)

all_pops = Reduce('+', seg_pars)


sub_pop_sizes <- list(    # internal subpopulation sizes 
  
  C_internal = 48,                         # population of a single CH
  S_internal = 48,                         # population of a single staff population # = C_internal (assuming 1:1 resident:staff member ratio per home)
  R_internal = 907580-2*48*seg_pars$C_pops # population of a single rest population # (assuming CH resident population is 0.57-0.6% of population)
)




N = sum(unlist(sub_pop_sizes)*unlist(seg_pars))      # total population
                                                     # "Mid-Year Population Estimates, Mid-2019" - National Records of Scotland Website 
                                                     # https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/population/population-estimates/mid-year-population-estimates/mid-2019
                                                     # NHS Lothian had population of 907580

# (sub_pop_sizes$C_internal*seg_pars$C_pops)/N # CH population as proportion of Lothian

# probDeath = c(0.32,0.011,0.011) # vector of mortality rates in order of C,S,R

probDeath = c(0.25,0.017,0.017) # vector of mortality rates in order of C,S,R

# CH death rate taken from Lothian data
# sum(Lothian_CH_death_data$`COVID-19 deaths 2020`)/(sum(Lothian_CH_7_day_positive_data$`7 day rolling average positive tests`)*100/68)

# Rest death rate
# sum(Weekly_Deaths_Lothian$WeeklyDeaths)/(sum(Lothian$DailyPositive)*100/7.7)

asymptomatic_proportions <- list(    
  
  # proportion of exposed that become asymptomatic (or go unreported) in each population
  # set all = 0 to get rid of asymptomatic compartment
  
  C_asymp = 1-0.53,               # proportion of exposed that become asymptomatic (or go unreported) in care homes -  https://www.gov.scot/publications/foi-202000067854/        https://www.careinspectorate.com/images/COVID-19_-_Letter_from_Cabinet_Secretary_for_Health_and_Sport_-_Social_care_guidance_-_13_March_2020.pdf https://web.archive.org/web/20200313234513/https://hpspubsrepo.blob.core.windows.net/hps-website/nss/2980/documents/1_COVID-19%20Guidance-for-Social-or-community-care-and-residentail-settings.pdf   
  S_asymp = 1-0.52,                # proportion of exposed that become asymptomatic (or go unreported) in staff/shielders
  R_asymp = 1-0.077                 # proportion of exposed that become asymptomatic (or go unreported) in rest  -  AK reference
)
 
# asymptomatic_proportions <- list(    
#   
#   # proportion of exposed that become asymptomatic (or go unreported) in each population
#   # set all = 0 to get rid of asymptomatic compartment
#   
#   C_asymp = 1-64/100,               # proportion of exposed that become asymptomatic (or go unreported) in care homes -  https://www.gov.scot/publications/foi-202000067854/        https://www.careinspectorate.com/images/COVID-19_-_Letter_from_Cabinet_Secretary_for_Health_and_Sport_-_Social_care_guidance_-_13_March_2020.pdf https://web.archive.org/web/20200313234513/https://hpspubsrepo.blob.core.windows.net/hps-website/nss/2980/documents/1_COVID-19%20Guidance-for-Social-or-community-care-and-residentail-settings.pdf   
#   S_asymp = 1-0.077,                # proportion of exposed that become asymptomatic (or go unreported) in staff/shielders
#   R_asymp = 1-0.077                 # proportion of exposed that become asymptomatic (or go unreported) in rest  -  AK reference
# )

num_homes_infected <- 4 # number of homes that have an outbreak intitally
num_S_infected <- 0 # number of staff compartments that have an outbreak intitally


Infect_init <- list(     # initial infected parameters
  
  # The initial total infected (reported + unreported, or symptomatic + asymptomatic)
  # in each subpopulation (e.g. R_total_infected) is split evenly between the 
  # subpopulations in that population that have an initial outbreak (e.g. R_subs_inf).

  C_subs_inf = C_initial_dist(seg_pars, num_homes_infected),    # C subpopulations that all initial C cases (C_total_infected) is split evenly between. Length must be less than or equal to seg_pars$C_pops. If no subpopulations have outbreak set as 0.
  S_subs_inf = C_initial_dist(seg_pars, num_S_infected),    # S subpopulations that all initial S cases (S_total_infected) is split evenly between. Length must be less than or equal to seg_pars$S_pops. If no subpopulations have outbreak set as 0.
  R_subs_inf = c(1),     # R subpopulations that all initial R cases (R_total_infected) is split evenly between. Length must be less than or equal to seg_pars$R_pops. If no subpopulations have outbreak set as 0.
  
  C_total_infected = num_homes_infected/(1-asymptomatic_proportions$C_asymp), # i.e. assuming 1 symptomatic case in "num_homes_infected" homes   # total initial infected across all care home sub-populations
  S_total_infected = num_S_infected/(1-asymptomatic_proportions$S_asymp),                     # total initial infected across all staff sub-populations
  R_total_infected = 110      # total initial infected across all rest sub-populations
  
  # Note: this requires (initial total infected + initital total exposed)/(number subpopulations with initial outbreak) = total infected or exposed in subpopulation initially < = capacity of that subpopulation
  #       for every subpopulation of C, S and R.
  
)


# Total infected to start is 
total_infect_init <- Infect_init$C_total_infected + Infect_init$S_total_infected + Infect_init$R_total_infected


Travel_pars <- list(    
                        # Parameters used to create the travel matrix, T. 
                        # This is the matrix whose [i,j] element is the
                        # the proportion of subpopulation i who travel to j
                        
                        # To create this matrix these get passed into create_TravelMatrix
  
                        # note that these proportions are "averaged out over a whole day"
  
  delta = 0.5,            # Proportion of staff who are at CH's instead of mixing with Rest. Controls how many staff are at each CH on average.
  # epsilon =           # Proportion of a CH's staff who work at other homes
  # gamma =             # Proportion of Rest that visit each CH. Total proportion who visit is seg_pars$C_pops*gamma. General form of gamma = x(sub_pop_sizes$C_internal/sub_pop_sizes$R_internal)(y), where x is the number of visitors per day per resident, and y is the proportion of a day a visitor spends at a care home
  
  # epsilon and gamma are time dependent (we model each of their trajectories using the logi function)
  # .....
  epsilon_pars =  list(start = -82,
                       end = 18+4,
                       rate = 0.5, 
                       low = 0.1, 
                       high = 0.1),
  
  gamma_pars =  list(start = -82,
                     end = 10,
                     rate = 3,
                     low = 0,
                     high = (sub_pop_sizes$C_internal/sub_pop_sizes$R_internal)*(2/24))
  

)


# 2a) Generate parameter lists/matrices to pass to functions ---------------------------

# generate initial SEIRD values per sub-population as a vector
SEIRD_vectors <- create_SEIRDVector(seg_pars, sub_pop_sizes, Infect_init, asymptomatic_proportions)

## FIRST WAVE matrices ######################################## 

# Create matrix of baseRLow values
baseRLow_matrix <- create_RMatrix(baseRLow_pars,seg_pars,between_pars)

# create matrix of baseRHigh Values
baseRHigh_matrix <- create_RMatrix(baseRHigh_pars, seg_pars,between_pars)

# create matrix of baseStart values
baseStart_matrix <- create_WaveMatrix(baseStart_pars, seg_pars)

# create matrix of baseEnd values
baseEnd_matrix <- create_WaveMatrix(baseEnd_pars, seg_pars)

# create matrix of baseRate values
baseRate_matrix <- create_WaveMatrix(baseRate_pars, seg_pars)



## SECOND WAVE matrices ######################################## 

# Create matrix of pulseRLow values
pulseRLow_matrix <- create_RMatrix(pulseRLow_pars,seg_pars,between_pars)

# create matrix of pulseRHigh Values
pulseRHigh_matrix <- create_RMatrix(pulseRHigh_pars, seg_pars,between_pars)

# create matrix of baseStart values
pulseStart_matrix <- create_WaveMatrix(pulseStart_pars, seg_pars)

# create matrix of baseEnd values
pulseEnd_matrix <- create_WaveMatrix(pulseEnd_pars, seg_pars)

# create matrix of pulseRate values
pulseRate_matrix <- create_WaveMatrix(pulseRate_pars, seg_pars)




## OTHER parameter lists ######################################## 

parms_passed <- list( # parameters to be passed to the SEIRD_multi function
                      
                      # infectious period appears to be longer than incubation period
                      #  https://theconversation.com/how-long-are-you-infectious-when-you-have-coronavirus-135295
  
  latency = 5.8,  # /incubation period
  Tau = 7,      # infectious period. 
  asymptomatic_ps = create_paramVector(seg_pars, asymptomatic_proportions),      # proportions of exposed that become asymptomatic in each subpopulation
  probDeath = create_paramVector(seg_pars, probDeath),     # create vector of probability of death values, length = all_pops
  seg_pars = seg_pars,     # number of sub-populations per segment
  all_pops = all_pops,
  
  # parameters decscribing travel between subpopulations - used to create travel matrix describing proportion of who travel from subpopulation i to j (has elements t_ij)
  Travel_pars = Travel_pars,  
  
  # matrices to generate the pulse(s) of transmission specific to sub-population interactions
  baseRHigh_matrix = baseRHigh_matrix,
  baseRLow_matrix = baseRLow_matrix,
  baseStart_matrix = baseStart_matrix,
  baseEnd_matrix = baseEnd_matrix,
  baseRate_matrix = baseRate_matrix,
  pulseRHigh_matrix = pulseRHigh_matrix,
  pulseRLow_matrix = pulseRLow_matrix,
  pulseStart_matrix = pulseStart_matrix,
  pulseEnd_matrix = pulseEnd_matrix,
  pulseRate_matrix = pulseRate_matrix,
  allowPulse = F
)



# 3) Data Generation ----------------------------------

# comment in below to get SEIRD curves for the parameters fixed above

# # create simulation data
# sims_multi_logi <- ode(y = SEIRD_vectors, times = simTime, func = seird_multi_logi, parms = parms_passed)
# 
# # array holds SEIRD values, e.g., S_c = Carehome_values[,,1], E_c = Carehome_values[,,2], etc.
# Carehome_values <- careHome_Aggregation(time = simTime, data = sims_multi_logi, seg_parms = seg_pars)
# 
# # array holds SEIRD values, e.g., S_s = Shielder_values[,,1], E_s = Shielder_values[,,2], etc.
# Shielder_values <- Shielder_Aggregation(time = simTime, data = sims_multi_logi, seg_parms = seg_pars)
# 
# # array holds SEIRD values, e.g., S_R = Shielder_values[,,1], E_R = Shielder_values[,,2], etc.
# Rest_values <- Rest_Aggregation(time = simTime, data = sims_multi_logi, seg_parms = seg_pars)
# 
# 
# # CH population deaths
# C_total<-rowSums(Carehome_values[,,6])[length(rowSums(Carehome_values[,,6]))]
# 
# # Staff population deaths
# S_total<-rowSums(Shielder_values[,,6])[length(rowSums(Shielder_values[,,6]))]
# 
# # Rest population deaths
# R_total<-Rest_values[,,6][length(Rest_values[,,6])]
# 
# data.frame("C_total" = C_total, "S_total" = S_total, "R_total" = R_total)
# 
