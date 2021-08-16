# Code Description ############################################################
# Date: 09_08_2021
# Written by: Paul McMenemy, Ewan McTaggart, Matthew Baister
# Script to model COVID dynamics using meta-populations of 
# VULNERABLE population segment - CARE HOME RESIDENTS (C)
# SHIELDING population segment - CARE HOME STAFF (S) and
# REST OF population segment - REST OF POPULATION (R)
# 
# This script contains the functions used in other scripts
#
# 1) Preamble - install required packages
# 2) Functions to create parameter matrices
# 3) Functions to create a subpopulation parameter vector e.g. create vector of mortality rates 
# 4) Single-pop diff. model - code to construct the basic metapopulation model
# 5) Multi-pop model - SEIRD model with M sub-populations AND incorporates 
#    smooth time-dependent reproduction/transmission rates
# 6) Population data tidy functions - extracts the C, S and R data into specific 3D arrays
# 7) AK's logi function
# 8) Create vectors of initial SEIRD values across all sub-populations.

# 1) Preamble -----------------------------------------------------


# if (!"deSolve" %in% installed.packages()){
#   install.packages("deSolve")}
library("deSolve")


start_timing <- Sys.time() # for timing how long code takes to run



# 2) Create Parameter matrix --------------------------------------

# general Parameter_matrix is in form:

#         C           S          R
#    - - - - - - - - - - - - - - - - -
#   |           |          |          |
# C |           |          |          |
#   |           |          |          |
#    - - - - - - - - - - - - - - - - -
#   |           |          |          |
# S |           |          |          |
#   |           |          |          |
#    - - - - - - - - - - - - - - - - -
#   |           |          |          |
# R |           |          |          |
#   |           |          |          |
#    - - - - - - - - - - - - - - - - -
# 

# function that generates matrices of parameters that determine:
# infection rates across time specific to the interactions between C,S,R

# i.e. this function creates matrix of parameters (high, low) to pass into 
# logi function that are affected by the "between parameters"

create_RMatrix <- function(R_pars, seg_parms, constant_pars){        # this version has has beta_rc and beta_cr as diagnal matrices. Therefore each C subpopulation is connected to one other R population. C_i is connected to (and will have visits only from) R_i.
  
  #test declarations
  #R_pars = baseRLow_pars

  
  with(as.list(c(R_pars,seg_parms,constant_pars)),{
    
    all_pops = Reduce('+',seg_parms) # total number of sub-populations
    
    # creates empty matrix to store transmission rates between segments and subcompartments of pop segments
    # create to allow flexibility in the number of subcompartments within seg_pars
    R_matrix <- matrix(NA, nrow = all_pops, ncol = all_pops)
    
    # care home to care home submatrix beta values
    
    cc_submatrix = diag(x = beta_c, nrow = C_pops, ncol = C_pops)
    cc_submatrix[lower.tri(cc_submatrix)] = beta_cc
    cc_submatrix[upper.tri(cc_submatrix)] = beta_cc
    
    R_matrix[1:C_pops, 1:C_pops] <- cc_submatrix
    
    # care home to shielders submatrix beta values
    
    cs_submatrix <- diag(x = beta_cs, nrow = C_pops, ncol = S_pops)
    cs_submatrix[lower.tri(cs_submatrix)] = beta_cs
    cs_submatrix[upper.tri(cs_submatrix)] = beta_cs
    
    R_matrix[1:C_pops,(C_pops+1):(C_pops+S_pops)] <- cs_submatrix
    
    # care home to rest of population submatrix beta values
    
    cr_submatrix = diag(x = beta_c, nrow = C_pops, ncol = R_pops)
    cr_submatrix[lower.tri(cr_submatrix)] = beta_c
    cr_submatrix[upper.tri(cr_submatrix)] = beta_c
    
    R_matrix[1:C_pops,(C_pops+S_pops+1):all_pops] <- cr_submatrix
    
    # shielder to care home submatrix beta values
    
    sc_submatrix <-  diag(x = beta_sc, nrow = S_pops, ncol = C_pops)
    sc_submatrix[lower.tri(sc_submatrix)] = beta_cs
    sc_submatrix[upper.tri(sc_submatrix)] = beta_cs
    
    R_matrix[(C_pops+1):(C_pops+S_pops),1:C_pops] <- sc_submatrix
    
    # shielder to shielder submatrix beta values
    
    ss_submatrix = diag(x = beta_s, nrow = S_pops, ncol = S_pops)
    ss_submatrix[lower.tri(ss_submatrix)] = beta_ss
    ss_submatrix[upper.tri(ss_submatrix)] = beta_ss
    
    R_matrix[(C_pops+1):(C_pops+S_pops),(C_pops+1):(C_pops+S_pops)] <- ss_submatrix
    
    # shielder to rest of population submatrix beta values
    
    R_matrix[(C_pops+1):(C_pops+S_pops),
             (C_pops+S_pops+1):(all_pops)] <- beta_rr*(1-alpha)
    
    # rest of population to care home submatrix beta values
    
    rc_submatrix = diag(x = beta_c, nrow = R_pops, ncol = C_pops)
    rc_submatrix[lower.tri(rc_submatrix)] = beta_c
    rc_submatrix[upper.tri(rc_submatrix)] = beta_c
    
    R_matrix[(C_pops+S_pops+1):(all_pops),1:C_pops] <- rc_submatrix
    
    # rest of population to shielder submatrix beta values
    
    R_matrix[(C_pops+S_pops+1):(all_pops),(C_pops+1):(C_pops+S_pops)] <- beta_rr*(1-alpha)
    
    # rest of population to rest of population beta values
    
    R_matrix[(C_pops+S_pops+1):(all_pops),
             (C_pops+S_pops+1):all_pops] <- beta_rr
    
    return(R_matrix)
  })
}


# function that generates matrices of parameters that determine:
# infection rates across time specific to the interactions between C,S,R

# i.e. this function creates matrix of parameters (start, end, rate) to pass into 
# logi function that are unaffected by the "between parameters"

create_WaveMatrix <- function(wave_pars, seg_parms){
  
  with(as.list(c(wave_pars, seg_parms)),{
    
    all_pops = Reduce('+',seg_parms) # total number of sub-populations
            
    # creates empty matrix to store infection wave parameters
    # create to allow flexibility in the number of subcompartments within seg_pars
    wave_matrix <- matrix(NA, nrow = all_pops, ncol = all_pops)
    
    # care home to care home submatrix wave values
    
    cc_submatrix = diag(x = c, nrow = C_pops, ncol = C_pops)
    cc_submatrix[lower.tri(cc_submatrix)] = cc
    cc_submatrix[upper.tri(cc_submatrix)] = cc
    
    wave_matrix[1:C_pops, 1:C_pops] <- cc_submatrix
    
    # care home to shielders submatrix wave values
    
    wave_matrix[1:C_pops,(C_pops+1):(C_pops+S_pops)] <- cs
    
    # care home to rest of population submatrix wave values
    
    wave_matrix[1:C_pops,(C_pops+S_pops+1):all_pops] <- rr
    
    # shielder to care home submatrix wave values
    
    wave_matrix[(C_pops+1):(C_pops+S_pops),1:C_pops] <- sc
    
    # shielder to shielder submatrix wave values
    
    ss_submatrix = diag(x = s, nrow = S_pops, ncol = S_pops)
    ss_submatrix[lower.tri(ss_submatrix)] = ss
    ss_submatrix[upper.tri(ss_submatrix)] = ss
    
    wave_matrix[(C_pops+1):(C_pops+S_pops),(C_pops+1):(C_pops+S_pops)] <- ss_submatrix
    
    # shielder to rest of population submatrix wave values
    
    wave_matrix[(C_pops+1):(C_pops+S_pops),
             (C_pops+S_pops+1):(all_pops)] <- rr
    
    # rest of population to care home submatrix wave values
    
    wave_matrix[(C_pops+S_pops+1):(all_pops),1:C_pops] <- rr
    
    # rest of population to shielder submatrix wave values
    
    wave_matrix[(C_pops+S_pops+1):(all_pops),(C_pops+1):(C_pops+S_pops)] <- rr
    
    # rest of population to rest of population wave values
    
    wave_matrix[(C_pops+S_pops+1):(all_pops),
             (C_pops+S_pops+1):all_pops] <- rr
    
    return(wave_matrix)
  })
}


# this function creates the travel matrix
# this is the matrix that describes travel between subpopulations
# entry [i,j] is the proportion of i that travel to j 
# i.e. [i,j] = [from, to]   
# note that the rows sum to 1

create_TravelMatrix <- function(Travel_pars, seg_pars, time){        
  
  with(as.list(c(Travel_pars, seg_pars, time)),{

    delta <- delta                                                 # Proportion of staff who are at CH's instead of mixing with Rest. Controls how many staff are at each CH on average.
    epsilon <- do.call(logi, as.list(c(t = time, epsilon_pars)))   # Proportion of a CH's staff who are at other homes
    gamma <- do.call(logi, as.list(c(t = time, gamma_pars)))       # Proportion of Rest that visit each CH
    
    all_pops = Reduce('+',seg_pars)         # total number of sub-populations
    
    # create empty matrix 
    T_matrix <- matrix(NA, nrow = all_pops, ncol = all_pops)
    
    # from care home to care home submatrix 
    
    cc_submatrix = diag(x = 1, nrow = C_pops, ncol = C_pops)
    
    T_matrix[1:C_pops, 1:C_pops] <- cc_submatrix
    
    # from care home to shielder submatrix 
    
    T_matrix[1:C_pops,(C_pops+1):(C_pops+S_pops)] <- 0
    
    # from care home to rest submatrix 
    
    T_matrix[1:C_pops,(C_pops+S_pops+1):all_pops] <- 0
    
    # from shielder to shielder submatrix 
    
    T_matrix[(C_pops+1):(C_pops+S_pops),(C_pops+1):(C_pops+S_pops)] <- 0
    
    # from shielder to care home submatrix 
    
    
    # option i) 
    # staff shared homogenously between two other homes in a circle 
    # communities of size 3 that overlap
    
    sc_submatrix <-  diag(x = (1-epsilon)*delta, nrow = S_pops, ncol = C_pops)   # diagonal elements
    sc_submatrix[abs(row(sc_submatrix) - col(sc_submatrix)) == 1]<-epsilon*delta/2  # elements next to diagonal
    sc_submatrix[1,C_pops]<-epsilon*delta/2  # top right corner
    sc_submatrix[S_pops,1]<-epsilon*delta/2  # bottom left corner

    T_matrix[(C_pops+1):(C_pops+S_pops),1:C_pops] <- sc_submatrix

    # alternative option ii) - not explored
    # staff shared homogenously between all homes
    
    # sc_submatrix <-  diag(x = (1-epsilon)*delta, nrow = S_pops, ncol = C_pops)
    # sc_submatrix[lower.tri(sc_submatrix)] = (epsilon/(C_pops-1))*delta
    # sc_submatrix[upper.tri(sc_submatrix)] = (epsilon/(C_pops-1))*delta
    # 
    # T_matrix[(C_pops+1):(C_pops+S_pops),1:C_pops] <- sc_submatrix
    
    
    # from shielder to rest submatrix 
    
    T_matrix[(C_pops+1):(C_pops+S_pops),
             (C_pops+S_pops+1):(all_pops)] <- 1-delta
    
    # from rest of population to care home submatrix 
    
    T_matrix[(C_pops+S_pops+1):(all_pops),1:C_pops] <- gamma
    
    # from rest of population to shielder submatrix 
    
    T_matrix[(C_pops+S_pops+1):(all_pops),(C_pops+1):(C_pops+S_pops)] <- 0
    
    # from rest of population to rest of population submatrix
    
    T_matrix[(C_pops+S_pops+1):(all_pops),
             (C_pops+S_pops+1):all_pops] <-1-C_pops*gamma 
    
    return(T_matrix)
  })
}


# 3) Create a subpopulation parameter vector -------------------

# param_values is a list/vector with three elements; each element
# is a fixed parameter value for each population (C, S and R). The parameter 
# is constant across each population i.e. across all CH's. An example is a 
# vector of mortality rates = c(C_prob_death, S_prob_death, R_prob_death).

# this function takes param_values and returns a vector of length = all_pops.
# The ith element of this vector is the fixed parameter value for the ith 
# subpopulation. For example, applying this function to the mortaility rates 
# vector mentioned would return a vector where the 1st element corresponds
# to the death rate for CH1.

# Note: the ith element of the returned vector matches up with the ith 
#       element of the state vector in the ode (and so with SEIRD_vectors)

create_paramVector <- function (seg_pars, param_values){
  
  C_param = unlist(param_values)[1]; S_param = unlist(param_values)[2]; R_param = unlist(param_values)[3]
  
  all_pops = Reduce('+', seg_pars)
  
  sub_pop_vector <- rep(NA, all_pops)
  
  sub_pop_vector[1:seg_pars$C_pops] = C_param
  sub_pop_vector[(seg_pars$C_pops+1):(seg_pars$C_pops+seg_pars$S_pops)] = S_param
  sub_pop_vector[(seg_pars$C_pops+seg_pars$S_pops+1): all_pops] = R_param
  
  # sub_pop_vector <- c(rep(C_param, seg_pars$C_pops),
  #                     rep(S_param, seg_pars$S_pops),
  #                     rep(R_param, seg_pars$R_pops)
  #                     )
  
  return(as.vector(sub_pop_vector))
}

# 4) Single-pop diff. model ------------------------------------

# needs a wrapper loop which redefines the parameters being passed to the function
# to take into account the split of the population into M (= all_pops) sub-populations
# seird_single <- function(t, y, parms){
#   
#   # length(y) must equal number of compartments in diff model
#   # Ensure to match y to order of SEIRD
#   S = y[1];  E = y[2];  I = y[3];  R = y[4];  D = y[5]
#   
#   with(as.list(parms),{
#     
#     Npop = S+E+I+R+D # total population of the sub-population
#     
#     # expressions for each of the derivatives in our diff. model
#     dS = -beta_value*I*S/Npop
#     dE = +beta_value*I*S/Npop - E/latency
#     dI = +E/latency - I/Tau
#     dR = +(1 - probDeath)*I/Tau
#     dD = +probDeath*I/Tau
#     return(list(c(dS,dE,dI,dR,dD)))
#   })
# }

# 5) Multi-pop model - works nicely now! -----------------------

seird_multi_logi <- function(time_t, y, parms){

  #test parameters - remove once function constructed/tested
  #y = SEIRD_vectors
  #parms = parms_passed
  
  # define vectors of initial SEIRD populations
  S = y[1:parms$all_pops];  E = y[(parms$all_pops+1):(2*parms$all_pops)];  IS = y[(2*parms$all_pops+1):(3*parms$all_pops)]; 
  IA = y[(3*parms$all_pops+1):(4*parms$all_pops)]; R = y[(4*parms$all_pops+1):(5*parms$all_pops)];
  D = y[(5*parms$all_pops+1):(6*parms$all_pops)]

  # create blank vectors to store derivative values
  Npop = rep(NA, parms$all_pops)
  dS = rep(NA, parms$all_pops); dE = rep(NA, parms$all_pops)
  dIS = rep(NA, parms$all_pops); dIA = rep(NA, parms$all_pops); dR = rep(NA, parms$all_pops); dD = rep(NA, parms$all_pops)
  logi_matrix <- matrix(NA, nrow = parms$all_pops, ncol = parms$all_pops)
  external_inf <- rep(NA, parms$all_pops)
  # for testing with AK
  # baseRLow_matrix1 <- diag(0.8,all_pops,all_pops)
  # baseRHigh_matrix1 <- diag(2.7,all_pops,all_pops)
  # pulseRLow_matrix1 <- diag(0.0,all_pops,all_pops)
  # pulseRHigh_matrix1 <- diag(0.8,all_pops,all_pops)
  # baseRLow_matrix1[1,2] <- 0.1
  # baseRLow_matrix1[2,1] <- 0.1
  
  with(as.list(parms),{
    
    # loops to calculate 
    # - logi matrix at time t
    # - external infection into subpopulation i at time t
    
    for (i in 1:all_pops){
      
      for (j in 1:all_pops){
        
        # entry [i,j] = transmission rate from i to j
        logi_matrix[i,j] = logi(time_t,baseStart_matrix[i,j],baseEnd_matrix[i,j],
                                baseRate_matrix[i,j],baseRLow_matrix[i,j],baseRHigh_matrix[i,j]) +
                            if(parms$allowPulse == T){logi(time_t,pulseStart_matrix[i,j],pulseEnd_matrix[i,j],
                                                      pulseRate_matrix[i,j],pulseRLow_matrix[i,j],pulseRHigh_matrix[i,j])} 
                            else {0}
      }
    }   
    
    
    Npop = S + E + IS + IA + R #+ D   # Comment out "+ D" as the dead do not come into contact with anyone
    
    # Note:
    # if Npop[i] = S[i] + E[i] + IS[i] + IA[i] + R[i] + D[i], 
    # then Npop[i] = internal size of subpopulation i
    
    
    # Equations inspired by paper
    # "Metapopulation Network Models for Understanding, Predicting, and Managing the Coronavirus Disease COVID-19"
    # https://doi.org/10.3389/fphy.2020.00261
      
      # Note: the following are scalars
      #
      #       beta_ji = logi_matrix[j,i]/parms$Tau
      #       N_j = Npop[j] 
      #       \hat{N}_k = Npop%*%Travel_Matrix[,k]
      
      # The travel matrix with elements t_ij is
      # Travel_Matrix = create_TravelMatrix(Travel_pars, seg_pars, time_t) 
      
      # we treat the FOI term as the sum of dot product of two vectors...
      
      # Firstly, 
      # sweep(Travel_Matrix,2,Npop%*%Travel_Matrix,`/`)  # = a matrix where [i,k]th element is t_ij[i,k]/(\hat{N}_k)
      # where we set the NaN's to zero for the later calculation
    
      Travel_Matrix <- create_TravelMatrix(Travel_pars, seg_pars, time_t)
    
      t_ik_w_columns_divided_by_hatNk<-matrix(mapply(function(x){ifelse(is.nan(x),0,x)},
                                                     sweep(Travel_Matrix,2,Npop%*%Travel_Matrix,`/`)),
                                              nrow = all_pops, ncol = all_pops)
      # so the kth element of t_ik_w_columns_divided_by_hatNk multiplied by S[i] 
      # is the kth term in the first sum in the equations    (*)
      

      # Secondly,
      # logi_matrix[,i]/parms$Tau # = column vector of beta_ji for a fixed i
      # (IS+IA)*t(logi_matrix[,i]/parms$Tau) # = a row vector where each element is (I_j+A_j)*beta_ji
      # Thus,
      # ((IS+IA)*t(logi_matrix[,i]/parms$Tau))%*%Travel_Matrix # row vector where kth element is the second sum in the equations for that fixed k. Note the matrix multiplication %*%.
      
      #  so for a fixed k this is the second sum in the equation
      # ((IS+IA)*t(logi_matrix[,i]/parms$Tau))%*%Travel_Matrix[,k]    (**)
      
      # Thus, to calculate the full FOI term we sum the dot product of (*) and (**)
      # and we must do this for every ith subpopulation
      
      # The following loop to calculates SEIRD for each sub-population
      
      for (i in 1:all_pops){

      dS[i] = -sum((t_ik_w_columns_divided_by_hatNk[i,]*S[i])*(((IS+IA)*t(logi_matrix[,i]/parms$Tau))%*%Travel_Matrix)) 
      dE[i] = sum((t_ik_w_columns_divided_by_hatNk[i,]*S[i])*(((IS+IA)*t(logi_matrix[,i]/parms$Tau))%*%Travel_Matrix)) - E[i]/parms$latency 
      dIS[i] = (1-asymptomatic_ps[i])*E[i]/parms$latency - IS[i]/parms$Tau 
      dIA[i] = asymptomatic_ps[i]*E[i]/parms$latency - IA[i]/parms$Tau
      dR[i] = (1 - probDeath[i])*(IS[i]+IA[i])/parms$Tau
      dD[i] = probDeath[i]*(IS[i]+IA[i])/parms$Tau

    }
    
    # return list of 5*all_pops values
    return(list(c(dS, dE, dIS, dIA, dR, dD)))
  })
}

# Inflexible Version - need to adapt this to allow variable values of m_i \in M (subpops \in all_pops)

# seird_multi <- function(t, y, parms){
# 
#   #test parameters - remove once function constructed/tested
#   # y = SEIRD_vectors
#   # parms = parms_passed
# 
#   # define vectors of SEIRD populations each of length(all_pops)
#   S = y[1:all_pops];  E = y[(all_pops+1):(2*all_pops)];  I = y[(2*all_pops+1):(3*all_pops)]
#   R = y[(3*all_pops+1):(4*all_pops)];  D = y[(4*all_pops+1):(5*all_pops)]
# 
#   # create blank vectors to store derivative values
#   Npop = rep(NA, all_pops); dS = rep(NA, all_pops); dE = rep(NA, all_pops)
#   dI = rep(NA, all_pops); dR = rep(NA, all_pops); dD = rep(NA, all_pops)
# 
#   with(as.list(parms),{
# 
#     for (i in 1:all_pops){
# 
#       Npop[i] = S[i]+E[i]+I[i]+R[i]+D[i] #(is fixed but may vary so worth checking)
#       # expressions for each of the derivatives in our diff. model
#       # extended to M (= all_pops) sub-populations
#       dS[i] = -Reduce('+',beta_matrix[i,]*I)*S[i]/Npop[i]
#       dE[i] = +Reduce('+',beta_matrix[i,]*I)*S[i]/Npop[i] - E[i]/latency
#       dI[i] = +E[i]/latency - I[i]/Tau
#       dR[i] = +(1 - probDeath[i])*I[i]/Tau
#       dD[i] = +probDeath[i]*I[i]/Tau
#     }
# 
#     # return list of 5*all_pops values
#     return(list(c(dS,dE,dI,dR,dD)))
#   })
# }

# 6) Population Data Tidy Functions - could be adapted into one function -------------------

careHome_Aggregation <-  function(time, data, seg_parms){
  
  with(as.list(seg_parms),{
    # empty array to store Care home sub-pop values
    # dims = (time, no. of Care home sub-pops, SEIRD-values)
    total_C = array(NA, dim = c(length(simTime), C_pops, 6))
    
    all_pops = Reduce('+',seg_parms)
    
    for (i in 1:6){
      total_C[,1:C_pops,i] = data[,(2+(i-1)*all_pops):(1+(i-1)*all_pops+C_pops)]
    }
    return(total_C)
  })
}

Shielder_Aggregation <-  function(time, data, seg_parms){
  
  with(as.list(seg_parms),{
    # empty array to store Shielder/Staff sub-pop values
    # dims = (time, no. of Shielder sub-pops, SEIRD-values)
    total_S = array(NA, dim = c(length(simTime), S_pops, 6 ))
    
    all_pops = Reduce('+',seg_parms)
    
    for (i in 1:6){
      total_S[,1:S_pops,i] = data[,(2+C_pops+(i-1)*all_pops):(1+C_pops+(i-1)*all_pops+S_pops)]
    }
    return(total_S)
  })
}

Rest_Aggregation <-  function(time, data, seg_parms){
  

  with(as.list(seg_parms),{
    # empty array to store Rest sub-pop values
    # dims = (time, no. of Rest sub-pops, SEIRD-values)
    total_R = array(NA, dim = c(length(simTime), R_pops, 6))
    
    all_pops = Reduce('+',seg_parms)
    
    for (i in 1:6){
      total_R[,1:R_pops,i] = data[,(2+C_pops+S_pops+(i-1)*all_pops):(1+C_pops+S_pops+(i-1)*all_pops+R_pops)]
    }
    return(total_R)
  })
}

# 7) logi function ----------------------------

# We need to call logi twice, for the first wave and for the second wave; it might be better to put these two calls together
# into one function which will have five (possibly 6) parameters:
#   
# 1. high value of R in the first wave
# 2. length of the first wave
# 3. low value after the first wave
# 4. start of the second wave
# 5. value in the second wave
# 6. possibly duration/end of the second wave

logi <- function(t,start,end,rate,low,high){
  
  l1 = 1/(1+exp(rate*(t-end)))
  l2 = 1/(1+exp(-rate*(t-start)))
  
  return((high-low)*l1*l2+low)
  
}


# 8) Functions to create SEIRD vector -------------


# Create SEIRD vectors for each sub-populations initial values by seeding equally sized 
# clusters of infection within each population.
#
# The total infected in each subpopulation (e.g. C_total_infected) is split
# evenly between the subpopulations in that population that have 
# an initial outbreak (e.g. C_subs_inf).
#
# Therefore, if C_subs_inf = c(1,2) then the total C infected to 
# start (C_total_infected) is seeded evenly between C_1 and C_2. If C_subs_inf = c(1)
# then the outbreak starts as a single cluster of size C_total_infected
# in C_1.



create_SEIRDVector <- function(seg_pars, sub_pop_sizes, Infect_init, asymptomatic_proportions){
  
  # current assumption is in any subpopulation the initial number of exposed
  # equals the total number infected in that subpopulation
  
  with(as.list(c(seg_pars, sub_pop_sizes, Infect_init, asymptomatic_proportions)),{
    
    all_pops = Reduce('+', seg_pars)
    
    C_infected<-rep(0, C_pops)     # initial number infected in each C subpopulation
    S_infected<-rep(0, S_pops)     # initial number infected in each S subpopulation         # infected = total infected i.e. symptomatic + asymptomatic or reported + unreported
    R_infected<-rep(0, R_pops)     # initial number infected in each R subpopulation
    
    
    # if(C_total_infected == 0/(1-C_asymp)){ 
    #   C_subs_inf<-c(1) 
    # }
    # 
    # if(C_total_infected == 1/(1-C_asymp)){ 
    #   C_subs_inf<-c(1)
    # }
    # 
    # if(C_total_infected == 2/(1-C_asymp)){ 
    #   C_subs_inf<-c(1,round(C_pops/2,0))
    # }
    
    
    C_infected[C_subs_inf]<-C_total_infected/(length(C_subs_inf))     # the total infected across each population (e.g. C_total_infected across CH's)
    S_infected[S_subs_inf]<-S_total_infected/(length(S_subs_inf))     # is split evenly between the subpopulations that have intial outbreak (e.g. C_subs_inf)
    R_infected[R_subs_inf]<-R_total_infected/(length(R_subs_inf))
    
    
    S_init<-c(c(rep(C_internal, C_pops)-1*C_infected),
              c(rep(S_internal, S_pops)-2*S_infected),
              c(rep(R_internal, R_pops)-2*R_infected)
    )
    E_init = c(1*C_infected,         
               1*S_infected,         
               1*R_infected
    )
    IS_init = c(0*(1-C_asymp)*C_infected,
                (1-S_asymp)*S_infected,
                (1-R_asymp)*R_infected
    )
    IA_init = c(0*C_asymp*C_infected,
                S_asymp*S_infected,
                R_asymp*R_infected
    )
    R_init = rep(0, all_pops)
    D_init = rep(0, all_pops)
    
    return(c(S_init,E_init,IS_init,IA_init,R_init,D_init))
  })
}


# This function below is used to create the vector for C_subs_inf

C_initial_dist<- function(seg_pars, number_homes_with_outbreak){   # number_homes_with_outbreak is the number of homes with an intitial outbreak
  
  # Assume all the care homes were arranged in a circle.
  # This function returns the indices of care homes
  # which have an outbreak and ensures that the homes in the vector
  # are equally spaced apart 
  
  if(number_homes_with_outbreak==0){   # if no homes have an outbreak return 0
    indices_homes_with_outbreak <- 0
  }
  
  else{
    indices_homes_with_outbreak <- seq(1, length.out = number_homes_with_outbreak, by = floor(seg_pars$C_pops/number_homes_with_outbreak))
  }
  
  
  
  return(indices_homes_with_outbreak)
}




# END -------------------