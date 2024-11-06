rm(list=ls())

library(apollo)
library(parallel)
nCores <- detectCores()


data=read.csv('all_data_norm_8.csv')

apollo_initialise()

apollo_control = list(
  modelName       = "MNL_base",
  indivID         = "id", 
  mixing          = FALSE,
  nCores          = 1,
  outputDirectory = "output_full_MNL"
)

# ################################################################# #
#### LOAD DATA AND APPLY ANY TRANSFORMATIONS                     ####
# ################################################################# #

database = data[!(data$action==4 & data$interseccion_semaforo==1),] ### 

# ################################################################# #
#### DEFINE MODEL PARAMETERS                                     ####
# ################################################################# #


### Vector of parameters, including any that are kept fixed in estimation
apollo_beta = c(w10=0, ## Utility parameters
                w11=0,
                w12=0,
                w13=0,
                #w14=0,
                #w15=0,
                w17=0,
                w1_dlow=0,
                w1_slow=0,
                w1_dhigh=0,
                w1_shigh=0,
                #w_sem=0,
                w_ceda_pare=0,
                w_sem=0,
                w20=0,
                w21=0,
                w212=0,
                w22=0,
                w23=0,
                #w24=0,
                #w25=0,
                #w28=0,
                w27=0,
                w2_dlow=0,
                w2_slow=0,
                w2_dhigh=0,
                w2_shigh=0,
                w30=0,
                w31=0,
                w32=0,
                w33=0,
                #w34=0,
                w37=0,
                #w38=0,
                w3_dlow=0,
                w3_slow=0,
                w3_dhigh=0,
                w3_shigh=0,
                w40=0,
                w41=0,
                #w42=0,
                #w44=0,
                #w45=0,
                #w46=0,
                w47=0,
                w4_dlow=0,
                w4_dhigh=0
                
)


### Vector with names (in quotes) of parameters to be kept fixed at their starting value in apollo_beta, use apollo_beta_fixed = c() if none
apollo_fixed = c() ## fix variances


# ################################################################# #
#### GROUP AND VALIDATE INPUTS                                   ####
# ################################################################# #

apollo_inputs = apollo_validateInputs()

# ################################################################# #
#### DEFINE MODEL AND LIKELIHOOD FUNCTION                        ####
# ################################################################# #

apollo_probabilities=function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  ### Attach inputs and detach after function exit
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ### Create list of probabilities P
  P = list()
  
  # ## Likelihood of choices
  # ## List of utilities: these must use the same names as in mnl_settings, order is irrelevant
  ## Define the utilities for speed changes: three alternatives
  V = list()
  V[["speed0"]] = 0 # mantener
  V[["speed1"]] = (w10+(w1_dlow*distance_to_intersection_low+w1_dhigh*distance_to_intersection_high+w11)*distance_to_intersection+(w1_slow*speed_kmph_low+w1_shigh*speed_kmph_high+w12)*speed_kmph+(w13 + w_ceda_pare*(interseccion_ceda_el_paso|interseccion_pare) + w_sem*(interseccion_semaforo))*intersección +w17*n_lanes) # frenar
  V[["speed2"]]= (w20+((w2_dlow*distance_to_intersection_low+w2_dhigh*distance_to_intersection_high+w21)+w212*conoce)*distance_to_intersection+(w2_slow*speed_kmph_low+w2_shigh*speed_kmph_high+w22)*speed_kmph+w23*interseccion_semaforo +w27*n_lanes) # acelerar
  V[["speed3"]] = (w30+(w3_dlow*distance_to_intersection_low+w3_dhigh*distance_to_intersection_high+w31)*distance_to_intersection+(w3_slow*speed_kmph_low+w3_shigh*speed_kmph_high+w32)*speed_kmph+w33*interseccion_semaforo +w37*n_lanes) # desacelera
  V[["speed4"]]= (w40+(w4_dlow*distance_to_intersection_low+w4_dhigh*distance_to_intersection_high+w41)*distance_to_intersection+w47*n_lanes) # esperar
  
  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives  = c(speed0=1, speed1=0,speed2=2, speed3=3,speed4=4),
    avail         = list(speed0=1, speed1=1, speed2=1, speed3=1,speed4=1), # availability of alternatives (=1 means always available)
    choiceVar     = action,
    V     = V,
    componentName = "Choice"
  )
  
  ### Compute probabilities for MNL model component
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  
  ### Likelihood of the whole model
  
  ### Take product across observation for same individual
  P = apollo_panelProd(P, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# ################################################################# #
#### What does apollo_probabilities return?                      ####
# ################################################################# #

#show initial probabilities
apollo_probabilities(apollo_beta, apollo_inputs, functionality="estimate")

#show initial log-likelihood
apollo_llCalc(apollo_beta, apollo_probabilities, apollo_inputs)

# Load estimated model 
# M=apollo_loadModel('modelo_ejemplo')

  #Estimate the model
model = apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs, 
                        estimate_settings = list(maxIterations=1000, estimationRoutine="bfgs"))#

#
#in case you want to estimate the model starting from previous parameters
#model=apollo_estimate(M$estimate,M$apollo_fixed, apollo_probabilities, apollo_inputs)

# ################################################################# #
#### MODEL OUTPUTS                                               ####
# ################################################################# #

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO SCREEN)                               ----
# ----------------------------------------------------------------- #
apollo_modelOutput(model)
# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #
apollo_saveOutput(model)
