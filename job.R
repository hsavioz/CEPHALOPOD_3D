#' ============================= MASTER SCRIPT =================================
#' For the Blue-Cloud 2026 project workbench
#' A. Schickele 2024
#' =============================================================================

# --- 0. Start up and load functions
# All will be called in the config file later
rm(list=ls())
closeAllConnections()
setwd("/net/sea/work/hsavioz/phd/cephalopod/CEPHALOPOD_3D")
source(file = "./code/00_config.R")
run_name <- "MEMENTO_2026-06-18_000_cephalopod_setup"
MAX_CLUSTER = 20

# --- 1. List the available species
# Within the user defined selection criteria
list_bio <- list_bio_wrapper(FOLDER_NAME = run_name,
                             DATA_SOURCE = "/net/sea/work/hsavioz/phd/cephalopod/input/TARGETS/file.csv", # occurrence ; abundance ; biomass ; MAG; or path to a .csv file
                             SAMPLE_SELECT = list(MIN_SAMPLE = 50, TARGET_MIN_DEPTH = 0, TARGET_MAX_DEPTH = 200, START_YEAR = 1970, STOP_YEAR = 2014))

# ------------------------------------------------------------------------------
# --- USER INPUT: Define the list of species to consider
# To extract all species available in a .csv file
sp_list <- list_bio %>%
  dplyr::select(worms_id) %>% 
  unique() %>% pull() %>% .[!grepl("No match", .)]

# ------------------------------------------------------------------------------

# --- 2. Create the output folder, initialize parallelisation and parameters
# (1) Create an output folder containing all species-level runs, (2) Stores the 
# global parameters in an object, (3) Builds a local list of monthly raster
subfolder_list <- run_init(FOLDER_NAME = run_name,
                           SP_SELECT = sp_list,
                           WORMS_CHECK = FALSE,
                           FAST = TRUE,
                           LOAD_FROM = NULL,
                           DATA_TYPE = "continuous", # presence_only ; continuous ; proportions
                           ENV_VAR = c("climatology_M_0_0","climatology_i_0_50","climatology_t_0_50","climatology_p_0_50","climatology_n_0_50","climatology_omega_ca_SODA","climatology_A_PAR_regridded"),
                           ENV_PATH = NULL, # replace by local path to environmental predictors : https://data.d4science.net/m9WC
                           METHOD_PA = "density",
                           PER_RANDOM = 0,
                           PA_ENV_STRATA = TRUE,
                           OUTLIER = FALSE,
                           RFE = TRUE,
                           ENV_COR = 0.8,
                           NFOLD = 3,
                           FOLD_METHOD = "lon",
                           MODEL_LIST = c("GLM","MLP","BRT","GAM","SVM","RF"), # light version
                           LEVELS = 3,
                           TARGET_TRANSFORMATION = NULL,
                           ENSEMBLE = TRUE,
                           N_BOOTSTRAP = 10,
                           CUT = 0)

# --- 3. Query biological data
# Get the biological data of the species we wish to model
mcmapply(FUN = query_bio_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 4. Query environmental data
# This functions returns an updated subfolder_list object to avoid computing
# species with less than the user defined minimum occurrence number
subfolder_list <- mcmapply(FUN = query_env,
                           FOLDER_NAME = run_name,
                           SUBFOLDER_NAME = subfolder_list,
                           mc.cores = min(length(subfolder_list), MAX_CLUSTERS), mc.preschedule = FALSE) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  .[grep("Error", ., invert = TRUE)] %>% # to exclude any API error or else
  as.vector()

# --- 5. Generate pseudo-absences if necessary
mcmapply(FUN = pseudo_abs,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 6. Outliers, Environmental predictor and MESS check 
# This functions returns an updated subfolder_list with meaningful feature set
subfolder_list <- mcmapply(FUN = query_check,
                           FOLDER_NAME = run_name,
                           SUBFOLDER_NAME = subfolder_list,
                           mc.cores = min(length(subfolder_list), MAX_CLUSTERS), mc.preschedule = FALSE) %>% 
  unlist() %>% 
  na.omit(subfolder_list) %>% 
  as.vector()

# --- 7. Generate split and re sampling folds
mcmapply(FUN = folds,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 8. Hyper parameters to train
hyperparameter(FOLDER_NAME = run_name)

# --- 9. Model fit
mcmapply(FUN = model_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 10. Model evaluation
# Performance metric and variable importance
mcmapply(FUN = eval_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# ---11. Model projections
mcmapply(FUN = proj_wrapper,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 12. Output plots
# --- 12.1. Standard maps per algorithms
mcmapply(FUN = standard_maps,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 12.2. Partial Dependence Plots (PDP)
mcmapply(FUN = pdp,
         FOLDER_NAME = run_name,
         SUBFOLDER_NAME = subfolder_list,
         mc.cores = min(length(subfolder_list), MAX_CLUSTERS), USE.NAMES = FALSE, mc.preschedule = FALSE)

# --- 12.4 User synthesis
user_synthesis(FOLDER_NAME = run_name)



#===============================================================================
# ADDED SECTION FOR DOWNSTREAM TREATMENTS AND PLOTS

library(rhdf5)
library(terra)
library(feather)

#---- 0. Create output directory outside CEPHALOPOD source directory
dir_out <- paste0("/net/sea/work/hsavioz/phd/cephalopod/output/", run_name)
dir.create(dir_out)

#---- 1. Load CALL object, save structure as .txt
load(paste0("./output/", run_name, "/CALL.RData"))
writeLines(capture.output(str(CALL, max.level = Inf)), 
           paste0(dir_out, "/CALL.txt"))

#---- 2. Loop over species and models, save outputs required for downstream workflow
for (sp in CALL$SP_SELECT) {
  
  # create subdirectory for each species
  dir_sp <- paste0(dir_out, "/", sp)
  dir.create(dir_sp)
  
  # load MODEL and QUERY objects
  load(paste0("./output/", run_name, "/", sp, "/MODEL.RData"))
  writeLines(capture.output(str(MODEL, max.level = Inf)),
             paste0(dir_out, "/MODEL.txt"))
  
  load(paste0("./output/", run_name, "/", sp, "/QUERY.RData"))
  writeLines(capture.output(str(QUERY, max.level = Inf)),
             paste0(dir_out, "/QUERY.txt"))
  
  # save MESS projection as .h5
  mess_raster <- rast(QUERY$MESS)
  crs(mess_raster) <- "EPSG:4326"
  mess_array <- as.array(mess_raster)
  fn_projmess <- paste0(dir_sp, "proj_mess.h5")
  h5createFile(fn_projmess)
  h5write(mess_array, fn_projmess, "proj_mess")
  
  # loop over MODEL_LIST
  for (model in CALL$MODEL_LIST) {
    
    # create subsubdirectory for each model
    dir_md <- paste0(dir_sp, "/", model)
    dir.create(dir_md)
    
    # save target projections as .h5
    fn_projyhat <- paste0(dir_md, "proj_yhat.h5")
    h5createFile(fn_projyhat)
    h5write(MODEL$model$proj$y_hat, fn_projyhat, "proj_yhat")
    
    # save tandard deviation projections as .h5
    fn_projsd <- paste0(dir_md, "proj_sd.h5")
    h5createFile(fn_projsd)
    h5write(MODEL$model$proj$r_sd, fn_projsd, "proj_sd")
    
    # save partial dependence plots as .feather
    fn_pdp <- paste0(dir_md, "pdp.feather")
    write_feather(MODEL$model$pdp_all[[sp]], fn_pdp)
    
    # save variable importance as .feather
    fn_vip <- paste0(dir_md, "vip.feather")
    write_feather(as.data.frame(MODEL$model$vip), fn_vip)
    
    # save input and estimated targets as .feather
    # MERGE IN ONE DATAFRAME
    fn_y <- paste0(dir_md, "y.feather")
    write_feather(as.data.frame(MODEL$model$res$y_all), fn_y)
    fn_yhat <- paste0(dir_md, "y_hat.feather")
    write_feather(as.data.frame(MODEL$model$res$y_hat_all), fn_yhat)
    
  }
  
}


# copy job.R to corresponding output directory


#===============================================================================



# --- END --- 