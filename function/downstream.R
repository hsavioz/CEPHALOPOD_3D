#' =============================================================================
#' @name downstream
#' @description Save CEPHALOPOD outputs for downstream treatment.
#'
#' @param run_name Name of the CEPHALOPOD run.
#' =============================================================================


save_outputs <- function(run_name) {
  
  # --- 0. Configuration
  
  dir_ceph <- file.path(project_wd, "output", run_name)
  
  dir_out <- file.path("/net/sea/work/hsavioz/phd/cephalopod/output", run_name)
  dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)
  
  # --- 1. Load CALL
  
  call_file <- file.path(dir_ceph, "CALL.RData")
  load(call_file)
  
  # This downstream workflow is currently intended for continuous data.
  if (!identical(CALL$DATA_TYPE, "continuous")) {
    stop("save_outputs() currently supports only DATA_TYPE = 'continuous'.")
  }
  
  
  # --- 2. Reconstruct the horizontal output grid
  #
  # In the future 3-D implementation, this section can be replaced by an
  # explicit output-grid definition without changing the NetCDF writer.
  
  CALL$ENV_DATA <- lapply(CALL$ENV_DATA, function(x) terra::rast(x))
  r0 <- CALL$ENV_DATA[[1]][[1]]
  
  ext_r <- terra::ext(r0)
  res_r <- terra::res(r0)
  
  nlon <- terra::ncol(r0)
  nlat <- terra::nrow(r0)
  
  lon <- seq(
    ext_r$xmin + res_r[1] / 2,
    ext_r$xmax - res_r[1] / 2,
    length.out = nlon
  )
  
  lat <- seq(
    ext_r$ymin + res_r[2] / 2,
    ext_r$ymax - res_r[2] / 2,
    length.out = nlat
  )
  
  
  # --- 3. Coordinates common to the present CEPHALOPOD workflow
  
  month <- seq_len(length(CALL$ENV_DATA))
  bootstrap <- seq_len(CALL$N_BOOTSTRAP)
  
  
  # --- 4. Determine the 2-D depth interval
  
  if (!is.null(CALL$SAMPLE_SELECT$TARGET_MIN_DEPTH) && !is.null(CALL$SAMPLE_SELECT$TARGET_MAX_DEPTH)) {
    layer_depth_range <- c(CALL$SAMPLE_SELECT$TARGET_MIN_DEPTH, CALL$SAMPLE_SELECT$TARGET_MAX_DEPTH)
  } 
  else {
    stop("CALL$SAMPLE_SELECT$TARGET_MIN_DEPTH and TARGET_MAX_DEPTH are required for the 2-D downstream output.")
  }
  

  # --- 5. Loop over species
  
  for (sp in CALL$SP_SELECT) {
    
    message("============================================================")
    message("DOWNSTREAM: ", sp)
    message("============================================================")
    
    
    dir_sp_ceph <- file.path(dir_ceph, sp)
    
    dir_sp <- file.path(dir_out, sp)
    dir.create(dir_sp, recursive=TRUE, showWarnings=FALSE)
    
    
    # --- 5.1 Load MODEL
    
    model_file <- file.path(dir_sp_ceph, "MODEL.RData")
    load(model_file)
    
    # --- 5.2 Load QUERY
    
    query_file <- file.path(dir_sp_ceph, "QUERY.RData")
    load(query_file)
    
    
    # --- 6. Models to save
    
    models_to_save <- CALL$MODEL_LIST
    
    # MODEL_LIST can have been reduced during the QC steps. Only models that
    # actually survived QC should therefore be exported.
    models_to_save <- models_to_save[
      models_to_save %in% names(MODEL)
    ]
    
    
    # --- 7. Save individual models
    
    for (model in models_to_save) {
      
      message("--- DOWNSTREAM: saving ", model)
      
      dir_md <- file.path(dir_sp, model)
      dir.create(dir_md, recursive=TRUE, showWarnings=FALSE)
      
      
      # --- 7.1 Projection
        
      yhat_proj <- MODEL[[model]]$proj$y_hat
      
      fn_projyhat <- file.path(
        dir_md,
        "proj_yhat.nc"
      )
      
      if (isFALSE(CALL$MODEL_3D)) {
        
        # Current CEPHALOPOD format:
        # [cell, bootstrap, month]
        write_projection_netcdf(
          array_in = yhat_proj,
          fn_out = fn_projyhat,
          lon = lon,
          lat = lat,
          month = month,
          bootstrap = bootstrap,
          layer_depth_range = layer_depth_range,
          varname = "proj_yhat",
          units = units_target,
          longname = paste("CEPHALOPOD", model, "projection")
        )
        
      } else {
        
        # Future CEPHALOPOD 3-D format:
        # [cell, depth, bootstrap, month]
        #
        # The depth coordinate must be provided by the 3-D projection
        # pipeline.
        
        depth <- attr(yhat_proj, "depth")
        
        write_projection_netcdf(
          array_in = yhat_proj,
          fn_out = fn_projyhat,
          lon = lon,
          lat = lat,
          depth = depth,
          month = month,
          bootstrap = bootstrap,
          varname = "proj_yhat",
          units = units_target,
          longname = paste(
            "CEPHALOPOD",
            model,
            "3-D projection"
          )
        )
      }
      
      
      # --- 7.2 Variable importance
        
        vip <- as.data.frame(MODEL[[model]]$vip)
        
        if ("variable" %in% names(vip)) {
          vip$variable <- as.character(vip$variable)
        }
        
        fn_vip <- file.path(dir_md, "vip.feather")
        feather::write_feather(vip, fn_vip)
        
      
      # --- 7.3 Evaluation / residuals
      
      eval <- build_eval_table(MODEL=MODEL, QUERY=QUERY, model=model)
      fn_eval <- file.path(dir_md,"eval.feather")
      feather::write_feather(eval, fn_eval)
      
    }
    
    
    # --- 8. Build and save ensemble
    
    if (isTRUE(CALL$ENSEMBLE) && length(models_to_save) > 1) {
      
      message("--- DOWNSTREAM: saving ENSEMBLE")
      
      dir_ens <- file.path(dir_sp, "ENSEMBLE")
      
      dir.create(dir_ens, recursive=TRUE, showWarnings=FALSE)
      
      
      # --- 8.1 Ensemble projection
      #
      # Every retained model uses the same bootstrap split object in
      # proj_continuous(), so we calculate the ensemble as the mean of the
      # model predictions for each bootstrap/month.
      #
      # Current individual-model shape:
      #
      #   cell x bootstrap x month
      #
      # Ensemble shape:
      #
      #   cell x bootstrap x month
      #
      
      yhat_list <- lapply(models_to_save, function(model) MODEL[[model]]$proj$y_hat)
      
      same_dim <- vapply(yhat_list, function(x) identical(dim(x), dim(yhat_list[[1]])), logical(1))
      
      if (!all(same_dim)) {
        stop(
          "The retained model projections do not have identical dimensions; ",
          "cannot construct the ensemble projection."
        )
      }
      
      yhat_ensemble <- Reduce(
        `+`,
        yhat_list
      ) / length(yhat_list)
      
      fn_projyhat <- file.path(dir_ens, "proj_yhat.nc")
      
      if (isFALSE(CALL$MODEL_3D)) {
        
        write_projection_netcdf(
          array_in = yhat_ensemble,
          fn_out = fn_projyhat,
          lon = lon,
          lat = lat,
          month = month,
          bootstrap = bootstrap,
          layer_depth_range = layer_depth_range,
          varname = "proj_yhat",
          units = units_target,
          longname = "CEPHALOPOD ensemble projection"
        )
        
      } else {
        
        if (is.null(attr(yhat_ensemble, "depth"))) {
          stop(
            "The ensemble 3-D projection does not contain a depth coordinate."
          )
        }
        
        write_projection_netcdf(
          array_in = yhat_ensemble,
          fn_out = fn_projyhat,
          lon = lon,
          lat = lat,
          depth = attr(yhat_ensemble, "depth"),
          month = month,
          bootstrap = bootstrap,
          varname = "proj_yhat",
          units = units_target,
          longname = "CEPHALOPOD ensemble 3-D projection"
        )
      }
      
      
      # --- 8.2 Ensemble VIP
        
      vip_ensemble <- as.data.frame(MODEL$ENSEMBLE$vip)
      
      if ("variable" %in% names(vip_ensemble)) {
        vip_ensemble$variable <- as.character(vip_ensemble$variable)
      }
      
      feather::write_feather(vip_ensemble, file.path(dir_ens, "vip.feather"))
      
      
      # --- 8.3 Ensemble evaluation
      #
      # The ensemble prediction at an observation is calculated as the mean
      # across the retained model predictions for the same CV fold and row.
      #
      
      eval_ensemble <- build_ensemble_eval_table(MODEL=MODEL, QUERY=QUERY, models=models_to_save)
      feather::write_feather(eval_ensemble,file.path(dir_ens,"eval.feather"))
    }
    
  }
  
  message("DOWNSTREAM: finished run ", run_name)
  
  invisible(NULL)
}



# =============================================================================
# Evaluation table
# =============================================================================

build_eval_table <- function(MODEL, QUERY, model) {
  
  final_fit_list <- MODEL[[model]]$final_fit
  
  if (is.null(final_fit_list) ||
      length(final_fit_list) == 0) {
    
    warning(
      "No final_fit found for model ",
      model
    )
    
    return(NULL)
  }
  
  
  # Observation metadata
  
  obs <- QUERY$S
  
  obs$observation_id <- seq_len(nrow(obs))
  obs$longitude <- as.numeric(obs$decimallongitude)
  obs$latitude <- as.numeric(obs$decimallatitude)
  
  # Depth is present in the current occurrence / abundance queries.
  # For future 3-D workflows it should remain a numeric column in QUERY$S.
  if ("depth" %in% names(obs)) {
    obs$depth <- suppressWarnings(
      as.numeric(obs$depth)
    )
  } else {
    obs$depth <- NA_real_
  }
  
  if ("month" %in% names(obs)) {
    obs$time <- obs$month
  } else {
    obs$time <- NA_real_
  }
  
  
  # Extract predictions for every CV fold
  
  eval_list <- lapply(
    seq_along(final_fit_list),
    function(fold) {
      
      pred <- final_fit_list[[fold]] %>%
        collect_predictions() %>%
        as.data.frame()
      
      # last_fit() normally supplies .row. If it is absent, fall back to the
      # row ordering of the predictions.
      if (".row" %in% names(pred)) {
        row_id <- pred$.row
      } else {
        row_id <- seq_len(nrow(pred))
      }
      
      tibble::tibble(
        observation_id = row_id,
        fold = fold,
        observed = pred$measurementvalue,
        predicted = pred$.pred
      )
    }
  ) %>%
    dplyr::bind_rows()
  
  
  # Join observation metadata
  eval <- eval_list %>%
    dplyr::left_join(
      obs %>%
        dplyr::select(
          observation_id,
          longitude,
          latitude,
          depth,
          time
        ),
      by = "observation_id"
    ) %>%
    dplyr::mutate(
      residual = observed - predicted
    ) %>%
    dplyr::select(
      longitude,
      latitude,
      depth,
      time,
      observed,
      predicted,
      residual,
      fold
    )
  
  return(eval)
}



# =============================================================================
# Ensemble evaluation table
# =============================================================================

build_ensemble_eval_table <- function(MODEL, QUERY, models) {
  
  if (length(models) == 0) {
    return(NULL)
  }
  
  
  eval_models <- lapply(
    models,
    function(model) {
      
      final_fit_list <- MODEL[[model]]$final_fit
      
      if (is.null(final_fit_list)) {
        return(NULL)
      }
      
      lapply(
        seq_along(final_fit_list),
        function(fold) {
          
          pred <- final_fit_list[[fold]] %>%
            collect_predictions() %>%
            as.data.frame()
          
          if (".row" %in% names(pred)) {
            row_id <- pred$.row
          } else {
            row_id <- seq_len(nrow(pred))
          }
          
          tibble::tibble(
            model = model,
            observation_id = row_id,
            fold = fold,
            predicted = pred$.pred
          )
        }
      ) %>%
        dplyr::bind_rows()
    }
  ) %>%
    dplyr::bind_rows()
  
  
  if (nrow(eval_models) == 0) {
    return(NULL)
  }
  
  
  # Average model predictions within each observation and CV fold
  
  ensemble_pred <- eval_models %>%
    dplyr::group_by(
      observation_id,
      fold
    ) %>%
    dplyr::summarise(
      predicted = mean(
        predicted,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  
  # Observation metadata
  
  obs <- QUERY$S
  
  obs$observation_id <- seq_len(
    nrow(obs)
  )
  
  obs$longitude <- as.numeric(
    obs$decimallongitude
  )
  
  obs$latitude <- as.numeric(
    obs$decimallatitude
  )
  
  if ("depth" %in% names(obs)) {
    obs$depth <- suppressWarnings(
      as.numeric(obs$depth)
    )
  } else {
    obs$depth <- NA_real_
  }
  
  if ("month" %in% names(obs)) {
    obs$time <- obs$month
  } else {
    obs$time <- NA_real_
  }
  
  
  # Observed values
  
  observed <- QUERY$Y$measurementvalue
  
  
  eval <- ensemble_pred %>%
    dplyr::mutate(
      observed = observed[observation_id],
      residual = observed - predicted
    ) %>%
    dplyr::left_join(
      obs %>%
        dplyr::select(
          observation_id,
          longitude,
          latitude,
          depth,
          time
        ),
      by = "observation_id"
    ) %>%
    dplyr::select(
      longitude,
      latitude,
      depth,
      time,
      observed,
      predicted,
      residual,
      fold
    )
  
  return(eval)
}