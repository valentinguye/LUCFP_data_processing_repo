make_APEs_discrete_interaction_noBeta2 <- function(#res_data, 
  interaction_dummy,
  K=1, 
  controls_pe = TRUE, 
  #SE = "cluster", 
  CLUSTER = "reachable", # paste0("reachable_",interaction_dummy),#"subdistrict",
  stddev = TRUE,
  rel_price_change = 0.01, 
  abs_price_change = 1, 
  rounding = 2){
  
  # store APEs and their deltaMethod statistics in this list 
  dM_ape_roi_list <- list()
  
  ## Redefine changes in regressor to one standard deviation if asked 
  if(stddev){
    # remove fixed effect variations from the regressor
    reg_sd <- fixest::feols(fml = as.formula(paste0(
      names(reg_res$coefficients)[1],
      " ~ 1 | ", 
      paste0(reg_res$fixef_vars, collapse = "+"))),
      data = d_clean)
    
    # and take the remaining standard deviation
    rel_price_change <- sd(reg_sd$residuals)
    abs_price_change <- sd(reg_sd$residuals)
  }
  
  # identify the nature of different variables
  coeff_names <- names(coef(reg_res))
  # define actual interaction terms (possibly lagged)
  interaction_effects <- coeff_names[grepl(pattern = "X", names(coef(reg_res)))]
  others <- coeff_names[!(grepl(pattern = "X", coeff_names))]
  interaction_terms <- others[paste0(others,"X",coeff_names[1]) %in% interaction_effects]
  
  # data POINTS that deltaMethod will grab 
  # average of the regressor of interest
  reg_bar <- mean(d_clean[,coeff_names[1]]) 
  
  # averages of the interaction terms -the controls that we interacted with the regressor of interest) 
  
  # averages of dep. variable are not necessary
  
  ## repeat the following for the K regressors of interest 
  linear_ape_fml_list <- list() # to store some results, see at end of loop
  
  ### COMPUTE PARTIAL EFFECTS OF INTERACTION TERMS 
  dM_ape_roi_list[[1]] <- list()
  
  # EFFECT OF A RELATIVE PRICE CHANGE ON THE PERCENTAGE CHANGE ACROSS THE INTERACTION DUMMY CONTRAST (ILLEGAL VS NOT)
  
  # WHEN RELATIVE PRICE CHANGE IS THE PARTIAL DERIVATIVE OF PRICE APPROXIMATELY TAKEN TO PERCENTAGE POINT (by multiplying the denominator by 100) 
  interaction_name = paste0(interaction_dummy, "X", coeff_names[1])
  # ape_fml_dit = paste0(interaction_name, "*exp(",interaction_dummy,"+",interaction_name,"*",reg_bar,")*100/100")
  
  # WHEN RELATIVE PRICE CHANGE IS A FULL 1% PRICE INCREASE 
  # for this formula we need to recover the average of the regressor in level 
  Pbar = exp(d_clean$ln_wa_cpo_price_imp1_4ya_lag1) %>% mean()
  
  ape_fml_dit = paste0("100*((",1+rel_price_change,"*",Pbar,")^(",interaction_name,") - ",Pbar,"^(",interaction_name,"))")
  # ape_fml_dit = paste0("100*((",1+rel_price_change,")^(",interaction_name,") - 1)")
  # This is exp(beta2)*Pbar^(beta3)*(1.01^(beta3)-1)*100
  
  # This is exactly equal to the unsimplified formula: 
  # ape_fml_dit = paste0("(exp(",interaction_dummy,"+",interaction_name,"*log(1.01*",reg_bar,")) - exp(",interaction_dummy,"+",interaction_name,"*log(",reg_bar,")))*100")
  
  dM <- deltaMethod(object = coef(reg_res), 
                    vcov. = vcov(reg_res, cluster = CLUSTER), #se = SE,
                    g. = ape_fml_dit, 
                    rhs = 0)
  
  row.names(dM) <- NULL
  dM <- as.matrix(dM)
  dM <- dM[,c("Estimate","2.5 %","97.5 %")]
  dM_ape_roi_list[[1]][[1]] <- dM # IN FIRST POSITION
  
  # make a one column matrix with all computed APEs' estimates, LB and HB values. 
  mat <- matrix(ncol = 1, 
                nrow = length(unlist(dM_ape_roi_list)), 
                data = unlist(dM_ape_roi_list))  
  
  mat <- round(mat, digits = rounding)
  
  k  <- 1
  while(k < nrow(mat)){
    mat[k+1,] <- paste0("[",mat[k+1,],"; ",mat[k+2,],"]")
    k <- k + 3
  } 
  row.names(mat) <- c(rep(c("Estimate","CI","delete"), nrow(mat)/3))
  mat <- mat[row.names(mat)!="delete",] %>% as.matrix()
  
  # add a row to indicate controls
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "Extra controls"
  mat[row.names(mat)=="Extra controls",] <- if_else(length(interaction_terms)>1, "X", " ")
  
  # add a row to indicate FE 1 
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "Set of reachable mills"
  mat[row.names(mat)=="Set of reachable mills",] <- if_else("reachable" %in% reg_res$fixef_vars, "X", " ")
  
  # add a row to indicate FE 2
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "Plantation site"
  mat[row.names(mat)=="Plantation site",] <- if_else("lonlat" %in% reg_res$fixef_vars, "X", " ")
  
  # add a row to indicate FE 3 
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "District-year"
  mat[row.names(mat)=="District-year",] <- if_else("district_year" %in% reg_res$fixef_vars, "X", " ")
  
  # add a row for p-value of equal variance test
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "Residual variance equality"
  mat[row.names(mat)=="Residual variance equality",] <- NA
  
  # add a row with the number of observations
  mat <- rbind(mat, reg_res$nobs)
  row.names(mat)[nrow(mat)] <- "Observations"
  mat[row.names(mat)=="Observations",] <- mat[row.names(mat)=="Observations",] %>% formatC(digits = 0, format = "f")
  
  # add a row with the number of clusters
  mat <- rbind(mat, length(unique(d_clean[,CLUSTER])))
  row.names(mat)[nrow(mat)] <- "Clusters"
  mat[row.names(mat)=="Clusters",] <- mat[row.names(mat)=="Clusters",] %>% formatC(digits = 0, format = "f")
  
  rm(coeff_names, interaction_effects, others, interaction_terms, d_clean, int_term_avg, 
     ape_fml_roi, dM_ape_roi, ape_fml_it, dM_ape_roi_list, linear_ape_fml_list)
  return(mat)
}
