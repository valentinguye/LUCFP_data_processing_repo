

# OLD CODE FOR HETEROGENEITY ANALYSIS ---------------------------------------------------------------
# compare_APEs_across_groups <- function(group1, group2, m0 = 0, alternative = "two.sided") { 
#   # important that it's selecting the "1" (resp. 2) row, and not the "Estimate" (resp. "SE") row in cases where there are several 
#   # rows with the same names, i.e. if APEs of interaction effects have been computed as well. 
#   ape1 <- ape_mat[1,group1] 
#   ape2 <- ape_mat[1,group2]
#   # # n1 <- ape_mat["Observations","Sumatra industrial"]
#   # # n2 <- ape_mat["Observations","Sumatra smallholders"]
#   sigma1 <- ape_mat[2,group1]^2
#   sigma2 <- ape_mat[2,group2]^2
#   
#   statistic <- (ape1 - ape2 - m0) / sqrt(sigma1 + sigma2)
#   
#   pval <- if (alternative == "two.sided") { 
#     2 * pnorm(abs(statistic), lower.tail = FALSE) 
#   } else if (alternative == "less") { 
#     pnorm(statistic, lower.tail = TRUE) 
#   } else { 
#     pnorm(statistic, lower.tail = FALSE) 
#   } 
#   # LCL <- (M1 - M2 - S * qnorm(1 - alpha / 2)) UCL <- (M1 - M2 + S * qnorm(1 - alpha / 2)) value <- list(mean1 = M1, mean2 = M2, m0 = m0, sigma1 = sigma1, sigma2 = sigma2, S = S, statistic = statistic, p.value = p, LCL = LCL, UCL = UCL, alternative = alternative) 
#   # print(sprintf("P-value = %g",p)) # print(sprintf("Lower %.2f%% Confidence Limit = %g", 
#   # alpha, LCL)) # print(sprintf("Upper %.2f%% Confidence Limit = %g", # alpha, UCL)) return(value) } test <- t.test_knownvar(dat1$sample1, dat1$sample2, V1 = 1, V2 = 1 )
#   return(pval)
# }

compare_coeffs_across_groups <- function(group1, group2, m0 = 0, alternative = "two.sided") { 
  # important that it's selecting the "1" (resp. 2) row, and not the "Estimate" (resp. "SE") row in cases where there are several 
  # rows with the same names, i.e. if APEs of interaction effects have been computed as well. 
  reg_sum1 <- res_data_list_full[[group1]][[1]] %>% summary()
  reg_sum2 <- res_data_list_full[[group2]][[1]] %>% summary()
  
  coeff1 <- reg_sum1$coefficients[1]
  coeff2 <- reg_sum2$coefficients[1]
  
  sigma1 <- reg_sum1$se[1]^2
  sigma2 <- reg_sum2$se[1]^2
  
  statistic <- (coeff1 - coeff2 - m0) / sqrt(sigma1 + sigma2)
  
  pval <- if (alternative == "two.sided") { 
    2 * pnorm(abs(statistic), lower.tail = FALSE) 
  } else if (alternative == "less") { 
    pnorm(statistic, lower.tail = TRUE) 
  } else { 
    pnorm(statistic, lower.tail = FALSE) 
  } 
  # LCL <- (M1 - M2 - S * qnorm(1 - alpha / 2)) UCL <- (M1 - M2 + S * qnorm(1 - alpha / 2)) value <- list(mean1 = M1, mean2 = M2, m0 = m0, sigma1 = sigma1, sigma2 = sigma2, S = S, statistic = statistic, p.value = p, LCL = LCL, UCL = UCL, alternative = alternative) 
  # print(sprintf("P-value = %g",p)) # print(sprintf("Lower %.2f%% Confidence Limit = %g", 
  # alpha, LCL)) # print(sprintf("Upper %.2f%% Confidence Limit = %g", # alpha, UCL)) return(value) } test <- t.test_knownvar(dat1$sample1, dat1$sample2, V1 = 1, V2 = 1 )
  return(pval)
}
# fill the matrix of p values of equality tests BY ROW

# comp_ape_mat["industrial = smallholders",paste0("both_a_all")] 
i_vs_sm <- compare_APEs_across_groups(group1 = paste0("both_i_all"), 
                                      group2 = paste0("both_sm_all")) 

# notill_vs_ill <- compare_APEs_across_groups(group1 = paste0("both_i_no_ill2"), 
#                                             group2 = paste0("both_i_ill2")) 

leg_vs_ill <- compare_APEs_across_groups(group1 = paste0("both_i_in_concession"), 
                                         group2 = paste0("both_i_ill2")) 


illindus_vs_sm <- compare_APEs_across_groups(group1 = paste0("both_i_ill2"), 
                                             group2 = paste0("both_sm_all")) 

# Old code to produce former table of p-values. 
# groups <- c("both_i_no_ill2", "both_i_ill2",
#             "both_i_all",
#             "both_sm_all")
# 
# comparisons <- c("industrial = smallholders",
#                  "legal indus. = illegal indus.")#,"primary forest = broadly def. forest"
# 
# comp_ape_mat <- matrix(ncol = length(groups), nrow = length(comparisons), data = NA)
# colnames(comp_ape_mat) <- groups
# row.names(comp_ape_mat) <- comparisons
# 
# for(SIZE in size_list){
#   comp_ape_mat["legal = illegal",paste0("both_",SIZE,"_all")] <- compare_APEs_across_groups(group1 = paste0("both_",SIZE,"_no_ill2"), 
#                                                                                             group2 = paste0("both_",SIZE,"_ill2")) 
# }
# 
# comp_ape_mat <- comp_ape_mat %>% formatC(digits = 4, format = "f")
# comp_ape_mat[comp_ape_mat=="   NA"] <- "" 
# comp_ape_mat
# colnames(comp_ape_mat) <- NULL
# 
# options(knitr.table.format = "latex")
# kable(comp_ape_mat, booktabs = T, align = "c",
#       caption = "p-values from equality tests of price elasticities") %>% #of 1 percentage change in medium-run price signal
#   kable_styling(latex_options = c("scale_down", "hold_position")) %>%
#   add_header_above(c("Ho" = 1,
#                      "Legal" = 1,
#                      "Illegal" = 1,
#                      "All" = 1,
#                      " " = 1,
#                      " " = 1),
#                    bold = F,
#                    align = "c") %>%  
#   add_header_above(c("Deforestation for:" = 1,
#                      "All plantations" = 3,
#                      "Industrial plantations" = 1,
#                      "Smallholder plantations" = 1),
#                    bold = F,
#                    align = "c") %>%
#   column_spec(column = 1,
#               width = "13em",
#               latex_valign = "m") %>% 
#   column_spec(column = c(2:(ncol(comp_ape_mat))),
#               width = "5em",
#               latex_valign = "m")



# ## ## ## # ## ## ## # ## ## ## # ## ## ## # ## ## ## # ## ## ## # ## ## ## # ## ## ## # ## ## ## 

## infrastructure to store results
ape_mat_list <- list()
elm <- 1

# make the APE. Note the difference: Standard errors are returnd with make_APEs_1regr, not CI95. 
rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_full)){
  res_data <- res_data_list_full[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs_1regr(res_data = res_data_list_full[[REGELM]]) # and for the same reason, this cannot be wrapped in other functions
  rm(d_clean, reg_res, res_data)
}
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()

rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_full)){
  res_data <- res_data_list_full[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat[[REGELM]] <- make_APEs() # and for the same reason, this cannot be wrapped in other functions
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
ape_mat <- bind_cols(ape_mat)  %>% as.matrix()

row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
colnames(ape_mat) <- names(res_data_list_full)
# keep ape_mat like this for later comparisons between APEs

# Code to extract regression coefficients
table_names = c("Estimate", 
                "95% CI", 
                "Estimate Interaction",
                "95% CI Interaction",
                "Estimate # reachable mills",
                "95% CI # reachable mills",
                "Estimate Initial forest cover (%)",
                "95% CI Initial forest cover (%)",
                "Estimate Mill ownership dom. private (%)",
                "95% CI Mill ownership dom. private (%)",
                "Estimate Mill ownership foreign (%)",
                "95% CI Mill ownership foreign (%)",
                "Estimate Mill CPO exports (%)",
                "95% CI Mill CPO exports (%)",
                "p-value",
                "Observations", 
                "Clusters")
tab_df = matrix(ncol = length(res_data_list_byplantation), 
                nrow = length(table_names)) 
row.names(tab_df) <- table_names

REGELM <- elm -1 

cluster_level <-  "reachable^illegal2" #  "reachable" #

for(REGELM in 1:length(res_data_list_byplantation)){
  reg_name <- res_data_list_byplantation[REGELM] %>% names() %>% print()
  if(cluster_level == "reachable^illegal2" & grepl("ill_or_concession", reg_name)){
    cluster_level <- "reachable^ill_or_concession"
  }
  d_clean = res_data_list_byplantation[[REGELM]][[2]]
  reg_res = res_data_list_byplantation[[REGELM]][[1]] #%>% summary(cluster = cluster_level) 
  
  CI_table =  etable(reg_res,
                     coefstat = "confint",
                     tex = F, 
                     se.below = TRUE,
                     digits = "s3",
                     signif.code=NA,
                     depvar = FALSE)
  
  row.names(CI_table) = CI_table[,1]
  # APE elasticity 
  ape_mat1 <- make_APEs(rounding = 3)
  # Works but useless.. 
  # if(nrow(ape_mat1)==6){row.names(ape_mat1)  = table_names[c(1:4, 10,11)]}
  # if(nrow(ape_mat1)==8){row.names(ape_mat1)  = table_names[c(1:6, 10,11)]}
  # if(nrow(ape_mat1)==10){row.names(ape_mat1) = table_names[c(1:8, 10,11)]}
  
  # Main price elasticity estimate and sample sizes
  tab_df[1:2,REGELM] <- ape_mat1[1:2]
  tab_df[(nrow(tab_df)-1):nrow(tab_df),REGELM] <- ape_mat1[(nrow(ape_mat1)-1):nrow(ape_mat1)] 
  
  
  # identify illegal var used 
  illegal_vars <- grep("ill", names(reg_res$coefficients), value = TRUE)
  illegal_var <- illegal_vars[1]
  illint_var <- illegal_vars[2]
  # Coefficients
  tab_df["Estimate Illegality", REGELM] <- reg_res$coefficients[illint_var] %>% round(3)
  tab_df["95% CI Illegality",   REGELM] <- CI_table[match(illint_var, row.names(CI_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate Illegality", REGELM] <- reg_res$coefficients[illint_var] %>% round(3)
  tab_df["95% CI Illegality",   REGELM] <- CI_table[match(illint_var, row.names(CI_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate # reachable mills", REGELM] <- reg_res$coefficients["n_reachable_umlXln_wa_cpo_price_imp1_4ya_lag1"] %>% round(3)
  tab_df["95% CI # reachable mills",   REGELM] <- CI_table[match("n_reachable_umlXln_wa_cpo_price_imp1_4ya_lag1", row.names(CI_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate Initial forest cover (%)", REGELM] <- reg_res$coefficients["pct_pfc2000_totalXln_wa_cpo_price_imp1_4ya_lag1"] %>% round(3)
  tab_df["95% CI Initial forest cover (%)",   REGELM] <- CI_table[match("pct_pfc2000_totalXln_wa_cpo_price_imp1_4ya_lag1", row.names(CI_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate Mill ownership dom. private (%)", REGELM] <- reg_res$coefficients["wa_pct_own_nat_priv_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1"] %>% round(3)
  tab_df["95% CI Mill ownership dom. private (%)",   REGELM] <- CI_table[match("wa_pct_own_nat_priv_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1", row.names(CI_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate Mill ownership foreign (%)", REGELM] <- reg_res$coefficients["wa_pct_own_for_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1"] %>% round(3)
  tab_df["95% CI Mill ownership foreign (%)",   REGELM] <- CI_table[match("wa_pct_own_for_imp_lag1Xln_wa_cpo_price_imp1_4ya_lag1", row.names(CI_table))+1,"reg_res"] %>% str_squish()
  
  tab_df["Estimate Mill CPO exports (%)", REGELM] <- reg_res$coefficients["wa_prex_cpo_imp1Xln_wa_cpo_price_imp1_4ya_lag1"] %>% round(3)
  tab_df["95% CI Mill CPO exports (%)",   REGELM] <- CI_table[match("wa_prex_cpo_imp1Xln_wa_cpo_price_imp1_4ya_lag1", row.names(CI_table))+1,"reg_res"] %>% str_squish()
  
  # Equal variance tests
  d_clean$resid = reg_res$residuals
  # unclustered test
  res_variance <- d_clean %>%
    group_by(lonlat, year, !!as.symbol(illegal_var)) %>% 
    summarise(var_resid = mean(resid^2), .groups = "drop")
  t.test(as.formula(paste0("var_resid ~ ",illegal_var)), data = res_variance) %>% print()
  
  # clustered test
  res_variance <- d_clean %>%
    group_by(reachable, !!as.symbol(illegal_var)) %>% # reachable, 
    summarise(var_resid = mean(resid^2), .groups = "drop")
  eqvars_clust_test <- t.test(as.formula(paste0("var_resid ~ ",illegal_var)), data = res_variance) %>% print()
  
  # save p-value of clustered test 
  tab_df["p-value", REGELM] =  eqvars_clust_test$p.value %>% round(3)
}

tab_df




REGELM <- elm - 1
res_data <- res_data_list_byplantation[[REGELM]]
reg_res <- res_data[[1]]
d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
rm(res_data)
ape_mat1 <- make_APEs(rounding = 3) # and for the same reason, this cannot be wrapped in other functions
reg_res %>% summary(cluster = "reachable")

res_data_list_byplantation[[REGELM]][[1]]$coeftable
## PARTIAL EFFECTS
# rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
ape_mat = list()
reg_res_list = list()
for(REGELM in 1:length(res_data_list_byplantation)){
  # get estimation results and data
  # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
  res_data <- res_data_list_byplantation[[REGELM]]
  reg_res <- res_data[[1]]
  d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
  rm(res_data)
  ape_mat1 <- make_APEs(rounding = 3) # and for the same reason, this cannot be wrapped in other functions
  reg_res %>% summary(cluster = "reachable")
  # arrange depending on interaction set - must make 12 rows in sum 
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byINDUS"){ 
    # This has 6 rows. Add 6 rows after 4th one. 
    ape_mat1 <- rbind(matrix(ape_mat1[1:4,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 6, data = ""), 
                      matrix(ape_mat1[5:6,]))
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byINDUSandIFC"){ 
    # This has 8 rows. Add 4 rows after 4th one (4 last rows are IFC and N)
    ape_mat1 <- rbind(matrix(ape_mat1[1:4,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 4, data = ""), 
                      matrix(ape_mat1[5:8,]))
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byINDUSandILL"){ 
    # This has 8 rows. Add 4 rows after 6th one (2 last rows are N)
    ape_mat1 <- rbind(matrix(ape_mat1[1:6,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 4, data = ""), 
                      matrix(ape_mat1[7:8,]))
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byINDUSandILLandIFC"){ 
    # This has 10 rows. Add 2 rows after 6th one (4 last rows are IFC and N)
    ape_mat1 <- rbind(matrix(ape_mat1[1:6,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 2, data = ""), 
                      matrix(ape_mat1[7:10,]))
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byILLINDUS"){ 
    # This has 6 rows. Add 4 rows after 2nd one, and 2 rows after the 4th one.  
    ape_mat1 <- rbind(matrix(ape_mat1[1:2,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 4, data = ""), 
                      matrix(ape_mat1[3:4,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 2, data = ""),
                      matrix(ape_mat1[5:6,])
    )
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byILLINDUSandIFC"){ 
    # This has 8 rows. Add 4 rows after 2nd one
    ape_mat1 <- rbind(matrix(ape_mat1[1:2,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 4, data = ""), 
                      matrix(ape_mat1[3:8,])
    )
  }
  if(names(res_data_list_byplantation[REGELM]) == "both_a_byIFC"){ 
    # This has 6 rows. Add 6 rows after 2nd one
    ape_mat1 <- rbind(matrix(ape_mat1[1:2,]), 
                      matrix(ncol = ncol(ape_mat1), nrow = 6, data = ""), 
                      matrix(ape_mat1[3:6,])
    )
  }
  ape_mat[[REGELM]] <- ape_mat1
  reg_res_list[[REGELM]] <- reg_res
  rm(d_clean, reg_res)
}
lapply(reg_res_list, FUN = function(x) x$convStatus)



rm(res_data_list_byplantation)


# 1. by INDUS 
res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
                                                   outcome_variable = paste0("lucpfap_pixelcount"),
                                                   interaction_terms =           c("share_indus"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                   controls = c("n_reachable_uml", "share_indus"),
                                                   n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                   offset = FALSE)
names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byINDUS")
elm <- elm + 1

# 2. by INDUS and IFC
res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
                                                   outcome_variable = paste0("lucpfap_pixelcount"),
                                                   interaction_terms =           c("share_indus", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
                                                   controls = c("n_reachable_uml", "share_indus", "pct_pfc2000_total"),
                                                   n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
                                                   offset = FALSE)
names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byINDUSandIFC")
elm <- elm + 1


# # 3. by INDUSxILL
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfip_pixelcount"),
#                                                    interaction_terms =           c("share_illindus"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_illindus"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byILLINDUS_inINDUS")
# elm <- elm + 1
# # 4. by INDUSxILL and IFC, controlling for indus
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfip_pixelcount"),
#                                                    interaction_terms =           c("share_illindus", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_illindus", "pct_pfc2000_total", "share_indus"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byILLINDUSandIFC_ctrlINDUS")
# elm <- elm + 1

# # 3. by INDUS and ILL
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfap_pixelcount"),
#                                                    interaction_terms =           c("share_indus", "illegal2"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_indus", "illegal2"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byINDUSandILL")
# elm <- elm + 1
# 
# # 4. by INDUS and ILL and IFC
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfap_pixelcount"),
#                                                    interaction_terms =           c("share_indus", "illegal2", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_indus", "illegal2", "pct_pfc2000_total"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byINDUSandILLandIFC")
# elm <- elm + 1
# 
# # 5. by INDUSxILL
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfap_pixelcount"),
#                                                    interaction_terms =           c("share_illindus"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_illindus"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byILLINDUS")
# elm <- elm + 1
# 
# # 6. by INDUSxILL and IFC
# res_data_list_byplantation[[elm]] <- make_base_reg(island = ISL, illegal = ILL,
#                                                    outcome_variable = paste0("lucpfap_pixelcount"),
#                                                    interaction_terms =           c("share_illindus", "pct_pfc2000_total"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "share_illindus", "pct_pfc2000_total"),
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
# names(res_data_list_byplantation)[elm] <- paste0(ISL,"_a_byILLINDUSandIFC")
# elm <- elm + 1


# # OR 
# ### For each plantation type 
# res_data_list_interact_ifc <- list()
# elm <- 1
# 
# size_list <- list("i","sm", "unr", "a")
# 
# # legality definition
# ill_def <- 2
# ill_status <- c("in_concession", paste0("ill",ill_def), "all")
# 
# for(SIZE in size_list){
#   if(SIZE == "i"){
#     # Industrial by illegal status
#     for(ILL in ill_status){
#       res_data_list_interact_ifc[[elm]] <- make_base_reg(island = ISL,
#                                                      outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                      interaction_terms = c("pfc2000_total_pixelcount"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                      controls = c("n_reachable_uml", "pfc2000_total_pixelcount"),
#                                                      illegal = ILL,
#                                                      n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                      offset = FALSE)
#       names(res_data_list_interact_ifc)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
#       elm <- elm + 1
#     }
#   } else {
#     # If size == sm or a, we do want to include all legal statuses.
#     # If size == unr, the selection of illegal indus is handled in make_base_reg 
#     ILL <- "all"
#     res_data_list_interact_ifc[[elm]] <- make_base_reg(island = ISL,
#                                                    outcome_variable = paste0("lucpf",SIZE,"p_pixelcount"),
#                                                    interaction_terms = c("pfc2000_total_pixelcount"),# "wa_pct_own_nat_priv_imp","wa_pct_own_for_imp",
#                                                    controls = c("n_reachable_uml", "pfc2000_total_pixelcount"),
#                                                    illegal = ILL,
#                                                    n_iter_glm = 2000,  # this is going to be long, but it's necessary. 
#                                                    offset = FALSE)
#     names(res_data_list_interact_ifc)[elm] <- paste0(ISL,"_",SIZE, "_",ILL)
#     elm <- elm + 1
#   }}
# 
# ## PARTIAL EFFECTS
# rm(ape_mat, d_clean) # it's necessary that no object called d_clean be in memory at this point, for vcov.fixest to fetch the correct data. 
# ape_mat = list()
# reg_res_list = list()
# for(REGELM in 1:length(res_data_list_interact_ifc)){
#   # get estimation results and data
#   # this is done outside the function now (R >= 4.5), otherwise not working (an environment problem)
#   res_data <- res_data_list_interact_ifc[[REGELM]]
#   reg_res <- res_data[[1]]
#   d_clean <- res_data[[2]] # it's necessary that the object is named equally to the data that was used in estimation.
#   rm(res_data)
#   ape_mat[[REGELM]] <- make_APEs(reg_elm = REGELM, rounding = 7) # and for the same reason, this cannot be wrapped in other functions
#   reg_res_list[[REGELM]] <- reg_res
#   rm(d_clean, reg_res)
# }
# # ape_mat <- lapply(res_data_list_interact_ifc, FUN = make_APEs) # this was the way to do it until it stopped working, on R 4.5 
# lapply(reg_res_list, FUN = function(x) x$convStatus)
# 
# ape_mat <- bind_cols(ape_mat)  %>% as.matrix()
# row.names(ape_mat) <- c(rep(c("Estimate","95% CI"), ((nrow(ape_mat)/2)-1)), "Observations", "Clusters") 
# ape_mat
# colnames(ape_mat) <- NULL
