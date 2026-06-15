# Title: Survival analysis in the young-onset lung cancer cohort (Figure S2)
# Description: Code used to generate the six Kaplan-Meier survival curves shown in Figure S2.
# Input: combined metadata file containing survival time, survival status,
#        age at diagnosis, sex, stage, driver, and fusion status.
# Output: Kaplan-Meier plots for survival by sex, stage, fusion status,
#         driver, ALK status, and EGFR status.

library(survival)
library(survminer)
library(ggplot2)

input_file <- "combined_data_ALL_snp_ancestry.tsv"
data <- read.delim(input_file, sep = "\t")

#convert time variable from character to numeric 
data$Survival_years <- as.numeric(data$Survival_years)

# Create the survival object
# Survival_years as time, Survival_status as event (0=alive or lost to follow-up, 1=death)
surv_obj <- with(data, Surv(Survival_years, Survival_status))

# Kaplan-Meier Survival Estimate, plots data without adjustments (eg, without using Cox PH which is required for statistics)
fit_km <- survfit(surv_obj ~ 1, data = data) #to look at data, everyone in single curve
fit_km <- survfit(surv_obj ~ sex, data = data) #to look at sex, no restrictions on stage etc -- p<0.0001 for Kaplan-Meier

# Plot the Kaplan-Meier Survival Curve
ggsurvplot(fit_km, 
           data = data,
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Curve",
           surv.median.line = "hv",    # Add median survival line
           pval = TRUE) 

#since Sex was significant (p<0.0001 in favor of F) then ran cox PH and included stage_dx as covariate, excluding SCLC patients with non-TNM staging
#subset the data to exclude "Ext SCLC" and "SCLC"
data_subset <- data[data$stage_dx_whole %in% c(1, 2, 3, 4) & 
                      !data$stage_dx_whole %in% c("Ext SCLC", "SCLC"), ]
surv_obj_subset <- with(data_subset, Surv(Survival_years, Survival_status))
fit_cox <- coxph(surv_obj_subset ~ + sex + age_dx + stage_dx_whole, data = data_subset)
summary(fit_cox)

#subset the data to exclude "Ext SCLC" and "SCLC", and include age_dx 45 and under only
data_subset <- data[data$stage_dx_whole %in% c(1, 2, 3, 4) & 
                      data$age_dx <= 45 & 
                      !data$stage_dx_whole %in% c("Ext SCLC", "SCLC"), ]
surv_obj_subset <- with(data_subset, Surv(Survival_years, Survival_status))
fit_cox <- coxph(surv_obj_subset ~ + sex + age_dx + stage_dx_whole, data = data_subset)
summary(fit_cox)

# Plot the Kaplan-Meier Survival Curve
ggsurvplot(fit_km, 
           data = data,
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Curve",
           surv.median.line = "hv",    # Add median survival line
           pval = TRUE) 

#to look at only stages 1,2,3,4
# Subset the data to include only the desired stages (1, 2, 3, 4) and age_dx 45 and under
data_subset <- data[data$stage_dx_whole %in% c(1, 2, 3, 4) & data$age_dx <= 45, ]
surv_obj_subset <- with(data_subset, Surv(Survival_years, Survival_status))
fit_km <- survfit(surv_obj_subset ~ stage_dx_whole, data = data_subset)
ggsurvplot(fit_km, 
           data = data_subset,
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Curve for Age 45 and Under",
           surv.median.line = "hv",    # Add median survival line
           pval = TRUE)

fit_km <- survfit(surv_obj_subset ~ sex, data = data_subset)
ggsurvplot(fit_km, 
           data = data_subset,
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Curve for Age 45 and Under",
           surv.median.line = "hv",    # Add median survival line
           pval = TRUE)

data_subset <- data[data$stage_dx_whole == 4 & data$age_dx <= 45, ]
surv_obj_subset <- with(data_subset, Surv(Survival_years, Survival_status))
fit_km <- survfit(surv_obj_subset ~ sex, data = data_subset)
ggsurvplot(fit_km, 
           data = data_subset,
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Curve for Age 45 and Under, Stage 4 only",
           surv.median.line = "hv",    # Add median survival line
           pval = TRUE)

fit_cox <- coxph(surv_obj_subset ~ + sex + age_dx, data = data_subset)
summary(fit_cox)

# ---- Restrict to stage IV only and perform stage IV analyses ----
stage_4_data <- data[data$stage_dx == "4", ]
surv_obj_stage_4 <- with(stage_4_data, Surv(Survival_years, Survival_status))
fit_km_sex <- survfit(surv_obj_stage_4 ~ sex, data = stage_4_data)

# Plot the Kaplan-Meier Survival Curve by Sex within Stage 4
ggsurvplot(fit_km_sex, 
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve by Sex within Stage 4",
           surv.median.line = "hv",
           pval = TRUE)

fit_cox <- coxph(surv_obj_stage_4 ~ + sex + age_dx, data = stage_4_data)
summary(fit_cox)

fit_km_fusion <- survfit(surv_obj_stage_4 ~ Fusion, data = stage_4_data)

# Plot the Kaplan-Meier Survival Curve by Fusion within Stage 4
ggsurvplot(fit_km_fusion, 
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve by Fusion within Stage 4",
           surv.median.line = "hv",
           pval = TRUE)

# ---- Subset to stage IV and exclude unknown fusion status ----
#restrict to </= 45? no need if only looking for stage 4, all NVS were stages 1-3
stage_4_data <- data[data$stage_dx == "4", ]
stage_4_data_no_unknown_fusion <- stage_4_data[stage_4_data$Fusion != "Unknown", ]
surv_obj_stage_4_no_unknown_fusion <- with(stage_4_data_no_unknown_fusion, Surv(Survival_years, Survival_status))

fit_km_fusion_no_unknown <- survfit(surv_obj_stage_4_no_unknown_fusion ~ Fusion, data = stage_4_data_no_unknown_fusion)
#p=0.0046 Fusion vs Not_fusion

# Plot the Kaplan-Meier Survival Curve by Fusion within Stage 4, excluding "Unknown"
ggsurvplot(fit_km_fusion_no_unknown, 
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve by Fusion within Stage 4",
           surv.median.line = "hv",
           pval = TRUE) 

#to get cox ph and hazard ratio
fit_cox <- coxph(surv_obj_stage_4_no_unknown_fusion ~ + sex + age_dx + Fusion, data = stage_4_data_no_unknown_fusion)
summary(fit_cox)

fit_cox <- coxph(surv_obj_stage_4_no_unknown_fusion ~ + sex + age_dx + histology + Fusion, data = stage_4_data_no_unknown_fusion)
summary(fit_cox)

# Convert the Fusion variable to a factor and relevel it (to get HR of fusion vs not_fusion); HR 0.51 (0.33-0.77)
stage_4_data_no_unknown_fusion$Fusion <- factor(stage_4_data_no_unknown_fusion$Fusion, levels = c("Not_fusion", "Fusion"))
fit_cox_relevel <- coxph(surv_obj_stage_4_no_unknown_fusion ~ sex + age_dx + Fusion, data = stage_4_data_no_unknown_fusion)
summary(fit_cox_relevel)
#adjust for histology, but >90% are LUAD
fit_cox_relevel <- coxph(surv_obj_stage_4_no_unknown_fusion ~ sex + age_dx + histology + Fusion, data = stage_4_data_no_unknown_fusion)
summary(fit_cox_relevel)

#to look at YL only (</= 45) and all stages; p=0.022
data_subset <- data[data$stage_dx_whole %in% c(1, 2, 3, 4) & data$age_dx <= 45, ]
data_subset_no_unknown_fusion <- data_subset[data_subset$Fusion != "Unknown", ]
surv_obj_subset <- with(data_subset_no_unknown_fusion, Surv(Survival_years, Survival_status))
fit_km <- survfit(surv_obj_subset ~ Fusion, data = data_subset_no_unknown_fusion)
ggsurvplot(fit_km, 
           data = data_subset_no_unknown_fusion,
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Curve for Age 45 and Under, all stages",
           surv.median.line = "hv",    # Add median survival line
           pval = TRUE)

fit_cox <- coxph(surv_obj_subset ~ + sex + age_dx + Fusion, data = data_subset_no_unknown_fusion)
summary(fit_cox)

# Convert the Fusion variable to a factor and relevel it (to get HR of fusion vs not_fusion); HR 0.51 (0.33-0.77)
data_subset_no_unknown_fusion$Fusion <- factor(data_subset_no_unknown_fusion$Fusion, levels = c("Not_fusion", "Fusion"))
fit_cox_relevel <- coxph(surv_obj_subset ~ sex + age_dx + Fusion, data = data_subset_no_unknown_fusion)
summary(fit_cox_relevel)

# ---- Subset to stage IV and compare drivers: ALK, ROS1, RET, EGFR, KRAS, HER2 ----
specific_drivers <- c("ALK", "ROS1", "RET", "EGFR", "KRAS", "HER2", "Unknown")
stage_4_data <- data[data$stage_dx == "4", ]
stage_4_data_specific_drivers <- stage_4_data[stage_4_data$reported_driver_gene %in% specific_drivers, ]
surv_obj_stage_4_specific_drivers <- with(stage_4_data_specific_drivers, Surv(Survival_years, Survival_status))
fit_km_driver_specific <- survfit(surv_obj_stage_4_specific_drivers ~ reported_driver_gene, data = stage_4_data_specific_drivers)

sample_sizes <- table(stage_4_data_specific_drivers$reported_driver_gene)
print(sample_sizes)

plot <- ggsurvplot(fit_km_driver_specific, 
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve by Driver within Stage 4",
           pval = TRUE) 

# Customize the legend labels
plot$plot <- plot$plot + 
  labs(color = "Driver") +  # Set legend title
  scale_color_discrete(labels = specific_drivers)  # Set legend labels
print(plot)


# ---- Subset to stage IV and compare EGFR versus non-EGFR ----
stage_4_data <- data[data$stage_dx == "4", ]
# Create a new binary variable 'EGFR_Status'
stage_4_data$EGFR_Status <- ifelse(stage_4_data$reported_driver_gene == "EGFR", "EGFR", "Non-EGFR")
surv_obj_egfr_vs_non_egfr <- with(stage_4_data, Surv(Survival_years, Survival_status))
fit_km_egfr_vs_non_egfr <- survfit(surv_obj_egfr_vs_non_egfr ~ EGFR_Status, data = stage_4_data)

# Plot the Kaplan-Meier Survival Curve comparing EGFR vs Non-EGFR
ggsurvplot(fit_km_egfr_vs_non_egfr, 
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve: EGFR vs. Non-EGFR within Stage 4",
           pval = TRUE)

# ---- Subset to stage IV and compare ALK versus non-ALK ----
stage_4_data <- data[data$stage_dx == "4", ]
# Create a new binary variable 'ALK_Status'
stage_4_data$ALK_Status <- ifelse(stage_4_data$reported_driver_gene == "ALK", "ALK", "Non-ALK")
surv_obj_alk_vs_non_alk <- with(stage_4_data, Surv(Survival_years, Survival_status))
fit_km_alk_vs_non_alk <- survfit(surv_obj_alk_vs_non_alk ~ ALK_Status, data = stage_4_data)

# Plot the Kaplan-Meier Survival Curve comparing ALK vs Non-ALK
ggsurvplot(fit_km_alk_vs_non_alk, 
           xlab = "Years", 
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve: ALK vs. Non-ALK within Stage 4",
           pval = TRUE)

#to get cox ph and hazard ratio
fit_cox <- coxph(surv_obj_alk_vs_non_alk ~ + sex + age_dx + ALK_Status, data = stage_4_data)
summary(fit_cox)

#to relevel HR - convert ALK_Status to a factor and set 'Non-ALK' as the reference level
stage_4_data$ALK_Status <- factor(stage_4_data$ALK_Status, levels = c("Non-ALK", "ALK"))
fit_cox_relevel <- coxph(surv_obj_alk_vs_non_alk ~ sex + age_dx + ALK_Status, data = stage_4_data)
summary(fit_cox_relevel)
