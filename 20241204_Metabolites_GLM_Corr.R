library('tidyverse')
library('ggplot2')
library('dplyr')
library('ggpubr')
library('see')
library('cowplot')
library('rstatix')
library('R.utils')
library('psych')
library('nlme')
library("car")
library("lme4")
library('introdataviz')
library('multcomp')

#############################################################################################################################################################
#
#
#
# PART I ################################ DATA LOADING AND PREPROCESSING ####################################################################################

# Create folder for analysis outputs 
directory_name = 'Metabolites_GLM_Corr' # Should change for each experiment type 
directory_name_agecorr = 'Metabolites_GLM_Corr_AgeCorrected' # Should change for each experiment type 
dir.create(directory_name)
dir.create(directory_name_agecorr)

# Read in experiment data file 
setwd('E://Iben_Daniel//R_Analysis_GLM_Corr_w_Correlations//')
data_all <- read.csv(file = 'E://Iben_Daniel//R_Analysis_GLM_Corr_w_Correlations//20241206_Human_CSF_ROC_analysis_Master.csv', header=TRUE, sep=',')
data_all$Age = as.numeric(data_all$Age); 
data_all$Group = as.factor(data_all$Group); 
data_all$Sex = as.factor(data_all$Sex); 

num_cols = dim(data_all)[2]; 
col_names = colnames(data_all); 

data_all_agecorr <- data_all 

age_beta_cutoff = 0.1; 
num_metabs = num_cols - 10

#############################################################################################################################################################
#
#
#
# PART II ################################ PREPROCESSED DATA PLOTTING : RELATIONSHIP TO AGE ################################################################# 

log_age_regression_stats = data.frame(matrix(NA,num_metabs,5))
colnames(log_age_regression_stats) <- c("Metabolite", "Control_B0", "SEM", "t", "p"); 

for(ii in 11:num_cols){
  metab_name = col_names[ii];
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
    
  data_regress = data.frame(Group=data_all[,4], Age=data_all[, 3], vartoplot = data_all[,ii])
  data_regress_control = subset(data_regress, Group=='Control')
  
  regress_age_cont <- lm(log(vartoplot) ~ Age, data = data_regress_control) 
  regress_age_cont_summary <- summary(regress_age_cont)
  
  age_beta = regress_age_cont[["coefficients"]][["Age"]]
  age_beta_p = regress_age_cont_summary[["coefficients"]][[8]]
  
  if (age_beta_p  < age_beta_cutoff){
    data_all_agecorr[, ii] = exp(log(data_all[, ii]) - age_beta * data_all$Age); 
  }
  
  filename_txt = paste0(directory_name,'//', metab_name,'_Regression_with_Age.txt') 
  filename_txt = file(filename_txt, 'w')

  cat("\n\n Regression with Age in CONTROL\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(regress_age_cont_summary), file = filename_txt, append = TRUE)
  
  close(filename_txt)
  
  log_age_regression_stats$Metabolite[ii-10] = metab_name_clean
  log_age_regression_stats$Control_B0[ii-10] = age_beta
  log_age_regression_stats$SEM[ii-10] = regress_age_cont_summary[["coefficients"]][[4]]
  log_age_regression_stats$t[ii-10] = regress_age_cont_summary[["coefficients"]][[6]]
  log_age_regression_stats$p[ii-10] = age_beta_p
}

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_AgeCorrected.csv') 
write.csv(data_all_agecorr, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_AgeRegression_ModelStats.csv') 
write.csv(log_age_regression_stats, filename_csv); 

##############################################################################################################################################################
#
#
#
# PART IIIA ################################ PREPROCESSED DATA PLOTTING, LINEAR MODELING, AND CORRELATIONS: AGE - UNCORRECTED ################################ 

# Calculate descriptive statistics by GROUP 
data_means <- data_all %>%
  group_by(Group) %>%
  get_summary_stats(type = "mean_sd")

# Calculate descriptive statistics by GROUP and SEX
data_means_sex <- data_all %>%
  group_by(Group, Sex) %>%
  get_summary_stats(type = "mean_sd")

# Prepare lists to store plots for final grid figures
plot_storage_violin = list(); 
plot_storage = list(); 
plot_storage_MoCA = list(); 
plot_storage_BA = list(); 
plot_storage_totaltau = list(); 
plot_storage_phosphotau = list(); 
plot_storage_AB42v40 = list(); 

# Prepare data frames to store regression and correlation statistics by metabolite
linear_age_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_age_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_age_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_age_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_MoCA_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_MoCA_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_MoCA_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_MoCA_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_BA_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_BA_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_BA_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_BA_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_totaltau_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_totaltau_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_totaltau_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_totaltau_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_phosphotau_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_phosphotau_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_phosphotau_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_phosphotau_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_AB42v40_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_AB42v40_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_AB42v40_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_AB42v40_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

for(ii in 11:num_cols){
  
  metab_name = col_names[ii];
  data_to_plot = data.frame(ID = data_all[, 1],  Group=data_all[, 4], Sex=data_all[, 2], Age=data_all[, 3], vartoplot = data_all[,ii], 
                            MoCA=data_all[, 5], Beta_amyloid = data_all[, 7], Total_tau=data_all[, 8], Phospho_tau=data_all[, 9], 
                            A42_A40_Ratio=data_all[, 10]); 
  colnames(data_to_plot) <- c("ID", "Group", "Sex", "Age", "vartoplot", "MoCA", "Beta_amyloid", "Total_tau", "Phospho_tau", "A42_A40_Ratio"); 
  
  data_series_var <- subset(data_means, variable==metab_name);
  data_series_sex_var <- subset(data_means_sex, variable==metab_name);
  
  # Plot variable relationship to group and sex 
  nrows = dim(data_to_plot)[1]
  
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 3]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2]; 
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'by group and sex')
  
  # Adapted from https://psyteachr.github.io/introdataviz/advanced-plots.html
  viol <-ggplot(datascatter, aes(x = as.factor(x), y = y, fill = Group)) +
    introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
    geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
    geom_point(aes(fill = Group), size = 5, shape = 21, position = position_jitterdodge()) +
    stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
                 position = position_dodge(.175)) +
    scale_x_discrete(name = "Sex", labels = c("F", "M")) +
    scale_y_continuous(name = metab_name_clean) +
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    ggtitle(title) +
    theme_minimal(); 
  
  plot_storage_violin[[ii-10]] = viol
  
  filename_png = paste0(directory_name, '//', metab_name, '_by_Sex_and_Group.png') 
  ggsave(filename_png, viol, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name, '//', metab_name, '_by_Sex_and_Group.pdf') 
  ggsave(filename_eps, plot = print(viol), device = pdf, width = 18, height = 6, dpi = 600)
  
  # General linear model 
  # Summary statistics by GROUP: mean and SD
  data_means_sum <- data_to_plot %>%
    group_by(Group) %>%
    get_summary_stats(vartoplot, type = "mean_sd")
  
  # Check for outliers by GROUP
  data_outliers_sum <- data_to_plot %>%
    group_by(Group) %>%
    identify_outliers(vartoplot)
  
  # Check for normality by GROUP
  data_sw_sum  = 0
  data_sw_sum <- try(data_to_plot %>%
                       group_by(Group) %>%
                       shapiro_test(vartoplot))
  rmaov_sum  = 0
  
  # Summary statistics by GROUP and SEX: mean and SD
  data_means_sex_sum <- data_to_plot %>%
    group_by(Group, Sex) %>%
    get_summary_stats(vartoplot, type = "mean_sd")
  
  # Check for outliers by GROUP and SEX
  data_outliers_sex_sum <- data_to_plot %>%
    group_by(Group, Sex) %>%
    identify_outliers(vartoplot)
  
  # Check for normality by GROUP and SEX
  data_sw_sex_sum  = 0
  data_sw_sex_sum <- try(data_to_plot %>%
                       group_by(Group, Sex) %>%
                       shapiro_test(vartoplot))
  rmaov_sum  = 0
  
  ############################ Run ANOVA on general linear model accounting for GROUP and SEX #########################
  modlm_sum <- try(glm(vartoplot ~ Group * Sex, data = data_to_plot))
  
  if (sum(residuals(modlm_sum)) != 0) {
    rmaov_sum <- try(car::Anova(modlm_sum, REML=TRUE, test="F"))
  }
  
  # Check linear model residuals for normality 
  rmaov_sum_sw = 0
  if (sum(residuals(modlm_sum)) != 0) {
    rmaov_sum_sw <- try(shapiro_test(residuals(modlm_sum)))
  }
  
  # Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  pwc_sum_group <- data_to_plot %>%
    group_by(Sex) %>%
    pairwise_t_test(
      vartoplot ~ Group, paired = FALSE, 
      p.adjust.method = "BH"
    )
  
  # GROUP: Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  glht_posthoc_group <-glht(modlm_sum, linfct = mcp(Group ="Sequen"))
  glht_posthoc_group_reported <- summary(glht_posthoc_group, test = adjusted("BH")) 
  
  posthoc_nonparametric_group <- pairwise.wilcox.test(data_to_plot$vartoplot, data_to_plot$Group, paired = FALSE, p.adjust.method = "BH")
  
  # SEX: Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  glht_posthoc_sex <-glht(modlm_sum, linfct = mcp(Sex ="Sequen"))
  glht_posthoc_sex_reported <- summary(glht_posthoc_sex, test = adjusted("BH")) 
  
  posthoc_nonparametric_sex <- pairwise.wilcox.test(data_to_plot$vartoplot, data_to_plot$Sex, paired = FALSE, p.adjust.method = "BH")
  
  
  ############################  Run ANOVA on general linear model accounting for GROUP and SEX and AGE ############################ 
  modlm_sum_age <- try(glm(vartoplot ~ Group * Sex * Age, data = data_to_plot))
  
  if (sum(residuals(modlm_sum_age)) != 0) {
    rmaov_sum_age <- try(car::Anova(modlm_sum_age, REML=TRUE, test="F"))
  }
  
  # Check linear model residuals for normality 
  rmaov_sum_sw_age = 0
  if (sum(residuals(modlm_sum_age)) != 0) {
    rmaov_sum_sw_age <- try(shapiro_test(residuals(modlm_sum_age)))
  }
  
  # Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  pwc_sum_group_age <- data_to_plot %>%
    group_by(Sex) %>%
    pairwise_t_test(
      vartoplot ~ Group, paired = FALSE, 
      p.adjust.method = "BH"
    )
  
  # GROUP: Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  glht_posthoc_group_age <-glht(modlm_sum_age, linfct = mcp(Group ="Sequen"))
  glht_posthoc_group_reported_age <- summary(glht_posthoc_group_age, test = adjusted("BH")) 
  
  posthoc_nonparametric_group_age <- pairwise.wilcox.test(data_to_plot$vartoplot, data_to_plot$Group, paired = FALSE, p.adjust.method = "BH")
  
  # SEX: Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  glht_posthoc_sex_age <-glht(modlm_sum_age, linfct = mcp(Sex ="Sequen"))
  glht_posthoc_sex_reported_age <- summary(glht_posthoc_sex_age, test = adjusted("BH")) 
  
  posthoc_nonparametric_sex_age <- pairwise.wilcox.test(data_to_plot$vartoplot, data_to_plot$Sex, paired = FALSE, p.adjust.method = "BH")
 
  
  
  ############################  Correlation and linear regression with age ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 4]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. age by group')

  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Age (years)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name, '//', metab_name, '_Scatter_vs_Age.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name, '//', metab_name, '_Scatter_vs_Age.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_age = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_age_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_age_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_age <- try(lm(y ~ x, data = datascatter)); 
  regress_age_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_age_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_age_cont_summary <- summary(regress_age_cont)
  regress_age_AD_summary <- summary(regress_age_AD)
  
  linear_age_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_age_regression_stats$Control_B0[ii-10] = regress_age_cont_summary[["coefficients"]][[2]]
  linear_age_regression_stats$Control_SEM[ii-10] = regress_age_cont_summary[["coefficients"]][[4]]
  linear_age_regression_stats$Control_t[ii-10] = regress_age_cont_summary[["coefficients"]][[6]]
  linear_age_regression_stats$Control_p[ii-10] = regress_age_cont_summary[["coefficients"]][[8]]
  linear_age_regression_stats$AD_B0[ii-10] = regress_age_AD_summary[["coefficients"]][[2]]
  linear_age_regression_stats$AD_SEM[ii-10] = regress_age_AD_summary[["coefficients"]][[4]]
  linear_age_regression_stats$AD_t[ii-10] = regress_age_AD_summary[["coefficients"]][[6]]
  linear_age_regression_stats$AD_p[ii-10] = regress_age_AD_summary[["coefficients"]][[8]]
  
  linear_age_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_age_correlation_stats$All_SpearmanRho[ii-10] = corr_age$r
  linear_age_correlation_stats$All_N[ii-10] = corr_age$n
  linear_age_correlation_stats$All_p[ii-10] = corr_age$p
  linear_age_correlation_stats$Control_SpearmanRho[ii-10] = corr_age_cont$r
  linear_age_correlation_stats$Control_N[ii-10] = corr_age_cont$n
  linear_age_correlation_stats$Control_p[ii-10] = corr_age_cont$p
  linear_age_correlation_stats$AD_SpearmanRho[ii-10] = corr_age_AD$r
  linear_age_correlation_stats$AD_N[ii-10] = corr_age_AD$n
  linear_age_correlation_stats$AD_p[ii-10] = corr_age_AD$p
  
  ############################  Correlation and linear regression with MoCA ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 6]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. Montreal Cognitive Assessment by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Montreal Cognitive Assessment (proportion of 30)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_MoCA[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name, '//', metab_name, '_Scatter_vs_MoCA.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name, '//', metab_name, '_Scatter_vs_MoCA.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_MoCA = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_MoCA_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_MoCA_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_MoCA <- try(lm(y ~ x, data = datascatter)); 
  regress_MoCA_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_MoCA_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_MoCA_cont_summary <- summary(regress_MoCA_cont)
  regress_MoCA_AD_summary <- summary(regress_MoCA_AD)
  
  linear_MoCA_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_MoCA_regression_stats$Control_B0[ii-10] = regress_MoCA_cont_summary[["coefficients"]][[2]]
  linear_MoCA_regression_stats$Control_SEM[ii-10] = regress_MoCA_cont_summary[["coefficients"]][[4]]
  linear_MoCA_regression_stats$Control_t[ii-10] = regress_MoCA_cont_summary[["coefficients"]][[6]]
  linear_MoCA_regression_stats$Control_p[ii-10] = regress_MoCA_cont_summary[["coefficients"]][[8]]
  linear_MoCA_regression_stats$AD_B0[ii-10] = regress_MoCA_AD_summary[["coefficients"]][[2]]
  linear_MoCA_regression_stats$AD_SEM[ii-10] = regress_MoCA_AD_summary[["coefficients"]][[4]]
  linear_MoCA_regression_stats$AD_t[ii-10] = regress_MoCA_AD_summary[["coefficients"]][[6]]
  linear_MoCA_regression_stats$AD_p[ii-10] = regress_MoCA_AD_summary[["coefficients"]][[8]]
  
  linear_MoCA_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_MoCA_correlation_stats$All_SpearmanRho[ii-10] = corr_MoCA$r
  linear_MoCA_correlation_stats$All_N[ii-10] = corr_MoCA$n
  linear_MoCA_correlation_stats$All_p[ii-10] = corr_MoCA$p
  linear_MoCA_correlation_stats$Control_SpearmanRho[ii-10] = corr_MoCA_cont$r
  linear_MoCA_correlation_stats$Control_N[ii-10] = corr_MoCA_cont$n
  linear_MoCA_correlation_stats$Control_p[ii-10] = corr_MoCA_cont$p
  linear_MoCA_correlation_stats$AD_SpearmanRho[ii-10] = corr_MoCA_AD$r
  linear_MoCA_correlation_stats$AD_N[ii-10] = corr_MoCA_AD$n
  linear_MoCA_correlation_stats$AD_p[ii-10] = corr_MoCA_AD$p
  
  ############################  Correlation and linear regression with beta-amyloid ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 7]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. beta-amyloid by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Beta-amyloid (pg/mL)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_BA[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name, '//', metab_name, '_Scatter_vs_Beta_Amyloid.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name, '//', metab_name, '_Scatter_vs_Beta_Amyloid.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_BA = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_BA_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_BA_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_BA <- try(lm(y ~ x, data = datascatter)); 
  regress_BA_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_BA_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_BA_cont_summary <- summary(regress_BA_cont)
  regress_BA_AD_summary <- summary(regress_BA_AD)
  
  linear_BA_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_BA_regression_stats$Control_B0[ii-10] = regress_BA_cont_summary[["coefficients"]][[2]]
  linear_BA_regression_stats$Control_SEM[ii-10] = regress_BA_cont_summary[["coefficients"]][[4]]
  linear_BA_regression_stats$Control_t[ii-10] = regress_BA_cont_summary[["coefficients"]][[6]]
  linear_BA_regression_stats$Control_p[ii-10] = regress_BA_cont_summary[["coefficients"]][[8]]
  linear_BA_regression_stats$AD_B0[ii-10] = regress_BA_AD_summary[["coefficients"]][[2]]
  linear_BA_regression_stats$AD_SEM[ii-10] = regress_BA_AD_summary[["coefficients"]][[4]]
  linear_BA_regression_stats$AD_t[ii-10] = regress_BA_AD_summary[["coefficients"]][[6]]
  linear_BA_regression_stats$AD_p[ii-10] = regress_BA_AD_summary[["coefficients"]][[8]]
  
  linear_BA_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_BA_correlation_stats$All_SpearmanRho[ii-10] = corr_BA$r
  linear_BA_correlation_stats$All_N[ii-10] = corr_BA$n
  linear_BA_correlation_stats$All_p[ii-10] = corr_BA$p
  linear_BA_correlation_stats$Control_SpearmanRho[ii-10] = corr_BA_cont$r
  linear_BA_correlation_stats$Control_N[ii-10] = corr_BA_cont$n
  linear_BA_correlation_stats$Control_p[ii-10] = corr_BA_cont$p
  linear_BA_correlation_stats$AD_SpearmanRho[ii-10] = corr_BA_AD$r
  linear_BA_correlation_stats$AD_N[ii-10] = corr_BA_AD$n
  linear_BA_correlation_stats$AD_p[ii-10] = corr_BA_AD$p
  
  ############################  Correlation and linear regression with total tau ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 8]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. total tau by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Total tau (pg/mL)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_totaltau[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name, '//', metab_name, '_Scatter_vs_Total_Tau.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name, '//', metab_name, '_Scatter_vs_Total_Tau.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_totaltau = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_totaltau_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_totaltau_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_totaltau <- try(lm(y ~ x, data = datascatter)); 
  regress_totaltau_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_totaltau_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_totaltau_cont_summary <- summary(regress_totaltau_cont)
  regress_totaltau_AD_summary <- summary(regress_totaltau_AD)
  
  linear_totaltau_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_totaltau_regression_stats$Control_B0[ii-10] = regress_totaltau_cont_summary[["coefficients"]][[2]]
  linear_totaltau_regression_stats$Control_SEM[ii-10] = regress_totaltau_cont_summary[["coefficients"]][[4]]
  linear_totaltau_regression_stats$Control_t[ii-10] = regress_totaltau_cont_summary[["coefficients"]][[6]]
  linear_totaltau_regression_stats$Control_p[ii-10] = regress_totaltau_cont_summary[["coefficients"]][[8]]
  linear_totaltau_regression_stats$AD_B0[ii-10] = regress_totaltau_AD_summary[["coefficients"]][[2]]
  linear_totaltau_regression_stats$AD_SEM[ii-10] = regress_totaltau_AD_summary[["coefficients"]][[4]]
  linear_totaltau_regression_stats$AD_t[ii-10] = regress_totaltau_AD_summary[["coefficients"]][[6]]
  linear_totaltau_regression_stats$AD_p[ii-10] = regress_totaltau_AD_summary[["coefficients"]][[8]]
  
  linear_totaltau_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_totaltau_correlation_stats$All_SpearmanRho[ii-10] = corr_totaltau$r
  linear_totaltau_correlation_stats$All_N[ii-10] = corr_totaltau$n
  linear_totaltau_correlation_stats$All_p[ii-10] = corr_totaltau$p
  linear_totaltau_correlation_stats$Control_SpearmanRho[ii-10] = corr_totaltau_cont$r
  linear_totaltau_correlation_stats$Control_N[ii-10] = corr_totaltau_cont$n
  linear_totaltau_correlation_stats$Control_p[ii-10] = corr_totaltau_cont$p
  linear_totaltau_correlation_stats$AD_SpearmanRho[ii-10] = corr_totaltau_AD$r
  linear_totaltau_correlation_stats$AD_N[ii-10] = corr_totaltau_AD$n
  linear_totaltau_correlation_stats$AD_p[ii-10] = corr_totaltau_AD$p
  
  ############################  Correlation and linear regression with phospho-tau ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 9]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. phosphorylated tau by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Phosphorylated tau (pg/mL)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_phosphotau[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name, '//', metab_name, '_Scatter_vs_Phospho_Tau.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name, '//', metab_name, '_Scatter_vs_Phospho_Tau.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_phosphotau = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_phosphotau_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_phosphotau_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_phosphotau <- try(lm(y ~ x, data = datascatter)); 
  regress_phosphotau_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_phosphotau_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_phosphotau_cont_summary <- summary(regress_phosphotau_cont)
  regress_phosphotau_AD_summary <- summary(regress_phosphotau_AD)
  
  linear_phosphotau_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_phosphotau_regression_stats$Control_B0[ii-10] = regress_phosphotau_cont_summary[["coefficients"]][[2]]
  linear_phosphotau_regression_stats$Control_SEM[ii-10] = regress_phosphotau_cont_summary[["coefficients"]][[4]]
  linear_phosphotau_regression_stats$Control_t[ii-10] = regress_phosphotau_cont_summary[["coefficients"]][[6]]
  linear_phosphotau_regression_stats$Control_p[ii-10] = regress_phosphotau_cont_summary[["coefficients"]][[8]]
  linear_phosphotau_regression_stats$AD_B0[ii-10] = regress_phosphotau_AD_summary[["coefficients"]][[2]]
  linear_phosphotau_regression_stats$AD_SEM[ii-10] = regress_phosphotau_AD_summary[["coefficients"]][[4]]
  linear_phosphotau_regression_stats$AD_t[ii-10] = regress_phosphotau_AD_summary[["coefficients"]][[6]]
  linear_phosphotau_regression_stats$AD_p[ii-10] = regress_phosphotau_AD_summary[["coefficients"]][[8]]
  
  linear_phosphotau_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_phosphotau_correlation_stats$All_SpearmanRho[ii-10] = corr_phosphotau$r
  linear_phosphotau_correlation_stats$All_N[ii-10] = corr_phosphotau$n
  linear_phosphotau_correlation_stats$All_p[ii-10] = corr_phosphotau$p
  linear_phosphotau_correlation_stats$Control_SpearmanRho[ii-10] = corr_phosphotau_cont$r
  linear_phosphotau_correlation_stats$Control_N[ii-10] = corr_phosphotau_cont$n
  linear_phosphotau_correlation_stats$Control_p[ii-10] = corr_phosphotau_cont$p
  linear_phosphotau_correlation_stats$AD_SpearmanRho[ii-10] = corr_phosphotau_AD$r
  linear_phosphotau_correlation_stats$AD_N[ii-10] = corr_phosphotau_AD$n
  linear_phosphotau_correlation_stats$AD_p[ii-10] = corr_phosphotau_AD$p
  
  ############################  Correlation and linear regression with AB42/AB40 ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 10]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. AB42 : AB40 ratio by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'AB42 : AB40 ratio',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_AB42v40[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name, '//', metab_name, '_Scatter_vs_AB42v40.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name, '//', metab_name, '_Scatter_vs_AB42v40.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_AB42v40 = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_AB42v40_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_AB42v40_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_AB42v40 <- try(lm(y ~ x, data = datascatter)); 
  regress_AB42v40_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_AB42v40_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_AB42v40_cont_summary <- summary(regress_AB42v40_cont)
  regress_AB42v40_AD_summary <- summary(regress_AB42v40_AD)
  
  linear_AB42v40_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_AB42v40_regression_stats$Control_B0[ii-10] = regress_AB42v40_cont_summary[["coefficients"]][[2]]
  linear_AB42v40_regression_stats$Control_SEM[ii-10] = regress_AB42v40_cont_summary[["coefficients"]][[4]]
  linear_AB42v40_regression_stats$Control_t[ii-10] = regress_AB42v40_cont_summary[["coefficients"]][[6]]
  linear_AB42v40_regression_stats$Control_p[ii-10] = regress_AB42v40_cont_summary[["coefficients"]][[8]]
  linear_AB42v40_regression_stats$AD_B0[ii-10] = regress_AB42v40_AD_summary[["coefficients"]][[2]]
  linear_AB42v40_regression_stats$AD_SEM[ii-10] = regress_AB42v40_AD_summary[["coefficients"]][[4]]
  linear_AB42v40_regression_stats$AD_t[ii-10] = regress_AB42v40_AD_summary[["coefficients"]][[6]]
  linear_AB42v40_regression_stats$AD_p[ii-10] = regress_AB42v40_AD_summary[["coefficients"]][[8]]
  
  linear_AB42v40_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_AB42v40_correlation_stats$All_SpearmanRho[ii-10] = corr_AB42v40$r
  linear_AB42v40_correlation_stats$All_N[ii-10] = corr_AB42v40$n
  linear_AB42v40_correlation_stats$All_p[ii-10] = corr_AB42v40$p
  linear_AB42v40_correlation_stats$Control_SpearmanRho[ii-10] = corr_AB42v40_cont$r
  linear_AB42v40_correlation_stats$Control_N[ii-10] = corr_AB42v40_cont$n
  linear_AB42v40_correlation_stats$Control_p[ii-10] = corr_AB42v40_cont$p
  linear_AB42v40_correlation_stats$AD_SpearmanRho[ii-10] = corr_AB42v40_AD$r
  linear_AB42v40_correlation_stats$AD_N[ii-10] = corr_AB42v40_AD$n
  linear_AB42v40_correlation_stats$AD_p[ii-10] = corr_AB42v40_AD$p
  
  ############################################## PRINT ALL FINDINGS #############################################
  #Define filename
  filename_txt = paste0(directory_name,'//', metab_name,'_Final_Inferential_Statistics.txt') 
  filename_txt = file(filename_txt, 'w')
  # 
  # cat("############################################# BETWEEN-GROUP STATISTICS #############################################\n", file = filename_txt, append=TRUE)
  cat("\n\n Variable Descriptives by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(data_means_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Outliers by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(data_outliers_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Shapiro-Wilk Test by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(data_sw_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Variable Descriptives by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(data_means_sex_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Outliers by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(data_outliers_sex_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Shapiro-Wilk Test by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(data_sw_sex_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n GLM by GROUP and SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(summary(modlm_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n ANOVA on GLM by GROUP and SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(get_anova_table(rmaov_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Shapiro-Wilk Test for Linear Model Residuals by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(rmaov_sum_sw, file = filename_txt, append = TRUE)
  
  cat("\n\n Pairwise Comparisons by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(pwc_sum_group), file = filename_txt, append = TRUE)
  cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
  
  cat("\n\n Post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by GROUP \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_group_reported, file = filename_txt, append = TRUE)
  
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric_group, file = filename_txt, append = TRUE)
  
  cat("\n\n Post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_sex_reported, file = filename_txt, append = TRUE)
  
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric_sex, file = filename_txt, append = TRUE)
  
  cat("\n\n GLM by GROUP and SEX and AGE \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(summary(modlm_sum_age), file = filename_txt, append = TRUE)
  
  cat("\n\n ANOVA on GLM by GROUP and SEX and AGE \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(get_anova_table(rmaov_sum_age), file = filename_txt, append = TRUE)
  
  cat("\n\n Shapiro-Wilk Test for Linear Model Residuals by GROUP and SEX and AGE \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(rmaov_sum_sw_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Pairwise Comparisons by GROUP and SEX and AGE \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(pwc_sum_group_age), file = filename_txt, append = TRUE)
  cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
  
  cat("\n\n Post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by GROUP \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_group_reported_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric_group_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_sex_reported_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric_sex_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Age\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_age, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Age: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_age_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Age: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_age_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Age\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_age), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Age: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_age_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Age: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_age_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with MoCA\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_MoCA, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with MoCA: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_MoCA_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with MoCA: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_MoCA_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with MoCA\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_MoCA), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with MoCA: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_MoCA_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with MoCA: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_MoCA_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Beta Amyloid\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_BA, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Beta Amyloid: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_BA_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Beta Amyloid: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_BA_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Beta Amyloid\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_BA), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Beta Amyloid: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_BA_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Beta Amyloid: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_BA_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Total Tau\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_totaltau, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Total Tau: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_totaltau_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Total Tau: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_totaltau_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Total Tau\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_totaltau), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Total Tau: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_totaltau_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Total Tau: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_totaltau_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Phosphorylated Tau\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_phosphotau, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Phosphorylated Tau: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_phosphotau_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Phosphorylated Tau: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_phosphotau_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Phosphorylated Tau\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_phosphotau), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Phosphorylated Tau: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_phosphotau_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Phosphorylated Tau: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_phosphotau_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with AB42:AB40 Ratio\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_AB42v40, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with AB42:AB40 Ratio: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_AB42v40_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with AB42:AB40 Ratio: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_AB42v40_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with AB42:AB40 Ratio\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_AB42v40), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with AB42:AB40 Ratio: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_AB42v40_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with AB42:AB40 Ratio: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_AB42v40_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  close(filename_txt)

}

# Correct p-values in additional models 
linear_MoCA_regression_stats$Control_p_adj = p.adjust(linear_MoCA_regression_stats[, 5], method = 'BH', n = length(linear_MoCA_regression_stats[, 5]))
linear_MoCA_regression_stats$AD_p_adj = p.adjust(linear_MoCA_regression_stats[, 9], method = 'BH', n = length(linear_MoCA_regression_stats[, 9]))

linear_MoCA_correlation_stats$All_p_adj = p.adjust(linear_MoCA_correlation_stats[, 4], method = 'BH', n = length(linear_MoCA_correlation_stats[, 4]))
linear_MoCA_correlation_stats$Control_p_adj = p.adjust(linear_MoCA_correlation_stats[, 7], method = 'BH', n = length(linear_MoCA_correlation_stats[, 7]))
linear_MoCA_correlation_stats$AD_p_adj = p.adjust(linear_MoCA_correlation_stats[, 10], method = 'BH', n = length(linear_MoCA_correlation_stats[, 10]))

linear_BA_regression_stats$Control_p_adj = p.adjust(linear_BA_regression_stats[, 5], method = 'BH', n = length(linear_BA_regression_stats[, 5]))
linear_BA_regression_stats$AD_p_adj = p.adjust(linear_BA_regression_stats[, 9], method = 'BH', n = length(linear_BA_regression_stats[, 9]))

linear_BA_correlation_stats$All_p_adj = p.adjust(linear_BA_correlation_stats[, 4], method = 'BH', n = length(linear_BA_correlation_stats[, 4]))
linear_BA_correlation_stats$Control_p_adj = p.adjust(linear_BA_correlation_stats[, 7], method = 'BH', n = length(linear_BA_correlation_stats[, 7]))
linear_BA_correlation_stats$AD_p_adj = p.adjust(linear_BA_correlation_stats[, 10], method = 'BH', n = length(linear_BA_correlation_stats[, 10]))

linear_totaltau_regression_stats$Control_p_adj = p.adjust(linear_totaltau_regression_stats[, 5], method = 'BH', n = length(linear_totaltau_regression_stats[, 5]))
linear_totaltau_regression_stats$AD_p_adj = p.adjust(linear_totaltau_regression_stats[, 9], method = 'BH', n = length(linear_totaltau_regression_stats[, 9]))

linear_totaltau_correlation_stats$All_p_adj = p.adjust(linear_totaltau_correlation_stats[, 4], method = 'BH', n = length(linear_totaltau_correlation_stats[, 4]))
linear_totaltau_correlation_stats$Control_p_adj = p.adjust(linear_totaltau_correlation_stats[, 7], method = 'BH', n = length(linear_totaltau_correlation_stats[, 7]))
linear_totaltau_correlation_stats$AD_p_adj = p.adjust(linear_totaltau_correlation_stats[, 10], method = 'BH', n = length(linear_totaltau_correlation_stats[, 10]))

linear_phosphotau_regression_stats$Control_p_adj = p.adjust(linear_phosphotau_regression_stats[, 5], method = 'BH', n = length(linear_phosphotau_regression_stats[, 5]))
linear_phosphotau_regression_stats$AD_p_adj = p.adjust(linear_phosphotau_regression_stats[, 9], method = 'BH', n = length(linear_phosphotau_regression_stats[, 9]))

linear_phosphotau_correlation_stats$All_p_adj = p.adjust(linear_phosphotau_correlation_stats[, 4], method = 'BH', n = length(linear_phosphotau_correlation_stats[, 4]))
linear_phosphotau_correlation_stats$Control_p_adj = p.adjust(linear_phosphotau_correlation_stats[, 7], method = 'BH', n = length(linear_phosphotau_correlation_stats[, 7]))
linear_phosphotau_correlation_stats$AD_p_adj = p.adjust(linear_phosphotau_correlation_stats[, 10], method = 'BH', n = length(linear_phosphotau_correlation_stats[, 10]))

linear_AB42v40_regression_stats$Control_p_adj = p.adjust(linear_AB42v40_regression_stats[, 5], method = 'BH', n = length(linear_AB42v40_regression_stats[, 5]))
linear_AB42v40_regression_stats$AD_p_adj = p.adjust(linear_AB42v40_regression_stats[, 9], method = 'BH', n = length(linear_AB42v40_regression_stats[, 9]))

linear_AB42v40_correlation_stats$All_p_adj = p.adjust(linear_AB42v40_correlation_stats[, 4], method = 'BH', n = length(linear_AB42v40_correlation_stats[, 4]))
linear_AB42v40_correlation_stats$Control_p_adj = p.adjust(linear_AB42v40_correlation_stats[, 7], method = 'BH', n = length(linear_AB42v40_correlation_stats[, 7]))
linear_AB42v40_correlation_stats$AD_p_adj = p.adjust(linear_AB42v40_correlation_stats[, 10], method = 'BH', n = length(linear_AB42v40_correlation_stats[, 10]))

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_LinearAgeRegression_ModelStats.csv') 
write.csv(linear_age_regression_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanAgeCorrelation_ModelStats.csv') 
write.csv(linear_age_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_LinearMoCARegression_ModelStats.csv') 
write.csv(linear_MoCA_regression_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanMoCACorrelation_ModelStats.csv') 
write.csv(linear_MoCA_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_LinearBARegression_ModelStats.csv') 
write.csv(linear_BA_regression_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanBACorrelation_ModelStats.csv') 
write.csv(linear_BA_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_LinearTotalTauRegression_ModelStats.csv') 
write.csv(linear_totaltau_regression_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanTotalTauCorrelation_ModelStats.csv') 
write.csv(linear_totaltau_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_LinearPhosphoTauRegression_ModelStats.csv') 
write.csv(linear_phosphotau_regression_stats, filename_csv);

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanPhosphoTauCorrelation_ModelStats.csv') 
write.csv(linear_phosphotau_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_LinearAB40v42Regression_ModelStats.csv') 
write.csv(linear_AB42v40_regression_stats, filename_csv); 

filename_csv = paste0(directory_name,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanAB40v42Correlation_ModelStats.csv') 
write.csv(linear_AB42v40_correlation_stats, filename_csv); 

############################################# MASTER GRID PLOT: AGE-UNCORRECTED ############################################# 

combined_plot_grid_violin <- plot_grid(plot_storage_violin[[1]], plot_storage_violin[[2]], plot_storage_violin[[3]], 
                                           plot_storage_violin[[4]], plot_storage_violin[[5]], plot_storage_violin[[6]], 
                                           plot_storage_violin[[7]], plot_storage_violin[[8]], plot_storage_violin[[9]], 
                                           plot_storage_violin[[10]], plot_storage_violin[[11]], plot_storage_violin[[12]], 
                                           plot_storage_violin[[13]], plot_storage_violin[[14]], plot_storage_violin[[15]], 
                                           plot_storage_violin[[16]], plot_storage_violin[[17]], plot_storage_violin[[18]], 
                                           ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name,'//Violin_Master_Grid_Figure.png') 
ggsave(filename_png, combined_plot_grid_violin, device = 'png', width = 36, height = 18, dpi = 600)
filename_eps = paste0(directory_name,'//Violin_Master_Grid_Figure.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_violin), device = pdf, width = 36, height = 18)

############################################# MASTER GRID PLOT vs AGE: AGE-UNCORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage[[1]], plot_storage[[2]], plot_storage[[3]], 
                                           plot_storage[[4]], plot_storage[[5]], plot_storage[[6]], 
                                           plot_storage[[7]], plot_storage[[8]], plot_storage[[9]], 
                                           plot_storage[[10]], plot_storage[[11]], plot_storage[[12]], 
                                           plot_storage[[13]], plot_storage[[14]], plot_storage[[15]], 
                                           plot_storage[[16]], plot_storage[[17]], plot_storage[[18]], 
                                           ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name,'//Scatter_Master_Grid_Figure.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name,'//Scatter_Master_Grid_Figure.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs MoCA: AGE-UNCORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_MoCA[[1]], plot_storage_MoCA[[2]], plot_storage_MoCA[[3]], 
                                        plot_storage_MoCA[[4]], plot_storage_MoCA[[5]], plot_storage_MoCA[[6]], 
                                        plot_storage_MoCA[[7]], plot_storage_MoCA[[8]], plot_storage_MoCA[[9]], 
                                        plot_storage_MoCA[[10]], plot_storage_MoCA[[11]], plot_storage_MoCA[[12]], 
                                        plot_storage_MoCA[[13]], plot_storage_MoCA[[14]], plot_storage_MoCA[[15]], 
                                        plot_storage_MoCA[[16]], plot_storage_MoCA[[17]], plot_storage_MoCA[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name,'//Scatter_MoCA_Master_Grid_Figure.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name,'//Scatter_MoCA_Master_Grid_Figure.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs AB: AGE-UNCORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_BA[[1]], plot_storage_BA[[2]], plot_storage_BA[[3]], 
                                        plot_storage_BA[[4]], plot_storage_BA[[5]], plot_storage_BA[[6]], 
                                        plot_storage_BA[[7]], plot_storage_BA[[8]], plot_storage_BA[[9]], 
                                        plot_storage_BA[[10]], plot_storage_BA[[11]], plot_storage_BA[[12]], 
                                        plot_storage_BA[[13]], plot_storage_BA[[14]], plot_storage_BA[[15]], 
                                        plot_storage_BA[[16]], plot_storage_BA[[17]], plot_storage_BA[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name,'//Scatter_AB_Master_Grid_Figure.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name,'//Scatter_AB_Master_Grid_Figure.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs TOTAL TAU: AGE-UNCORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_totaltau[[1]], plot_storage_totaltau[[2]], plot_storage_totaltau[[3]], 
                                        plot_storage_totaltau[[4]], plot_storage_totaltau[[5]], plot_storage_totaltau[[6]], 
                                        plot_storage_totaltau[[7]], plot_storage_totaltau[[8]], plot_storage_totaltau[[9]], 
                                        plot_storage_totaltau[[10]], plot_storage_totaltau[[11]], plot_storage_totaltau[[12]], 
                                        plot_storage_totaltau[[13]], plot_storage_totaltau[[14]], plot_storage_totaltau[[15]], 
                                        plot_storage_totaltau[[16]], plot_storage_totaltau[[17]], plot_storage_totaltau[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name,'//Scatter_totaltau_Master_Grid_Figure.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name,'//Scatter_totaltau_Master_Grid_Figure.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs PHOSPHORYLATED TAU: AGE-UNCORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_phosphotau[[1]], plot_storage_phosphotau[[2]], plot_storage_phosphotau[[3]], 
                                        plot_storage_phosphotau[[4]], plot_storage_phosphotau[[5]], plot_storage_phosphotau[[6]], 
                                        plot_storage_phosphotau[[7]], plot_storage_phosphotau[[8]], plot_storage_phosphotau[[9]], 
                                        plot_storage_phosphotau[[10]], plot_storage_phosphotau[[11]], plot_storage_phosphotau[[12]], 
                                        plot_storage_phosphotau[[13]], plot_storage_phosphotau[[14]], plot_storage_phosphotau[[15]], 
                                        plot_storage_phosphotau[[16]], plot_storage_phosphotau[[17]], plot_storage_phosphotau[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name,'//Scatter_phosphotau_Master_Grid_Figure.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name,'//Scatter_phosphotau_Master_Grid_Figure.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs AB42 : AB40: AGE-UNCORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_AB42v40[[1]], plot_storage_AB42v40[[2]], plot_storage_AB42v40[[3]], 
                                        plot_storage_AB42v40[[4]], plot_storage_AB42v40[[5]], plot_storage_AB42v40[[6]], 
                                        plot_storage_AB42v40[[7]], plot_storage_AB42v40[[8]], plot_storage_AB42v40[[9]], 
                                        plot_storage_AB42v40[[10]], plot_storage_AB42v40[[11]], plot_storage_AB42v40[[12]], 
                                        plot_storage_AB42v40[[13]], plot_storage_AB42v40[[14]], plot_storage_AB42v40[[15]], 
                                        plot_storage_AB42v40[[16]], plot_storage_AB42v40[[17]], plot_storage_AB42v40[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name,'//Scatter_AB42v40_Master_Grid_Figure.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name,'//Scatter_AB42v40_Master_Grid_Figure.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)


#############################################################################################################################################################
#
#
#
# PART IIIB ################################ PREPROCESSED DATA PLOTTING, LINEAR MODELING, AND CORRELATIONS: AGE - CORRECTED ################################# 

# Calculate descriptive statistics by GROUP 
data_means_agecorr <- data_all_agecorr %>%
  group_by(Group) %>%
  get_summary_stats(type = "mean_sd")

# Calculate descriptive statistics by GROUP and SEX
data_means_sex_agecorr <- data_all_agecorr %>%
  group_by(Group, Sex) %>%
  get_summary_stats(type = "mean_sd")

# Prepare lists to store plots for final grid figures
plot_storage_violin_agecorr = list(); 
plot_storage_agecorr = list(); 
plot_storage_MoCA_agecorr = list(); 
plot_storage_BA_agecorr = list(); 
plot_storage_totaltau_agecorr = list(); 
plot_storage_phosphotau_agecorr = list(); 
plot_storage_AB42v40_agecorr = list(); 

# Prepare data frames to store descriptive statistics and linear model parameters on age-corrected metabolite values for paper 
group_sex_desc_stats = data.frame(matrix(NA,num_metabs,9))
colnames(group_sex_desc_stats) = c("Metabolite", "AD_F_Mean", "AD_F_SD", "Cont_F_Mean", "Cont_F_SD", "AD_M_Mean", "AD_M_SD", 
                                   "Cont_M_Mean", "Cont_M_SD")

glm_group_sex_stats = data.frame(matrix(NA,num_metabs,13))
colnames(glm_group_sex_stats) <- c("Metabolite", "Group_B", "Group_SEM", "Group_t", "Group_p", "Sex_B", "Sex_SEM", "Sex_t", "Sex_p", 
                                   "Group_Sex_B", "Group_Sex_SEM", "Group_Sex_t", "Group_Sex_p"); 

# Prepare data frames to store regression and correlation statistics by metabolite
linear_age_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_age_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_age_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_age_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_MoCA_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_MoCA_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_MoCA_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_MoCA_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_BA_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_BA_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_BA_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_BA_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_totaltau_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_totaltau_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_totaltau_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_totaltau_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_phosphotau_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_phosphotau_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_phosphotau_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_phosphotau_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

linear_AB42v40_regression_stats = data.frame(matrix(NA,num_metabs,9))
colnames(linear_AB42v40_regression_stats) <- c("Metabolite", "Control_B0", "Control_SEM", "Control_t", "Control_p", "AD_B0", "AD_SEM", "AD_t", "AD_p"); 

linear_AB42v40_correlation_stats = data.frame(matrix(NA,num_metabs,10))
colnames(linear_AB42v40_correlation_stats) <- c("Metabolite", "All_SpearmanRho", "All_N", "All_p", "Control_SpearmanRho", "Control_N", "Control_p", "AD_SpearmanRho", "AD_N", "AD_p"); 

for(ii in 11:num_cols){
  
  metab_name = col_names[ii];
  data_to_plot = data.frame(ID = data_all_agecorr[, 1],  Group=data_all_agecorr[, 4], Sex=data_all_agecorr[, 2], Age=data_all_agecorr[, 3], 
                            vartoplot = data_all_agecorr[,ii], MoCA=data_all_agecorr[, 5], Beta_amyloid = data_all_agecorr[, 7], 
                            Total_tau=data_all_agecorr[, 8], Phospho_tau=data_all_agecorr[, 9], A42_A40_Ratio=data_all_agecorr[, 10]); 
  colnames(data_to_plot) <- c("ID", "Group", "Sex", "Age", "vartoplot", "MoCA", "Beta_amyloid", "Total_tau", "Phospho_tau", "A42_A40_Ratio");
  
  data_series_var <- subset(data_means, variable==metab_name);
  data_series_sex_var <- subset(data_means_sex, variable==metab_name);
  
  # Plot variable relationship to group and sex 
  nrows = dim(data_to_plot)[1]
  
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 3]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2]; 
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'by group and sex')
  
  # Adapted from https://psyteachr.github.io/introdataviz/advanced-plots.html
  viol <-ggplot(datascatter, aes(x = as.factor(x), y = y, fill = Group)) +
    introdataviz::geom_split_violin(alpha = .4, trim = FALSE) +
    geom_boxplot(width = .2, alpha = .6, fatten = NULL, show.legend = FALSE) +
    geom_point(aes(fill = Group), size = 5, shape = 21, position = position_jitterdodge()) +
    stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F, 
                 position = position_dodge(.175)) +
    scale_x_discrete(name = "Sex", labels = c("F", "M")) +
    scale_y_continuous(name = metab_name_clean) +
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    ggtitle(title) +
    theme_minimal(); 
  
  plot_storage_violin_agecorr[[ii-10]] = viol
  
  filename_png = paste0(directory_name_agecorr, '//', metab_name, '_by_Sex_and_Group_AgeCorrected.png') 
  ggsave(filename_png, viol, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name_agecorr, '//', metab_name, '_by_Sex_and_Group_AgeCorrected.pdf') 
  ggsave(filename_eps, plot = print(viol), device = pdf, width = 18, height = 6, dpi = 600)
  
  # General linear model 
  # Summary statistics by GROUP: mean and SD
  data_means_sum <- data_to_plot %>%
    group_by(Group) %>%
    get_summary_stats(vartoplot, type = "mean_sd")
  
  # Check for outliers by GROUP
  data_outliers_sum <- data_to_plot %>%
    group_by(Group) %>%
    identify_outliers(vartoplot)
  
  # Check for normality by GROUP
  data_sw_sum  = 0
  data_sw_sum <- try(data_to_plot %>%
                       group_by(Group) %>%
                       shapiro_test(vartoplot))
  rmaov_sum  = 0
  
  # Summary statistics by GROUP and SEX: mean and SD
  data_means_sex_sum <- data_to_plot %>%
    group_by(Group, Sex) %>%
    get_summary_stats(vartoplot, type = "mean_sd")
  
  group_sex_desc_stats$Metabolite[ii-10] = metab_name_clean 
  group_sex_desc_stats$AD_F_Mean[ii-10] = data_means_sex_sum[[1,5]]
  group_sex_desc_stats$AD_F_SD[ii-10] = data_means_sex_sum[[1,6]]
  group_sex_desc_stats$Cont_F_Mean[ii-10] = data_means_sex_sum[[3,5]]
  group_sex_desc_stats$Cont_F_SD[ii-10] = data_means_sex_sum[[3,6]]
  group_sex_desc_stats$AD_M_Mean[ii-10] = data_means_sex_sum[[2,5]]
  group_sex_desc_stats$AD_M_SD[ii-10] = data_means_sex_sum[[2,6]]
  group_sex_desc_stats$Cont_M_Mean[ii-10] = data_means_sex_sum[[4,5]]
  group_sex_desc_stats$Cont_M_SD[ii-10] = data_means_sex_sum[[4,6]]
  
  # Check for outliers by GROUP and SEX
  data_outliers_sex_sum <- data_to_plot %>%
    group_by(Group, Sex) %>%
    identify_outliers(vartoplot)
  
  # Check for normality by GROUP and SEX
  data_sw_sex_sum  = 0
  data_sw_sex_sum <- try(data_to_plot %>%
                           group_by(Group, Sex) %>%
                           shapiro_test(vartoplot))
  rmaov_sum  = 0
  
  ############################ Run ANOVA on general linear model accounting for GROUP and SEX #########################
  modlm_sum <- try(glm(vartoplot ~ Group * Sex, data = data_to_plot))
  modlm_sum_summary <- summary(modlm_sum)
  
  if (sum(residuals(modlm_sum)) != 0) {
    rmaov_sum <- try(car::Anova(modlm_sum, REML=TRUE, test="F"))
  }
  
  # Check linear model residuals for normality 
  rmaov_sum_sw = 0
  if (sum(residuals(modlm_sum)) != 0) {
    rmaov_sum_sw <- try(shapiro_test(residuals(modlm_sum)))
  }
  
  # Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  pwc_sum_group <- data_to_plot %>%
    group_by(Sex) %>%
    pairwise_t_test(
      vartoplot ~ Group, paired = FALSE, 
      p.adjust.method = "BH"
    )
  
  # GROUP: Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  glht_posthoc_group <-glht(modlm_sum, linfct = mcp(Group ="Sequen"))
  glht_posthoc_group_reported <- summary(glht_posthoc_group, test = adjusted("BH")) 
  
  posthoc_nonparametric_group <- pairwise.wilcox.test(data_to_plot$vartoplot, data_to_plot$Group, paired = FALSE, p.adjust.method = "BH")
  
  # SEX: Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  glht_posthoc_sex <-glht(modlm_sum, linfct = mcp(Sex ="Sequen"))
  glht_posthoc_sex_reported <- summary(glht_posthoc_sex, test = adjusted("BH")) 
  
  posthoc_nonparametric_sex <- pairwise.wilcox.test(data_to_plot$vartoplot, data_to_plot$Sex, paired = FALSE, p.adjust.method = "BH")
  
  glm_group_sex_stats$Metabolite[ii-10] = metab_name_clean
  glm_group_sex_stats$Group_B[ii-10] = modlm_sum_summary[["coefficients"]][[2]]
  glm_group_sex_stats$Group_SEM[ii-10] = modlm_sum_summary[["coefficients"]][[6]]
  glm_group_sex_stats$Group_t[ii-10] = modlm_sum_summary[["coefficients"]][[10]]
  glm_group_sex_stats$Group_p[ii-10] = modlm_sum_summary[["coefficients"]][[14]]
  glm_group_sex_stats$Sex_B[ii-10] = modlm_sum_summary[["coefficients"]][[3]]
  glm_group_sex_stats$Sex_SEM[ii-10] = modlm_sum_summary[["coefficients"]][[7]]
  glm_group_sex_stats$Sex_t[ii-10] = modlm_sum_summary[["coefficients"]][[11]]
  glm_group_sex_stats$Sex_p[ii-10] = modlm_sum_summary[["coefficients"]][[15]]
  glm_group_sex_stats$Group_Sex_B[ii-10] = modlm_sum_summary[["coefficients"]][[4]]
  glm_group_sex_stats$Group_Sex_SEM[ii-10] = modlm_sum_summary[["coefficients"]][[8]]
  glm_group_sex_stats$Group_Sex_t[ii-10] = modlm_sum_summary[["coefficients"]][[12]]
  glm_group_sex_stats$Group_Sex_p[ii-10] = modlm_sum_summary[["coefficients"]][[16]]
  
  ############################  Run ANOVA on general linear model accounting for GROUP and SEX and AGE ############################ 
  modlm_sum_age <- try(glm(vartoplot ~ Group * Sex * Age, data = data_to_plot))
  
  if (sum(residuals(modlm_sum_age)) != 0) {
    rmaov_sum_age <- try(car::Anova(modlm_sum_age, REML=TRUE, test="F"))
  }
  
  # Check linear model residuals for normality 
  rmaov_sum_sw_age = 0
  if (sum(residuals(modlm_sum_age)) != 0) {
    rmaov_sum_sw_age <- try(shapiro_test(residuals(modlm_sum_age)))
  }
  
  # Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  pwc_sum_group_age <- data_to_plot %>%
    group_by(Sex) %>%
    pairwise_t_test(
      vartoplot ~ Group, paired = FALSE, 
      p.adjust.method = "BH"
    )
  
  # GROUP: Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  glht_posthoc_group_age <-glht(modlm_sum_age, linfct = mcp(Group ="Sequen"))
  glht_posthoc_group_reported_age <- summary(glht_posthoc_group_age, test = adjusted("BH")) 
  
  posthoc_nonparametric_group_age <- pairwise.wilcox.test(data_to_plot$vartoplot, data_to_plot$Group, paired = FALSE, p.adjust.method = "BH")
  
  # SEX: Pairwise comparisons post hoc with Benjamini-Hochberg correction 
  glht_posthoc_sex_age <-glht(modlm_sum_age, linfct = mcp(Sex ="Sequen"))
  glht_posthoc_sex_reported_age <- summary(glht_posthoc_sex_age, test = adjusted("BH")) 
  
  posthoc_nonparametric_sex_age <- pairwise.wilcox.test(data_to_plot$vartoplot, data_to_plot$Sex, paired = FALSE, p.adjust.method = "BH")
  
  
  
  ############################  Correlation and linear regression with age ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 4]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  # cut_scatter <- ggscatter(datascatter, x='x', y='y', conf.int = TRUE, 
  #                          color = "Group",
  #                          ellipse = TRUE, mean.point = TRUE,
  #                          star.plot = TRUE)   + 
  #   scale_fill_brewer(palette = "Dark2", name = "Group") +
  #   theme_minimal() 
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. age by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Age (years)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_agecorr[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_Age_AgeCorrected.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_Age_AgeCorrected.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_age = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_age_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_age_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_age <- try(lm(y ~ x, data = datascatter)); 
  regress_age_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_age_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  
  linear_age_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_age_regression_stats$Control_B0[ii-10] = regress_age_cont_summary[["coefficients"]][[2]]
  linear_age_regression_stats$Control_SEM[ii-10] = regress_age_cont_summary[["coefficients"]][[4]]
  linear_age_regression_stats$Control_t[ii-10] = regress_age_cont_summary[["coefficients"]][[6]]
  linear_age_regression_stats$Control_p[ii-10] = regress_age_cont_summary[["coefficients"]][[8]]
  linear_age_regression_stats$AD_B0[ii-10] = regress_age_AD_summary[["coefficients"]][[2]]
  linear_age_regression_stats$AD_SEM[ii-10] = regress_age_AD_summary[["coefficients"]][[4]]
  linear_age_regression_stats$AD_t[ii-10] = regress_age_AD_summary[["coefficients"]][[6]]
  linear_age_regression_stats$AD_p[ii-10] = regress_age_AD_summary[["coefficients"]][[8]]
  
  linear_age_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_age_correlation_stats$All_SpearmanRho[ii-10] = corr_age$r
  linear_age_correlation_stats$All_N[ii-10] = corr_age$n
  linear_age_correlation_stats$All_p[ii-10] = corr_age$p
  linear_age_correlation_stats$Control_SpearmanRho[ii-10] = corr_age_cont$r
  linear_age_correlation_stats$Control_N[ii-10] = corr_age_cont$n
  linear_age_correlation_stats$Control_p[ii-10] = corr_age_cont$p
  linear_age_correlation_stats$AD_SpearmanRho[ii-10] = corr_age_AD$r
  linear_age_correlation_stats$AD_N[ii-10] = corr_age_AD$n
  linear_age_correlation_stats$AD_p[ii-10] = corr_age_AD$p
  
  
  
  ############################  Correlation and linear regression with MoCA ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 6]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. Montreal Cognitive Assessment by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Montreal Cognitive Assessment (proportion of 30)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_MoCA_agecorr[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_MoCA.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_MoCA.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_MoCA = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_MoCA_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_MoCA_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_MoCA <- try(lm(y ~ x, data = datascatter)); 
  regress_MoCA_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_MoCA_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_MoCA_cont_summary <- summary(regress_MoCA_cont)
  regress_MoCA_AD_summary <- summary(regress_MoCA_AD)
  
  linear_MoCA_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_MoCA_regression_stats$Control_B0[ii-10] = regress_MoCA_cont_summary[["coefficients"]][[2]]
  linear_MoCA_regression_stats$Control_SEM[ii-10] = regress_MoCA_cont_summary[["coefficients"]][[4]]
  linear_MoCA_regression_stats$Control_t[ii-10] = regress_MoCA_cont_summary[["coefficients"]][[6]]
  linear_MoCA_regression_stats$Control_p[ii-10] = regress_MoCA_cont_summary[["coefficients"]][[8]]
  linear_MoCA_regression_stats$AD_B0[ii-10] = regress_MoCA_AD_summary[["coefficients"]][[2]]
  linear_MoCA_regression_stats$AD_SEM[ii-10] = regress_MoCA_AD_summary[["coefficients"]][[4]]
  linear_MoCA_regression_stats$AD_t[ii-10] = regress_MoCA_AD_summary[["coefficients"]][[6]]
  linear_MoCA_regression_stats$AD_p[ii-10] = regress_MoCA_AD_summary[["coefficients"]][[8]]
  
  linear_MoCA_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_MoCA_correlation_stats$All_SpearmanRho[ii-10] = corr_MoCA$r
  linear_MoCA_correlation_stats$All_N[ii-10] = corr_MoCA$n
  linear_MoCA_correlation_stats$All_p[ii-10] = corr_MoCA$p
  linear_MoCA_correlation_stats$Control_SpearmanRho[ii-10] = corr_MoCA_cont$r
  linear_MoCA_correlation_stats$Control_N[ii-10] = corr_MoCA_cont$n
  linear_MoCA_correlation_stats$Control_p[ii-10] = corr_MoCA_cont$p
  linear_MoCA_correlation_stats$AD_SpearmanRho[ii-10] = corr_MoCA_AD$r
  linear_MoCA_correlation_stats$AD_N[ii-10] = corr_MoCA_AD$n
  linear_MoCA_correlation_stats$AD_p[ii-10] = corr_MoCA_AD$p
  
  
  
  ############################  Correlation and linear regression with beta-amyloid ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 7]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. beta-amyloid by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Beta-amyloid (pg/mL)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_BA_agecorr[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_Beta_Amyloid.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_Beta_Amyloid.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_BA = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_BA_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_BA_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_BA <- try(lm(y ~ x, data = datascatter)); 
  regress_BA_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_BA_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_BA_cont_summary <- summary(regress_BA_cont)
  regress_BA_AD_summary <- summary(regress_BA_AD)
  
  linear_BA_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_BA_regression_stats$Control_B0[ii-10] = regress_BA_cont_summary[["coefficients"]][[2]]
  linear_BA_regression_stats$Control_SEM[ii-10] = regress_BA_cont_summary[["coefficients"]][[4]]
  linear_BA_regression_stats$Control_t[ii-10] = regress_BA_cont_summary[["coefficients"]][[6]]
  linear_BA_regression_stats$Control_p[ii-10] = regress_BA_cont_summary[["coefficients"]][[8]]
  linear_BA_regression_stats$AD_B0[ii-10] = regress_BA_AD_summary[["coefficients"]][[2]]
  linear_BA_regression_stats$AD_SEM[ii-10] = regress_BA_AD_summary[["coefficients"]][[4]]
  linear_BA_regression_stats$AD_t[ii-10] = regress_BA_AD_summary[["coefficients"]][[6]]
  linear_BA_regression_stats$AD_p[ii-10] = regress_BA_AD_summary[["coefficients"]][[8]]
  
  linear_BA_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_BA_correlation_stats$All_SpearmanRho[ii-10] = corr_BA$r
  linear_BA_correlation_stats$All_N[ii-10] = corr_BA$n
  linear_BA_correlation_stats$All_p[ii-10] = corr_BA$p
  linear_BA_correlation_stats$Control_SpearmanRho[ii-10] = corr_BA_cont$r
  linear_BA_correlation_stats$Control_N[ii-10] = corr_BA_cont$n
  linear_BA_correlation_stats$Control_p[ii-10] = corr_BA_cont$p
  linear_BA_correlation_stats$AD_SpearmanRho[ii-10] = corr_BA_AD$r
  linear_BA_correlation_stats$AD_N[ii-10] = corr_BA_AD$n
  linear_BA_correlation_stats$AD_p[ii-10] = corr_BA_AD$p
  
  
  
  ############################  Correlation and linear regression with total tau ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 8]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. total tau by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Total tau (pg/mL)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_totaltau_agecorr[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_Total_Tau.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_Total_Tau.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_totaltau = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_totaltau_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_totaltau_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_totaltau <- try(lm(y ~ x, data = datascatter)); 
  regress_totaltau_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_totaltau_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_totaltau_cont_summary <- summary(regress_totaltau_cont)
  regress_totaltau_AD_summary <- summary(regress_totaltau_AD)
  
  linear_totaltau_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_totaltau_regression_stats$Control_B0[ii-10] = regress_totaltau_cont_summary[["coefficients"]][[2]]
  linear_totaltau_regression_stats$Control_SEM[ii-10] = regress_totaltau_cont_summary[["coefficients"]][[4]]
  linear_totaltau_regression_stats$Control_t[ii-10] = regress_totaltau_cont_summary[["coefficients"]][[6]]
  linear_totaltau_regression_stats$Control_p[ii-10] = regress_totaltau_cont_summary[["coefficients"]][[8]]
  linear_totaltau_regression_stats$AD_B0[ii-10] = regress_totaltau_AD_summary[["coefficients"]][[2]]
  linear_totaltau_regression_stats$AD_SEM[ii-10] = regress_totaltau_AD_summary[["coefficients"]][[4]]
  linear_totaltau_regression_stats$AD_t[ii-10] = regress_totaltau_AD_summary[["coefficients"]][[6]]
  linear_totaltau_regression_stats$AD_p[ii-10] = regress_totaltau_AD_summary[["coefficients"]][[8]]
  
  linear_totaltau_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_totaltau_correlation_stats$All_SpearmanRho[ii-10] = corr_totaltau$r
  linear_totaltau_correlation_stats$All_N[ii-10] = corr_totaltau$n
  linear_totaltau_correlation_stats$All_p[ii-10] = corr_totaltau$p
  linear_totaltau_correlation_stats$Control_SpearmanRho[ii-10] = corr_totaltau_cont$r
  linear_totaltau_correlation_stats$Control_N[ii-10] = corr_totaltau_cont$n
  linear_totaltau_correlation_stats$Control_p[ii-10] = corr_totaltau_cont$p
  linear_totaltau_correlation_stats$AD_SpearmanRho[ii-10] = corr_totaltau_AD$r
  linear_totaltau_correlation_stats$AD_N[ii-10] = corr_totaltau_AD$n
  linear_totaltau_correlation_stats$AD_p[ii-10] = corr_totaltau_AD$p
  
  
  
  ############################  Correlation and linear regression with phospho-tau ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 9]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. phosphorylated tau by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'Phosphorylated tau (pg/mL)',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_phosphotau_agecorr[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_Phospho_Tau.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_Phospho_Tau.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_phosphotau = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_phosphotau_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_phosphotau_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_phosphotau <- try(lm(y ~ x, data = datascatter)); 
  regress_phosphotau_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_phosphotau_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_phosphotau_cont_summary <- summary(regress_phosphotau_cont)
  regress_phosphotau_AD_summary <- summary(regress_phosphotau_AD)
  
  linear_phosphotau_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_phosphotau_regression_stats$Control_B0[ii-10] = regress_phosphotau_cont_summary[["coefficients"]][[2]]
  linear_phosphotau_regression_stats$Control_SEM[ii-10] = regress_phosphotau_cont_summary[["coefficients"]][[4]]
  linear_phosphotau_regression_stats$Control_t[ii-10] = regress_phosphotau_cont_summary[["coefficients"]][[6]]
  linear_phosphotau_regression_stats$Control_p[ii-10] = regress_phosphotau_cont_summary[["coefficients"]][[8]]
  linear_phosphotau_regression_stats$AD_B0[ii-10] = regress_phosphotau_AD_summary[["coefficients"]][[2]]
  linear_phosphotau_regression_stats$AD_SEM[ii-10] = regress_phosphotau_AD_summary[["coefficients"]][[4]]
  linear_phosphotau_regression_stats$AD_t[ii-10] = regress_phosphotau_AD_summary[["coefficients"]][[6]]
  linear_phosphotau_regression_stats$AD_p[ii-10] = regress_phosphotau_AD_summary[["coefficients"]][[8]]
  
  linear_phosphotau_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_phosphotau_correlation_stats$All_SpearmanRho[ii-10] = corr_phosphotau$r
  linear_phosphotau_correlation_stats$All_N[ii-10] = corr_phosphotau$n
  linear_phosphotau_correlation_stats$All_p[ii-10] = corr_phosphotau$p
  linear_phosphotau_correlation_stats$Control_SpearmanRho[ii-10] = corr_phosphotau_cont$r
  linear_phosphotau_correlation_stats$Control_N[ii-10] = corr_phosphotau_cont$n
  linear_phosphotau_correlation_stats$Control_p[ii-10] = corr_phosphotau_cont$p
  linear_phosphotau_correlation_stats$AD_SpearmanRho[ii-10] = corr_phosphotau_AD$r
  linear_phosphotau_correlation_stats$AD_N[ii-10] = corr_phosphotau_AD$n
  linear_phosphotau_correlation_stats$AD_p[ii-10] = corr_phosphotau_AD$p
  
  
  
  ############################  Correlation and linear regression with AB42/AB40 ############################ 
  datascatter = data.frame(matrix(NA, nrow = nrows, ncol = 3)) 
  datascatter[,1] = data_to_plot[, 10]; 
  datascatter[,2] = data_to_plot[, 5]; 
  datascatter[,3] = data_to_plot[, 2];
  colnames(datascatter) <- c('x', 'y', 'Group')
  
  datascatter_control = subset(datascatter, Group == 'Control')
  datascatter_AD = subset(datascatter, Group == 'AD')
  
  metab_name_clean = gsub('.', ' ', metab_name, fixed=TRUE)
  title = paste(metab_name_clean, 'vs. AB42 : AB40 ratio by group')
  
  cut_scatter <- ggscatter(datascatter, x='x', y='y', 
                           color = "Group", add = "reg.line", conf.int = FALSE, 
                           ellipse = FALSE, mean.point = TRUE,
                           star.plot = TRUE,  star.plot.lty = 'dashed',
                           title = title,
                           xlab = 'AB42 : AB40 ratio',
                           ylab = metab_name_clean)   + 
    scale_fill_brewer(palette = "Dark2", name = "Group") +
    theme_lucid() + 
    theme(legend.position = 'top')
  
  plot_storage_AB42v40_agecorr[[ii-10]] = cut_scatter
  
  filename_png = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_AB42v40.png') 
  ggsave(filename_png, cut_scatter, device = 'png', width = 18, height = 6, dpi = 600)
  filename_eps = paste0(directory_name_agecorr, '//', metab_name, '_Scatter_vs_AB42v40.pdf') 
  ggsave(filename_eps, plot = print(cut_scatter), device = pdf, width = 18, height = 6, dpi = 600)
  
  corr_AB42v40 = corr.test(datascatter[,1], datascatter[,2], method="spearman"); 
  corr_AB42v40_cont = corr.test(datascatter_control[,1], datascatter_control[,2], method="spearman"); 
  corr_AB42v40_AD = corr.test(datascatter_AD[,1], datascatter_AD[,2], method="spearman"); 
  
  regress_AB42v40 <- try(lm(y ~ x, data = datascatter)); 
  regress_AB42v40_cont <- try(lm(y ~ x, data = datascatter_control)); 
  regress_AB42v40_AD <- try(lm(y ~ x, data = datascatter_AD)); 
  regress_AB42v40_cont_summary <- summary(regress_AB42v40_cont)
  regress_AB42v40_AD_summary <- summary(regress_AB42v40_AD)
  
  linear_AB42v40_regression_stats$Metabolite[ii-10] = metab_name_clean
  linear_AB42v40_regression_stats$Control_B0[ii-10] = regress_AB42v40_cont_summary[["coefficients"]][[2]]
  linear_AB42v40_regression_stats$Control_SEM[ii-10] = regress_AB42v40_cont_summary[["coefficients"]][[4]]
  linear_AB42v40_regression_stats$Control_t[ii-10] = regress_AB42v40_cont_summary[["coefficients"]][[6]]
  linear_AB42v40_regression_stats$Control_p[ii-10] = regress_AB42v40_cont_summary[["coefficients"]][[8]]
  linear_AB42v40_regression_stats$AD_B0[ii-10] = regress_AB42v40_AD_summary[["coefficients"]][[2]]
  linear_AB42v40_regression_stats$AD_SEM[ii-10] = regress_AB42v40_AD_summary[["coefficients"]][[4]]
  linear_AB42v40_regression_stats$AD_t[ii-10] = regress_AB42v40_AD_summary[["coefficients"]][[6]]
  linear_AB42v40_regression_stats$AD_p[ii-10] = regress_AB42v40_AD_summary[["coefficients"]][[8]]
  
  linear_AB42v40_correlation_stats$Metabolite[ii-10] = metab_name_clean
  linear_AB42v40_correlation_stats$All_SpearmanRho[ii-10] = corr_AB42v40$r
  linear_AB42v40_correlation_stats$All_N[ii-10] = corr_AB42v40$n
  linear_AB42v40_correlation_stats$All_p[ii-10] = corr_AB42v40$p
  linear_AB42v40_correlation_stats$Control_SpearmanRho[ii-10] = corr_AB42v40_cont$r
  linear_AB42v40_correlation_stats$Control_N[ii-10] = corr_AB42v40_cont$n
  linear_AB42v40_correlation_stats$Control_p[ii-10] = corr_AB42v40_cont$p
  linear_AB42v40_correlation_stats$AD_SpearmanRho[ii-10] = corr_AB42v40_AD$r
  linear_AB42v40_correlation_stats$AD_N[ii-10] = corr_AB42v40_AD$n
  linear_AB42v40_correlation_stats$AD_p[ii-10] = corr_AB42v40_AD$p
  
  
  
  ############################################## PRINT ALL FINDINGS #############################################
  #Define filename
  filename_txt = paste0(directory_name_agecorr,'//', metab_name,'_Final_Inferential_Statistics_AgeCorrected.txt') 
  filename_txt = file(filename_txt, 'w')
  # 
  # cat("############################################# BETWEEN-GROUP STATISTICS #############################################\n", file = filename_txt, append=TRUE)
  cat("\n\n Variable Descriptives by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(data_means_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Outliers by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(data_outliers_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Shapiro-Wilk Test by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(data_sw_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Variable Descriptives by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(data_means_sex_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Outliers by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(data_outliers_sex_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Shapiro-Wilk Test by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(data_sw_sex_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n GLM by GROUP and SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(summary(modlm_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n ANOVA on GLM by GROUP and SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(get_anova_table(rmaov_sum), file = filename_txt, append = TRUE)
  
  cat("\n\n Shapiro-Wilk Test for Linear Model Residuals by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(rmaov_sum_sw, file = filename_txt, append = TRUE)
  
  cat("\n\n Pairwise Comparisons by GROUP and SEX\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(pwc_sum_group), file = filename_txt, append = TRUE)
  cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
  
  cat("\n\n Post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by GROUP \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_group_reported, file = filename_txt, append = TRUE)
  
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric_group, file = filename_txt, append = TRUE)
  
  cat("\n\n Post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_sex_reported, file = filename_txt, append = TRUE)
  
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric_sex, file = filename_txt, append = TRUE)
  
  cat("\n\n GLM by GROUP and SEX and AGE \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(summary(modlm_sum_age), file = filename_txt, append = TRUE)
  
  cat("\n\n ANOVA on GLM by GROUP and SEX and AGE \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(get_anova_table(rmaov_sum_age), file = filename_txt, append = TRUE)
  
  cat("\n\n Shapiro-Wilk Test for Linear Model Residuals by GROUP and SEX and AGE \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(rmaov_sum_sw_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Pairwise Comparisons by GROUP and SEX and AGE \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(as.data.frame(pwc_sum_group_age), file = filename_txt, append = TRUE)
  cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
  
  cat("\n\n Post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by GROUP \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_group_reported_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by GROUP\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric_group_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_sex_reported_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values) by SEX \n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric_sex_age, file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Age\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_age, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Age: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_age_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Age: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_age_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Age\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_age), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Age: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_age_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Age: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_age_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with MoCA\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_MoCA, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with MoCA: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_MoCA_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with MoCA: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_MoCA_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with MoCA\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_MoCA), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with MoCA: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_MoCA_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with MoCA: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_MoCA_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Beta Amyloid\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_BA, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Beta Amyloid: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_BA_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Beta Amyloid: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_BA_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Beta Amyloid\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_BA), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Beta Amyloid: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_BA_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Beta Amyloid: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_BA_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Total Tau\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_totaltau, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Total Tau: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_totaltau_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Total Tau: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_totaltau_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Total Tau\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_totaltau), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Total Tau: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_totaltau_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Total Tau: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_totaltau_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Phosphorylated Tau\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_phosphotau, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Phosphorylated Tau: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_phosphotau_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with Phosphorylated Tau: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_phosphotau_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Phosphorylated Tau\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_phosphotau), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Phosphorylated Tau: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_phosphotau_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with Phosphorylated Tau: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_phosphotau_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with AB42:AB40 Ratio\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_AB42v40, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with AB42:AB40 Ratio: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_AB42v40_AD, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Correlation with AB42:AB40 Ratio: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(corr_AB42v40_cont, short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with AB42:AB40 Ratio\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_AB42v40), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with AB42:AB40 Ratio: AD\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_AB42v40_AD), short=FALSE), file = filename_txt, append = TRUE)
  
  cat("\n\n Regression with AB42:AB40 Ratio: Control\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(print(summary(regress_AB42v40_cont), short=FALSE), file = filename_txt, append = TRUE)
  
  close(filename_txt)
}


filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_GLM_Group_Sex_AgeCorr_ModelStats.csv') 
write.csv(glm_group_sex_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_Descriptives_Group_Sex_AgeCorr_ModelStats.csv') 
write.csv(group_sex_desc_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_LinearAgeRegression_AgeCorr_ModelStats.csv') 
write.csv(linear_age_regression_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanAgeCorrelation_AgeCorr_ModelStats.csv') 
write.csv(linear_age_correlation_stats, filename_csv); 

# Correct p-values in additional models 
linear_MoCA_regression_stats$Control_p_adj = p.adjust(linear_MoCA_regression_stats[, 5], method = 'BH', n = length(linear_MoCA_regression_stats[, 5]))
linear_MoCA_regression_stats$AD_p_adj = p.adjust(linear_MoCA_regression_stats[, 9], method = 'BH', n = length(linear_MoCA_regression_stats[, 9]))

linear_MoCA_correlation_stats$All_p_adj = p.adjust(linear_MoCA_correlation_stats[, 4], method = 'BH', n = length(linear_MoCA_correlation_stats[, 4]))
linear_MoCA_correlation_stats$Control_p_adj = p.adjust(linear_MoCA_correlation_stats[, 7], method = 'BH', n = length(linear_MoCA_correlation_stats[, 7]))
linear_MoCA_correlation_stats$AD_p_adj = p.adjust(linear_MoCA_correlation_stats[, 10], method = 'BH', n = length(linear_MoCA_correlation_stats[, 10]))

linear_BA_regression_stats$Control_p_adj = p.adjust(linear_BA_regression_stats[, 5], method = 'BH', n = length(linear_BA_regression_stats[, 5]))
linear_BA_regression_stats$AD_p_adj = p.adjust(linear_BA_regression_stats[, 9], method = 'BH', n = length(linear_BA_regression_stats[, 9]))

linear_BA_correlation_stats$All_p_adj = p.adjust(linear_BA_correlation_stats[, 4], method = 'BH', n = length(linear_BA_correlation_stats[, 4]))
linear_BA_correlation_stats$Control_p_adj = p.adjust(linear_BA_correlation_stats[, 7], method = 'BH', n = length(linear_BA_correlation_stats[, 7]))
linear_BA_correlation_stats$AD_p_adj = p.adjust(linear_BA_correlation_stats[, 10], method = 'BH', n = length(linear_BA_correlation_stats[, 10]))

linear_totaltau_regression_stats$Control_p_adj = p.adjust(linear_totaltau_regression_stats[, 5], method = 'BH', n = length(linear_totaltau_regression_stats[, 5]))
linear_totaltau_regression_stats$AD_p_adj = p.adjust(linear_totaltau_regression_stats[, 9], method = 'BH', n = length(linear_totaltau_regression_stats[, 9]))

linear_totaltau_correlation_stats$All_p_adj = p.adjust(linear_totaltau_correlation_stats[, 4], method = 'BH', n = length(linear_totaltau_correlation_stats[, 4]))
linear_totaltau_correlation_stats$Control_p_adj = p.adjust(linear_totaltau_correlation_stats[, 7], method = 'BH', n = length(linear_totaltau_correlation_stats[, 7]))
linear_totaltau_correlation_stats$AD_p_adj = p.adjust(linear_totaltau_correlation_stats[, 10], method = 'BH', n = length(linear_totaltau_correlation_stats[, 10]))

linear_phosphotau_regression_stats$Control_p_adj = p.adjust(linear_phosphotau_regression_stats[, 5], method = 'BH', n = length(linear_phosphotau_regression_stats[, 5]))
linear_phosphotau_regression_stats$AD_p_adj = p.adjust(linear_phosphotau_regression_stats[, 9], method = 'BH', n = length(linear_phosphotau_regression_stats[, 9]))

linear_phosphotau_correlation_stats$All_p_adj = p.adjust(linear_phosphotau_correlation_stats[, 4], method = 'BH', n = length(linear_phosphotau_correlation_stats[, 4]))
linear_phosphotau_correlation_stats$Control_p_adj = p.adjust(linear_phosphotau_correlation_stats[, 7], method = 'BH', n = length(linear_phosphotau_correlation_stats[, 7]))
linear_phosphotau_correlation_stats$AD_p_adj = p.adjust(linear_phosphotau_correlation_stats[, 10], method = 'BH', n = length(linear_phosphotau_correlation_stats[, 10]))

linear_AB42v40_regression_stats$Control_p_adj = p.adjust(linear_AB42v40_regression_stats[, 5], method = 'BH', n = length(linear_AB42v40_regression_stats[, 5]))
linear_AB42v40_regression_stats$AD_p_adj = p.adjust(linear_AB42v40_regression_stats[, 9], method = 'BH', n = length(linear_AB42v40_regression_stats[, 9]))

linear_AB42v40_correlation_stats$All_p_adj = p.adjust(linear_AB42v40_correlation_stats[, 4], method = 'BH', n = length(linear_AB42v40_correlation_stats[, 4]))
linear_AB42v40_correlation_stats$Control_p_adj = p.adjust(linear_AB42v40_correlation_stats[, 7], method = 'BH', n = length(linear_AB42v40_correlation_stats[, 7]))
linear_AB42v40_correlation_stats$AD_p_adj = p.adjust(linear_AB42v40_correlation_stats[, 10], method = 'BH', n = length(linear_AB42v40_correlation_stats[, 10]))

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_LinearMoCARegression_AgeCorr_ModelStats.csv') 
write.csv(linear_MoCA_regression_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanMoCACorrelation_AgeCorr_ModelStats.csv') 
write.csv(linear_MoCA_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_LinearBARegression_AgeCorr_ModelStats.csv') 
write.csv(linear_BA_regression_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanBACorrelation_AgeCorr_ModelStats.csv') 
write.csv(linear_BA_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_LinearTotalTauRegression_AgeCorr_ModelStats.csv') 
write.csv(linear_totaltau_regression_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanTotalTauCorrelation_AgeCorr_ModelStats.csv') 
write.csv(linear_totaltau_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_LinearPhosphoTauRegression_AgeCorr_ModelStats.csv') 
write.csv(linear_phosphotau_regression_stats, filename_csv);

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanPhosphoTauCorrelation_AgeCorr_ModelStats.csv') 
write.csv(linear_phosphotau_correlation_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_LinearAB40v42Regression_AgeCorr_ModelStats.csv') 
write.csv(linear_AB42v40_regression_stats, filename_csv); 

filename_csv = paste0(directory_name_agecorr,'//','20241206_Human_CSF_ROC_analysis_Master_SpearmanAB40v42Correlation_AgeCorr_ModelStats.csv') 
write.csv(linear_AB42v40_correlation_stats, filename_csv); 

############################################# MASTER GRID PLOT: AGE-CORRECTED ############################################# 

combined_plot_grid_violin <- plot_grid(plot_storage_violin_agecorr[[1]], plot_storage_violin_agecorr[[2]], plot_storage_violin_agecorr[[3]], 
                                           plot_storage_violin_agecorr[[4]], plot_storage_violin_agecorr[[5]], plot_storage_violin_agecorr[[6]], 
                                           plot_storage_violin_agecorr[[7]], plot_storage_violin_agecorr[[8]], plot_storage_violin_agecorr[[9]], 
                                           plot_storage_violin_agecorr[[10]], plot_storage_violin_agecorr[[11]], plot_storage_violin_agecorr[[12]], 
                                           plot_storage_violin_agecorr[[13]], plot_storage_violin_agecorr[[14]], plot_storage_violin_agecorr[[15]], 
                                           plot_storage_violin_agecorr[[16]], plot_storage_violin_agecorr[[17]], plot_storage_violin_agecorr[[18]], 
                                           ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name_agecorr,'//Violin_Master_Grid_Figure_Age-Corrected.png') 
ggsave(filename_png, combined_plot_grid_violin, device = 'png', width = 36, height = 18, dpi = 600)
filename_eps = paste0(directory_name_agecorr,'//Violin_Master_Grid_Figure_Age-Corrected.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_violin), device = pdf, width = 36, height = 18)

############################################# MASTER GRID PLOT: AGE-CORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_agecorr[[1]], plot_storage_agecorr[[2]], plot_storage_agecorr[[3]], 
                                           plot_storage_agecorr[[4]], plot_storage_agecorr[[5]], plot_storage_agecorr[[6]], 
                                           plot_storage_agecorr[[7]], plot_storage_agecorr[[8]], plot_storage_agecorr[[9]], 
                                           plot_storage_agecorr[[10]], plot_storage_agecorr[[11]], plot_storage_agecorr[[12]], 
                                           plot_storage_agecorr[[13]], plot_storage_agecorr[[14]], plot_storage_agecorr[[15]], 
                                           plot_storage_agecorr[[16]], plot_storage_agecorr[[17]], plot_storage_agecorr[[18]], 
                                           ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name_agecorr,'//Scatter_Master_Grid_Figure_Age-Corrected.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name_agecorr,'//Scatter_Master_Grid_Figure_Age-Corrected.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs MoCA: AGE-CORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_MoCA_agecorr[[1]], plot_storage_MoCA_agecorr[[2]], plot_storage_MoCA_agecorr[[3]], 
                                        plot_storage_MoCA_agecorr[[4]], plot_storage_MoCA_agecorr[[5]], plot_storage_MoCA_agecorr[[6]], 
                                        plot_storage_MoCA_agecorr[[7]], plot_storage_MoCA_agecorr[[8]], plot_storage_MoCA_agecorr[[9]], 
                                        plot_storage_MoCA_agecorr[[10]], plot_storage_MoCA_agecorr[[11]], plot_storage_MoCA_agecorr[[12]], 
                                        plot_storage_MoCA_agecorr[[13]], plot_storage_MoCA_agecorr[[14]], plot_storage_MoCA_agecorr[[15]], 
                                        plot_storage_MoCA_agecorr[[16]], plot_storage_MoCA_agecorr[[17]], plot_storage_MoCA_agecorr[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name_agecorr,'//Scatter_MoCA_Master_Grid_Figure_Age-Corrected.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name_agecorr,'//Scatter_MoCA_Master_Grid_Figure_Age-Corrected.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs AB: AGE-CORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_BA_agecorr[[1]], plot_storage_BA_agecorr[[2]], plot_storage_BA_agecorr[[3]], 
                                        plot_storage_BA_agecorr[[4]], plot_storage_BA_agecorr[[5]], plot_storage_BA_agecorr[[6]], 
                                        plot_storage_BA_agecorr[[7]], plot_storage_BA_agecorr[[8]], plot_storage_BA_agecorr[[9]], 
                                        plot_storage_BA_agecorr[[10]], plot_storage_BA_agecorr[[11]], plot_storage_BA_agecorr[[12]], 
                                        plot_storage_BA_agecorr[[13]], plot_storage_BA_agecorr[[14]], plot_storage_BA_agecorr[[15]], 
                                        plot_storage_BA_agecorr[[16]], plot_storage_BA_agecorr[[17]], plot_storage_BA_agecorr[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name_agecorr,'//Scatter_AB_Master_Grid_Figure_Age-Corrected.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name_agecorr,'//Scatter_AB_Master_Grid_Figure_Age-Corrected.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs TOTAL TAU: AGE-CORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_totaltau_agecorr[[1]], plot_storage_totaltau_agecorr[[2]], plot_storage_totaltau_agecorr[[3]], 
                                        plot_storage_totaltau_agecorr[[4]], plot_storage_totaltau_agecorr[[5]], plot_storage_totaltau_agecorr[[6]], 
                                        plot_storage_totaltau_agecorr[[7]], plot_storage_totaltau_agecorr[[8]], plot_storage_totaltau_agecorr[[9]], 
                                        plot_storage_totaltau_agecorr[[10]], plot_storage_totaltau_agecorr[[11]], plot_storage_totaltau_agecorr[[12]], 
                                        plot_storage_totaltau_agecorr[[13]], plot_storage_totaltau_agecorr[[14]], plot_storage_totaltau_agecorr[[15]], 
                                        plot_storage_totaltau_agecorr[[16]], plot_storage_totaltau_agecorr[[17]], plot_storage_totaltau_agecorr[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name_agecorr,'//Scatter_totaltau_Master_Grid_Figure_Age-Corrected.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name_agecorr,'//Scatter_totaltau_Master_Grid_Figure_Age-Corrected.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs PHOSPHORYLATED TAU: AGE-CORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_phosphotau_agecorr[[1]], plot_storage_phosphotau_agecorr[[2]], plot_storage_phosphotau_agecorr[[3]], 
                                        plot_storage_phosphotau_agecorr[[4]], plot_storage_phosphotau_agecorr[[5]], plot_storage_phosphotau_agecorr[[6]], 
                                        plot_storage_phosphotau_agecorr[[7]], plot_storage_phosphotau_agecorr[[8]], plot_storage_phosphotau_agecorr[[9]], 
                                        plot_storage_phosphotau_agecorr[[10]], plot_storage_phosphotau_agecorr[[11]], plot_storage_phosphotau_agecorr[[12]], 
                                        plot_storage_phosphotau_agecorr[[13]], plot_storage_phosphotau_agecorr[[14]], plot_storage_phosphotau_agecorr[[15]], 
                                        plot_storage_phosphotau_agecorr[[16]], plot_storage_phosphotau_agecorr[[17]], plot_storage_phosphotau_agecorr[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name_agecorr,'//Scatter_phosphotau_Master_Grid_Figure_Age-Corrected.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name_agecorr,'//Scatter_phosphotau_Master_Grid_Figure_Age-Corrected.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)

############################################# MASTER GRID PLOT vs AB42 : AB40: AGE-CORRECTED ############################################# 

combined_plot_grid_scatter <- plot_grid(plot_storage_AB42v40_agecorr[[1]], plot_storage_AB42v40_agecorr[[2]], plot_storage_AB42v40_agecorr[[3]], 
                                        plot_storage_AB42v40_agecorr[[4]], plot_storage_AB42v40_agecorr[[5]], plot_storage_AB42v40_agecorr[[6]], 
                                        plot_storage_AB42v40_agecorr[[7]], plot_storage_AB42v40_agecorr[[8]], plot_storage_AB42v40_agecorr[[9]], 
                                        plot_storage_AB42v40_agecorr[[10]], plot_storage_AB42v40_agecorr[[11]], plot_storage_AB42v40_agecorr[[12]], 
                                        plot_storage_AB42v40_agecorr[[13]], plot_storage_AB42v40_agecorr[[14]], plot_storage_AB42v40_agecorr[[15]], 
                                        plot_storage_AB42v40_agecorr[[16]], plot_storage_AB42v40_agecorr[[17]], plot_storage_AB42v40_agecorr[[18]], 
                                        ncol = 6, align = "v", rel_widths = c(1, 1, 1, 1, 1, 1), scale = 0.85)
filename_png = paste0(directory_name_agecorr,'//Scatter_AB42v40_Master_Grid_Figure_Age-Corrected.png') 
ggsave(filename_png, combined_plot_grid_scatter, device = 'png', width = 28, height = 14, dpi = 600)
filename_eps = paste0(directory_name_agecorr,'//Scatter_AB42v40_Master_Grid_Figure_Age-Corrected.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_scatter), device = pdf, width = 28, height = 14)