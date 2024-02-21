library(patchwork)
library(afex)
library(tidyverse)
library(MASS)
library(emmeans)

random_matrix <- function(k, epsilon_target, max_iter=1000) {
  
  for (i6 in 1:max_iter) {
    cov_matrix <- trialr::rlkjcorr(1, 3)
    cov_matrix <- abs(cov_matrix)
    
    
    a <- mean(diag(cov_matrix)) # mean of diagonal
    b <- mean(cov_matrix) # grand mean
    c <- sum(cov_matrix^2) # sum of squares
    d <- sum(rowMeans(cov_matrix)^2) # sum of the mean of each row squared
    
    
    epsilon <- k^2 * (a - b)^2 / ((k - 1) * (c - 2 * k * d + k^2 * b^2))
    
    if (epsilon >= epsilon_lower && epsilon <= epsilon_upper) {
      return(list(cov_matrix = cov_matrix, epsilon = epsilon))
    }
  }
}


NSIM <- 10

# A list to store the dataframes of pairwise comparisons and anova
uni_list <- vector("list", NSIM)
multi_list <- vector("list", NSIM)
anova_list <- vector("list", NSIM)
n_tmp_list_uni <- vector("list", 3)
n_tmp_list_multi <- vector("list", 3)
n_tmp_list_anova <- vector("list", 3)
n_tmp_list_uni2 <- vector("list", 3)
n_tmp_list_multi2 <- vector("list", 3)
n_tmp_list_anova2 <- vector("list", 3)

n_loop <- c(5, 10, 20, 30, 60)
matrix_list <- vector("list", NSIM * length(n_loop) * 3)
matrix_counter <- 1

for (i5 in 1:3) {
  
  for (i4 in seq_along(n_loop)) {
    n <- n_loop[i4]
    #For Loop
    for (i in seq_len(NSIM)) {
      k <- 3
      if (i5 == 1){
        j <- runif(1, 0.2, 0.6)
        cov_matrix_100 <- matrix(0.5, k, k)
        diag(cov_matrix_100) <- 1
        cov_matrix <- cov_matrix_100
        epsilon <- 1
        matrix_list[[matrix_counter]] <- list(matrix=cov_matrix_100, epsilon=epsilon, n=n)
        matrix_counter <- matrix_counter + 1
      }
      
      
      if (i5 == 2){
        epsilon_upper <- 0.76
        epsilon_lower <- 0.74
        max_iter <- 10000 
        cov_matrix_75_all <- random_matrix(k, epsilon_target, max_iter)
        cov_matrix_75 <- cov_matrix_75_all$cov_matrix
        cov_matrix <- cov_matrix_75
        epsilon <- cov_matrix_75_all$epsilon
        matrix_list[[matrix_counter]] <- list(matrix=cov_matrix_75, epsilon=epsilon, n=n)
        matrix_counter <- matrix_counter + 1
      }
      
      if (i5 == 3){
        epsilon_upper <- 0.51
        epsilon_lower <- 0.49
        max_iter <- 10000
        cov_matrix_50_all <- random_matrix(k, epsilon_target, max_iter)
        cov_matrix_50 <- cov_matrix_50_all$cov_matrix
        cov_matrix <- cov_matrix_50
        epsilon <- cov_matrix_50_all$epsilon
        matrix_list[[matrix_counter]] <- list(matrix=cov_matrix_50, epsilon=epsilon, n=n)
        matrix_counter <- matrix_counter + 1
      }
      
      # Generate multivariate normal data
      simdat <- mvrnorm(n, mu = rep(0, k), Sigma = cov_matrix)
      
      # Convert to data frame
      dat <- as.data.frame(simdat)
      dat <- dat %>% rename(k1 = V1, k2 = V2, k3 = V3)
      
      #Convert to long dataframe
      dat1 <- dat %>% mutate(ID = 1:n)
      long <- dat1 %>% pivot_longer(cols = k1:k3, names_to = "Conditions", values_to = "Values")
      
      #Conduct ANOVA
      anova <- aov_ez(id = "ID", dv = "Values", data = long, within = "Conditions", include_aov = TRUE)
      summary(anova)
      #Extract ANOVA, convert to dataframe and store in list
      anova_tab <- as.data.frame(anova$anova_table)
      anova_list [[i]] <- anova_tab
      anova_list [[i]] <- anova_list[[i]] %>% mutate(iteration = i) #Create new column 'iteration' and assigned value of i
      anova_list [[i]] <- anova_list[[i]] %>% mutate(n = n_loop[i4])
      #Conduct pairwise comparisons
      multi <- emmeans(anova, "Conditions", model = "multivariate")
      
      uni <- emmeans(anova, "Conditions", model = "univariate")
      
      #Inner loop for multivariate corrections
      
      corrections <- c("bonferroni", "mvt", "holm", "hochberg")
      correction_list_multi <- vector("list", length(corrections))
      for (i2 in seq_along(corrections)) {
        correction_list_multi[[i2]] <- as.data.frame(summary(pairs(multi), adjust = corrections[i2])) %>% 
          mutate(correction = corrections[i2])
      }
      multi_list[[i]] <- do.call(rbind, correction_list_multi) 
      
      #Inner loop for univariate corrections
      
      correction_list_uni <- vector("list", length(corrections))
      for (i3 in seq_along(corrections)) {
        correction_list_uni[[i3]] <- as.data.frame(summary(pairs(uni), adjust = corrections[i3])) %>%
          mutate(correction = corrections[i3])
      }
      
      uni_list[[i]] <- do.call(rbind, correction_list_uni)
      
      # Iteration identifier to dataframe
      multi_list[[i]] <- multi_list[[i]] %>% mutate(iteration = i)
      uni_list[[i]] <- uni_list[[i]] %>% mutate(iteration = i)
      multi_list[[i]] <- multi_list[[i]] %>% mutate(n = n_loop[i4])
      uni_list[[i]] <- uni_list[[i]] %>% mutate(n = n_loop[i4])
      multi_list[[i]] <- multi_list[[i]] %>% mutate(epsilon = i5)
      uni_list[[i]] <- uni_list[[i]] %>% mutate(epsilon = i5)
      
      
    }
    n_tmp_list_uni[[i4]] <- do.call(rbind, uni_list)
    n_tmp_list_multi[[i4]] <- do.call(rbind, multi_list)
    n_tmp_list_anova[[i4]] <- do.call(rbind, anova_list)
    
  }

  n_tmp_list_uni2[[i5]] <- do.call(rbind, n_tmp_list_uni)
  n_tmp_list_multi2[[i5]] <- do.call(rbind, n_tmp_list_multi)
  n_tmp_list_anova2[[i5]] <- do.call(rbind, n_tmp_list_anova)
  
}

all_anova <- do.call(rbind, n_tmp_list_anova2)
all_uni <- do.call(rbind, n_tmp_list_uni2)
all_multi <- do.call(rbind, n_tmp_list_multi2)


#save(all_multi, all_uni, all_anova, matrix_list, file = "results3.rda")

#Load in one of the RDA files first so the list can calculate length of dataframes

file_list <- c("results1.rda", "results2.rda", "results3.rda", "results4.rda", "results5.rda", 
               "results6.rda", "results7.rda", "results8.rda", "results9.rda", "results10.rda")
for(i in seq_along(file_list)) {
  load(file_list[i])
  print(i)
  print(nrow(all_multi))
  if (i == 1) {
    combined_anova <- all_anova
    combined_multi <- all_multi
    combined_uni <- all_uni
  } else {
    all_anova$iteration <- all_anova$iteration + (i-1)*1000
    all_multi$iteration <- all_multi$iteration + (i-1)*1000
    all_uni$iteration <- all_uni$iteration + (i-1)*1000
    combined_anova <- rbind(combined_anova, all_anova)
    combined_multi <- rbind(combined_multi, all_multi)
    combined_uni <- rbind(combined_uni, all_uni)
  }
}
#save(combined_anova, combined_multi, combined_uni, file = "all_results_ws.rda")

all_combined <- bind_rows(
  mutate(combined_multi, method = "multi"),
  mutate(combined_uni, method = "uni")
)

all_combined <- all_combined %>%
  mutate(epsilon = ifelse(epsilon == 2, 0.75, 
                     ifelse(epsilon == 3, 0.5, epsilon)))

all_combined_sum <- all_combined %>%
  group_by(method, n, correction, epsilon, iteration) %>%
  summarise(p_signif = any(p.value < .05),
            num_signif = sum(p.value < .05)) #number of significant results in each group
all_combined_sum

combined_out <- all_combined_sum %>%
  summarise(n_signif = sum(p_signif),  #number of groups with at least one significant p-value
            n_signif_l1 = sum(num_signif > 1),#number of groups with more than one significant p-value
            total_n = n(),
            prop_signif = mean(p_signif)) %>%
  mutate(lower = binom::binom.agresti.coull(x = n_signif, n = total_n)$lower,
         upper = binom::binom.agresti.coull(x = n_signif, n = total_n)$upper)
combined_out

# Combined Plot
ws_plot <- combined_out %>%
  ggplot(aes(x = n, group = method, colour = method)) +
  geom_hline(yintercept = 0.05) +
  geom_pointrange(aes(y = prop_signif, ymin = lower, ymax = upper),
                 position = position_dodge(2)) +
  facet_grid(rows = vars(epsilon), cols = vars(correction)) + 
  ylim(0.02, 0.1) +
  ggtitle("Within-Subjects Design") +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
# Individual Plot
# ws_plot <- combined_out %>%
#   ggplot(aes(x = n, group = method, colour = method)) +
#   geom_hline(yintercept = 0.05) +
#   geom_pointrange(aes(y = prop_signif, ymin = lower, ymax = upper),
#                   position = position_dodge(2)) +
#   facet_grid(rows = vars(epsilon), cols = vars(correction)) +
#   ylim(0.02, 0.1) +
#   ggtitle("Within-Subjects Design") +
#   labs(x = "Sample Size (n)", y = "Type I Error Rate")




# Extension: Only significant ANOVA

combined_anova$epsilon <- rep(rep(1:3, each = 5000), times = 10)


combined_anova1 <- combined_anova %>% dplyr::select(6:9)
comb_anova_multi <- inner_join(combined_anova1, combined_multi, by = c("n", "iteration", "epsilon"))
comb_anova_uni <- inner_join(combined_anova1, combined_uni, by = c("iteration", "n", "epsilon"))
multi_only_sig <- comb_anova_multi %>% filter(`Pr(>F)` < 0.05)
uni_only_sig <- comb_anova_uni %>% filter(`Pr(>F)` < 0.05)

all_combined_only_sig <- bind_rows(
  mutate(multi_only_sig, method = "multi"),
  mutate(uni_only_sig, method = "uni")
)

all_combined_only_sig <- all_combined_only_sig %>%
  mutate(epsilon = ifelse(epsilon == 2, 0.75, 
                          ifelse(epsilon == 3, 0.5, epsilon)))

all_combined_only_sig_sum <- all_combined_only_sig %>%
  group_by(method, n, correction, epsilon, iteration) %>%
  summarise(p_signif = any(p.value < .05),
            num_signif = sum(p.value < .05)) #number of significant results in each group
all_combined_only_sig_sum

combined_out_only_sig <- all_combined_only_sig_sum %>%
  summarise(n_signif = sum(p_signif),  #number of groups with at least one significant p-value
            n_signif_l1 = sum(num_signif > 1),#number of groups with more than one significant p-value
            total_n = 10000,
            prop_signif = sum(p_signif)/10000) %>%
  mutate(lower = binom::binom.agresti.coull(x = n_signif, n = total_n)$lower,
         upper = binom::binom.agresti.coull(x = n_signif, n = total_n)$upper)
combined_out_only_sig

ws_plot_only_sig <- combined_out_only_sig %>%
  ggplot(aes(x = n, group = method, colour = method)) +
  geom_hline(yintercept = 0.05) +
  geom_pointrange(aes(y = prop_signif, ymin = lower, ymax = upper),
                  position = position_dodge(2)) +
  facet_grid(rows = vars(epsilon), cols = vars(correction)) + 
  ylim(0, 0.1) +
  ggtitle("Within-Subjects Design") +
  ylab("Conditional Type I Error Rate") +
  theme(axis.title.x = element_blank(), legend.position = "none")

# Individual Plot
# ws_plot_only_sig <- combined_out_only_sig %>%
#   ggplot(aes(x = n, group = method, colour = method)) +
#   geom_hline(yintercept = 0.05) +
#   geom_pointrange(aes(y = prop_signif, ymin = lower, ymax = upper),
#                   position = position_dodge(2)) +
#   facet_grid(rows = vars(epsilon), cols = vars(correction)) +
#   ylim(0, 0.1) +
#   ggtitle("Within-Subjects Design") +
#   labs(y = "Conditional Type I Error Rate", x = "Sample Size (n)")

# ANOVA TYPE I ERROR

anova_sum <- combined_anova %>% mutate(epsilon = ifelse(epsilon == 2, 0.75, 
                                                         ifelse(epsilon == 3, 0.5, epsilon))) %>% 
  group_by(n, epsilon, iteration) %>%
  summarise(p_signif = any(`Pr(>F)` < .05),
            num_signif = sum(`Pr(>F)` < .05)) #number of significant results in each group

combined_out_anova <- anova_sum %>%
  summarise(n_signif = sum(p_signif),  #number of groups with at least one significant p-value
            n_signif_l1 = sum(num_signif > 1),#number of groups with more than one significant p-value
            total_n = n(),
            prop_signif = mean(p_signif)) %>%
  mutate(lower = binom::binom.agresti.coull(x = n_signif, n = total_n)$lower,
         upper = binom::binom.agresti.coull(x = n_signif, n = total_n)$upper)

ws_plot_anova <- combined_out_anova %>%
  ggplot(aes(x = n)) +
  geom_hline(yintercept = 0.05) +
  geom_pointrange(aes(y = prop_signif, ymin = lower, ymax = upper),
                  position = position_dodge(2)) +
  facet_grid(rows = vars(epsilon)) + 
  ylim(0.02, 0.1) +
  ggtitle("Within-Subjects Design") +
  ylab("ANOVA Type I Error Rate") +
  theme(axis.title.x = element_blank(), plot.title = element_text(size = 10))



# Follow-up procedure combined plot 
combined_plot <- equal_plot + unequal_plot + ws_plot
print(combined_plot)
#1800x700

# Follow-up procedure only significant anova combined plot
combined_plot_only_sig <- ws_plot_only_sig + bs_same_plot_only_sig + bs_diff_plot_only_sig
print(combined_plot_only_sig)
#1800x800

# ANOVA combined plot
combined_plot_anova <- ws_plot_anova + bs_same_plot_anova + bs_diff_plot_anova
print(combined_plot_anova)
#1000x700
