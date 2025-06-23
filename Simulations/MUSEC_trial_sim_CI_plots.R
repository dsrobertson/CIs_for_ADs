### Code for Section 4 (Figure 1): Plots of realised CIs for first N simulation
### replicates.

################################################################################


library(ggplot2)
library(dplyr)

source('MUSEC_trial.R')
load('sim_results.Rdata')


#### Unconditional plots

N = 5 # Number of simulation replicates to plot

lower = upper = Procedure = x = rep(NA, 9*N)

CI_df = data.frame(lower, upper, Procedure, x)


x = 1:N

CI_df$lower[x] = standard_CI[x,1]
CI_df$upper[x] = standard_CI[x,2]
CI_df$Procedure[x] = 'Wald test (standard)'

CI_df$lower[x+N] = exact_CI[x,1]
CI_df$upper[x+N] = exact_CI[x,2]
CI_df$Procedure[x+N] = 'Unconditional final'

CI_df$lower[x+2*N] = RCI[x,1]
CI_df$upper[x+2*N] = RCI[x,2]
CI_df$Procedure[x+2*N] = 'RCI'

CI_df$lower[x+3*N] = adjusted_asym_CI[x,1]
CI_df$upper[x+3*N] = adjusted_asym_CI[x,2]
CI_df$Procedure[x+3*N] = 'Adjusted asymptotic'

CI_df$lower[x+4*N] = bootstrap_uncond_CI[x,1]
CI_df$upper[x+4*N] = bootstrap_uncond_CI[x,2]
CI_df$Procedure[x+4*N] = 'Bootstrap'

CI_df$lower[x+5*N] = exact_cond_CI[x,1]
CI_df$upper[x+5*N] = exact_cond_CI[x,2]
CI_df$Procedure[x+5*N] = 'Conditional final'

CI_df$lower[x+6*N] = r_cond_CI[x,1]
CI_df$upper[x+6*N] = r_cond_CI[x,2]
CI_df$Procedure[x+6*N] = 'Restricted conditional'

CI_df$lower[x+7*N] = clik_CI[x,1]
CI_df$upper[x+7*N] = clik_CI[x,2]
CI_df$Procedure[x+7*N] = 'Conditional likelihood'

CI_df$lower[x+8*N] = plik_CI[x,1]
CI_df$upper[x+8*N] = plik_CI[x,2]
CI_df$Procedure[x+8*N] = 'Penalised likelihood'

CI_df$x = rep(1:N,9)

# Define the color palette
# Colorblind-friendly palette
cbbPalette = c("#000000", "#D55E00", "#009E73",
               "#56B4E9", "#A020F0", "#CC79A7", "#CC79A7",
               "#F0E442", "#999999")


CI_names_original = c('Wald test (standard)', 'Unconditional final',
                      'RCI', 'Adjusted asymptotic', 'Bootstrap',
                      'Conditional final', 'Restricted conditional',
                      'Conditional likelihood', 'Penalised likelihood')


letters_prefix <- LETTERS[1:length(CI_names_original)]

CI_names_with_letters <- paste0(letters_prefix, ". ", CI_names_original)

CI_df$Procedure = factor(CI_df$Procedure, levels = CI_names_original)

label_map <- setNames(letters_prefix, CI_names_original)


CI_df <- CI_df %>%
  mutate(plot_label = label_map[Procedure])

linetypes_vector_names = c(rep("solid", 6), "dotdash", rep("solid", 2))

names(linetypes_vector_names) <- CI_names_original


# Create the ggplot object
ggplot(CI_df, aes(x = as.factor(x), group = Procedure)) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper, colour = Procedure,
                    linetype = Procedure),
                position = position_dodge(width = 0.8), 
                width = 0.25) + 

  geom_text(aes(y = upper, label = plot_label, colour = Procedure),
            position = position_dodge(width = 0.8), 
            vjust = -1, 
            size = 3, 
            show.legend = FALSE) + 
  
  geom_hline(yintercept=MLE, linetype="dashed", color = "red", linewidth = 1) +
  
  scale_colour_manual(name = "CI Procedure",
                      values=cbbPalette,
                      labels = CI_names_with_letters) + 
  
  scale_linetype_manual(name = "CI Procedure",
                        values = linetypes_vector_names, 
                        labels = CI_names_with_letters) + 

  guides(colour = guide_legend(nrow = 3, override.aes = list(shape = NA)),
         linetype = guide_legend(nrow = 3, override.aes = list(shape = NA))) +
  
  coord_cartesian(ylim = c(-2, 0.5),
                  clip = "on") + 
  xlab('Simulation replicate') +
  ylab('Confidence Interval') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position = "bottom",
        legend.box = "vertical",
        plot.margin = margin(t = 15, r = 5.5, b = 5.5, l = 5.5, unit = "pt")) 
