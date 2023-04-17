# Creates Table S4, which: ####
# reports the results of ANOVAs tests for whether diversity differs among tree and
# needle age classes more than expected by chance.

# Creates Fig. 4, which: ####
# (A) shows how OTU richness varies among needle age classes,
# (B) models the relationship between estimated richness and closure 
# (C) how the Shannon diversity index differs among exposure groups, and
# (D) how the Shannon diversity index differs among trees.

# Creates Fig. S3, which: ####
# models the relationship between the estimated Shannon index and closure

# Load packages ####
library(patchwork)
library(cowplot)
library(ggthemes)
library(gt)
library(iNEXT)
library(phyloseq)
library(nlme)
library(car)
library(tidyverse)
library(magrittr)
library(here)
# Consult "code/functions.R" for descriptions of function applications, inputs, and outputs
source(here('code', 'functions.R')
       )

# Set out-paths ####
figure.out <- here('output', 'figures')
dir.create(figure.out, recursive = T)
table.out <- here('output', 'tables')
dir.create(table.out, recursive = T)

# Read in files ####
perf.n <- readRDS(here('data', 'compile', 'perf.n.clean.rds')
                  )

# Estimate richness and diversity ####
perf.n.div <- inext.div(perf.n$counts$all, cores = 6, level = 1000)

rich <- perf.n.div %>% filter(metric == 'q0')

shannon <- perf.n.div %>% filter(metric == 'q1') %>% mutate(estimate = log2(estimate))

# Assess analysis assumptions ####
# Levene test for homogeneous variance
# Richness
leveneTest(estimate ~ age*tree, data = rich)

# Shannon's index
leveneTest(estimate ~ age*tree, data = shannon)

# Shapiro-Wilk test
# Richness
rich %>% group_by(age) %>% rstatix::shapiro_test(estimate)
rich %>% group_by(tree) %>% rstatix::shapiro_test(estimate)
# Shannon's index
shannon %>% group_by(age) %>% rstatix::shapiro_test(estimate)
shannon %>% group_by(tree) %>% rstatix::shapiro_test(estimate)

# Plots examining homoscedasticity and normality
# Richness
rich.aov <- aov(estimate ~ age*tree, data = rich)
rich.aov %>% plot(ask = F)

#Shannon's index
shannon.aov <- aov(estimate ~ age*tree, data = shannon)
shannon.aov %>% plot(ask = F)

# Richness
# Needle age class
rich.age <- ggplot(rich,
                   aes(x = age, y = estimate)
                   ) +
  geom_violin(fill = 'white', color = 'black', draw_quantiles = c(0.25, 0.5, 0.75)
              ) +
  ggbeeswarm::geom_quasirandom(aes(fill = tree),
                               width = 0.4, shape = 21, size = 3) +
  scale_fill_colorblind() +
  xlab('') +
  ylab('Estimated richness\n') +
  labs(fill = 'Tree') +
  theme_cowplot() +
  theme(axis.text.x = element_text(face = 'bold', size = 7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold', size = 7),
        legend.text = element_text(size = 7)
        )
# Tree
rich.tree <- ggplot(rich,
                    aes(x = tree, y = estimate)
                    ) +
  geom_violin(fill = 'white', color = 'black', draw_quantiles = c(0.25, 0.5, 0.75)
              ) +
  ggbeeswarm::geom_quasirandom(aes(fill = age),
                               width = 0.4, shape = 21, size = 3) +
  scale_fill_colorblind() +
  xlab('') +
  ylab('Estimated richness\n') +
  labs(fill = 'Age') +
  theme_cowplot() +
  theme(axis.text.x = element_text(face = 'bold', size = 7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold', size = 7),
        legend.text = element_text(size = 7)
  )

rich.age / rich.tree

# Shannon's index
# Needle age class
shannon.age <- ggplot(shannon,
                      aes(x = age, y = estimate)
                      ) +
  geom_violin(fill = 'white', color = 'black', draw_quantiles = c(0.25, 0.5, 0.75),
              width = 0.4) +
  ggbeeswarm::geom_quasirandom(aes(fill = tree),
                               width = 0.2, shape = 21, size = 2) +
  scale_fill_colorblind() +
  xlab('') +
  ylab('Estimated Shannon index\n') +
  labs(fill = 'Tree') +
  theme_cowplot() +
  theme(axis.text.x = element_text(face = 'bold', size = 7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold', size = 7),
        legend.text = element_text(size = 7)
        )
# Tree
 shannon.tree <- ggplot(shannon,
                        aes(x = tree, y = estimate)
                        ) +
   geom_violin(fill = 'white', color = 'black', draw_quantiles = c(0.25, 0.5, 0.75)
               ) +
   ggbeeswarm::geom_quasirandom(aes(fill = age),
                                width = 0.4, shape = 21, size = 2) +
   scale_fill_colorblind() +
   xlab('') +
   ylab('Estimated Shannon index\n') +
   labs(fill = 'Age') +
   theme_cowplot() +
   theme(axis.text.x = element_text(face = 'bold', size = 7),
         axis.ticks.x = element_blank(),
         axis.title.y = element_text(face = 'bold', size = 7),
         axis.text.y = element_text(size = 7),
         legend.title = element_text(face = 'bold', size = 7),
         legend.text = element_text(size = 7)
   )

shannon.age / shannon.tree

# ANOVA ####
# Tree and age ####
# Homoscedasticity and semi-normality suggest that ANOVA can be applied
rich.f.test <- Anova(rich.aov,
              white.adjust = F) %>% data.frame() %>%
  rownames_to_column(var = 'Term') %>%
  select(Term, Df, `F` = F.value, P = Pr..F.) %>% bind_rows(data.frame(Term = 'Total', Df = sum(.$Df)
                                                                       )
                                                            )
rich.f.test %<>% mutate(P = round(P, digits = 3),
                           `F` = round(`F`, digits = 3)
                        )
rich.f.test$P[rich.f.test$P < 0.001] <- '< 0.001'
rich.f.test$metric <- 'Richness'

# Levene's test suggests the use of a White adjustment for this ANOVA
shannon.f.test <- Anova(shannon.aov,
                     white.adjust = T) %>% data.frame() %>%
  rownames_to_column(var = 'Term') %>%
  select(Term, Df, `F`, P = Pr..F.) %>% bind_rows(data.frame(Term = 'Total', Df = sum(.$Df)
                                                                       )
                                                            )
shannon.f.test %<>% mutate(P = round(P, digits = 3),
                           `F` = round(`F`, digits = 3)
                           )
shannon.f.test$P[shannon.f.test$P < 0.001] <- '< 0.001'
shannon.f.test$metric <- 'Shannon index'

div.f.test <- rbind(rich.f.test, shannon.f.test)

rbind(div.f.test) %>%
  gt(groupname_col = 'metric') %>%
  tab_header(title = 'Table S4. ANOVA testing variation in alpha diversity among trees and needle age classes.') %>%
  tab_options(
    table.font.size = px(15),
    heading.title.font.size = px(18),
    heading.align = 'left',
    heading.border.bottom.color = 'white',
    table.border.top.color = 'white',
    table.border.bottom.color = 'white',
    table_body.border.top.color = 'white',
    table_body.border.bottom.color = 'white',
    column_labels.font.weight = 'bold',
    column_labels.border.bottom.color = 'white',
    row_group.font.weight = 'bold',
    row_group.border.top.color = 'white',
    row_group.border.bottom.color = 'white',
    table_body.hlines.color = 'white',
    source_notes.font.size = px(12)
  ) %>%
  tab_style(
    style = cell_text(style = 'italic'),
    locations = cells_body(
      columns = 'Term',
      rows = grepl('^[[:lower:]]', Term)
      )
    ) %>%
  tab_style(
    style = cell_text(color = 'white'),
    locations = cells_column_labels(columns = 'Term')
  ) %>%
  tab_style(
    style = cell_borders(
      sides = 'top',
      color = 'black',
      weight = px(1.5),
      style = 'solid'),
    locations = cells_body(
      columns = everything(),
      rows = Term == 'Total')
  ) %>%
  fmt_missing(columns = everything(), missing_text = '') %>%
  tab_source_note(source_note = 'Estimates of OTU richness and the Shannon index of diversity were interpolated or extrapolated to a depth of 1000 reads.') %>%
  opt_table_font(font = google_font('Crimson Text')
                 ) %>%
  cols_align(align = 'center') %>%
  cols_width(everything()~px(100)
             ) %>%
  gtsave(filename = here(table.out, 'table.s4.png')
         )

# Tukey ####
# Age ####
tukey.age <- perf.n.div %>%
  select(sample, age, tree, group, metric, estimate) %>%
  group_by(metric) %>%
  group_map(~ tukey(.x, form = 'age'), .keep = T) %>%
  bind_rows()
tukey.age$estimate <- ifelse(tukey.age$metric == 'q0', 22, 8.5)

# Tree ####
tukey.tree <- perf.n.div %>%
  select(sample, age, tree, group, metric, estimate) %>%
  group_by(metric) %>%
  group_map(~ tukey(.x, form = 'tree'), .keep = T) %>%
  bind_rows()
tukey.tree$estimate <- ifelse(tukey.tree$metric == 'q0', 21.5, 3.5)

# Bootstrapped means and confidence intervals ####
# Age ####
boot.age <- perf.n.div %>%
  select(age, tree, group, metric, estimate) %>%
  group_by(metric, age) %>%
  group_map(~ mean.boot(.x, .y, var = 'age')
            ) %>%
  bind_rows() %>%
  group_by(metric, factor) %>%
  summarize(mean = mean(estimate),
            lci = coxed::bca(estimate)[1],
            uci = coxed::bca(estimate)[2]
            ) %>%
  select(age = factor, everything()
                )

# Tree ####
boot.tree <- perf.n.div %>%
  select(age, tree, group, metric, estimate) %>%
  group_by(metric, tree) %>%
  group_map(~ mean.boot(.x, .y, var = 'tree')
            ) %>%
  bind_rows %>%
  group_by(metric, factor) %>%
  summarize(mean = mean(estimate),
            lci = coxed::bca(estimate)[1],
            uci = coxed::bca(estimate)[2]) %>%
  select(tree = factor, everything()
  )

# Age and tree plots ####
# Richness, age, pairwise ####
age.rich.pair <- ggplot(rich,
                      aes(x = age, y = estimate)
                      ) +
  geom_violin(fill = 'white',
              color = 'black'
              ) +
  geom_text(aes(label = Letters),
            filter(tukey.age, metric == 'q0'),
            fontface = 'bold',
            size = 2.5
            ) +
  geom_pointrange(aes(y = mean,
                      ymin = lci, ymax = uci),
                  filter(boot.age, metric == 'q0'),
                  size = 0.5,
                  stroke = 0.5,
                  shape = 3,
                  color = 'black'
                  ) +
  xlab('') +
  ylab('Estimated richness\n') +
  labs(fill = 'Age') +
  scale_fill_colorblind() +
  theme_cowplot() +
  theme(axis.text.x = element_text(face = 'bold', size = 7),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7)
        )
age.rich.pair

# Shannon, tree, pairwise ####
tree.shannon.pair <- ggplot(shannon,
                       aes(x = tree, y = estimate)
                       ) +
  geom_violin(aes(fill = tree),
              color = 'white'
              ) +
  geom_text(aes(label = Letters, y = estimate),
            filter(tukey.tree, metric == 'q1'),
            fontface = 'bold',
            size = 2.5
            ) +
  geom_pointrange(aes(y = mean,
                      ymin = lci, ymax = uci),
                  filter(boot.tree, metric == 'q1'),
                  size = 0.5,
                  stroke = 0.5,
                  shape = 3,
                  color = 'white') +
  xlab('') +
  ylab('Estimated Shannon index\n') +
  labs(fill = 'Tree') +
  scale_fill_colorblind() +
  theme_cowplot() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(face = 'bold', size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.75, 'line')
        )
tree.shannon.pair

# Linear mixed-effect modeling ####
# Richness vs crown variables ####
# Full
rich.lme.full <- lme(estimate ~ age + height + closure,
                     random = ~ 1|tree,
                     data = rich,
                     method = 'ML')

# Alternate
rich.lme.alt <- lme(estimate ~ age + height,
                    random = ~ 1|tree,
                    data = rich,
                    method = 'ML')

# Reduced
rich.lme.red <- lme(estimate ~ age + closure,
                    random = ~ 1|tree,
                    data = rich,
                    method = 'ML')

# Random intercept only
rich.lme.base <- lme(estimate ~ age,
                     random = ~1|tree,
                     data = rich,
                     method = 'ML')
# Check assumptions
# Homoscedasticity
rich.lme.red %>% plot()
# Normally distributed residuals
rich.lme.red %>% residuals() %>% qqnorm()
rich.lme.red %>% residuals() %>% qqline()
# Normally distributed random effects
# One tree appears to be somewhat divergent
ranef(rich.lme.red)$`(Intercept)` %>% qqnorm()
ranef(rich.lme.red)$`(Intercept)` %>% qqline()

# Overview
summary(rich.lme.red)
# Type II ANOVA
car::Anova(rich.lme.red)
# Compare full, reduced, and alternate models
AICcmodavg::aictab(list(rich.lme.full, rich.lme.red, rich.lme.alt, rich.lme.base)
                   )
AICcmodavg::bictab(list(rich.lme.full, rich.lme.red, rich.lme.alt, rich.lme.base)
                   )
# Compare reduced and base model to obtain overall model significance
rich.lme.p <- anova(rich.lme.red, rich.lme.base)$`p-value`[[2]] %>% round(digits = 3)
rich.lme.p[rich.lme.p < 0.001] <- '< 0.001'
# Obtain an R2 for the model
rich.lme.r2 <-  piecewiseSEM::rsquared(rich.lme.red)$Marginal %>% round(digits = 3)

# Diversity vs crown variables ####
# Full
shannon.lme.full <- lme(estimate ~ height + closure,
                        random = ~1|tree,
                        data = shannon,
                        method = 'ML')

# Alternate
shannon.lme.alt <- lme(estimate ~ height,
                       random = ~1|tree,
                       data = shannon,
                       method = 'ML')

# Reduced
shannon.lme.red <- lme(estimate ~ closure,
                       random = ~1|tree,
                       data = shannon,
                       method = 'ML')
# Random intercept only
shannon.lme.base <- lme(estimate ~ 1,
                        random = ~1|tree,
                        data = shannon,
                        method = 'ML')
# Check residuals
# Homoscedasticity
# Somewhat funnel-shaped residuals observed vs fitted
shannon.lme.red %>% plot()
# Normally distributed residuals
shannon.lme.red %>% residuals() %>% qqnorm()
shannon.lme.red %>% residuals() %>% qqline()
# Normally distributed random effects
ranef(shannon.lme.red)$`(Intercept)` %>% qqnorm()
ranef(shannon.lme.red)$`(Intercept)` %>% qqline()

# Overview
summary(shannon.lme.red)
# Compare full, reduced, and alternate models
AICcmodavg::aictab(list(shannon.lme.full, shannon.lme.red, shannon.lme.alt, shannon.lme.base)
)
AICcmodavg::bictab(list(shannon.lme.full, shannon.lme.red, shannon.lme.alt, shannon.lme.base)
)
# Compare reduced and base model to obtain overall model significance
shannon.lme.p <- anova(shannon.lme.red, shannon.lme.base)$`p-value`[[2]] %>% round(digits = 3)
shannon.lme.p[shannon.lme.p < 0.001] <- '< 0.001'
# Obtain an R2 for the model
shannon.lme.r2 <- piecewiseSEM::rsquared(shannon.lme.red)$Marginal %>% round(digits = 3)

# Linear mixed-model plots ####
# Richness vs closure plot ####
rich$lme.fixed <- rich.lme.red$fitted %>% data.frame() %>% .$fixed
rich.lab <- substitute(P~p*','~R^2~"="~r2, list(p = rich.lme.p, r2 = rich.lme.r2)
                       ) %>% as.expression()

age.cols <- c('#0172e2', '#fd8a1a', '#dcbaf0', '#b3003b')

closure.rich.lme <- ggplot(rich, aes(x = closure, y = estimate)
                           ) +
  geom_point(size = 1.5, shape = 21, aes(fill = tree), show.legend = F) +
  geom_line(aes(y = lme.fixed, color = age), size = 1) +
  xlab('\nCrown closure index') +
  ylab('Estimated richness\n') +
  labs(color = 'Age', fill = 'Tree') +
  annotate('text', x = 10, y = 20, label = rich.lab, parse = T, size = 2.5) +
  scale_fill_colorblind() +
  scale_color_manual(values = age.cols) +
  coord_equal() +
  theme_cowplot() +
  theme(legend.title = element_text(size = 7, face = 'bold'),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7)
        )
closure.rich.lme

# Figure output
(age.rich.pair + closure.rich.lme) / (tree.shannon.pair + plot_spacer()
                                      )  +
  plot_annotation(tag_levels = list(c('A', 'B', 'C', 'D')
                                    )
                  ) &
  theme(plot.tag = element_text(size = 10, face = 'bold')
        )
ggsave(here(figure.out, 'fig.4.tiff'), units = 'mm', width = 190, height = 120,
       dpi = 300, compression = 'lzw')

# Diversity vs closure plot ####
shannon$lme.fixed <- shannon.lme.red$fitted %>% data.frame() %>% .$fixed
shannon.lab <- substitute(P~p*','~R^2~"="~r2, list(p = shannon.lme.p, r2 = shannon.lme.r2)
                          ) %>% as.expression()

closure.shannon.lme <- ggplot(shannon, aes(x = closure, y = estimate)
                                   ) +
  geom_point(size = 1.5, shape = 21, aes(fill = tree)
             ) +
  geom_line(aes(y = lme.fixed), size = 1.5) +
  xlab('\nCrown closure index') +
  ylab('Estimated Shannon index\n') +
  labs(color = 'Age', fill = 'Tree') +
  annotate('text', x = 10, y = 3, label = shannon.lab, parse = T, size = 2.5) +
  scale_fill_colorblind() +
  theme_cowplot() +
  theme(legend.title = element_text(size = 7, face = 'bold'),
        legend.text = element_text(size = 7),
        axis.title.x = element_text(face = 'bold', size = 7),
        axis.text.x = element_text(size = 7),
        axis.title.y = element_text(face = 'bold', size = 7),
        axis.text.y = element_text(size = 7),
        plot.background = element_rect(fill = 'white', color = 'white')
        )

ggsave(here(figure.out, 'fig.s3.png'), units = 'mm', width = 140,
       dpi = 300)

# Get session info ####
session.path <- here('output', 'sessions')
dir.create(session.path, recursive = T)
sessionInfo() %>% saveRDS(here(session.path, 'div.sesh.rds')
                          )
