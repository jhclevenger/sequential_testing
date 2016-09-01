# Load packages
library("BayesFactor")
library("dplyr")
library("tidyr")

##############################
# Simulations
##############################

# Number of simulations
num_sims = 5000
# Effect size (difference between groups)
effect_size = c(.2, .5, .8, 1)
# Initial sample size of each group
n = c(10)
# bf
bayes_factor = c(5)
# max n for each group
max_n = 100

# which variable are we testing? Also change hack index
variable = effect_size

bf = matrix(nrow = num_sims , ncol = 2*length(variable))
bf_subs = matrix(nrow = num_sims , ncol = 2*length(variable))
bf_cohens_d = matrix(nrow = num_sims , ncol = 2*length(variable))
bf_hedges_g = matrix(nrow = num_sims , ncol = 2*length(variable))

### Run tests ###
for (i_bf in 1:length(bayes_factor)) {
  for (i_n in 1:length(n)) {
    for (i_effect in 1:length(effect_size)) {
      for (i_sims in 1:num_sims) {
        
        # BF simulations
        # Generate initial data points
        group_01 <- rnorm(n[i_n])
        group_02 <- rnorm(n[i_n],mean=effect_size[i_effect])
        
        hack_index = (i_effect - 1) * 2 + 1
        no_hack_index = (i_effect - 1) * 2 + 2
        
        bf[i_sims,hack_index] = as.numeric(as.vector(ttestBF(group_01,group_02)))
        
        # If BF > bayes_factor or < 1 / bayes_factor, gather more data and re-test
        while (bf[i_sims,hack_index] < bayes_factor[i_bf] & bf[i_sims,hack_index] > 1/bayes_factor[i_bf]) {
          group_01 <- c(group_01,rnorm(1))
          group_02 <- c(group_02,rnorm(1,mean=effect_size[i_effect]))
          bf[i_sims,hack_index] <- as.numeric(as.vector(ttestBF(group_01,group_02)))
          if (length(group_01) >= max_n) {
            break
          }
        }
        
        # Sequential n and observed effect size
        bf_subs[i_sims,hack_index] = length(group_01)
        t_value = as.numeric(t.test(group_02,group_01)$statistic)
        bf_cohens_d[i_sims,hack_index] = t_value * sqrt(1/bf_subs[i_sims,hack_index] + 1/bf_subs[i_sims,hack_index])
        bf_hedges_g[i_sims,hack_index] = bf_cohens_d[i_sims,hack_index] * ( 1 - (3 / (4 * (bf_subs[i_sims,hack_index] + bf_subs[i_sims,hack_index]) - 9)))
        
        # Single sample
        group_01 <- rnorm(bf_subs[i_sims,hack_index])
        group_02 <- rnorm(bf_subs[i_sims,hack_index],mean=effect_size[i_effect])
        bf[i_sims,no_hack_index] = as.numeric(as.vector(ttestBF(group_01,group_02)))
        
        # Single sample n and observed effect size
        bf_subs[i_sims,no_hack_index] = length(group_01)
        t_value = as.numeric(t.test(group_02,group_01)$statistic)
        bf_cohens_d[i_sims,no_hack_index] = t_value * sqrt((1/(bf_subs[i_sims,no_hack_index]) + (1/bf_subs[i_sims,no_hack_index])))
        bf_hedges_g[i_sims,no_hack_index] = bf_cohens_d[i_sims,no_hack_index] * ( 1 - (3 / (4 * (bf_subs[i_sims,no_hack_index] + bf_subs[i_sims,no_hack_index]) - 9)))
      }
    }
  }
}

bf = as.data.frame(bf)
cohens_d= as.data.frame(bf_cohens_d)
hedges_g = as.data.frame(bf_hedges_g)
bf_subs = as.data.frame(bf_subs)
cols = c("two_seq", "two", "five_seq", "five", "eight_seq", "eight", "ten_seq", "ten")
colnames(bf) = cols
colnames(bf_cohens_d) = cols
colnames(bf_hedges_g) = cols
colnames(bf_subs) = cols

# Convert bf to long format
bf_long = gather(bf, condition, bf, two_seq:ten)

# Convert cohen's d to long format
cohens_d_long = gather(cohens_d, condition, cohens_d, two_seq:ten)

# Convert hedges g to long format
hedges_g_long = gather(hedges_g, condition, hedges_g, two_seq:ten)

# Summarize effect size
cohens_sum =
  cohens_d_long %>%
  group_by(condition) %>%
  summarise(
    num_sims = length(cohens_d),
    mean_d = mean(cohens_d),
    sd   = sd(cohens_d),
    se   = sd / sqrt(num_sims)
  ) %>%
  arrange(condition)

# Summarize effect size
hedges_sum =
  hedges_g_long %>%
  group_by(condition) %>%
  summarise(
    num_sims = length(hedges_g),
    mean_d = mean(hedges_g),
    sd   = sd(hedges_g),
    se   = sd / sqrt(num_sims)
  ) %>%
  arrange(condition)

# Summarize bf
bf_sum =
  bf_long %>%
  group_by(condition) %>%
  summarise(
    num_sims = length(bf),
    mean_bf = mean(bf),
    sd   = sd(bf),
    se   = sd / sqrt(num_sims)
  ) %>%
  arrange(condition)

# Summarize cor(n,d)
cor_sum =
  bf_long %>%
  group_by(condition) %>%
  summarise(
    num_sims = length(bf),
    mean_bf = mean(bf),
    sd   = sd(bf),
    se   = sd / sqrt(num_sims)
  ) %>%
  arrange(condition)

plot(bf_subs$five,cohens_d$five)
plot(bf_subs$five_seq,cohens_d$five_seq)

plot(bf_subs$two,cohens_d$two)
plot(bf_subs$two_seq,cohens_d$two_seq)


plot(bf_subs$eight_seq,cohens_d$eight_seq)
plot(bf_subs$eight,cohens_d$eight)

# plot
theme_set(theme_gray(base_size = 20))
theme = theme_update(legend.position="top")
#theme = theme_update(legend.position="top", panel.grid.minor.x=element_blank())
#theme = theme_update(panel.grid.minor.x=element_blank())

ggplot(bf_subs, aes(two_seq)) + geom_point(aes(y = cohens_d$two_seq), color='#F8766D', alpha = .4) + 
  xlab("N at time of stop") + ylab("Observed Cohen's d") +
  scale_x_continuous(limits = c(10,100)) +
  scale_y_continuous(breaks = c(-2,-1,0,1,2)) + ggtitle("Sequential Samples, d = .2")


