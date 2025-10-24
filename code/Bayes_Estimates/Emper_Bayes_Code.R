####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####                     IMPORT FILES                     ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

#[] Libraries ####
# Load all libraries to be used for this analysis
library(multilevel)
library(nlme)
library(lme4)
library(lattice)
library(dplyr)

#[] Functions ####
# Copy to clipboard
copy.table <- function(obj, size = 4096) {
  clip <- paste('clipboard-', size, sep = '')
  f <- file(description = clip, open = 'w')
  write.table(obj, f, row.names = FALSE, sep = '\t')
  close(f)  
}

####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####                   Background Info                    ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# Extracting Bayesian Estimates from Random Coefficient Growth Models for 
# Momentum Analysis

# --- Introduction --- #

# Momentum in longitudinal data refers to the systematic improvement or decline 
# in a variable over time. In other words, it's the trajectory or rate of change
# (slope) of a repeated measure for each individual. Prior research (e.g., Chen 
# et al., 2011) demonstrated that such momentum can predict outcomes above and
# beyond what static or average levels of the variable explain. For example, an 
# employee whose job satisfaction is steadily improving may have different turnover
# intentions than one with the same average satisfaction but a declining trend.


# --- Goal --- # 
# This guide will show how to extract Bayesian estimates (empirical Bayes or 
# BLUPs) of individual intercepts and slopes from a random coefficient growth 
# model and use them to model momentum effects. The Bayesian estimates of random
# effects are essentially the model-based adjusted values for each individual's 
# intercept and slope, combining individual data and overall trends. These provide 
# a more reliable estimate of an individual's trajectory than using raw values alone, 
# due to shrinkage (extreme individual estimates are pulled toward the grand mean
# based on data reliability).


# --- Procedure Overview --- # 
# We will go through the following steps to demonstrate momentum modeling:
#   Simulate Longitudinal Data: Create a dataset of individuals measured across 
#       time (e.g., an "Attitude" measured over 5 time points), and simulate an 
#       outcome variable influenced by the trajectory.
#   Fit a Growth Model: Use a linear mixed-effects model to capture each individual's 
#       intercept (initial level) and slope (rate of change) as random effects.
#   Extract Bayesian Estimates: Obtain each individual's estimated intercept and
#       slope (empirical Bayes estimates from the model).
#   Momentum in Downstream Model: Use the extracted slope (momentum) in a subsequent 
#       analysis to predict an outcome, controlling for the intercept or average 
#       level, and compare approaches.
#
# Throughout, we will use clear R code with annotations so that researchers can 
# adapt it to their own data. We use the lme4 package for mixed modeling (for 
# simplicity), but other packages (e.g., nlme) could also be used.


####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####             Simulating Longitudinal Data             ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# First, let's simulate a longitudinal dataset. Suppose we have 100 individuals, 
# each with an "Attitude" score measured at 5 time points (say, Time 0 through 4). 
# We'll also simulate a single outcome per individual (e.g., a future outcome 
# influenced by their attitude momentum). For realism, individuals will vary in 
# their starting Attitude (random intercepts) and in their change over time (random 
# slopes). We set up the simulation so that the outcome is related to both the 
# individual's intercept and slope:
#   Intercept effect: Individuals with higher initial Attitude tend to have a lower 
#       outcome (e.g., if outcome is a negative behavior like turnover intention, 
#       higher satisfaction -> lower intention).
#   Slope effect (Momentum): Individuals whose Attitude is increasing (positive 
#       slope) will have a lower outcome (e.g., improving attitude -> less turnover 
#       intention), whereas those with declining attitudes have higher outcome.


# --- R code to simulate the data --- #
set.seed(123)  # for reproducibility

# Parameters for simulation
N <- 100                         # number of individuals
time_points <- 0:4               # five measurements (coded 0 to 4 for convenience)
n_obs <- length(time_points) * N # total observations

# Random true intercepts and slopes for each individual
true_intercepts <- rnorm(N, mean = 5, sd = 1)    # baseline Attitude around 5 (e.g., on a 1-7 scale)
true_slopes    <- rnorm(N, mean = 0, sd = 0.3)   # on average no change, some up or down

# Residual variance for Attitude measurements
resid_sd <- 0.5

# Construct the long-format data frame for repeated measures
PersonID <- rep(1:N, each = length(time_points))
Time     <- rep(time_points, times = N)

# True attitude without noise = intercept + slope * time
true_attitude <- true_intercepts[PersonID] + true_slopes[PersonID] * Time

# Observed Attitude with noise
Attitude <- rnorm(n_obs, mean = true_attitude, sd = resid_sd)

data_long <- data.frame(PersonID, Time, Attitude)
head(data_long)   # preview the first few rows


# In this simulated data:
#   PersonID identifies the individual.
#   Time is the measurement occasion (0 through 4).
#   Attitude is the observed score, which fluctuates around a true linear trajectory 
#      plus some noise.

# Next, we simulate an outcome for each person that depends on their intercept and
# slope. For example, think of Outcome as a metric like "Turnover Intention" 
# (higher values = more likely to quit).

# We set a model: Outcome = β0 + β1*(Intercept) + β2*(Slope) + error. We choose 
# β1 negative (so higher initial Attitude leads to lower turnover intention) and 
# β2 negative (so a positive Attitude slope leads to lower turnover intention). 
# We also add some random noise to Outcome.

# Parameters for outcome simulation
B0 <- 5    # baseline outcome level
B1 <- -0.5 # effect of true intercept on outcome
B2 <- -1.0 # effect of true slope (momentum) on outcome
outcome_sd <- 1  # noise for outcome

# Compute outcome for each individual
# (Each individual gets one outcome value, so use the true intercept and slope per PersonID)
Outcome <- B0 + B1 * true_intercepts + B2 * true_slopes + rnorm(N, mean = 0, sd = outcome_sd)

# Combine into a person-level data frame
person_data <- data.frame(PersonID = 1:N, 
                          Intercept_true = true_intercepts, 
                          Slope_true = true_slopes, 
                          Outcome = Outcome)

# Also compute each person's average observed Attitude across time (for later comparison)
person_data <- data_long %>% 
  group_by(PersonID) %>% 
  summarize(Average_Attitude = mean(Attitude, na.rm = TRUE)) %>%
  right_join(person_data, by = "PersonID")

head(person_data)  # preview person-level data

# Note: In a real study, you would start from an observed dataset. Here we simulate 
# data to know the "truth" (true intercepts/slopes) and demonstrate the method. 
# The Average_Attitude we calculate is each person's mean observed attitude across 
# the five time points – this represents that individual's overall level (which 
# we will compare with model-estimated intercept later).


####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####                    Fitting a RCGM                    ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# Now we fit a linear mixed-effects model (random coefficient model) to the 
# longitudinal Attitude data. This model will estimate:
#   A fixed effect for Time (the average linear trend across all individuals)
#   A random effects for each person’s intercept and slope (allowing each individual 
#     to have their own starting level and rate of change).

# We'll use lmer from lme4, specifying Attitude ~ Time + (Time | PersonID). This 
# notation means Time is a predictor, and (Time | PersonID) allows the intercept 
# and Time slope to vary by person (random intercept and random slope).

# --- Fit a linear mixed-effects model with random intercepts and slopes for Time --- #
model <- lmer(Attitude ~ Time + (Time | PersonID), data = data_long)
summary(model)

# The summary will show the fixed effects (overall intercept and slope) and the 
# variance of random effects. For example, the fixed intercept might be around 5
# (overall average Attitude at Time 0) and the fixed slope around 0 (since we 
# simulated mean no change). The random effects' standard deviations indicate 
# variability between individuals in intercepts and slopes.

# Tip: The significance of random effects can be tested by model comparison (or 
# examining the variance estimates). In our simulation, we know there is true 
# variance in slopes, so including (Time|PersonID) is appropriate. In practice, 
# you might compare with a random-intercept-only model to see if adding random 
# slopes improves fit.


####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####         Extracting Empirical Bayes Estimates         ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# A key benefit of the mixed model is that it produces empirical Bayes estimates 
# (also called BLUPs: Best Linear Unbiased Predictors) for each person's intercept 
# and slope. These are sometimes called "Bayesian estimates" of the random effects. 
# They combine the individual's observed data with the group trends to form a 
# shrunken estimate (extreme values are pulled toward the overall mean depending 
# on data reliability).

# We can extract these estimates using the coef() or ranef() functions:
#   ranef(model) gives the deviation of each individual's intercept and slope from 
#       the fixed effects.
#   coef(model) gives the individual-specific coefficients (i.e., fixed + random 
#       = the estimated intercept and slope for each person).

# We'll use coef(model)$PersonID to get a data frame of each Person's intercept 
# and slope. Then we'll merge that with our person-level data.

# Extract the empirical Bayes estimates of intercept and slope for each person
coef_df <- coef(model)$PersonID
colnames(coef_df) <- c("Intercept_est", "Slope_est")  # rename for clarity
coef_df$PersonID <- as.numeric(rownames(coef_df))     # ensure PersonID is a numeric column

# Merge these estimates into the person_data frame
person_data <- merge(person_data, coef_df, by = "PersonID")
head(person_data[, c("PersonID","Average_Attitude","Intercept_est","Slope_est","Outcome")])

# After this, person_data contains:
#   Average_Attitude: the average of observed Attitude for each person (a raw calculation).
#   Intercept_est: the model-estimated intercept for each person (their predicted Attitude 
#       at Time 0, borrowing strength from the overall model).
#   Slope_est: the model-estimated slope (trajectory) for each person.


####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####     Comparing Bayesian Intecerpt to Raw Average      ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# These two measures of "level" (Average_Attitude vs. Intercept_est) will be similar
# but not identical:
#   Average_Attitude is a simple mean of that person’s data. If an individual has 
#       an increasing or decreasing trend, this mean is influenced by their slope.
#   Intercept_est is the estimated starting point (Time 0) from the model, adjusted 
#       for noise via shrinkage. For example, if someone had an unusually high first
#       measurement due to noise, the model's intercept estimate may shrink it 
#       toward the overall mean, using information from the other time points.

# You can examine the correlation between Intercept_est and Average_Attitude to 
# see the difference. In our simulation it will be very high (close to 1) because 
# we have many time points and low noise, but in real data with fewer waves or 
# more noise, the shrinkage can be more substantial.

cor(person_data$Intercept_est, person_data$Average_Attitude)

####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####        Using Momentum in a Downstream Model          ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# Now we will use the extracted momentum (slope) in a second-stage analysis. 
# Typically, to test momentum effects on an outcome, we include both the level of 
# the predictor and its momentum as predictors of the outcome. The level can be 
# represented by either the individual's intercept (initial status) or their 
# average value over time (overall level). The momentum is represented by the 
# individual's slope (rate of change).

# In our simulated scenario, the outcome is influenced by both the true intercept 
# and slope (by construction). We'll see how well using the estimated intercept 
# and slope recovers that relationship.

# Let's fit two regression models at the person level:
#   Model 1: Outcome ~ Intercept_est + Slope_est (using the Bayesian intercept 
#        from the model as the level control).
#   Model 2: Outcome ~ Average_Attitude + Slope_est (using the observed average 
#       as the level control).


# Model 1: Using model-based intercept estimate as control for level
model1 <- lm(Outcome ~ Intercept_est + Slope_est, data = person_data)
summary(model1)$coefficients

# Model 2: Using observed average Attitude as control for level
model2 <- lm(Outcome ~ Average_Attitude + Slope_est, data = person_data)
summary(model2)$coefficients


# Both models include a term for momentum (Slope_est) and a term for overall level 
# (either Intercept_est or Average_Attitude). The coefficients will tell us:
#   The effect of the slope (Slope_est) on the outcome, holding the level constant.
#   The effect of the level (intercept or average) on the outcome.

# ---- Interpretation ---- #
# In our simulation, we expect:
#   A negative coefficient for Slope_est (momentum) – since a higher slope 
#       (improvement in Attitude) should reduce the outcome (e.g., lower turnover
#       intention).
#   A negative coefficient for the level (intercept or average) – since higher 
#       overall Attitude should reduce the outcome.

# After running the above, you might see results (for example):
#   Model 1: Outcome = β0 + (-0.6)*Intercept_est + (-0.8)*Slope_est + ...
#   Model 2: Outcome = β0 + (-0.5)*Average_Attitude + (-1.0)*Slope_est + ...

# Both models should find significant negative effects for Slope_est. This indicates 
# a momentum effect: even after controlling for a person's general Attitude level, 
# the direction and rate of change in Attitude has a unique impact on the outcome. 
# In practical terms, two individuals with the same average Attitude can have 
# different outcomes if one’s attitude was rising and the other’s was falling.


# --- Bayesian Intercept vs. Average as Level Controls --- #

# There are subtle differences between using the model-based intercept vs. 
# the raw average as the level control:
#   Bayesian Intercept (Initial level): This focuses on where each individual 
#       started (Time 0), which is useful if baseline status is theoretically 
#       important. It is less influenced by later measurements due to shrinkage. 
#       Controlling for intercept means you're comparing individuals who started 
#       at the same level but had different trajectories afterward.
#   Average Level: This represents the typical level over the whole period. 
#       Controlling for the average means you're comparing individuals with the 
#       same overall average value but different patterns of change. This approach
#       was used by Chen et al. (2011), who refer to it as absolute (average) level.

# In many cases, both approaches lead to qualitatively similar conclusions (as 
# they did in our simulation). However, they answer slightly different questions:
#   Using the average level asks: Given the same overall level of the predictor, 
#       does momentum (change) make a difference? This treats the person's average
#       as the reference point.
#   Using the initial level asks: Given the same starting point, does a higher 
#       or lower trajectory lead to a different outcome? This emphasizes baseline 
#       as a reference.

# When the trajectory is approximately linear, the average and intercept are highly 
# correlated. If there is an interest in the very beginning vs. overall exposure, 
# choose accordingly. In practice, you should also consider theory and measurement:
#   If the first measurement is special (e.g., a true baseline before an intervention), 
#       using the intercept might be logical.
#   If the variable fluctuates and no single time point is privileged, using the 
#       person's mean level might be more stable.

# In either case, the slope (momentum) is the key variable capturing the dynamic
# change. Our analysis shows how to extract that slope as a Bayesian estimate and
# include it in an outcome model. By leveraging the random coefficient model, 
# we obtain more reliable individual slopes that account for measurement error 
# and unequal observations, rather than simple two-point changes. This momentum 
# can then be interpreted as the velocity of change that has predictive power 
# beyond the static level of the variable.


####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####                      Conclusion                      ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

# We demonstrated how to:
#   Fit a multilevel growth model and obtain each individual's intercept and slope 
#       (empirical Bayes estimates).
#   Use these estimates to test momentum effects on an outcome, controlling for 
#       the variable's level (either via model-based intercept or raw average).
#   Interpret the results, recognizing that momentum (change over time) can have 
#       a unique influence on outcomes above and beyond average levels.

# Researchers can apply this approach to their own longitudinal data. By doing so, 
# they can capture momentum – the dynamic aspect of change – which enriches our 
# understanding of processes that unfold over time. This approach is general and 
# not tied to any specific variables, so you can use it for any scenario where the 
# trajectory of a predictor is theorized to impact a future outcome. Just remember 
# to center your time variable appropriately (e.g., Time 0 = start) based on the 
# interpretation you want for the intercept, and ensure your model fits well before 
# trusting the extracted random effect estimates.

# Lastly, keep in mind that the term "Bayesian estimates" here refers to the empirical
# Bayes/BLUP estimates from the mixed model – these are obtained via maximum likelihood
# estimation with the benefits of Bayesian shrinkage. This is a convenient approach 
# that does not require a full Bayesian sampling procedure, yet yields the 
# individualized estimates needed for momentum analysis. Using these in secondary 
# analyses (sometimes called two-stage modeling) is a powerful way to link growth 
# trajectories to later outcomes.
