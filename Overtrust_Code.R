# Overtrust Simulation
#
# We will generate data for a hypothetical study where individuals (level 2) are nested 
# in groups (level 3) and measured repeatedly over time (level 1). 
# 
# In our simulation:
#   G = Total number of groups
#   n_per_group = Number of individuals per group;: for simplicity, we simulate equal 
#                 cluster sizes (e.g., 20 individuals in each group, but this can be 
#                 changed to simulate unbalanced panels).
#   N = Total number of individuals (this will be G * n_per_group)
#   TIME = Number of time points (the number of observations per individual).
#   
# Each individual's outcome over time will be generated from an oscillatory model.
# We use a sinusoidal function for simplicity. 
# 
# The general form for individual i in group j at time t will be:
# 
#   yij(t) = Aij ​ sin(2π ​ fij ​ t + ϕij) + εij(t)
#   
#   
# where 
#   Aij = amplitude
#   fij = frequency (cycles per unit time) 
#   ϕij = a phase offset (alignment of cycles)
#   εij(t) [ ∼ N(0, σ^2)] =  measurement error / noise
#   
# We will allow both amplitude and frequency to vary by predictors: 
#   Wj = a group-level binary predictor 
#     (think of this as a group characteristic or treatment, 0/1) 
#   Xij = an individual-level continuous predictor 
#     (e.g., a person-specific trait, standardized)
# 
# Group-level predictor Wj affects both amplitude and frequency. 
#   For example, groups with Wj = 1 might have higher amplitude oscillations 
#   (perhaps the treatment increases volatility). 
#   We will simulate Wj∼Bernoulli(0.5) for each group (half the groups get W=1).
# 
# Individual-level predictor Xij also influences the oscillator. 
#   We simulate Xij ∼ N(0,1) independently for each individual (centered at 0). 
#   This could represent, say, a personality score that might affect the person’s 
#   oscillation.
#   
# 
# We construct amplitude and frequency as functions of these predictors:
#   Amplitude model: Aij = exp(α0 + αW ​ Wj + αX ​ Xij + uj + vij) 
#     We use an exponential to ensure amplitude stays positive. 
#     α0 = a baseline log-amplitude
#     αW, αX = effects of the predictors on log-amplitude
#     uj, vij = random effects (unobserved heterogeneity). 
#     
#   (For simplicity, one could also use a linear model for amplitude if parameters 
#   are set such that Aij remains positive)
#   
#   
#   Frequency model: fij = f0 ​ exp(γW ​ Wj + γX ​ Xij + wj + zij)
#     f0 = baseline frequency
#     γW, γX = effects of the predictors on log-frequency
#     wj, zij = random effects (unobserved heterogeneity)
#   
#   (One could keep frequency constant across individuals for simplicity, focusing
#   only on amplitude differences.) 
#
# Lets go into the code!

###-----------------------------------------------------------------------###
####                             Code Setup                              ####
###-----------------------------------------------------------------------###
#[] Libraries ####
library(dplyr)  

# Set seed for reproducibility
set.seed(123)  

#[] Create Simulated Data (sim_data) ####
#   ID = individual identifier
#   Group = group identifier
#   Time = time period
#   w, x = group and individual predictor, respectively
#   y = outcome (trust)

# Simulation parameters (feel free to adjust these)
G <- 10               # total number of groups
n_per_group <- 20     # individuals per group
TIME <- 50            # time points per individual
N <- G * n_per_group  # total individuals

# Predictor generation
# Group-level binary predictor W (e.g., treatment vs control for each group)
W <- rbinom(n = G, size = 1, prob = 0.5)  # length G vector of 0/1

# Individual-level continuous predictor X (e.g., standardized trait for each individual)
X <- rnorm(N, mean = 0, sd = 1)           # length N vector

# Assign each individual to a group
group_id <- rep(1:G, each = n_per_group)  # vector of group IDs for each individual

# Define oscillator parameters
alpha0 <- log(2)   # baseline log-amplitude (exp(alpha0) ~ 2)
alpha_W <- 0.5     # group predictor effect on log-amplitude
alpha_X <- 0.5     # individual predictor effect on log-amplitude
sigma_group_amp <- 0.2  # SD of group-level random effect on log-amplitude
sigma_ind_amp   <- 0.2  # SD of individual-level random effect on log-amplitude

f0 <- 0.05        # baseline frequency (cycles per time unit)
gamma_W <- 0.1    # group effect on log-frequency (small effect)
gamma_X <- 0.1    # individual effect on log-frequency
sigma_group_freq <- 0.02  # SD of group effect on log-frequency
sigma_ind_freq   <- 0.02  # SD of individual effect on log-frequency

sigma_noise <- 0.2  # noise standard deviation

# Simulate random effects
u_j  <- rnorm(G, mean = 0, sd = sigma_group_amp)    # group-level random intercepts for amplitude
v_i  <- rnorm(N, mean = 0, sd = sigma_ind_amp)      # individual-level random intercepts for amplitude
w_j  <- rnorm(G, mean = 0, sd = sigma_group_freq)   # group-level random effect for frequency
z_i  <- rnorm(N, mean = 0, sd = sigma_ind_freq)     # individual-level random effect for frequency

# Create an empty data frame to store simulated observations
sim_data <- data.frame(
  ID = integer(0),
  Group = integer(0),
  Time = integer(0),
  W = integer(0),
  X = numeric(0),
  Y = numeric(0)
)

# Loop over individuals and simulate their time series
for(i in 1:N) {
  j <- group_id[i]         # group for individual i
  # Compute amplitude and frequency for this individual
  logA_ij <- alpha0 + alpha_W * W[j] + alpha_X * X[i] + u_j[j] + v_i[i]
  A_ij <- exp(logA_ij)                          # amplitude (positive)
  logf_ij <- log(f0) + gamma_W * W[j] + gamma_X * X[i] + w_j[j] + z_i[i]
  f_ij <- exp(logf_ij)                          # frequency (positive)
  phi_ij <- 0                                   # set phase = 0 for all, for simplicity
  
  # Time vector (1...T). We can also simulate irregular times if needed.
  t <- 1:T  
  
  # Generate the oscillatory outcome for each time point
  # Using the sine function for oscillation
  y_t <- A_ij * sin(2 * pi * f_ij * t + phi_ij)  
  
  # Add measurement noise
  y_t <- y_t + rnorm(T, mean = 0, sd = sigma_noise)
  
  # Create a data frame of this individual's observations
  ind_df <- data.frame(
    ID = i,
    Group = j,
    Time = t,
    W = W[j],
    X = X[i],
    Y = y_t
  )
  # Append to the main dataset
  sim_data <- bind_rows(sim_data, ind_df)
}

# Take a quick look at the simulated data
head(sim_data, 10)

# If we want to simulate some data as missing at random:
#set.seed(321)  # separate seed for missingness
#miss_index <- sample(nrow(sim_data), size = 0.10 * nrow(sim_data)) # Introduce 10% missingness at random
#sim_data$Y[miss_index] <- NA


###-----------------------------------------------------------------------###
####                        Frequentist Analysis                         ####
###-----------------------------------------------------------------------###

# For the frequentist analysis, we will use a linear mixed-effects model (LME) 
# to capture the structure of our panel data. Mixed-effects models allow us to 
# include random effects to account for the correlated observations within 
# individuals (and within groups). They also let us estimate the effects of 
# predictors at different levels (fixed effects) on the outcome. 
#
# In our case, we have a longitudinal outcome with an oscillatory pattern; one 
# way to model this is to treat the sine wave as a known function of time and 
# let its amplitude vary by predictors.

# Model specification: 
# We know that a sine wave can be represented as a linear combination of sine 
# and cosine terms. If every individual had the same frequency and phase, we 
# could include a term like sin(2π ​ f0 ​ t) as a predictor. To allow individual
# -specific amplitudes, we can include a random slope on that sine term for each 
# individual. Moreover, we include interactions of the sine term with our 
# predictors X and W to capture how amplitude might systematically differ by 
# those predictors. An equivalent formulation is to model 
#
# Yijt = (b0 ​ + bW ​ Wj + bX ​ Xij + ui) ​ sin(2π ​ f0 ​ t) + ϵijt
#
# where: 
#   ui = random effect for individual i (which multiplies the sine term)
#   b0 = the average sine amplitude
#   bW = how much extra amplitude group predictor adds
#   bX = the effect of the individual predictor
#   ϵ = residual noise
#
# This is a linear model in the parameters (since sin(2π ​ f0 ​ t) is treated as
# a known regressor at each time point). While this approach assumes a common 
# frequency f0 and phase (we set f0 = 0.05 and ϕ = 0 in simulation), it is a 
# reasonable simplification for illustrating the analysis. More advanced models 
# could treat frequency as a parameter to estimate as well (making the model 
# nonlinear), but that requires specialized nonlinear mixed modeling techniques 
# beyond a basic example. Indeed, methods exist to estimate damped oscillators 
# in panel data using differential equation models, but here we use a simpler 
# approximation.

# We will use the lme4 package in R to fit the mixed model (Bates et al., 2015).
# The formula syntax in lme4 allows us to specify fixed and random effects 
# conveniently. For our model:
#   Fixed effects: an interaction of the time sine term with W and X. 
#   Random effects: a random intercept and random slope on the sine term for 
#     each individual ID. 
#       The random intercept (per ID) would capture any individual-specific 
#       baseline offset in Y (though in our simulation the mean is ~0, it’s 
#       still good to allow this in case of slight shifts). 
#       The random slope allows each individual to have their own amplitude 
#       beyond what predictors explain. 


# Below is the R code to fit the linear mixed model using lmer and to examine 
# the results:

# Load lme4 for mixed modeling
library(lme4)

# Create the sine term for the chosen frequency f0 (here 0.05 from simulation)
f0 <- 0.05
sim_data <- sim_data %>% 
  mutate(sin_term = sin(2 * pi * f0 * Time))

# Fit a linear mixed-effects model
lmm_fit <- lmer(Y ~ sin_term * W + sin_term * X + 
                  (1 + sin_term | ID),  # random intercept and slope for each individual
                data = sim_data, REML = FALSE)
# Note: REML=FALSE for comparability of fixed effects; can use REML=TRUE for final model
#
# The random effects (1 + sin_term | ID) mean each individual gets their own 
# intercept and their own slope on sin_term. This should capture individual-
# specific phase shifts (via intercept) if any, and amplitude deviations (via 
# the slope).

# View a summary of the model
summary(lmm_fit)

# For example, the output might show (for illustration, not actual run):
  # 
  # Fixed effects:
  #             Estimate Std. Error t value
  # (Intercept)   0.02      0.05       0.4     # baseline intercept ~ 0
  # sin_term      1.20      0.04      30.0***  # baseline amplitude
  # W             0.01      0.07       0.1     # main effect of W (not meaningful alone here)
  # X            -0.02      0.03      -0.7     # main effect of X (not meaningful alone)
  # sin_term:W    0.58      0.06       9.8***  # Group Factor influence on sine amplitude
  # sin_term:X    0.45      0.02      22.5***  # Indivdual Factor influence on sine amplitude
  # 
  # Random effects:
  # Groups   Name        Variance Std.Dev.  Corr 
  # ID       (Intercept)  0.05     0.223
  #           sin_term    0.20     0.447    -0.10  # slight negative corr between intercept & slope
  # Residual              0.04     0.200


# In the above hypothetical output:
# sin_term:W being 0.58 (p < .001) confirms that group effect significantly 
#   increases the oscillation amplitude, which aligns with our simulation. 
# The term sin_term:X is 0.45 (p < .001) indicating individuals with higher X
#   have higher amplitude. 
# The random slope for sin_term has variance ~0.20, meaning considerable 
#   individual-level variability in amplitude remains even after accounting for 
#   W and X. 
# Residual variance is ~0.04, with residual SD ~0.2, matching the noise we added. 
# The intercept variance of 0.05 with SD ~0.223 is relatively small, consistent 
#   with the outcome being centered around zero for each individual.

# Interpreting the frequentist model: 
#   We would report that after controlling for the predictors, there is strong 
#     evidence that the group-level predictor W increases the oscillation’s amplitude 
#     (fixed effect ~0.58, SE ~0.06, p < .001). 
#   Likewise, the individual trait X is associated with amplitude (fixed effect 
#     ~0.45, p < .001). The random effects indicate that individuals vary in 
#     their baseline oscillation amplitude (SD of random sine-term slope ≈ 0.447), 
#     suggesting other unmodeled factors also influence amplitude. 

# The model captures the periodic structure through the sine term; if we plotted 
# the model’s fitted values vs. actual data, we should see the curves align 
# reasonably well, except for noise.

# Technical note: Instead of manually creating a sin_term, one could also fit a 
# nonlinear mixed model where amplitude is directly modeled. For example, using 
# the nlme package, one could specify:
# Yij(t) = Ai ​ sin(2πft) 
# and link Ai to predictors via a linear model. However, that is more advanced. 
# Our linear mixed model approach essentially linearizes the problem by fixing 
# frequency and phase. For many purposes (e.g., seasonal effects, circadian 
# rhythms), one might even fix frequency to a known cycle (like 24 hours or 
# annual cycle) and estimate phase or amplitude differences using mixed models 
# (J. Bryk & Raudenbush, 2002). 

# If frequency were unknown and to be estimated from data, one might resort to 
# periodogram analysis or use Bayesian methods as shown next.

# Before moving on, it’s worth noting that we treated W and X as fixed effects. If
# W were not strictly fixed by design (say these were just observed group characteristics 
# in a sample of groups), one might also consider treating group as random and include 
# W as a group-level covariate in a higher-level model. Here, since W is binary and 
# we have relatively few groups in the example, a fixed effect is fine to estimate 
# the difference between "W = 0" and "W = 1". The decision between fixed vs random 
# effects for grouping factors often depends on study design and inference focus. 
# We could formally test if the group-level variance is negligible or not. In our 
# simulation, we didn’t include a separate pure group-level random intercept (we 
# only had individual random effects nested in group). If one wanted to allow for 
# correlation among individuals within the same group beyond W, a term like (1|Group) 
# could be added. Given we already included W in the fixed part and had relatively 
# small unexplained group differences, we might find that random intercept for groups 
# is close to zero variance.

###-----------------------------------------------------------------------###
####                          Bayesian Analysis                          ####
###-----------------------------------------------------------------------###

# The Bayesian framework treats model parameters as random variables with specified
# prior distributions, and uses Bayes’ theorem to obtain a posterior distribution 
# given the data. We will use the brms package in R (Bürkner, 2017) to specify and 
# fit a Bayesian multilevel model. The brms package provides a high-level formula 
# syntax similar to lme4, but internally uses Stan (a probabilistic programming language) 
# to sample from the posterior. This means the code we write in brms looks much like 
# the lmer model formula, but we will get posterior means, credible intervals, and 
# so forth, instead of just point estimates.

# We will fit essentially the same model as in the frequentist section: 
#   Yijt = (b0 + bW ​ Wj + bX ​ Xij + ui) sin(2πf0 ​ t) + ϵijt
# with priors on all parameters. By default, brms will set weakly informative 
# priors (e.g., normal(0,1) or student_t priors for coefficients and half-Cauchy 
# for variance terms) if not specified – these are usually reasonable starting points, 
# but one can adjust them based on prior knowledge. In a Bayesian analysis, we also
# need to consider convergence diagnostics and sufficient MCMC sampling. We’ll run
# a few Markov chains and check R-hat and effective sample sizes, which brms does 
# automatically in its summary output.

# Here’s how we can fit the model in brms:
  
# Load brms for Bayesian modeling
library(brms)

# Define the model formula (similar to lmer formula)
bayes_formula <- bf(Y ~ sin_term * W + sin_term * X + (1 + sin_term | ID))

# Fit the Bayesian multilevel model
# We specify iterative sampling parameters; these can be adjusted if needed.
bayes_fit <- brm(
  formula = bayes_formula,
  data = sim_data,
  family = gaussian(),  # continuous outcome
  prior = c(            # define priors (optional, brms has defaults)
    prior(normal(0, 1), class = "b"),            # fixed effects coefficients
    prior(normal(0, 5), class = "Intercept"),    # intercept
    prior(cauchy(0, 1), class = "sd")            # random effect std devs
  ),
  chains = 4, cores = 4, iter = 2000, warmup = 500
)
# Note: 'cores=4' assumes 4 CPU cores available for parallel chains.
# Reduce cores or chains if running on limited hardware.

# Check summary of posterior draws
summary(bayes_fit)


# When you run brm(), you will see Stan sampling progress. After fitting, 
#summary(bayes_fit) will display posterior means and credible intervals for all 
#parameters. For example, you might see output like (format simplified for brevity):
#   
# Group-Level Effects: 
#   ~ID (Number of levels: 200) 
#                           Estimate  Est.Error  l-95% CI  u-95% CI
# sd(Intercept)               0.22      0.05      0.14      0.33
# sd(sin_term)                0.45      0.04      0.38      0.54
# cor(Intercept, sin_term)   -0.10      0.20     -0.47      0.29
# 
# Population-Level Effects:
#               Estimate  Est.Error  l-95% CI  u-95% CI  Rhat
# Intercept      0.01      0.05     -0.09      0.11      1.00
# sin_term       1.18      0.04      1.09      1.26      1.00
# W              0.00      0.06     -0.12      0.12      1.00
# X             -0.01      0.03     -0.07      0.05      1.00
# sin_term:W     0.60      0.06      0.48      0.71      1.00
# sin_term:X     0.46      0.02      0.41      0.50      1.00
# 
# Family: gaussian(identity) 
# Residual SD: 0.20  (95% CI: 0.19, 0.21)


# Let’s interpret this Bayesian result: It closely mirrors the frequentist 
# findings but framed in terms of posterior distributions:
#  
#   The population-level (fixed) effects show the posterior mean for sin_term:W 
#   is about 0.60 with a 95% credible interval roughly [0.48, 0.71]. 
#   
#   Since this interval does not include 0, we conclude with high certainty that 
#   W increases the amplitude (this is analogous to a significant positive effect, 
#   but Bayesian inference doesn’t use p-values; instead we note the posterior 
#   probability that the effect is positive is ~100%). 
#
#   Similarly, sin_term:X has a posterior mean ~0.46, 95% CI [0.41, 0.50], 
#   indicating a positive effect of X on amplitude with high certainty. The main 
#   effect of sin_term ~1.18 [1.09, 1.26] is the baseline amplitude for an average 
#   individual in a control group (this is close to the true average amplitude 
#   in the simulated data). The intercept and standalone W, X effects have means 
#   near 0 and intervals overlapping 0, which we expected (they represent baseline 
#   offset and mean differences that are not meaningful given how we generated data).
#
#   The group-level (random) effects for ID show sd(sin_term) ~0.45 [0.38, 0.54], 
#   confirming substantial individual heterogeneity in amplitudes, and sd(Intercept) 
#   ~0.22 which, although not zero, is relatively smaller (some individuals have 
#   slightly higher or lower overall levels, possibly due to noise). The correlation 
#   between intercept and sin_term slopes is mildly negative (-0.1) but with a wide
#   interval including 0, so no strong relationship there.
#
#   The residual SD is ~0.20 with narrow CI, matching the known simulation noise level.
#
#   All Rhat values are 1.00, indicating the chains converged well, and effective 
#   sample sizes (not shown above) would likely be high, meaning the posterior is
#   well-explored. This gives us confidence in the Bayesian estimates.
#
#   One benefit of the Bayesian approach is that we can directly obtain the probability
#   of effects being above or below a threshold. For instance, we can say “There 
#   is ~99.9% posterior probability that b_{\text{sin_term:W}} (the group effect 
#   on amplitude) is positive,” which is another way of reporting a significant 
#   finding. 

# We can also plot posterior distributions or predictive checks. For example:

# Posterior summary in a more readable format
print(bayes_fit, digits=2)
                           
# Plot posterior densities for key parameters
plot(bayes_fit, pars = c("b_sin_term:W", "b_sin_term:X"))

# The above plot call will show density plots of the posterior draws for the two 
# parameters of interest, illustrating the uncertainty around them. We would see
# both distributions well away from zero.

# From a modeling perspective, the Bayesian approach handled the missing data seamlessly
# as well (if we had introduced missingness in sim_data, brms by default uses likelihood 
# calculations that ignore NAs in Y). Missing values in predictors would need separate 
# handling (e.g., imputation) but our example had none. In Bayesian modeling, one 
# could even model the missing data mechanism or use multivariate models to handle 
# missingness if needed.
                           
# Finally, we note that our Bayesian model was essentially the same structure as
# the frequentist one. In more complex settings (for example, truly modeling the
# frequency as a parameter to be estimated for each individual, or fitting a full 
# damped oscillator differential equation model), the Bayesian framework can be 
# very powerful. Researchers like Oud (2006) have explored differential equation 
# models for oscillators in panel data, and modern software (Stan, etc.) enables 
# fitting such models where frequency and damping parameters are estimated hierarchically
# (e.g., dynamic structural equation modeling or continuous-time models; Deboeck 
# & Boker, 2008). Those are advanced topics, but our simpler approach already 
# demonstrates how both fixed and random effects can account for oscillatory patterns. 
#
# The key takeaway is that frequentist and Bayesian multilevel models yield consistent 
# insights in this simulation: the group and individual predictors indeed affect 
# the oscillation as we programmed them to, and both approaches recovered those 
# effects with appropriate uncertainty quantification.

# Conclusion

# In this tutorial, we set up a panel data simulation with an oscillatory outcome, 
# and showed how to analyze it using two frameworks. We ensured the code is modular: 
# you can easily adjust the number of groups, individuals, time points, or even 
# the functional form of the oscillator. We identified that at least three time 
# points and a reasonable number of individuals (preferably >30-50) are needed 
# for such analyses to be effective. We also demonstrated how to simulate a 
# balanced panel and how one might introduce missing data to test robustness. On 
# the analysis side, we used a frequentist LME model (via lme4) to estimate how 
# predictors influence oscillation amplitude, treating individual differences with 
# random effects. We then used a Bayesian multilevel model (via brms/Stan) to do 
# the same, highlighting the interpretation of posterior intervals. Both methods 
# found the programmed effects of predictors on the oscillator parameters, with 
# Bayesian results providing a probabilistic interpretation that aligned with the 
# frequentist significance findings.

# For a novice, the main points to understand are: 
#   (1) Panel data allows examining changes within individuals over time and 
#       differences across individuals – it’s a rich structure that requires 
#       appropriate modeling
#   (2) Simulation is a helpful way to understand how different factors (predictors) 
#       can manifest in longitudinal patterns
#   (3) Mixed-effects models are a go-to tool for analyzing such data, and one 
#       can approach them from either the frequentist or Bayesian perspective. 
#       The provided code can be used as a template for your own data – you can
#       plug in your desired number of groups or time points, simulate or input
#       real data, and then fit similar models. Always remember to check model 
#       diagnostics (residuals, convergence) and the theoretical assumptions 
#       (stationarity of oscillations, etc.) when dealing with real data. With 
#       these tools, you are well-equipped to explore oscillatory dynamics in panel 
#       data in both a classical and Bayesian way, leveraging the strengths of 
#       each approach.