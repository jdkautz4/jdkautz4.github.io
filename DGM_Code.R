####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####                     IMPORT FILES                     ######
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

#[] Libraries ####
# Load all libraries to be used for this analysis
library(multilevel)
library(nlme)
library(lattice)

#[] Functions ####
# Copy to clipboard
copy.table <- function(obj, size = 4096) {
  clip <- paste('clipboard-', size, sep = '')
  f <- file(description = clip, open = 'w')
  write.table(obj, f, row.names = FALSE, sep = '\t')
  close(f)  
}

LLAB.Out <- function(obj) {
  return(writeLines(paste(
    paste("-2 Log Lik = ", round(-2*logLik(obj),2), "\n"),
    paste("AIC       = ", round(AIC(obj),2), "\n"),
    paste("BIC       = ",round(BIC(obj),2)) 
    )))
}

#[] LME Optimizer ####
csetting.basic <- lmeControl(opt="optim") # From Bliese & Lang (2016)
csetting.advanced <- lmeControl(maxIter=3000,msMaxIter=3000,opt="optim",optimMethod="Nelder-Mead") # From Lang et al. (2018)
  
#[] Set Working Directory ####
  # Set the working directory ("wd" below) to be the place where the 
  # data file is saved
wd <- 'XXXXX'
setwd(wd) # Sets the working directory in R
getwd() # Check that the working directory was set correctly


####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####
#####                   Example Analyses                    #####
####>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>####

#### _______________________________________________________ ####

# US State GDP ####
  # Note: GDP is reported in millions of dollars
dat.GDP.W <- read.csv(paste(wd,"SMA PDW - Yearly GDP Example 1.csv",sep="/"))
dat.GDP <- make.univ(dat.GDP.W, dvs=dat.GDP.W[,c(2:28)], tname="TIME",outname="GDP")[,c("State","TIME","GDP")]

# For ease of interpretation, lets divide the GDP by 1k
dat.GDP$GDP1k <- dat.GDP$GDP / 1000

# Also, lets put some years in the data
dat.GDP$Year <- dat.GDP$TIME + 1997

#[] Exploring Trends ####
#_[] General trend ####
plot(
  x = c(1997:2023),
  y = as.numeric(gsub("[^0-9]", "", summary(dat.GDP.W[,c(2:28)], na.rm=T)[seq(from=4,to=162,by=6)]))/10000,
  type = "o", col = "blue", 
  xlab = "Years", ylab = "Mean GDP (1k)",
  main = "Average GDP growth (1997 to 2023)")

#_[] Individual state trends ####
xyplot(GDP1k ~ Year|as.factor(State),
       data = dat.GDP[dat.GDP$State %in% 
                        c("Wisconsin","Kentucky","Nevada",
                            "Oregon","South Dakota","Colorado",
                            "Iowa","Minnesota","Missouri"),],
       layout = c(3,3),
       type = c("p","g","r"), col="dark blue",
       col.line = "black",
       xlab = "Years",
       ylab= "GDP (1k Mil)")


#[] Selecting the data range we want to explore ####
# We see that there are two decreases in GDP, one at 2008 and one at 2020
# Did something happen that would cause these dips? (i.e., exogenous shocks?)

# Lets say we want to explore things around the 2008 recession. We would need to
#   refine the selection of the data because we do not want to leverage too much
#   data prior to the event and too much data after the event. Ideally, we would
#   have it balanced on both sides. In this data, we can! Lets grab 4 years on either 
#   side of the 2008 discontinuity 

dat1 <- dat.GDP[dat.GDP$Year >= 2004 & dat.GDP$Year <= 2012,]

#[] Create DGM Matrix ####
dat1$TIME <- dat1$TIME - 7                                 # Fix the TIME variable based on our subset
dat1$TIME.A <- ifelse(dat1$Year < 2008, dat1$TIME, 3)      # Create the absolute variant of TIME
dat1$EVENT <- ifelse(dat1$Year < 2008, 0, 1)               # Set the EVENT indicator after 2008
dat1$POST <- ifelse(dat1$Year < 2008, 0, dat1$Year - 2008) # Create a new count variable after EVENT
## --- Check DGM Matrix --- ##
dat1[1:9, c("Year","TIME","EVENT","POST")] 
dat1[1:9, c("Year","TIME.A","EVENT","POST")]


#[] Model Building Steps (Bliese & Ployhart 2002) ####

#_[] Step 0 - Check ICC ####
mod0 <- lme(GDP1k ~ 1, random=~1|State, data = dat1, na.action = na.omit)
gmeanrel(mod0)$ICC 
  # ICC = .996; 
  # Conclusion: Nearly all the variance is due to between-state factors
  # (Does this ICC level make sense?)
  
  # Hand Calculation: ICC1 = Intercept / Total
    # VarCorr(mod0)
    # as.numeric(VarCorr(mod0)[1]) / (as.numeric(VarCorr(mod0)[1]) + as.numeric(VarCorr(mod0)[2]))


#_[] Step 1 - Create basic OLS model ####
mod1 <- gls(GDP1k ~ TIME + EVENT + POST
           ,data=dat1)


#_[] Step 2 - Check for random intercepts ####
  # Create a random intercept model
mod2 <- lme(GDP1k ~ TIME + EVENT + POST
            ,random = ~1|State
            ,data = dat1
            ,control = csetting.basic)
  # Compare the random intercept model to the OLS model
anova(mod1, mod2)
  # L.Ratio = 2200.56, p < .01
  # Conclusion: Random Intercepts are useful

  # Note that estimates are nearly identical across models
  #   but the standard errors are different
round(summary(mod1)$tTable,2)
round(summary(mod2)$tTable,2)


#_[] Step 3 - Check for random slopes ####
mod3a <- update(mod2,random=~TIME|State)
mod3b <- update(mod2,random=~TIME+EVENT|State)
mod3c <- update(mod2,random=~TIME+EVENT+POST|State)
anova(mod2, mod3a, mod3b, mod3c)
  # Full random effects model (mod 3c) L.Ratio = 97.88, p < .01
  # Conclusion: A model that allows all three time covariates (TIME, EVENT, POST) to randomly vary
  #   across States is the best fitting model


#_[] Step 4 - Test for autocorrelation and heteroscedasticity ####
## --- Autocorrelation -- ##
mod4a <- update(mod3c, correlation=corAR1(form=~TIME|State))
anova(mod3c, mod4a) 
  # L.Ratio = 42.05, p < .01
  # Conclusion: Random intercept model with autocorrelation accounted provides a better fit

## --- Heteroscedasticity --- ##
mod4b <- update(mod4a, weights=varExp(form=~TIME))
anova(mod4a, mod4b)
  # L.Ratio = 13.45, p < .01
  # Conclusion: Heteroscedasticity exists in the data


#_[] Step 5 - Contrasting models ####
## --- TIME+EVENT -- ##
mod5a <- lme(GDP1k ~ TIME + EVENT, random = ~TIME + EVENT|State, dat1, correlation=corAR1(form=~TIME+EVENT|State),
             weights=varExp(form=~TIME), control=csetting.advanced)
round(summary(mod5a)$tTable,2)

## --- TIME+RECOV --- ##
mod5b <- lme(GDP1k ~ TIME + POST, random = ~TIME + POST|State, dat1, correlation=corAR1(form=~TIME+POST|State),
             weights=varExp(form=~TIME), control=csetting.advanced)
round(summary(mod5b)$tTable,2)

logLik(mod4b)
AIC(mod4b)
BIC(mod4b)

logLik(mod5a)
AIC(mod5a)
BIC(mod5a)

logLik(mod5b)
AIC(mod5b)
BIC(mod5b)

#_[] Final Working Model ####
mod <- lme(GDP1k ~ TIME + EVENT + POST
           ,random = ~TIME + EVENT + POST|State
           ,correlation=corAR1()
           ,weights=varExp(form=~TIME)
           ,data = dat1
           ,control = csetting.basic)
round(summary(mod)$tTable, dig=2)
VarCorr(mod)
LLAB.Out(mod)
  # We find that all three time covariates are significant
    # TIME  = 8.75   (t-value =  4.64; p < .01)
    # EVENT = -8.56  (t-value = -3.49; p < .01)
    # POST  = -6.58  (t-value = -6.58; p < .01)
  
  # What does this pattern suggest? 
    # (Remember: GDP numbers are "thousands of million dollars", and we divided by 1000 to make it more manageable)
    # TIME  -> Before the 2008 recession, average state GDP grew at about 8750 million dollars per year
    # EVENT -> Average state GDP dropped by about 8,560 million dollars in 2008, relative to where they were in 2007
    # POST  -> After the 2008 recession, average state GDP grew at about 2,170 million dollars per year
    #           Remember, this is relative coding: TIME (8.75) - POST (6.58) = Post-event growth (2.17)

  # We can check our interpretation of these results by comparing it to the Absolute coding variant
  # Following the same 4 steps outlined above, we find a final model that allows all time covariates
  # to randomly vary and accounts for autocorrelation as well as heteroscedasticity fits the model best
mod.absolute <- lme(GDP1k ~ TIME.A + EVENT + POST
                    ,random = ~TIME.A + EVENT + POST|State
                    ,correlation=corAR1()
                    ,weights=varExp(form=~TIME.A)
                    ,data = dat1
                    ,control = csetting.advanced)
round(summary(mod.absolute)$tTable, dig=2)
  
# Interestingly, we find now that EVENT is not significant at p < .05
    # TIME.A  = 9.37  (t-value = 4.82; p < .01)
    # EVENT =   -.26  (t-value = -.15; p = .88)
    # POST  =   2.90  (t-value = 2.45; p = .01)
  # What does this tell us? 
  # (Other than we might've messed up the coding of TIME.A or forgotten to update either the "random=" term or the 
  # "weights=" term and we want to double check that!)
  # Average state GDP grows over time. Prior to the 2008 recession, that growth was 8750 million dollars; after 
  # the recession it was 2,900 million dollars. However, the impact of the recession is not significant in absolute terms,
  # indicating that 2008 GDP was similar to 2007 GDP.
  # Again, we see the same pattern in the relative coding model, but with slightly different interpretations. See the slides
  # for how these interpretations change and what questions are best answered within each model.

  # Example of a final interpretation write-up: 
    #   Annual GDP grew, on average, at about $8.75k per year prior to the 2008 recession. However, the 
    #   2008 recession caused GDP to drop by $8.56k. While GDP continued to grow after the 2008 recession, 
    #   the rate was reduced by $6.58k annually relative to pre-recession levels.


  #_[] Visualize The Graph ####
# What if we want to see what the general trends look like... well, Paul knows how to do it in R...
# but thats way to complicated for me... I brute force it in Excel :P
# 
# First, pull the tTable:
copy.table(round(summary(mod.Pop)$tTable, dig=2))
# Next, copy the DGM matrix:
copy.table(dat2[1:9, c("Year","TIME","EVENT","POST")])
# Finally, calculate by hand using formulas in Excel....
# OR!!!
# Use the predict command to let R calcualte for you!
tdata <- data.frame(dat2[1:9, c("Year","TIME","EVENT","POST")])
tdata$predicted <- predict(mod4b, # The model you are using for analysis
                           tdata,   # the dataframe built above
                           level=0  # level=0 indicates the sample average.
)
# Copy the information out into excel to graph it!
copy.table(tdata)


#### _______________________________________________________ ####

# US State GDP - Example 2 ####
# What if we are interested in some factors that predict the changes in GDP?
# We can incorporate between-entity and time-varying effects within the RCDGM models
# Lets reimport our GDP data with some additional variables (Taxes paid on production and Population)
dat.GDP2.W <- read.csv(paste(wd,"SMA PDW - Yearly GDP Example 2.csv",sep="/"))

# Lets look at the data
names(dat.GDP2.W)
# Oh no! While we have GDP and Tax data for each year, we only have population data for each decade!
# So, we have to be creative in how we transition the data into long format. We can use make.mult.univ 
# from the multilevel package to handle GDP and Tax, but we can't do that for population. Luckily,
# make.mult.univ returns you the original dataframe with the new columns appended. So we can just grab the
# Population decade we want from there!
dat.GDP2 <- mult.make.univ(dat.GDP2.W
                           ,dvlist=list(
                            GDP = c("GDP.1997","GDP.1998","GDP.1999","GDP.2000","GDP.2001","GDP.2002","GDP.2003","GDP.2004",
                                    "GDP.2005","GDP.2006","GDP.2007","GDP.2008","GDP.2009","GDP.2010","GDP.2011","GDP.2012",
                                    "GDP.2013","GDP.2014","GDP.2015","GDP.2016","GDP.2017","GDP.2018","GDP.2019","GDP.2020",
                                    "GDP.2021","GDP.2022","GDP.2023"),
                            Tax = c("Tax.1997","Tax.1998","Tax.1999","Tax.2000","Tax.2001","Tax.2002","Tax.2003","Tax.2004",
                                    "Tax.2005","Tax.2006","Tax.2007","Tax.2008","Tax.2009","Tax.2010","Tax.2011","Tax.2012",
                                    "Tax.2013","Tax.2014","Tax.2015","Tax.2016","Tax.2017","Tax.2018","Tax.2019","Tax.2020",
                                    "Tax.2021","Tax.2022","Tax.2023")
                            ))[,c("State","TIME","GDP","Tax","Pop.2010")]

# For ease of interpretation, lets divide the GDP and Tax by 1k and grand-mean center Population (divided by 100,000)
dat.GDP2$GDP1k <- dat.GDP2$GDP / 1000
dat.GDP2$Tax1k <- dat.GDP2$Tax / 1000
dat.GDP2$Pop.GM <- (dat.GDP2$Pop.2010 / 100000) - mean(dat.GDP2$Pop.2010 / 100000)

# Also, lets put some years in the data
dat.GDP2$Year <- dat.GDP2$TIME + 1997

# Selecting the data range we want to explore ####
dat2 <- dat.GDP2[dat.GDP2$Year >= 2004 & dat.GDP2$Year <= 2012,]

# Create DGM Matrix ####
dat2$TIME <- dat2$TIME - 7                                 # Fix the TIME variable based on our subset
dat2$TIME.A <- ifelse(dat2$Year < 2008, dat2$TIME, 3)      # Create the absolute variant of TIME
dat2$EVENT <- ifelse(dat2$Year < 2008, 0, 1)               # Set the EVENT indicator after 2008
dat2$POST <- ifelse(dat2$Year < 2008, 0, dat2$Year - 2008) # Create a new count variable after EVENT
## --- Check DGM Matrix --- ##
dat2[1:9, c("Year","TIME","EVENT","POST")] 
dat2[1:9, c("Year","TIME.A","EVENT","POST")]


#[] Between-State factor ####
# We found that a model that allows random slopes for all time covariates, accounts for autocorrelation
# and heteroscedasticity was the best fitting model
mod.base <- lme(GDP1k ~ TIME + EVENT + POST
           ,random = ~TIME + EVENT + POST|State
           ,correlation=corAR1()
           ,weights=varExp(form=~TIME)
           ,data = dat2
           ,control = csetting.basic)
round(summary(mod.base)$tTable, dig=2)
VarCorr(mod.base)
LLAB.Out(mod.base)


#_[] Main effect (average or intercept change) ####
# We can now simply add in the between-state factor (Population of the 2010s) as a predictor
# the time covariates. (Also, note that the basic control settings lead to convergence errors... gotta switch)
# to the advanced settings!)
mod.Pop <- lme(GDP1k ~ TIME + EVENT + POST + Pop.GM
                ,random = ~TIME + EVENT + POST|State
                ,correlation=corAR1()
                ,weights=varExp(form=~TIME)
                ,data = dat2
                ,control = csetting.advanced)
round(summary(mod.Pop)$tTable, dig=2)
VarCorr(mod.Pop)
LLAB.Out(mod.Pop)
# Population is a significant predictor of GDP (5.41, p < .01), but this is predicting just the AVERAGE GDP
# That is, for every 100 thousand people in a state above the average state population, there is an average (intercept) 
# increase of $5,410 million


#_[] Interactive effect (individual time covariates) ####
# If we want to explore the impact on the GROWTH of GDP, we need to interact the between-state factor with
# the time covariates. Note that in order to interpret the interactions on EVENT and POST in the relative form
# We should include the interaction on TIME as well.
mod.Pop <- lme(GDP1k ~ TIME + EVENT + POST + Pop.GM
               + TIME*Pop.GM
               #+ EVENT*Pop.GM
               #+ POST*Pop.GM
               ,random = ~TIME + EVENT + POST|State
               ,correlation=corAR1()
               ,weights=varExp(form=~TIME)
               ,data = dat2
               ,control = csetting.advanced)
round(summary(mod.Pop)$tTable, dig=2)
VarCorr(mod.Pop)
LLAB.Out(mod.Pop)
  # Time*Pop.GM = .12, p < .01 -> This indicates that there is a significant effect of population level on 
  #                               the growth of GDP prior to the 2008 recession
mod.Pop <- lme(GDP1k ~ TIME + EVENT + POST + Pop.GM
               + TIME*Pop.GM
               + EVENT*Pop.GM
               #+ POST*Pop.GM
               ,random = ~TIME + EVENT + POST|State
               ,correlation=corAR1()
               ,weights=varExp(form=~TIME)
               ,data = dat2
               ,control = csetting.advanced)
round(summary(mod.Pop)$tTable, dig=2)
VarCorr(mod.Pop)
LLAB.Out(mod.Pop)
  # EVENT*Pop.GM = -.18, p < .01 -> This indicates that there is a significant effect of population level on 
  #                                 the immediate decrease in GDP assocaited with the 2008 recession, relative
  #                                 2007 levels.
mod.Pop <- lme(GDP1k ~ TIME + EVENT + POST + Pop.GM
               + TIME*Pop.GM
               #+ EVENT*Pop.GM
               + POST*Pop.GM
               ,random = ~TIME + EVENT + POST|State
               ,correlation=corAR1()
               #,weights=varExp(form=~TIME) # Convergence error when we account for this
               ,data = dat2
               ,control = csetting.advanced)
round(summary(mod.Pop)$tTable, dig=2)
VarCorr(mod.Pop)
LLAB.Out(mod.Pop)
  # POST*Pop.GM = -.07, p < .01 -> This indicates that there is a significant effect of population level on 
  #                                the relative growth in GDP following the 2008 recession
  #                                (but what does this mean, effect on the "relative growth"?)


#_[] Interactive effect (simultaneous) ####
# When looked at piece meal, we get a bit of an understanding... but what if we interact on all time covariates
# at once?
mod.Pop <- lme(GDP1k ~ TIME + EVENT + POST + Pop.GM
               + TIME*Pop.GM
               + EVENT*Pop.GM
               + POST*Pop.GM
               ,random = ~TIME + EVENT + POST|State
               ,correlation=corAR1()
               #,weights=varExp(form=~TIME) # Convergence error when we account for this
               ,data = dat2
               ,control = csetting.basic) # Also, for some reason, basic settings work best!
round(summary(mod.Pop)$tTable, dig=2)
VarCorr(mod.Pop)
LLAB.Out(mod.Pop)
  # Notice how the effect sizes have changed. Due to the multicolinarity of the time covariates, it is not uncommon
  # for the interaction terms with the time covariats to be highly inter-dependent. Conceptually, this also makes sense 
  # as it is often difficult to claim a factor will affect only one time covariate and not the other two (a counter example
  # may be something that is specifically impacting the intensity of the EVENT phenomenon ala Ballinger and Rockmann, 2010). 

  # We can check real quick to see if there is likely a true effect on POST by subsetting the data to only look at post-event
  # observations as indexed by EVENT==1.
    mod.Test <- lme(GDP1k ~ POST + Pop.GM
                   + POST*Pop.GM
                   ,random = ~POST|State
                   ,correlation=corAR1()
                   ,weights=varExp(form=~POST)
                   ,data = subset(dat2, EVENT==1)
                   ,control = csetting.advanced)
    round(summary(mod.Test)$tTable, dig=2)
  # We find that the POST*Pop interaction is now signficant, lending weight to the fact that it is likely a true effect
  # that we are finding in the full interaction model that explores the simulatensously intearctive effects on all time covariates

    
    
mod.Pop <- lme(GDP1k ~ TIME + EVENT + POST + Pop.GM
               + TIME*Pop.GM
               + EVENT*Pop.GM
               + POST*Pop.GM
               ,random = ~TIME + EVENT + POST|State
               ,correlation=corAR1()
               #,weights=varExp(form=~TIME) # Convergence error when we account for this
               ,data = dat2
               ,control = csetting.basic) # Also, for some reason, basic settings work best!
round(summary(mod.Pop)$tTable, dig=2)
  # Interpreting the effects is a bit tricky because the relative form of the DGM "reads" like an interaction already...
  # so adding in an interaction term makes the interpretation almost like a three-way interaction.
  # The effects we found were:
  #     TIME      =  8.39, p < .01
  #     EVENT     = -8.38, p < .01      
  #     POST      = -5.22, p < .01      
  #     TIME*Pop  =   .18, p < .01
  #     EVENT*Pop =  -.19, p < .01
  #     POST*Pop  =  -.10, p < .01
  # We can read this as: 
  #     "In general, there is a positive trend in state GDP growth (8,390M per year), and this positive trend is strengthened by population 
  #      such that for every 100k more people in the state compared to the average state population, GDP grows by an 
  #      additional $180M per year. The 2008 recession lead to an immediate decrease in state GDP, which was exacerbated
  #      by population such that states with greater population than average experienced a larger decline ($-190M per 
  #      each 100k people above the average state population), realtive to where we would expect them to be absent the 2008 recesion. 
  #      Following the 2008 recession, average state GDP grew by $3,170M (a decline of $5,140M compared to pre-recession growth). 
  #      This decline in the growth rate was exacerbated by population such that states with greater population than average experienced 
  #      a slower rate of growth following the 2008 recession (-100M per year each 100k people above the average state population)."
  # Or, to simplify it even further:
  #     "In general, there is a positive trend in state GDP growth, which is stronger in states that have greater population.
  #      The 2008 recession lead to an immediate decrease in state GDP relative to where we would expect them to be, which was 
  #      worse in states with greater population. Following the 2008 recession, state GDP growth did not return to pre-recession 
  #      levels and states with greater population experienced a greater slow down in their GDP growth rate."
    
    
#_[] Visualizing the interactive effect ####
# What if we want to see what the general trends look like... well, Paul knows how to do it in R...
# but thats way to complicated for me... I brute force it in Excel :P
# 
# First, pull the tTable:
copy.table(round(summary(mod.Pop)$tTable, dig=2))
# Next, copy the DGM matrix:
copy.table(dat2[1:9, c("Year","TIME","EVENT","POST")])
# Finally, calculate by hand using formulas in Excel....
# OR!!!
# Use the predict command to let R calcualte for you!
tdata <- data.frame(dat2[1:9, c("Year","TIME","EVENT","POST")])
# Choose the state group you want
tdata$Pop.GM <- mean(dat2$Pop.GM)
tdata$Pop.GM <- mean(dat2$Pop.GM)-(.5*sd(dat2$Pop.GM)) # Due to the skew, one SD takes us way below the lower bounds... so only use .5
tdata$Pop.GM <- mean(dat2$Pop.GM)+(.5*sd(dat2$Pop.GM))

# Once you choose the state group based on population you want, then use predict() to estimate the scores
tdata$predicted <- predict(mod.Pop, # The model you are using for analysis
                           tdata,   # the dataframe built above
                           level=1  # level=0 indicates the sample average.
)
# Copy the information out into excel
copy.table(tdata)


#### _______________________________________________________ ####

# US State GDP - Example 3 ####
# I have included a time-varying variable (taxes paid) in the data. Try for yourself to model the effects.

# I found the following results:
#     Tax (1k)  =   6.08, p < .01
#     TIME      =   4.25, p < .01
#     EVENT     =   2.03, p = .31      
#     POST      = -10.79, p < .01      
#     TIME*Tax  =   -.12, p < .01
#     EVENT*Tax =   -.25, p < .01
#     POST*Tax  =    .45, p < .01

