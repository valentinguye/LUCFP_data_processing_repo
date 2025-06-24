## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

set.seed(0)

if(requireNamespace("data.table", quietly = TRUE)) library(data.table)

require_DT_ON = function(){
  if(!requireNamespace("data.table", quietly = TRUE)){
    knitr::opts_chunk$set(eval = FALSE)
    cat("Evaluation of the next chunks requires 'data.table', which is not present.")
  }
}

require_DT_OFF = function(){
  knitr::opts_chunk$set(eval = TRUE)
}

library(fixest)
setFixest_nthreads(1)

## ----echo=TRUE----------------------------------------------------------------
library(fixest)
data(trade)


## ---- echo=FALSE, results='asis'----------------------------------------------
tab = head(trade)
knitr::kable(tab)

## -----------------------------------------------------------------------------
gravity_pois = fepois(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)

## -----------------------------------------------------------------------------
print(gravity_pois)

## -----------------------------------------------------------------------------
summary(gravity_pois, se = "twoway")

## ---- eval = FALSE------------------------------------------------------------
#  # Two ways to summon clustering on the Product variable
#  # - by reference:
#  summary(gravity_pois, cluster = "Product")
#  # - with a formula:
#  summary(gravity_pois, cluster = ~Product)

## ---- eval = TRUE-------------------------------------------------------------
summary(gravity_pois, cluster = ~Product)

## -----------------------------------------------------------------------------
gravity_simple = fepois(Euros ~ log(dist_km), trade)
# We use a formula to specify the variables used for two way clustering
# (note that the values of the variables are fetched directly in the original database)
summary(gravity_simple, cluster = ~Origin + Destination)

## -----------------------------------------------------------------------------
fepois(Euros ~ log(dist_km), trade, cluster = ~Product)

## -----------------------------------------------------------------------------
gravity_ols = feols(log(Euros) ~ log(dist_km) | Origin + Destination + Product + Year, trade)

## -----------------------------------------------------------------------------
gravity_negbin = fenegbin(Euros ~ log(dist_km) | Origin + Destination + Product + Year, trade)


## ---- eval=FALSE--------------------------------------------------------------
#  etable(gravity_pois, gravity_negbin, gravity_ols,
#           se = "twoway", subtitles = c("Poisson", "Negative Binomial", "Gaussian"))

## ---- echo=FALSE, results='asis'----------------------------------------------
tab = etable(gravity_pois, gravity_negbin, gravity_ols, se = "twoway", subtitles = c("Poisson", "Negative Binomial", "Gaussian"))
# problem to display the second empty line in markdown
knitr::kable(tab[-2, ])

## -----------------------------------------------------------------------------
gravity_subfe = list()
all_FEs = c("Year", "Destination", "Origin")
for(i in 0:3){
	gravity_subfe[[i+1]] = fepois(Euros ~ log(dist_km), trade, fixef = all_FEs[0:i])
}

## ---- eval=FALSE--------------------------------------------------------------
#  etable(gravity_subfe, cluster = ~Origin+Destination)

## ---- echo=FALSE, results='asis'----------------------------------------------
tab = etable(gravity_subfe, cluster = ~Origin+Destination)
knitr::kable(tab)

## -----------------------------------------------------------------------------
res_multi = fepois(Euros ~ log(dist_km) | csw0(Year, Destination, Origin), trade)

## -----------------------------------------------------------------------------
# with two-way clustered SEs
etable(res_multi, cluster = ~Origin+Destination, tex = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  # we set the dictionary once and for all
#  myDict = c("log(dist_km)" = "$\\ln (Distance)$", "(Intercept)" = "Constant")
#  # 1st export: we change the signif code and drop the intercept
#  etable(res_multi, signifCode = c("a" = 0.01, "b" = 0.05),
#         drop = "Const", dict = myDict, file = "Estimation Tables.tex",
#         replace = TRUE, title = "First export -- normal Standard-errors")
#  # 2nd export: clustered S-E + distance as the first coefficient
#  etable(res_multi, cluster = ~Product, order = "Dist",
#         dict = myDict, file = "Estimation Tables.tex",
#         title = "Second export -- clustered standard-errors (on Product variable)")
#  

## -----------------------------------------------------------------------------
fixedEffects = fixef(gravity_pois)
summary(fixedEffects)

## -----------------------------------------------------------------------------
fixedEffects$Year

## ---- fig.width=7-------------------------------------------------------------
plot(fixedEffects)

## -----------------------------------------------------------------------------
base = iris
names(base) = c("y", "x1", "x_endo_1", "x_inst_1", "fe")
set.seed(2)
base$x_inst_2 = 0.2 * base$y + 0.2 * base$x_endo_1 + rnorm(150, sd = 0.5)
base$x_endo_2 = 0.2 * base$y - 0.2 * base$x_inst_1 + rnorm(150, sd = 0.5)

est_iv = feols(y ~ x1 | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
est_iv

## -----------------------------------------------------------------------------
fitstat(est_iv, ~ ivf1 + ivwald1 + ivf2 + ivwald2, cluster = "fe")

## -----------------------------------------------------------------------------
setFixest_print(fitstat = ~ . + ivwald2)
est_iv

## -----------------------------------------------------------------------------
est_iv_fe = feols(y ~ x1 | fe | x_endo_1 + x_endo_2 ~ x_inst_1 + x_inst_2, base)
est_iv_fe

## -----------------------------------------------------------------------------
summary(est_iv_fe, stage = 1)

## -----------------------------------------------------------------------------
etable(summary(est_iv_fe, stage = 1:2), fitstat = ~ . + ivfall + ivwaldall.p)

## -----------------------------------------------------------------------------
# Our base data for this section
base = iris
names(base) = c("y", paste0("x", 1:3), "fe1")
# Create another "fixed-effect"
base$fe2 = rep(letters[1:5], 30)
head(base)

## -----------------------------------------------------------------------------
est_comb = feols(y ~ x1 | fe1^fe2, base)
est_comb

## -----------------------------------------------------------------------------
fixef(est_comb)[[1]]

## -----------------------------------------------------------------------------
head(base)

## -----------------------------------------------------------------------------
est_vs = feols(y ~ x1 | fe1[x2], base)
est_vs

## -----------------------------------------------------------------------------
summary(fixef(est_vs))

## -----------------------------------------------------------------------------
data(airquality)
res_i1 = feols(Ozone ~ Solar.R + i(Month), airquality)
res_i2 = feols(Ozone ~ Solar.R + i(Month, ref = 8), airquality)
res_i3 = feols(Ozone ~ Solar.R + i(Month, keep = 5:6), airquality)

etable(res_i1, res_i2, res_i3, dict = c("6" = "June", "Month::5" = "May"), 
       order = c("Int|May", "Mon"))

## ---- eval = TRUE-------------------------------------------------------------
# Sample data illustrating the DiD
data(base_did)
head(base_did)

## ---- eval = TRUE-------------------------------------------------------------
# Estimation of treatment × period effects
# We also add individual and period fixed-effects:
est_did = feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)
est_did

## ---- fig.width=7-------------------------------------------------------------
iplot(est_did)

## -----------------------------------------------------------------------------
data(base_stagg)

head(base_stagg)

## ---- echo=FALSE--------------------------------------------------------------
if(requireNamespace("ggplot2", quietly = TRUE)){
  library(ggplot2)
  ggplot(aggregate(base_stagg[, c('year_treated', 'treatment_effect_true')], 
                   by = list(year = base_stagg$year, group = to_integer(base_stagg$year_treated)), 
                   mean), 
         aes(year, group, fill = year>=year_treated, alpha = treatment_effect_true)) + 
    geom_tile(colour = "white", lwd = 1) + 
    scale_fill_brewer('Treated?', palette = 'Set1') +
    scale_alpha('Avg. treatment\neffect') +
    labs(x = 'Year', y = 'Group') +
    theme_minimal()
} else {
  print("This graph requires ggplot2 which is currently not available.")
}

## -----------------------------------------------------------------------------
# "Naive" TWFE DiD (note that the time to treatment for the never treated is -1000)
# (by using ref = c(-1, -1000) we exclude the period just before the treatment and 
# the never treated)
res_twfe = feols(y ~ x1 + i(time_to_treatment, treated, ref = c(-1, -1000)) | id + year, base_stagg)

# To implement the Sun and Abraham (2020) method,
# we use the sunab(cohort, period) function
res_sa20 = feols(y ~ x1 + sunab(year_treated, year) | id + year, base_stagg)

## -----------------------------------------------------------------------------
# Plot the two TWFE results
iplot(list(res_twfe, res_sa20), sep = 0.5)

# Add the true results
att_true = tapply(base_stagg$treatment_effect_true, base_stagg$time_to_treatment, mean)[-1]
points(-9:8, att_true, pch = 15, col = 4)

legend("topleft", col = c(1, 4, 2), pch = c(20, 17, 15), 
       legend = c("TWFE", "Truth", "Sun & Abraham (2020)"))

## -----------------------------------------------------------------------------
# The full ATT
summary(res_sa20, agg = "att")

## -----------------------------------------------------------------------------
# Full disaggregation (you could have used summary instead of etable)
head(etable(res_sa20, agg = FALSE), 20)

## -----------------------------------------------------------------------------
base = iris
names(base) = c("y", "x1", "x2", "x3", "species")
# Defining the macro variables
setFixest_fml(..ctrl = ~poly(x2, 2) + poly(x3, 2))
# Accessing them
xpd(y ~ x1 + ..ctrl)

# Definition at run time
vars = c("x2", "x2^2", "x3")
for(i in 1:3){
  print(xpd(y ~ x1 + ..ctrl, ..ctrl = vars[1:i]))
}

## -----------------------------------------------------------------------------
feols(y ~ x1 + ..ctrl, base)

## -----------------------------------------------------------------------------
data(longley)
xpd(Armed.Forces ~ Population + ..("GNP|ployed"), data = longley)

## -----------------------------------------------------------------------------
feols(Armed.Forces ~ Population + ..("GNP|ployed"), longley)

## -----------------------------------------------------------------------------
est1 = feols(y ~ l(x1, 0:1), base_did, panel.id = ~id+period)
est2 = feols(f(y) ~ l(x1, -1:1), base_did, panel.id = ~id+period)
est3 = feols(l(y) ~ l(x1, 0:3), base_did, panel.id = ~id+period)
etable(est1, est2, est3, order = "f", drop = "Int")

## -----------------------------------------------------------------------------
# setting up the panel
pdat = panel(base_did, ~id + period)
# Now the panel.id argument is not required
est1 = feols(y ~ l(x1, 0:1), pdat)
est2 = feols(f(y) ~ l(x1, -1:1), pdat)
# You can use sub selections of the panel data
est_sub = feols(y ~ l(x1, 0:1), pdat[!pdat$period %in% c(2, 4)])
etable(est1, est2, est_sub, order = "f", drop = "Int")

## ---- include = FALSE---------------------------------------------------------
require_DT_ON()

## -----------------------------------------------------------------------------
library(data.table)
pdat_dt = panel(as.data.table(base_did), ~id+period)
# we create a lagged value of the variable x1
pdat_dt[, x1_l1 := l(x1)]
# Now 
pdat_dt[, c("x1_l1_fill0", "y_f2") := .(l(x1, fill = 0), f(y, 2))]
head(pdat_dt)

## ---- include = FALSE---------------------------------------------------------
require_DT_OFF()

## -----------------------------------------------------------------------------
base_lag = base_did
# we create a lagged value of the variable x1
base_lag$x1.l1 = lag(x1 ~ id + period, 1, base_lag)
head(base_lag)

## ---- include = FALSE---------------------------------------------------------
require_DT_ON()

## -----------------------------------------------------------------------------
library(data.table)
base_lag_dt = as.data.table(base_did)
# we create a lagged value of the variable x1
base_lag_dt[, x1.l1 := lag(x1 ~ id + period, 1)]

## ---- include = FALSE---------------------------------------------------------
require_DT_OFF()

## -----------------------------------------------------------------------------
# Generating data:
n = 1000
# x and y: two positive random variables
x = rnorm(n, 1, 5)**2
y = rnorm(n, -1, 5)**2
# E(z) = 2*x + 3*y and some noise
z = rpois(n, 2*x + 3*y) + rpois(n, 1)
base = data.frame(x, y, z)

## -----------------------------------------------------------------------------
result_NL = feNmlm(z~0, base, NL.fml = ~ log(a*x + b*y), NL.start = list(a=1, b=1), lower = list(a=0, b=0))

## -----------------------------------------------------------------------------
print(result_NL)

## -----------------------------------------------------------------------------
# the class of each observation
id = sample(20, n, replace = TRUE)
base$id = id
# the vector of fixed-effects
gamma = rnorm(20)**2
# the new vector z_bis
z_bis = rpois(n, gamma[id] * (2*x + 3*y)) + rpois(n, 1)
base$z_bis = z_bis

## -----------------------------------------------------------------------------
# we add the fixed-effect in the formula
result_NL_fe = feNmlm(z_bis~0|id, base, NL.fml = ~ log(2*x + b*y), NL.start = list(b=1), lower = list(b=0))
# The coef should be around 3
coef(result_NL_fe)
# the gamma and the exponential of the fixed-effects should be similar
rbind(gamma, exp(fixef(result_NL_fe)$id[as.character(1:20)]))


## ---- eval = FALSE------------------------------------------------------------
#  # Sample of results:
#  # 1 nthreads: 3.13s
#  system.time(fenegbin(Euros ~ log(dist_km)|Origin+Destination+Product+Year, trade, nthreads = 1))
#  # 2 nthreads: 1.82s
#  system.time(fenegbin(Euros ~ log(dist_km)|Origin+Destination+Product+Year, trade, nthreads = 2))
#  # 4 nthreads: 1.17s
#  system.time(fenegbin(Euros ~ log(dist_km)|Origin+Destination+Product+Year, trade, nthreads = 4))

