# brglm2 0.9.4

## Improvements, updates and additions

* Implemented the `mdyplFit()` method for the `glm()` function, which
  estimates logistic regression models using maximum
  Diaconis-Ylvisaker prior penalized likelihood, and associated
  methods. Corrections to estimates, estimated standard errors, z
  statistics, and so on, can be performed using the developments
  Sterzinger & Kosmidis (2024, arXiv:2311.07419v2) using
  `hd_correction = TRUE` in the `summary()` and `confint()` methods
  for `"mdyplFit"` objects.
  
* Documentation updates

# brglm2 0.9.3

## Bug fixes

* Added `add1()` and `drop1()` methods for `brglmFit` objects, so that
  `step()` does not default to the methods for `glm` objects (which
  would return nonsense results); thanks to jamiahuswalton@github
  (https://github.com/ikosmidis/brglm2/issues/33) for reporting the
  issue.

* Sane staring values, graceful failing, updates to optimization
  procedure, and improvements to methods for `brnb()`; thanks to
  cperk@github (https://github.com/ikosmidis/brglm2/issues/31) for
  reporting issues with `brnb()` when infinite ML estimates are
  encountered.

## Improvements, updates and additions

* Documentation updates

# brglm2 0.9.2

## Improvements, updates and additions

* Convergence of the `brglm_fit` iterations is now determined if the L^Inf norm of the step size (rather than the L^1 as it was previously) of the quasi-Fisher scoring procedure is less than `epsilon` (see `?brglm_control` for the definition of `epsilon`). This is more natural as `epsilon` then determines directly the precision of the reported estimates and does not depend on their number.

* `brglm_control()` now checks that the supplied value of `max_step_factor` is numeric and greater or equal to `1`. If not, then it is set to the default value of `12`.

* Vignette updates


# brglm2 0.9.1

## Improvements, updates and additions

* Added the `enzymes` and `hepatitis` data sets (from the [*pmlr*](https://cran.r-project.org/package=pmlr)) to support examples and tests.

# brglm2 0.9.0

## New functionality

* The `expo()` method for `brglmFit` and `glm` objects estimates the exponential of parameters of generalized linear models with maximum likelihood or various mean and median bias reduction methods (see `?expo` for details). The `expo()` method is particularly useful for computing (corrected) estimates of the multiplicative impact of a unit increase on a covariate on the mean of a Poisson log-linear model (`family = poisson("log")` in `glm()`) while adjusting for other covariates, the odds ratio associated with a unit increase on a covariate in a logistic regression model (`family = binomial("logit")` in `glm()`) while adjusting for other covariates, the relative risk associated with a unit increase on a covariate in a relative risk regression model (`family = binomial("log")` in `glm()`) while adjusting for other covariates, among others.

## Bug fixes

* Fixed a bug where the dispersion in the resulting object would not be transformed even if `transformation != "identity"` when `type` is `ML` or `AS_median` or `AS_mixed`.

## Improvements, updates and additions

* Moved unit tests to [**tinytest**](https://cran.r-project.org/package=tinytest).

* Moved documentation to markdown markup through [**roxygen2**](https://cran.r-project.org/package=roxygen2).

* New vignette titled "Estimating the exponential of regression parameters using **brglm2**", to demonstrate the `expo()` method.

* Various documentation fixes.

# brglm2 0.8.2

## Improvements, updates and additions

* Housekeeping.
* Removed lpSolveAPI from imports.

# brglm2 0.8.1

## Bug fixes

* Fixed a bug when predicting from `bracl` objects with non-identifiable parameters.

## Improvements, updates and additions

* Work on output consistently from `print()` methods for `summary.XYZ`
  objects; estimator type is now printed and other fixes.
  
* Enriched warning when algorithm does not converge with more informative text.
  
* Documentation fixes and updates

# brglm2 0.8.0
 
## New functionality

* `brnb()` allows fitting negative binomial regression models using
  implicit and explicit bias reduction methods. See vignettes for a
  case study.
* `simulate()` method for objects of class `brmultinom` and `bracl`
* `ordinal_superiority()` method to estimate Agresti and Kateri
  (2017)'s ordinal superiority measures, and compute bias corrections
  for those.

## Bug fixes

* Fixed a bug that would return an error when `Wald.ratios = TRUE` in
  `summary.brmultinom`.
* Fixed bug in `vcov.bracl` that would return an error if the
  `"bracl"` object was computed using `bracl()` with `parallel = TRUE`
  and one covariate.
* Fixed a bug in `bracl()` related to the handling or zero weights
  that could result in hard-to-traceback errors.
* Fixed a bug in `bracl()` that could cause errors in fits with one
  covariate.
* `brglmFit()` iteration returns last estimates that worked if
  iteration fails.

## Improvements, updates and additions

* Documentation and example updates.

# brglm2 0.7.1

## Bug fixes

* Fixed bug where `confint()` was not returning anything when applied
  to objects of class `brmultinom`.
* Fixed bug where and error could result when the `control` `glm()`.
  argument was specified using the output from `brglmControl()` or
  `brglm_control()`.

## New functionality
* added `check_aliasing` option in `brglmControl()` to tell
  `brglm_fit()` to skip (`check_aliasing = TRUE`) or not
  (`check_aliasing = FALSE`) rank deficiency checks (through a QR
  decomposition of the model matrix), saving some computational effort.

## Improvements, updates and additions
* updated DOI links in documentation and some http -> https fixes.

# brglm2 0.7.0

## Bug fixes
* Fixed bug that resulted in `NA` coefficients when `brglmFit()` was
  called with a vector `x` or an `x` with no column names.

## New functionality
* `confint` method for `brmulitnom` objects

## Improvements, updates and additions
* Updated reference to [Kenne Pagui et al (2017)](https://doi.org/10.1093/biomet/asx046).q
* Updated reference to [Kosmidis and Firth (2020)](https://doi.org/10.1093/biomet/asaa052).
* Fixed issues with references.
* Updated documentation.

# brglm2 0.6.2

## Improvements, updates and additions
* `vcov.brglmFit()` now uses `vcov.summary.glm()` and supports the
  `complete` argument for controlling whether the variance covariance
  matrix should include rows and columns for aliased parameters.
* Deprecated `detect_sepration()` and `check_infinite_estimates()`, which
  will be removed from **brglm2** at version 0.8. New versions of
  `detect_sepration()` and `check_infinite_estimates()` are now maintained
  in the
  [**detectseparation**](https://cran.r-project.org/package=detectseparation)
  R package.
* Fixed typos in NEWS.

# brglm2 0.6.1
## Bug fixes
* Fixed bug in AIC reported by `print.summary()` for `brmultinom` and
  `bracl` objects.
* `detect_separation()` now handles one-column model matrices correctly.

## Improvements, updates and additions
* Documentation improvements and typo fixes.

# brglm2 0.6
## New functionality
* `brglmFit()` can now do maximum penalized likelihood with powers of
  the Jeffreys prior as penalty (`type = "MPL_Jeffreys`) for all
  supported generalized linear models. See the help files of
  `brglmControl()` and `brglmFit()` for details.

## Improvements, updates and additions
* Documentation updates and improvements.
* Updated vignettes to include maximum penalized likelihood with
  powers of the Jeffreys prior as penalty.
* New examples in `?brglmFit`.

# brglm2 0.5.2
## Bug fixes
* `print.brmultinom()` is now exported, so `bracl` and `brmultinom`
  objects print correctly.

## New functionality
* Added `response_adjustment` argument in `brglmControl()` to allow
  for more fine-tuning of the starting values when `brglmFit()` is
  called with `start = NULL`.

## Improvements, updates and additions
* Documentation updates and improvements.
* Added Kosmidis et al (2019) in the description file.
* Added tests for `brglmControl()`.

# brglm2 0.5.1

## Improvements, updates and additions
* Fixed typos in vignettes and documentation.
* Added ORCHID for Ioannis Kosmidis in DESCRIPTION.

# brglm2 0.5.0
## Bug fixes
* `brglmFit()` now works as expected with custom link functions (mean
  and median bias reduction).
* `brglmFit()` respects the specification of the transformation
  argument in `brglmControl()`.
* Fixed bug in the computation of the QR decomposition under aliasing
  in `brglmFit()`.
* Other minor bug fixes and performance improvements.
* Protection against use of `quasi()`, `quasibinomial()` and
  `quasibinomial()` families and documentation update.

## New functionality
* Added `bracl()` for fitting adjacent category logit models for ordinal
  responses using maximum likelihood, mean bias reduction, and median
  bias reduction and associated methods (`logLik`, `summary` and so
  on).
* Added `predict()` methods for `brmultinom` and `bracl` objects.
Added `residuals()` methods for `brmultinom` and `bracl` objects
(residuals of the equivalent Poisson log-linear model)
* Added the `mis()` link functions for accounting for
  misclassification in binomial response models (Neuhaus, 1999,
  Biometrika).

## Improvements, updates and additions
* Improved `summary()` method for `brmultinom` objects.
* Better starting values for null fits.
* Added references to [arxiv:1804.04085](https://arxiv.org/abs/1804.04085) in
  documentation.
* Updated reference to [Kenne Pagui et al (2017)](https://doi.org/10.1093/biomet/asx046).

# brglm2 0.1.8
## Improvements, updates and additions
* Improved documentation examples.
* Removed warning about observations with non-positive weights from brmultinom.
* Updated email address for Ioannis Kosmidis in brglmFit.

## Bug fixes
* brmultinom returns a fitted values matrix that respects the
  dimension of data.
* Fixed bug on condition for `NA` dispersion for models with `0` df
  resid.

# brglm2 0.1.7

## Improvements, updates and additions
* Eliminated errors from markdown chunks in multinomial vignette.

# brglm2 0.1.6

## Bug fixes
* Compatibility with new version of enrichwith.

## Improvements, updates and additions
* New email for Ioannis Kosmidis.

# brglm2 0.1.5

## Bug fixes

## New functionality
* Added `type = AS_mixed` as an option to use **mean-bias reducing
  score functions** for the regression parameters and **median-bias
  reducing score functions** for the dispersion in models with unknown
  dispersion.
* `check_infinite_estimates()` now accepts `brmultinom` objects.
* Added `singular.ok` argument to `brglmFit()` and
  `detect_separation()` methods in line with the update of
  `glm.fit()`.

## Improvements, updates and additions
* less strict tolerance in `brglm_control()`.
* Updates to help files.
* Fixed typos in iteration vignette.
* Added URL and bugreports in Description.
* Added new tests.

# brglm2 0.1.4

## Bug fixes
* `brglmControl()` is now exported.
* `slowit` did nothing; now included in iteration.

## New functionality
* The `detect_separation()` method for the `glm()` function can be used
  to check for separation in binomial response settings without
  fitting the model. This relies on a port of Kjell Konis'
  `safeBinaryRegression:::separator()` function (see ?detect_separation).
* **brglm2** provides estimation via **median-bias reducing score
  functions** with `type = "AS_median"`.
* **brglm2** provides camel and underscored aliases for basic methods
  (`brglmFit()`, `brglm_fit()`, `detectSeparation()`,
  `detect_separation()`, `brglm_control()`, `brglmControl()`,
  `detectSeparationControl()`, `detect_separation_control()`,
  `checkInfiniteEstimates()`, `check_infinite_estimates()`).

## Improvements, updates and additions
* Minor enhancements in the codebase.
* The inverse expected information matrix is computed internally using
  `cho2inv()`.
* Internal changes to have more meaningful variable names.
* Renamed detect_infinite* to check_infinite.

# brglm2 0.1.3

## Improvements, updates and additions
* Fixed typo in $f_{Y_i}(y)$ in iteration vignette (thanks to Eugene
  Clovis Kenne Pagui for spotting),

# brglm2 0.1.2

* First release.



