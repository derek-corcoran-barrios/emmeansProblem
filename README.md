
<!-- README.md is generated from README.Rmd. Please edit that file -->

# emmeansProblem

<!-- badges: start -->
<!-- badges: end -->

The goal of emmeansProblem is to generate a reproducible example to
solve the issues within our model averaging within the emmeans package

first we load the needed r packages and the dataset

``` r
library(emmeans)
library(lme4)
library(MuMIn)
library(doParallel)

Data <- readRDS("Data.rds")
```

Now we fit a general model:

``` r
Model <- glmer(richness ~ aspect + elevation + 
                     initial_habitat  +
                       I(abs(year - 1)) +
                       I((year-1)^2) +
                     slope +
                     treatment:initial_habitat +
                     year:initial_habitat +
                     year:treatment + 
                     year:treatment:initial_habitat + 
                     (1 | block_no), family = poisson, data = Data, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
```

And we make a model selection (skip this step since it takes a while,
the “SelectRichness.rds” file is in
[this](https://github.com/derek-corcoran-barrios/emmeansProblem) github)

``` r
options(na.action = "na.fail")

library(doParallel)
cl <- makeCluster(4) 
registerDoParallel(cl)

clusterEvalQ(cl, library(lme4))
#> [[1]]
#> [1] "lme4"      "Matrix"    "stats"     "graphics"  "grDevices" "utils"    
#> [7] "datasets"  "methods"   "base"     
#> 
#> [[2]]
#> [1] "lme4"      "Matrix"    "stats"     "graphics"  "grDevices" "utils"    
#> [7] "datasets"  "methods"   "base"     
#> 
#> [[3]]
#> [1] "lme4"      "Matrix"    "stats"     "graphics"  "grDevices" "utils"    
#> [7] "datasets"  "methods"   "base"     
#> 
#> [[4]]
#> [1] "lme4"      "Matrix"    "stats"     "graphics"  "grDevices" "utils"    
#> [7] "datasets"  "methods"   "base"
clusterExport(cl, "Data")

Select <- MuMIn::pdredge(Model, extra = list(R2m = function(x) r.squaredGLMM(x)[1, 1], R2c = function(x) r.squaredGLMM(x)[1, 2]),fixed = ~YEAR:Treatment, cluster = cl)

stopCluster(cl)

saveRDS(Select, "SelectRichness.rds")
```

And now we Select the best models, I will do this twice, since the
outcome of the `subset` function will be used to show how the best model
does not have issues and the averaged model from `get.models` which is
the result I need is not working

``` r
Select <- readRDS("SelectRichness.rds")
Selected <- subset(Select, delta < 2)
SelectedList <- get.models(Select, delta < 2)
```

## Working with the best model works

As specified above the goal is to find if the treatments do yieald
differences by year 4, based on the model. So first we will show this
with the best model

``` r
BestModel <- get.models(Selected, 1)[[1]]

noise.emm <- emmeans(BestModel, pairwise ~ year  + initial_habitat +  initial_habitat:year + year:treatment, at = list(year = 4), data = Data)

pairs(noise.emm, simple = "treatment") |> 
  as.data.frame() |>  dplyr::filter(p.value < 0.05) |> 
  dplyr::arrange(initial_habitat, estimate) |> 
  dplyr::select(-SE, -df, -z.ratio) |> 
  knitr::kable()
```

| contrast                     | year | initial_habitat |   estimate |  p.value |
|:-----------------------------|-----:|:----------------|-----------:|---------:|
| PermanentExclosure - Control |    4 | Forest          | -0.2193592 | 3.06e-05 |
| PermanentExclosure - Control |    4 | Meadow          | -0.2193592 | 3.06e-05 |
| PermanentExclosure - Control |    4 | Rangeland       | -0.2193592 | 3.06e-05 |

## Working with the average model does not works

This does not work

``` r
AV <- model.avg(SelectedList, fit = TRUE)

noise.emm_av <- emmeans(AV, pairwise ~ year  + initial_habitat +  initial_habitat:year + year:treatment, at = list(year = 4), data = Data)

pairs(noise.emm_av, simple = "treatment") |> 
  as.data.frame() |>  dplyr::filter(p.value < 0.05) |> 
  dplyr::arrange(initial_habitat, estimate) |> 
  dplyr::select(-SE, -df, -z.ratio) |> 
  knitr::kable()
```

Giving the following error

``` r
Error in (mth$objs[[1]])(object, trms, xlev, grid, ...) : 
  Unable to match model terms
```
