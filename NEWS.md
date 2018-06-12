# amen 1.4.2 

#### Changes 

* Added a function `rrho_fc` that updates the dyadic correlation `rho` from its full conditional distribution, leading to faster-mixing Markov chains. 

* Changed `rZ_fc_bin` to avoid numerical instabilities that lead to infinities. 

* Changed `ame_rep` so that the GOF statistics are calculated for each rep. 

* Changed the default priors. For fitting non-normal models to sparse 
data, the priors on the covariance matrices `Sab` and `Suv` have a smaller
scale. Also changed the prior on `beta` to be a g-prior on the non-intercept 
coefficients but a more diffuse prior on the intercept. 

* Changed the name of the covariance of `U` and `V` from  `Psi` to `Suv`. 

* Moved the update of `Suv` outside of the update for `U` and `V`. 

* Added an argument `offset` to most Gibbs sampling functions. The function documentation indicates what things should be subtracted off (offset) for each update. 

* Fixed a bug in plot.ame that created the wrong number of panels when `p=0`.

* Added a secondary plotting parameter to  `circplot`.




# amen 1.4.1

This update from 1.3 is a medium-sized step towards making 
`amen` more modular, so that individuals can build 
their own custom AME models. 
There are some medium-sized changes detailed below. 
The vignette sill works fine, but please let me know 
if this screws up something you've been working on, 
or if you have strong objections to some of the changes. 

#### Changes    
     
* Change to `rbeta_ab_fc`: The function now takes as optional 
  arguments a prior mean
  and a prior precision matrix for beta.    
  The default prior mean is zero and the default prior
  precision is much smaller than in 1.3.   

* Change to `rbeta_ab_fc`: The function no longer takes as arguments all the 
      items  computed from the design array. Instead, these items 
      are either attributes of the design array `X`, or if 
      they aren't, they are calculated in this function. Since
      in most applications these items don't change during the
      MCMC algorithm, you should precompute these items
      with the `design_array` command or the `precomputeX`
      command. 

* New function `precomputeX`: Precomputes various quantities 
   from `X` that will be repeatedly needed for the MCMC, and 
   returns a new `X` where the precomputed items are stored 
   as attributes. If you construct your own design array, you 
   should run `X<-precomputeX(X)` (unless you are using the 
   canned `ame` or `ame_rep` functions, as this do it 
   automatically). The precomputation is also done automatically
   if you construct your design array with the `design_array` 
   function.  

* Change to `design_array`: Derived quantities for the MCMC 
   are precomputed using the `precomputeX` function and stored 
   as attributes in the constructed design array.  

* Changes to `ame`: There is now an additional parameter 
   `prior`, which is a list of hyperparameters (empty by default). 
   Parameters for which priors may be set include    
   `beta` (see `rbeta_ab_fc`), 
   `s2` (see `rs2_fc`),
   `Sab` (see `rSab_fc`),
   `Suv` (see `rUV_fc`).  

* Modified `simZ`, just changing some notation. 

* Modified `rs2_fc` in two ways:     
  the diagonal of the error matrix is now part of the update;
  the function now takes optional prior parameter values.  
 
* Added a secondary triadic dependence gof statistic. Also 
   modified `ame` and `plot.ame` to plot the additional statistic. 

* Changed all `rZ` functions to appropriately update the diagonal. 

* Changed all Wishart prior degrees of freedom to be two plus the number of parameters. 

* In `ame`, use standard update for Sab whenever appropriate, and still use "specialty" updates in certain cases. 

* Created a new function `rSab_fc` to update `Sab` given `a` and `b`. 

* The `Sab` update functions, including `rSab_fc`, `raSab_bin_fc`,
    `raSab_cbin_fc` and `raSab_frn_fc` all take optional parameters 
    for the inverse-Wishart prior on `Sab`. 

* For the monk example in the vignette, the model is now fit without 
    an intercept. There is not really any information about the 
    intercept from these data. 

* The "model" parameter is now the "family" parameter and is now required  
    to prevent inadvertent fits of the normal model to binary data. 

* The order of the "model/family" and "R" parameter have been changed in the 
     `ame` and `ame_rep` functions. 

* The `rUV_fc` function now takes arguments for the prior distribution 
    over Psi, the covariance matrix of U and V. 

* The `zscores` function now takes an optional argument for dealing with ties. 

### To do list:

* Allow prior means for (`a`,`b`,`U`,`V`) and a prior for `rho`. 
   The former will facilitate hierarchical or longitudinal modeling. 

* Incorporate prior specification into the `ame_rep` wrapper function. 

* Incorporate more full prior specifications for the symmetric case. 

* Update `rho` using the marginal likelihood (integrating over the 
   additive row and column effects).   

* Allow for user-defined GOF stats.  

* Should remove the plotting from the `ame` wrapper, just call `plot.ame`. 

* The functions `rZ_nrm_fc` and `rZ_cbin_fc` don't handle the dyadic 
   correlation in missing data properly if data are missing 
   asymmetrically.

* Random reorderings can be avoided for FRN and RRL. 

* Simplify the update for `U` and `V`. I don't think the complicated 
    method is saving much time. 


