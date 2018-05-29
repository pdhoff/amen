#### The following changes have been made in version 1.4.1:

This update is a medium-sized step towards making 
`amen` more modular, so that individuals can build 
their own custom AME models. 

There are some medium-sized changes detailed below. 
The vignette sill works fine, but please let me know 
if this screws up something you've been working on, 
or if you have strong objections to some of the changes. 


1. Changes to `rbeta_ab_fc`:
  a. The function now takes as optional arguments a prior mean
     and a prior precision matrix for beta.    
  b. The default prior mean is zero and the default prior 
     precision is much smaller than in 1.3.   
  c. The function no longer takes as arguments all the items 
     computed from the design array. Instead, these items 
     are either attributes of the design array `X`, or if 
     they aren't, they are calculated in this function. Since 
     in most applications these items don't change during the
     MCMC algorithm, you should precompute these items 
     with the `design_array` command or the `precomputeX` 
     command. 

2. New function `precomputeX`: Precomputes various quantities 
   from `X` that will be repeatedly needed for the MCMC, and 
   returns a new `X` where the precomputed items are stored 
   as attributes. If you construct your own design array, you 
   should run `X<-precomputeX(X)` (unless you are using the 
   canned `ame` or `ame_rep` functions, as this do it 
   automatically). The precomputation is also done automatically
   if you construct your design array with the `design_array` 
   function. 

3. Changes to `design_array`: Derived quantities for the MCMC 
   are precomputed using the `precomputeX` function and stored 
   as attributes in the constructed design array.  

4. Changes to `ame`: There is now an additional parameter 
   `prior`, which is a list of hyperparameters (empty by default). 
   Parameters for which priors may be set include    
   `beta` (see `rbeta_ab_fc`), 
   `s2` (see `rs2_fc`)
   `Sab` (see `rSab_fc`)  
   `Suv` (see `rUV_fc`). 

5. Modified `simZ`, just changing some notation. 

6. Modified `rs2_fc` in two ways:
  a. The diagonal of the error matrix is now part of the update.  
  b. The function now takes optional prior parameter values.  
 
7. Added a secondary triadic dependence gof statistic. Also 
   modified `ame` and `plot.ame` to plot the additional statistic. 

8. Changed all `rZ` functions to appropriately update the diagonal. 

9. Changed all Wishart prior degrees of freedom to be two plus the number of parameters. 

10. In `ame`, use standard update for Sab whenever appropriate, and still use "specialty" updates in certain cases. 

11. Created a new function `rSab_fc` to update `Sab` given `a` and `b`. 

12. The `Sab` update functions, including `rSab_fc`, `raSab_bin_fc`,
    `raSab_cbin_fc` and `raSab_frn_fc` all take optional parameters 
    for the inverse-Wishart prior on `Sab`. 

13. For the monk example in the vignette, the model is now fit without 
    an intercept. There is not really any information about the 
    intercept from these data. 

14. The "model" parameter is now the "family" parameter and is now required  
    to prevent inadvertent fits of the normal model to binary data. 

15. The order of the "model/family" and "R" parameter have been changed in the 
     `ame` and `ame_rep` functions. 

16. The `rUV_fc` function now takes arguments for the prior distribution 
    over Psi, the covariance matrix of U and V. 

#### Things to do:

1. Allow prior means for (`a`,`b`,`U`,`V`) and a prior for `rho`. 
   The former will facilitate hierarchical or longitudinal modeling. 

2. Incorporate prior specification into the `ame_rep` wrapper function. 

3. Incorporate more full prior specifications for the symmetric case. 

4. Update `rho` using the marginal likelihood (integrating over the 
   additive row and column effects).   

5. Allow for user-defined GOF stats.  

6. Should remove the plotting from the `ame` wrapper, just call `plot.ame`. 

7. The functions `rZ_nrm_fc` and `rZ_cbin_fc` don't handle the dyadic 
   correlation in missing data properly if data are missing 
   asymmetrically.

8. Random reorderings can be avoided for FRN and RRL. 

9. Simplify the update for `U` and `V`. I don't think the complicated 
    method is saving much time. 


