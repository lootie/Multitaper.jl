
## Plots recipes for Multitaper.jl

For those using `Plots.jl`, we use `RecipesBase.jl` to quickly plot the output
structs given by the base functions. 

## Quick synopsis of recipes

One can plot univariate structs

* `MtSpec` structs, which are outputs of the `multispec` function. These are ordinary
  power spectrum estimates.

* `MtAcf`, `Mtacvf` structs, which are outputs of the `mtacf` and `mtacvf` function.
  These are autocorrelation or autocovariance estimates, respectively. 

* `MtCeps` structs, which are multitaper cepstra.

* `Demodulate` structs, which contain complex demodulates.

And multivariate structs

* `MtCoh` structs, which are outputs of the multivariate `multispec` function. These
  are coherence estimates. 

* Outputs of coherence calculations are often given as
  `Tuple{Array{MtSpec{C,T,P,1},Array{MtCoh,2},Nothing}`, and the corresponding recipe
will plot all of the spectra, coherences, and phases on a grid. 

* `MtCcvf` and `Mtccf` structs, which are outputs of the `mtccvf` and `mtccf` 
functions, respectively. These are plots of the multitaper cross-covariance function
or cross-correlation function.

### Univariate recipes

#### Plotting spectra

To plot a univariate spectrum estimate, use the recipe with signature

```
@recipe function plot(S::MtSpec; cross = true)
``` 

The following arguments are necessary

* `MtSpec` struct is the output of a call to `multispec` where only a single time
  series has been given, or two time series were given and the output desired is a 
  cross-spectrum.

With the vanilla call to the plot recipe, one will get the power spectrum plotted on
a logarithmic base 10 scale in the y-axis and appropriately labeled axes. If the
jackknifed variance field (jkvar) in the `MtSpec` struct contains values, then 95
percent confidence intervals will be shown as shaded bands around the spectrum. If
the phase field contains values, then the y-axis will say "Cross-Spectrum" instead of
"Spectrum".  

The optional keyword argument specifies

* `cross` shows a small cross-shaped mark in the corner of the plot having width
  twice the bandwidth, and height equal to the expected jackknifed variance. 

#### Plotting autocorrelations

To plot a multitaper estimate of autocorrelation, use the recipe with signature

```
@recipe function acfplt(A::MtAcf)
```

The following arguments are necessary

* `MtAcf` struct is the output of a call to `mtacf` with a single time series given. 

The resulting plot will be a stem plot labeled as lags on the x-axis and
autocorrelations on the y-axis. 

### Plotting autocovariances

To plot a multitaper estimate of autocovariance, use the recipe with signature

```
@recipe function acvfplt(A::MtAcvf)
```

The following arguments are necessary

* `MtAcvf` struct is the output of a call to `mtacvf` with a single time series given. 

The resulting plot will be a stem plot labeled as lags on the x-axis and
autocovariances on the y-axis. 

### Cepstrum plot

To plot the cepstrum, use the recipe with signature 

```
@recipe function cepsplt(A::MtCeps)
```

The following arguments are necessary

* `MtCeps` which is the output of a `mtceps` function call.

The resulting plot will show the cepstrum coefficient as a stem plot in terms of
lags. 


#### Demodulate plot

To plot the complex demodulate, use the recipe with signature 

```
@recipe function demodplt(cdm::Demodulate)
```

The following arguments are necessary

* `cdm` which is the output of a `Demodulate` function call.

The resulting plot will show the complex demoodulate in two stacked panels, one for
the magnitude and one for the phase of the complex quantity. 

### Multivariate recipes

#### Coherence plotting

To plot a multitaper estimate of the coherence, use

```
@recipe function plot(C::MtCoh; phase = false, sigMax = 0)
```

The following arguments are necessary

* `MtCoh` which is the output of a `multispec` call with more than one time series as
  input and coherence as output. 

The resulting plot shows the magnitude squared coherences on a scale from 0 to 1, if
no jackknifing has been done. If jackknifing has been done, the magnitude squared
coherence is shown on a transformed scale (inverse hyperbolic tangent) and 95 percent
confidence intervals will be shown as a shaded band. 

The following keyword arguments can be set

* `phase` when set to true produces a second subplot below the first which shows the
  (unwrapped) values of the phase estimate.

* `sigMax` is the number of significance lines to put on the plot. These lines are at
  (1.0 - 10^{1:sigMax}). For example, if one selects the default, namely `sigMax = 4`,
one obtains lines at 0.9, 0.99, 0.999, 0.9999. 

#### Comparing coherences

To plot multiple estimates of the coherence together on the same axes for comparison, 
use

```
@recipe function plot(C::Matrix{MtCoh}; phase = false, jk = false, sigMax = 0)
```

The following arguments are necessary

* `MtCoh` is as above.

The following keyword arguments can be set

* `phase` and `sigMax` are as above.

* `jk` plots confidence intervals when set to true.

#### Coherences and spectra

To plot the coherences and spectra from, say, a call to `multispec` where there are
two or more time series, one uses the recipe with the following signature:

```
@recipe function
specCohMatrix(Spec::Tuple{Array{MtSpec{C,T,P},1},Array{MtCoh,2},Nothing}; 
                                sigMax = 0) where C where T where P
```

One obtains a plot with $p\times p$ panels (if $p$ is the number of time series)
where spectra are shown in red on the diagonals, MSC are on the upper triagonal panels of the
plot, and phases are on the lower triagonal panels of the plot. 

The required input is 

* `Spec` which is of type `Tuple{Array{MtSpec{C,T,P},1},Array{MtCoh,2},Nothing}`;
  which is the result one gets when one calls `multispec` on two or more time series. 

The keywork argument

* `sigMax` is as above.


#### Note on coherence plotting

Note that at the time of this writing, it was impossible to not require Plots.jl as a
dependency when using an advanced layoout, so some bonus plotting routines are in
Multitaper.jl/src/PlotsRecipes/Coherence_Plot.jl. This code can be included on the
fly for an even more advanced layout. The layout one gets has three axes, one is the
transformed MSC, one is the MSC in units from 0 to 1, and the third axes is the
significance if the true MSC were zero. 

The main function in this piece of code has the
following signature:

```
@recipe function f(h::MtCoh; siglines = true, msclines = true, sigMax = 4, legtext = false, 
        force_ylims = nothing, mscaxis = true, sigaxis = true, jk = true)
```

The `MtCoh` type input is explanatory (it is as above), but one has the additional
keyword arguments as follows

* `siglines` allows one to put a selected number of horizontal lines on the plot
  (specifically, `sigMax` number of lines are added) 
  indicating the significance of a magnitude squared coherence value under the null
hypothesis that the two series have zero coherence at that frequency. 

* `msclines` allows one to put horizontal lines at various levels of the msc

* `sigMax` is the number of significance lines to put on the plot. These lines are at
  (1.0 - 10^{1:sigMax}). For example, if one selects the default, namely `sigMax = 4`,
one obtains lines at 0.9, 0.99, 0.999, 0.9999. 

* `legtext` puts text in the legend.

* `force_ylims` gives the y-limits for the transformed axis, if desired. 

* `mscaxis` when toggled to false will delete the MSC axis on the scale from 0 to 1. 

* `sigaxis` when toggled to false will delete the significance axis.

* `jk` indicates whether jackknife confidence intervals will be drawn. 
 
#### Cross-covariance and Cross-correlation plots

The following two recipe signatures:

```
@recipe function ccvfplot(A::MtCcvf)
```

and 

```
@recipe function ccfplot(A::MtCcf)
```

will show either the cross-covariance estimate (if the input is of type `MtCcvf`) or
the cross-correlation estimate (if the input is of type `MtCcf`) as a stem plot.

