<!-- badges: start -->
  [![R-CMD-check](https://github.com/mengyu-he/MIDASim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mengyu-he/MIDASim/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

# Overview

The **MIDASim** package (v2.0) provides a flexible simulator for microbiome sequencing data, preserving both the marginal distributions and taxon–taxon correlation structure of a user-supplied template dataset. You can optionally introduce treatment or covariate effects via model parameters. Version 2.0 further improves the simulation procedure when covariates are incorporated.

## Required packages 

MIDASim relies on the following CRAN packages:

- **psych**
- **MASS**
- **vegan**
- **pracma**
- **scam**

These packages are all available on CRAN.

## Installation

Install from the local file (the .tar.gz file), which you can download from https://github.com/mengyu-he/MIDASim.

```{r}
devtools::install_local("MIDASim_2.0.tar.gz", dependencies = TRUE)
```

A github installation is also available via

```{r}
devtools::install_github("mengyu-he/MIDASim")
```


## Open the Vignette in R

```{r}
browseVignettes("MIDASim")
```

or

```{r}
vignette("MIDASim_vignette", package = "MIDASim")
```


# MIDASim: Microbiome Data Simulator

There are three main functions in MIDASim package: 

| Function | Purpose |
|----------|---------|
| `MIDASim.setup()`  | Extracts taxon‑specific parameters (mean relative abundances, marginal probability of zeros, etc.) from a template count table and stores them in a simulation object. |
| `MIDASim.modify()` | Adds user‑specified effects (e.g., covariate shifts, changes in library sizes) to the stored parameters. |
| `MIDASim()`        | Generates realistic count tables under the current parameter configuration. |

## Template data description

To illustrate the workflow, we use a filtered inflammatory bowel disease (IBD) 16S dataset from the Human Microbiome Project 2 (HMP2) [1].
The processed otu table can be loaded directly: 

```{r}
data("count.ibd")
```

This HMP2 data were modified through package **HMP2Data** on Bioconductor by filtering out samples with library sizes smaller than 3000 and taxa that appear in less than 2 samples. The filtered data have 146 rows (samples) and 614 columns (taxa).

```{r, eval = FALSE}
library(HMP2Data)
IBD16S()

count.ibd <- t(IBD16S_mtx[, colSums(IBD16S_mtx)>3000] )
count.ibd <- count.ibd[, colSums(count.ibd>0)>1]

dim(count.ibd)
```


Additional examples datasets that come with the package include a vaginal microbiome dataset from HMP2 [1] and a upper-respiratory-tract microbiome dataset [2].

```{r}
data("count.vaginal")
data("throat.otu.tab")
```

## Setup the MIDASim model

The function `MIDASim.setup()` fits the underlying MIDASim model and extracts all quantities needed to simulate new count tables from a template dataset. The function offers two modeling options for relative abundances: nonparametric and parametric. Both approaches utilize a two-step procedure in which the first step models the presence-absence relationship between taxa and the second step models the non-zero values in the microbiome community. These two approaches differ by the way how the second step marginal distributions are modeled. The nonparametric approach generates data that resemble the template data better than the parametric approach. However, the parametric approach offers modifications of the template data in more controlled way.

### Non-parametric mode

In nonparametric mode (the default), a rank based approach is used to obtain the relative abundances for taxa that are considered to be non-zero in the first step. 

To fit MIDASim with nonparametric mode to the `count.ibd` data, use the following:

```{r}
otu.tab <- count.ibd
count.ibd.setup <- MIDASim.setup(otu.tab, 
                                 mode = 'nonparametric', 
                                 n.break.ties = 10)
```

The argument `n.break.ties` specifies the number of replicates used to break ties when ranking relative abundances for each taxon, and the default is `n.break.ties = 100`. 

### Parametric mode

In parametric mode, MIDASim fits a generalized gamma model to the relative abundance data of each taxon. Generalized gamma model is a three-parameter distribution in the location-scale family that was proposed for analyzing right-censored survival data. 

Denote the simulated relative abundance for \(i\)-th subject and \(j\)-th taxon as \(\tilde\pi_{ij}\). We define **survival time**:
$$
\tilde t_{ij} =
\begin{cases}
\dfrac{1}{\tilde\pi_{ij}}, & \text{if } \tilde\pi_{ij} > 0,\\[6pt]
N_i,                     & \text{if } \tilde\pi_{ij} = 0.
\end{cases}
$$
which corresponds to treating \(\tilde t_{ij}\) as right-censored when \(\pi_{ij}<\tfrac{1}{N_i}\). The generalized gamma model then assumes \(\tilde t_{ij}\) has the distribution
$$
\ln(\tilde t_{ij}) = -\mu_j + s_j\,\sigma_j\;\omega_{ij}\,.
$$
where $e^{\omega_{ij}}$ follows a gamma distribution with shape parameter $k_j = 1/|Q_j|$ and scale parameter 1, and $s_j=\text{sign}(Q_j)$}. 

We estimate the three parameters $\mu_j, \sigma_j, Q_j$ in `MIDASim.setup()` using the relative abundances $\pi_j$ in the template.

To fit MIDASim in parametric mode to the count.ibd data, use the following command:

```{r}
count.ibd.setup <- MIDASim.setup(otu.tab, mode = 'parametric')
```

The resulting list of *MIDASim.setup()* contains

| Value                          | Description                                         |
| ------------------------------ | --------------------------------------------------- |
| `tetra.corr`                   | tetrachoric correlation (presence/absence)          |
| `corr.rel.corrected`           | Pearson correlation among non‑zeros                 |
| `taxa.1.prop`                  | proportion of non‑zeros per taxon                   |
| `sample.1.prop`                | proportion of non‑zeros per sample                  |
| `mean.rel.abund`               | mean relative abundance per taxon                   |
| `rel.abund.1`                  | vector of non‑zero abundances (each taxon)          |
| `mu.est`, `sigma.est`, `Q.est` | generalized gamma parameters (parametric mode only) |


## Modify the fitted MIDASim model

The function `MIDASim.modify()` prepares (and optionally tweaks) the parameters returned by `MIDASim.setup()` before data generation. Note that this is a **mandatory** step even if no adjustment is made. The following arguments can be modified ($n$ = target number of samples, $p$ = the number of taxa):

| Argument        | What it controls | Expected type          |
|-----------------|------------------|------------------------|
| `lib.size`      | Target library sizes for each sample.  If `NULL`, original depths are reused. | vector of length $n$ |
| `mean.rel.abund`| Taxon‑specific mean relative abundances. | vector of length $p$ |
| `gengamma.mu`   | Location shifts \(\mu_j\) for the generalized‑gamma fit. | numeric vector of length $p$ |
| `sample.1.prop` | Desired proportion of non‑zero cells per sample. | numeric vector of length $n$ or single number |
| `taxa.1.prop`   | Desired proportion of non‑zero cells per taxon.  | numeric vector of length $p$ |
| `individual.rel.abund`   | Individual‑by‑taxon relative‑abundance matrix, overriding `mean.rel.abund`.  | numeric $n$-by-$p$ matrix |


### No modification at all

The default is no additional adjustment is wanted. The number of samples and their target library sizes for the simulated data are the same as that of the template. The following code is an example.

```{r}
count.ibd.modified <- MIDASim.modify(
  count.ibd.setup,
  lib.size               = NULL,   # keep original depths
  mean.rel.abund         = NULL,   # keep original means
  gengamma.mu            = NULL,   # keep original location parameters
  sample.1.prop          = NULL,   # keep original sparsity per sample
  taxa.1.prop            = NULL,   # keep original sparsity per taxon
  individual.rel.abund   = NULL,   # keep original individual relative abundances
)
```

Next, we will explain how to modify the data features to accommodate a wide range of simulation scenarios, separately for the nonparametric and parametric modes.

### Model modification in parametric mode

When `MIDASim.setup()` is called with `mode = "parametric"`, users can modify:

- **Sample size and library sizes**  
- **Taxon‑level location parameters** (\(\mu_j\))  
- **Mean or individual‑level relative abundances**


Changes in the parameters of the parametric model (including library sizes) imply coordinated changes in other quantities. After either modification of the parameters, we predict the probability of being non-zero of $i$-th subject and $j$-th taxon by 

$$
P\bigl(\widetilde{Z}_{ij} = 1\bigr)
\;=\;
F_{j}\bigl(N_{i}\;;\;\widehat{\mu}_{ij},\;\widehat{\sigma}_{j},\;\widehat{Q}_{j}\bigr)\,.
$$

to obtain the number of non-zero cells in each subject $Z_{i\cdot }$ and that in each taxon $Z_{\cdot j}$.

#### Modify the number of samples and their library sizes

To change the library sizes, set the `lib.size` argument equal to a vector with the target library sizes.  The length of this vector becomes the new number of samples.  The advantage of the parametric mode is that it handles coordinated changes automatically. As a result, when changing library sizes in parametric mode, there's no need for the user to provide the marginal totals of non-zero cells explicitly.

For example, to generate data having library sizes that are uniformly distributed between 1,000 and 10,000 while changing the number of samples to 700, we execute the code

```{r}
new.lib.size <- sample(1000:10000, size=700, replace=TRUE)

count.ibd.modified <- MIDASim.modify(
  count.ibd.setup, 
  lib.size = new.lib.size   # other arguments left NULL
  )
```


#### Modify taxa features by mean relative abundances or location parameters

Specifying new values for mean relative abundances (`mean.rel.abund`) or location parameters (`gengamma.mu`) can modify taxa features. 

You can alter taxon‑level characteristics in either of two equivalent ways:

- `mean.rel.abund` provide supply a new vector of target mean relative abundances. MIDASim will calculate the corresponding location parameters of the generalized‑gamma model, then updates all downstream probabilities and sparsity counts.
- `gengamma.mu` - provide the mean specification and shift $\mu_j$ directly. This is mathematically identical to changing the means on the log scale.

Because the two arguments modify the same quantity, specify one or the other, not both.

For example, if we want to increase $log(p_j)$ by $\beta=0.1$ for the first 10 taxa (e.g. modify these taxa according to a compositional model), the following is the code for generating a set of target relative abundances as described, 

```{r}
beta <- 0.1
new.mean.rel.abund <- count.ibd.setup$mean.rel.abund
new.mean.rel.abund[1:10] <- exp(beta) * new.mean.rel.abund[1:10]
new.mean.rel.abund <- new.mean.rel.abund / sum(new.mean.rel.abund)  # renormalize
```

then the target taxa relative abundances can be used in *MIDASim.modify()* as the following,

```{r}
count.ibd.modified <- MIDASim.modify(
  count.ibd.setup, 
  mean.rel.abund = new.mean.rel.abund
  )
```

Alternatively, the users can also directly change the location parameters $\mu_j$ by

$$\mu_j \rightarrow \mu_j + \beta $$

```{r}
count.ibd.modified <- MIDASim.modify(
  count.ibd.setup, 
  gengamma.mu = count.ibd.setup$mu.est + beta
  )
```


#### Modify taxa features by individual‑specific relative abundances

MIDASim v2.0 allows you to pass an individual-by-taxon matrix of expected relative abundances (individual.rel.abund). This is the most straightforward way to encode subject‑level covariate effects (binary or continuous, one or many) in a single step.

Illustration with two covariates:

- $Y1$ is binary: 50 controls ($Y=0$) and 50 cases ($Y=1$). Cases get a uniform 0.1 log-fold increase in the first 10 taxa.
- $Y2$ is continuous: 100 draws from a standard normal. Each unit increase in $Y2$ affects the second 10 taxa (11-20) by taxon-specific log-fold changes sampled from uniform distribution, Unif(−1, 1).


```{r}
Y1 <- rep(c(0, 1), each = 50)  # binary case-control
Y2 <- rnorm(100)   # continuous
Y.all <- cbind(Y1, Y2)

n.taxa <- length(count.ibd.setup$mean.rel.abund)
beta1 <- beta2 <- rep(0, n.taxa) # effect on each taxon
beta1[1:10] <- rep(0.1, 10)
beta2[11:20] <- runif(10, -1, 1)
beta.all <- cbind(beta1, beta2)

new.individual.rel.abund <- matrix(count.ibd.setup$mean.rel.abund, 
                                   nrow = 100, ncol = n.taxa, byrow = T)
new.individual.rel.abund <- exp(Y.all %*% t(beta.all))  * new.individual.rel.abund  # introduce signals
new.individual.rel.abund <- new.individual.rel.abund/rowSums(new.individual.rel.abund)  # normalize
```

Feed this matrix to *MIDASim.modify()* as the following,

```{r}
count.ibd.modified <- MIDASim.modify(
  count.ibd.setup, 
  individual.rel.abund = new.individual.rel.abund
  )
```

The resulting setup object now carries subject‑specific compositions that reflect the desired subject-level signals and is ready for data generation with `MIDASim()`.


### Model modification in nonparametric mode (Not recommended)

We provide options for modifying the model under the nonparametric mode. However, it's important to note that such modifications may not be as controlled as in parametric mode. For example, if we wish to change the library sizes of certain observations or the relative abundances of various taxa, it is not clear how the proportion of non-zero taxa should change in the nonparametric mode.

In nonparametric mode, users are able to modify the sample size, library sizes, taxa relative abundances, taxa proportion of non-zero cells. 

#### Modify the number of samples and their library sizes

To change the library sizes, set the `lib.size` argument equal to a vector with the target library sizes. The number of samples can be modified implicitly by specifying library sizes with a vector length equal to the desired sample size.

Note that, if we wish to change the library sizes, the proportion of non-zero taxa would also change. Typically, larger library sizes are accompanied by fewer zero cells. However, since the model is nonparametric, there is no simple parametric relationship between the changes in zero cells and relative abundances. *MIDASim.modify()* requires the user to specify the marginal totals of non-zeros together when a change in `lib.size` is made. This can be achieved by providing the marginal proportions using the `sample.1.prop` and `taxa.1.prop` arguments within the function."

An example way determine the appropriate number of non-zero cells for the new library size is to fit a Shape Constrained Additive model (SCAM) model 
$$ log10(Z_{i.}) = f\left(log_{10}(N_{i})\right)+\epsilon_i$$
where $f$ is a monotone smoothing spline, $Z_{i\cdot}$ is the number of non-zero cells in each sample, $\epsilon_i$ is the error term. We could adjust the number of non-zero cells for each taxon in a way that satisfies the constraint that the product of `sample.1.prop` and the number of taxa should be equal to the product of `taxa.1.prop` and the number of samples,

$$Z_{i.}\times J = Z_{.j}\times n$$

For example, if we want to generate data with library sizes uniformly distributed between 1,000 and 10,000 while changing the number of samples to 700 and allowing the proportion of non-zero cells to adjust according to the SCAM model, we can start by generating the marginal totals of non-zero cells and transforming the totals to proportions using the following:

```{r}
new.lib.size <- sample(1000:10000, size=700, replace=TRUE)

obs.sample.1.ct <- rowSums(count.ibd > 0)
xvar <- log10(rowSums(count.ibd))
scamfit.non0 <- scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
sample.1.ct <- 10^(predict(scamfit.non0,
                           newdata = data.frame(xvar = log10(new.lib.size) )) )

n.taxa <- ncol(count.ibd)
new.sample.1.prop <- sample.1.ct/n.taxa
new.taxa.1.prop <- fitted$taxa.1.prop * (sum(new.sample.1.prop) * n.taxa / 700) / sum(fitted$taxa.1.prop)
```

Then we can use *MIDASim.modify()* by

```{r}
count.ibd.modified <- MIDASim.modify(
  count.ibd.setup, 
  lib.size = new.lib.size,
  sample.1.prop = new.sample.1.prop,
  taxa.1.prop = new.taxa.1.prop
  )
```


#### Modify features of taxa by mean relative abundances or proportion of non-zero cells

If the mean relative abundances are modified, adjustments to the values of non-zero relative abundances may be necessary. Let $\hat{p}_j$ represent the mean relative abundances of taxa (mean.rel.abund), $\hat{p}^{(1)}_j$ denote the mean relative abundances of taxa among non-zero samples (mean.rel.abund.1), and $\hat{\delta}_j$ represent the proportion of non-zero cells for taxa (taxa.1.prop).

These three quantities are related by
$$\hat{p}_j = \hat{\delta}_j\hat{p}^{(1)}_j,$$

For example, if $\hat{p}_j$ is altered while keeping $\hat{\delta}j$ constant, MIDASim calculates the mean relative abundance of non-zero cells as $p^{(1)}_j=p_j/\delta_j$. Subsequently, it determines the value $\alpha_j$ for each taxon in a way that ensures that ${ \pi^\alpha{i,j} | \pi{ij}>0 }$ has a mean equal to $p^{(1)}_j$ for each taxon.

Suppose we want to increase $log(p_j)$ by $\beta=0.1$ for the first 10 taxa (e.g. modify these taxa according to a compositional model).  The following is the code for generating a set of target relative abundances as described, 

```{r}
beta <- 0.1
new.mean.rel.abund <- count.ibd.setup$mean.rel.abund
new.mean.rel.abund[1:10] <- exp(beta) * new.mean.rel.abund[1:10]
new.mean.rel.abund <- new.mean.rel.abund / sum(new.mean.rel.abund)
```

then the target taxa relative abundances can be used in *MIDASim.modify()* as the following,

```{r}
count.ibd.modified <- MIDASim.modify(count.ibd.setup,
                                     mean.rel.abund = new.mean.rel.abund)
```


Users also have the option to modify the proportion of non-zero cells for all taxa. However, because this modification affects the overall total count of non-zero cells throughout the entire dataset, users are required to specify the marginal proportions of non-zero cells for all samples. 

Suppose we wish to increase the proportion of non-zero cells for first 10 taxa by 0.05. We can achieve this by 

```{r}
new.taxa.1.prop <- count.ibd.setup$taxa.1.prop 
new.taxa.1.prop[1:10] <- new.taxa.1.prop[1:10] + 0.05
```

Then to ensure that appropriate proportions of non-zero cells are provided for all samples, here we adjust them accordingly by

```{r}
r <- sum(new.taxa.1.prop) / sum(count.ibd.setup$taxa.1.prop)
new.sample.1.prop <- count.ibd.setup$sample.1.prop * r
```

then the target proportion of non-zero cells for taxa can be used in *MIDASim.modify()* as the following,

```{r}
count.ibd.modified <- MIDASim.modify(
  count.ibd.setup,
  taxa.1.prop = new.taxa.1.prop,
  sample.1.prop = new.sample.1.prop
  )
```

Users can also simultaneously adjust both the mean relative abundances and the proportion of non-zero cells for taxa using the following code by

```{r}
count.ibd.modified <- MIDASim.modify(
  count.ibd.setup,
  mean.rel.abund = new.mean.rel.abund,
  taxa.1.prop = new.taxa.1.prop,
  sample.1.prop = new.sample.1.prop
  )
```


## Simulate microbiome data with the fitted MIDASim model

The function *MIDASim()* takes the quantities obtained from *MIDASim.modify()* and simulates microbiome data. The output is comprised of the simulated 0/1 (presence-absence) data, `sim_rel`, the table of relative abundances, and `sim_count`, the table of count data. 

```{r}
simulated.data <- MIDASim(count.ibd.modified)
summary(simulated.data)
```


## References

[1] Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N. *et al*. Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* 569, 655–662 (2019). https://doi.org/10.1038/s41586-019-1237-9.

[2] Charlson, E.S., Chen, J., Custers-Allen, R. *et al*. Disordered microbial communities in the upper respiratory tract of cigarette smokers. *PloS one*, 5(12), e15216 (2010). https://doi.org/10.1371/journal.pone.0015216

