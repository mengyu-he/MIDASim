<!-- badges: start -->
  [![R-CMD-check](https://github.com/mengyu-he/MIDASim/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mengyu-he/MIDASim/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

# MIDASim: a fast and simple simulator for realistic microbiome data


# Overview

The MIDASim package (v0.1.0) implements the **MI**crobiome **DA**ta **Sim**ulator functions that simulate realistic microbiome data, maintaining the distributional and taxon-taxon correlation of a template microbiome datase. Users may also change the parameter setup in the model to introduce ``effect'' to facilitate simulation.  

## Required packages 

Required packages for functions in MIDASim include: psych, MASS, vegan, pracma, scam. These packages are all available on CRAN.

## Installation

Install from the local file (the .tar.gz file), which you can download from https://github.com/mengyu-he/MIDASim.

```{r}
devtools::install_local("MIDASim_0.1.0.tar.gz", dependencies = TRUE)
```

A github installation is also available via

```{r}
install_github("mengyu-he/MIDASim")
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

There are three main functions in MIDASim package: *MIDASim.setup()*, *MIDASim.modify()*, *MIDASim()*.

## Template data description

We demonstrate the method of MIDASim and functions in this package using a filtered microbiome dataset of patients with IBD(Inflammatory Bowel Disease) in Human Microbiome Project 2 (HMP2) [1]. The data can be directly loaded through 

```{r}
data("count.ibd")
```

This HMP2 data were modified through package **HMP2Data** on Bioconductor by filtering out samples with library sizes smaller than 3000 and taxa that appear in less than 2 samples. The filtered data have 146 rows (samples) and 614 columns (taxa).

```{r, eval = FALSE}
library(HMP2Data)
IBD16S()

count.ibd = t(IBD16S_mtx[, colSums(IBD16S_mtx)>3000] )
count.ibd = count.ibd[, colSums(count.ibd>0)>1]

dim(count.ibd)
```


Other datasets come with the package include a vaginal microbiome dataset from HMP2 [1] and a upper-respiratory-tract microbiome dataset [2].

## Setup the MIDASim model

The function *MIDASim.setup()* fits the underlying MIDASim model and extracts information from the template data and returns the estimated parameters for microbiome simulation. MIDASim offers two modeling options for relative abundances: nonparametric and parametric. Both approaches utilize a two-stage procedure in which the first stage models the presence-absence relationship between taxa and the second stage models the non-zero values in the microbiome community. These two approaches differ by the way how the second stage marginal distributions are modeled. The nonparametric approach generates data that resemble the template data better than the parametric approach. However, the parametric approach offers modifications of the template data in more controlled way.

In nonparametric mode (the default), a rank based approach is used to obtain the relative abundances for taxa that are considered to be non-zero in the simulated data. 

To fit MIDASim with nonparametric mode to the *count.ibd* data, use the following command:

```{r}
count.ibd.setup = MIDASim.setup(otu.tab, mode = 'nonparametric', n.break.ties = 10)
```

The argument `n.break.ties` specifies the number of replicates used to break ties when ranking relative abundances for each taxon, and the default is `n.break.ties = 100`. 

In parametric mode, MIDASim fits the generalized gamma model, a three-parameter distribution in the location-scale family that was proposed for analyzing right-censored survival data to the relative abundance data of each taxon separately.

Denote the simulated relative abundance for $i$-th subject and $j$-th taxon as $\tilde{\pi}_{ij}$. We define ``survival time'' 
\begin{equation} \label{eq:survivalTime}
    \tilde{t}_{ij} =
    \begin{cases}
\frac{1}{\tilde{\pi}_{ij}} & \text{if } \tilde{\pi}_{ij} > 0
\\
N_i & \text{if } \tilde{\pi}_{ij} = 0 \\
\end{cases}
\end{equation}
which corresponds to treating $\tilde{t}_{ij}$ as right-censored when $\pi_{ij}<\frac{1}{N_i}$.  The generalized gamma model then assumes $\widetilde{t}_{ij}$ has the distribution specified by
\begin{equation}\label{eq:logLinear}
 \ln(\tilde{t}_{ij}) =- \mu_j + s_j\sigma_j \cdot \omega_{ij}~, 
\end{equation}
where $e^{\omega_{ij}}$ follows a gamma distribution with shape parameter $k_j = 1/|Q_j|$ and scale parameter 1, and $s_j=\text{sign}(Q_j)$}. 

We estimate the three parameters $\mu_j, \sigma_j, Q_j$ in *MIDASim.setup()* using the relative abundances $\pi_j$ in the template.

To fit MIDASim in parametric mode to the count.ibd data, use the following command:

```{r}
count.ibd.setup = MIDASim.setup(otu.tab, mode = 'parametric')
```

The returned values of *MIDASim.setup()* include the estimated tetrachoric correlation of observed presence-absence data (`tetra.corr`), the estimated Pearson correlation of observed relative abundances (`corr.rel.corrected`), the proportion of non-zero cells for each taxon (`taxa.1.prop`), the proportion of non-zero cells for each subject (`sample.1.prop`), the observed mean relative abundances of each taxon (`mean.rel.abund`), and the observed relative abundances among non-zero cells for each taxon (`rel.abund.1`). If mode = 'parametric', the function will also return the estimated parameters of the generalized gamma distribution (`mu.est`, `sigma.est`, and `Q.est`).

## Modify the fitted MIDASim model

The function *MIDASim.modify()* sets up quantities required to simulate data using MIDAS.  This function also allows users to modify library sizes, taxa relative abundances, the location parameters in the parametric model, taxa proportion of non-zero cells, number of samples. Note that this is a **required** step even if no adjustment is made.

### No modification at all

The default is no additional adjustment is wanted. The number of samples and their target library sizes for the simulated data are the same as that of the template. The following code is an example.

```{r}
count.ibd.modified = MIDASim.modify(count.ibd.setup, 
                                   lib.size = NULL,
                                   mean.rel.abund = NULL,
                                   gengamma.mu = NULL,
                                   sample.1.prop = NULL,
                                   taxa.1.prop = NULL)
```

Next, we will explain how to modify the data features to accommodate a wide range of simulation scenarios, separately for the nonparametric and parametric modes.

### Model modification in nonparametric mode

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
new.lib.size = sample(1000:10000, size=700, replace=TRUE)

obs.sample.1.ct = rowSums(count.ibd > 0)
xvar = log10(rowSums(count.ibd))
scamfit.non0 = scam::scam( log10(obs.sample.1.ct) ~  s( xvar, bs = "mpi" ))
sample.1.ct = 10^(predict(scamfit.non0,
                           newdata = data.frame(xvar = log10(new.lib.size) )) )

n.taxa = ncol(count.ibd)
new.sample.1.prop = sample.1.ct/n.taxa
new.taxa.1.prop = fitted$taxa.1.prop * (sum(new.sample.1.prop) * n.taxa / 700) / sum(fitted$taxa.1.prop)
```

Then we can use *MIDASim.modify()* by

```{r}
count.ibd.modified = MIDASim.modify(count.ibd.setup, 
                                    lib.size = new.lib.size,
                                    sample.1.prop = new.sample.1.prop,
                                    taxa.1.prop = new.taxa.1.prop)
```


#### Modify features of taxa by mean relative abundances or proportion of non-zero cells

If the mean relative abundances are modified, adjustments to the values of non-zero relative abundances may be necessary. Let $\hat{p}_j$ represent the mean relative abundances of taxa (mean.rel.abund), $\hat{p}^{(1)}_j$ denote the mean relative abundances of taxa among non-zero samples (mean.rel.abund.1), and $\hat{\delta}_j$ represent the proportion of non-zero cells for taxa (taxa.1.prop).

These three quantities are related by
$$\hat{p}_j = \hat{\delta}_j\hat{p}^{(1)}_j,$$

For example, if $\hat{p}_j$ is altered while keeping $\hat{\delta}j$ constant, MIDASim calculates the mean relative abundance of non-zero cells as $p^{(1)}_j=p_j/\delta_j$. Subsequently, it determines the value $\alpha_j$ for each taxon in a way that ensures that ${ \pi^\alpha{i,j} | \pi{ij}>0 }$ has a mean equal to $p^{(1)}_j$ for each taxon.

Suppose we want to increase $log(p_j)$ by $\beta=0.1$ for the first 10 taxa (e.g. modify these taxa according to a compositional model).  The following is the code for generating a set of target relative abundances as described, 

```{r}
beta = 0.1
new.mean.rel.abund = count.ibd.setup$mean.rel.abund
new.mean.rel.abund[1:10] = exp(beta) * new.mean.rel.abund[1:10]
new.mean.rel.abund = new.mean.rel.abund / sum(new.mean.rel.abund)
```

then the target taxa relative abundances can be used in *MIDASim.modify()* as the following,

```{r}
count.ibd.modified = MIDASim.modify(count.ibd.setup,
                                    mean.rel.abund = new.mean.rel.abund)
```


Users also have the option to modify the proportion of non-zero cells for all taxa. However, because this modification affects the overall total count of non-zero cells throughout the entire dataset, users are required to specify the marginal proportions of non-zero cells for all samples. 

Suppose we wish to increase the proportion of non-zero cells for first 10 taxa by 0.05. We can achieve this by 

```{r}
new.taxa.1.prop = count.ibd.setup$taxa.1.prop 
new.taxa.1.prop[1:10] = new.taxa.1.prop[1:10] + 0.05
```

Then to ensure that appropriate proportions of non-zero cells are provided for all samples, here we adjust them accordingly by

```{r}
r = sum(new.taxa.1.prop) / sum(count.ibd.setup$taxa.1.prop)
new.sample.1.prop = count.ibd.setup$sample.1.prop * r
```

then the target proportion of non-zero cells for taxa can be used in *MIDASim.modify()* as the following,

```{r}
count.ibd.modified = MIDASim.modify(count.ibd.setup,
                                    taxa.1.prop = new.taxa.1.prop,
                                    sample.1.prop = new.sample.1.prop)
```

Users can also simultaneously adjust both the mean relative abundances and the proportion of non-zero cells for taxa using the following code by

```{r}
count.ibd.modified = MIDASim.modify(count.ibd.setup,
                                    mean.rel.abund = new.mean.rel.abund,
                                    taxa.1.prop = new.taxa.1.prop,
                                    sample.1.prop = new.sample.1.prop)
```

### parametric mode

If the parametric mode is chosen in *MIDASim.setup()*, then the user are able to modify the sample size, library sizes, taxa relative abundances, the location parameters in the generalized gamma distribution. 

Changes in the parameters of the parametric model (including library sizes) imply coordinated changes in all other quantities. After either modification of the parameters, we predict the probability of being non-zero of $i$-th subject and $j$-th taxon by 

\begin{equation} \label{gengamma_pred01}
P(\widetilde{Z}_{ij} = 1) = F_j(N_i ~;~ \widehat{\mu}_j, \widehat{\sigma},_j \widehat{Q}_j),
\end{equation}
to obtain the number of non-zero cells in each subject $Z_{i\cdot }$ and that in each taxon $Z_{\cdot j}$.

#### Modify the number of samples and their library sizes

To change the library sizes, set the `lib.size` argument equal to a vector with the target library sizes.  The length of this vector becomes the new number of samples.  The advantage of the parametric mode is that it handles coordinated changes automatically. As a result, when changing library sizes in parametric mode, there's no need for the user to provide the marginal totals of non-zero cells explicitly.

For example, to generate data having library sizes that are uniformly distributed between 1,000 and 10,000 while changing the number of samples to 700, we execute the code

```{r}
new.lib.size = sample(1000:10000, size=700, replace=TRUE)
count.ibd.modified = MIDASim.modify(count.ibd.setup, 
                                    lib.size = new.lib.size)
```


### Modify taxa features by mean relative abundances or location parameters

Specifying new values for mean relative abundances (`mean.rel.abund`) or location parameters (`gengamma.mu`) can modify taxa features. 

When mean relative abundances are changed, it effectively adjusts the location parameters in the parametric model. In this case, when the user specifies mean relative abundances for taxa, MIDASim automatically estimates new location parameters for the generalized gamma distribution. These changes in location parameters influence the predicted probabilities of non-zero values for each subject and each taxon, subsequently leading to adjustments in the marginal totals of non-zero cells.

For example, if we want to increase $log(p_j)$ by $\beta=0.1$ for the first 10 taxa (e.g. modify these taxa according to a compositional model).  The following is the code for generating a set of target relative abundances as described, 

```{r}
beta = 0.1
new.mean.rel.abund = count.ibd.setup$mean.rel.abund
new.mean.rel.abund[1:10] = exp(beta) * new.mean.rel.abund[1:10]
new.mean.rel.abund = new.mean.rel.abund / sum(new.mean.rel.abund)
```

then the target taxa relative abundances can be used in *MIDASim.modify()* as the following,

```{r}
count.ibd.modified = MIDASim.modify(count.ibd.setup, 
                                    mean.rel.abund = new.mean.rel.abund)
```

Alternatively, the users can also directly change the location parameters $\mu_j$ by

$$\mu_j \rightarrow \mu_j + \beta $$

```{r}
count.ibd.modified = MIDASim.modify(count.ibd.setup, 
                                    gengamma.mu = count.ibd.setup$mu.est + beta)
```


## Simulate microbiome data with the fitted MIDAS model

The function *MIDASim()* takes the quantities obtained from *MIDASim.modify()* and simulates microbiome data. The output is comprised of the simulated 0/1 (presence-absence) data, `sim_rel`, the table of relative abundances, and `sim_count`, the table of count data. 

```{r}
simulated.data = MIDASim(count.ibd.modified)
summary(simulated.data)
```


## References

[1] Lloyd-Price, J., Arze, C., Ananthakrishnan, A.N. *et al*. Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases. *Nature* 569, 655–662 (2019). https://doi.org/10.1038/s41586-019-1237-9.

[2] Charlson, E.S., Chen, J., Custers-Allen, R. *et al*. Disordered microbial communities in the upper respiratory tract of cigarette smokers. *PloS one*, 5(12), e15216 (2010). https://doi.org/10.1371/journal.pone.0015216

