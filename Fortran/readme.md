# Practical in fortran for the processes in Climate change and terrestrial ecosystem modeling

## 2. Quantitative description of ecosystems

Terestial biosphere models characterize ecosystems by features that control biogeochemical cycles and energy, massn and momentum fluxes with the atmosphere.
These include:

- Leaf area index and its vertical profile in the canopy
- leaf angle distribution
- the vertical profile of leaf mass and leaf nitrogen in the canopy
- the profile of roots in the soil
- the size structure of plants
- the distribution of carbon within an ecosystem


This chapter defines these decriptors of ecosystems

1. Leaf area density

> Use beta distribution probability density function to calculate the leaf area density profile.
>
> - It calculate the profile for different values of $p$ and $q$ in the beta function (一种由正态分布演化来的函数) with $ \frac{L}{h_c}$

2. Leaf angle distribution

> Uses the beta distribution to calculate leaf angle distribution.
>
> - It uses the mean and standard deviations for the five leaf angel distribution in table 2.1, and 
> - compares the numverical solution obtained with nine leaf angle classes ($10^\circ$) to the analytical solutions.

## 5. Soil temperature

1. The code in supplemental program 5.2 calculates the diurnal cycle of soil temperature. Compare results after one simulation day, 10 simulations days, and 100 days. How does soil temperature differ between days. How does soil temperature differ between days? when does soil temperature reach equilibrium?

2. Figure 5.8a shows the diurnal cycle of soil temperature for sand. Contrast this when the soil is covered by 5 cm of organic meterial with $k = 0.5 W/mK$. What is the effect of an organic layer on the diurnal cycle?

3. Modify th code in supplemental program 5.2 to calculate the annual cycle of soil temperature with a minimum surface temperature of $-10^\circ C$ on day 15 and a maximum of $30^\circ C$ on dayt 197.5. In this case, $T_0(t) = \bar{T_0} + A_0sin[2\pi(t-106.25)/365] $ for t in days of the year (including fractions of the current day). Set the depth of snow to 10 cm during November through February; otherwise, there is no snow. Use $K_{snow} = 0.3W/mK$. Compare soil temperatures with and without the snow. Describe the effect of snow on soil temperature.

4. Investigate the effects of phase change on the annual cycle of soil temperature. compare simulations with and without phase change.

## 6. Turbulent fluxes and scalar profiles in the surface layer

## 7. Surface energy fluxes

## 8. Soil moisture

## 9. Hydrologic scaling and spatial heterogeneity

## 10. Leaf temperature and energy fluxes

## 11. Leaf photosynthesis

## 12. Stomatal conductance

## 13. Plant hydraulics

## 14. Radiative transfer

## 16. Plant canopies

## 17. Scalar canopy profiles

## 18. Soil biogeochemistry