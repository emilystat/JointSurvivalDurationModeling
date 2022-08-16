# JointSpatialModeling
This is a study of joint spatial modeling of time to infection (Weibull) and disease duration (Poisson) in risk assessment of botanical epidemics.<br/>
The random effect is simulated from four models: IndNorm, MvNorm, UniCAR and GMCAR. Each of the simulated data is then fitted to five models: IndNorm, MvNorm, UniCAR, GMCAR and GMCAR (reversed order).<br/>
AMSE, DIC and LOOIC are calculated and used to compare different models.

# Required R packages
`mvtnorm`, `rstan`, `loo`

# References
1. Banerjee, Sudipto, Melanie M. Wall, and Bradley P. Carlin. "[Frailty modeling for spatially correlated survival data, with application to infant mortality in Minnesota.](https://academic.oup.com/biostatistics/article/4/1/123/246234?login=true)" _Biostatistics_ 4.1 (2003): 123-142.
2. Jin, Xiaoping, Bradley P. Carlin, and Sudipto Banerjee. "[Generalized hierarchical multivariate CAR models for areal data.](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1541-0420.2005.00359.x)" _Biometrics_ 61.4 (2005): 950-961.

# Contributors
J. Zhang, E. L. Kang, P. S. Ojiambo
