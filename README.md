# GP-LV
Gaussian process local volatility model.

Code accompanies this [paper](https://arxiv.org/abs/1901.06021) by Martin Tegner and Stephen Roberts.

## Install

Clone and enter the repo.

```bash
git clone https://github.com/martnj/GP-LV
cd GP-LV
```

## Datasets

`data_gp_lv.Rdata` contains a synthetic dataset with price and implied volatility surfaces from a single date. 

![fig1](fig1.png)

Similarly, `data_sequence_gp_lv.Rdata` contains a series of surfaces from ten dates.

## Functions

Functions in `functions_gp_lv.R` include MCMC samplers, GP covariance functions and prediction equations, Black-Scholes and local volatility pricing functions.

```bash
Tgrid = c(0.5, 2, 3)
Kgrid = c(90, 95, 100, 105, 110, 115, 120)
LV = matrix(0.5, nrow=length(Tgrid), ncol=length(Kgrid))
IV = localVolCalls(S0=100, rf=0.015, q=0.02, LV=LV, Kgrid=Kgrid, Tgrid=Tgrid, impVol=TRUE)
```
will produce a 6x7 matrix `IV` with implied volatilities over Tgrid x Kgrid:

![fig1](fig2.png)

## Demo

`main_gp_lv.R` 

