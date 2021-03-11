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
Tgrid = c(0.1, 0.25,0.5, 1, 2, 3)
Kgrid = c(90, 95, 100, 105, 110, 115, 120)
LV = matrix(0.25, nrow=length(Tgrid), ncol=length(Kgrid))
Mext = seq(0.1, 4, by=0.2)
IV = localVolCalls(S0=100, rf=0.015, q=0.02, LV=LV, Kgrid=Kgrid, Tgrid=Tgrid, KflatExt=100*Mext, impVol=TRUE)
```
will produce a 6x7 matrix `IV` with implied volatilities over Tgrid x Kgrid:

```bash
> IV
          [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
[1,] 0.2098161 0.2214284 0.2642055 0.2252830 0.2334082 0.2474540 0.2583035
[2,] 0.2123303 0.2414574 0.2370104 0.2504015 0.2439533 0.2406356 0.2383015
[3,] 0.2265204 0.2329941 0.2438070 0.2409527 0.2449153 0.2443306 0.2413687
[4,] 0.2386151 0.2419249 0.2400347 0.2442403 0.2427800 0.2427263 0.2423808
[5,] 0.2464329 0.2453425 0.2472783 0.2447980 0.2456464 0.2450801 0.2442844
[6,] 0.2474535 0.2480237 0.2462238 0.2476969 0.2465076 0.2463319 0.2461990
```

![fig1](fig2.png)

## Demo

`main_gp_lv.R` 

