## Results

Parameters: S = 100, K = 100, r = 5%, σ = 20%, T = 1 year (European call), seed = 42

Black-Scholes closed-form price: **10.4506**

| Simulations (n) | MC price | Absolute error |
|----------------:|---------:|---------------:|
| 1,000           | 10.7783  | 0.3277         |
| 10,000          | 10.6298  | 0.1792         |
| 100,000         | 10.5120  | 0.0614         |
| 1,000,000       | 10.4723  | 0.0218         |
| 10,000,000      | 10.4501  | 0.0005         |

The Monte Carlo estimate converges toward the analytical price as n increases.
Error is expected to shrink at a rate of 1/√n. These numbers come from a single
seed, so individual errors are noisy; confidence intervals will be added next.
