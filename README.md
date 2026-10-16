### Delta

| Method                              | Delta  | Std error | 95% CI           |
|-------------------------------------|-------:|----------:|------------------|
| Black-Scholes (closed form)         | 0.6368 | -         | -                |
| MC finite difference (h = 1, 100k)  | 0.6370 | 0.0018    | [0.6335, 0.6406] |
| MC pathwise (100k)                  | 0.6370 | 0.0018    | [0.6334, 0.6405] |

Both estimators use the same seed, so their results are nearly identical.
For a call option the finite-difference and pathwise estimators have similar
variance. Pathwise has no bias from the bump size h and needs only one
payoff evaluation per path.
### Coverage test

To check that the 95% confidence interval is calibrated, I ran 100 independent
simulations (seeds 1 to 100, n = 10,000 each) and counted how often the interval
contained the Black-Scholes price.

Result: **96 / 100** intervals contained the true price.

For a correctly calibrated 95% interval, the expected count is about 95, with a
standard deviation of roughly 2.2, so this result is consistent with theory
Plain MC      (100000 paths): price 10.5120, std err 0.0469
Antithetic MC (50000 pairs): price 10.4520, std err 0.0328
Std err ratio (antithetic / plain): 0.7002
Average std err ratio over 20 seeds: 0.7070
BS Delta: 0.6368
MC Delta (finite diff, h = 1): 0.6370, std err 0.0018, 95% CI [0.6335, 0.6406], BS Delta in CI? yes
MC Delta (pathwise): 0.6370, std err 0.0018, 95% CI [0.6334, 0.6405], BS Delta in CI? yes