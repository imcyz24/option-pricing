| n          | MC price | Std error | 95% CI low | 95% CI high | BS in CI? |
|-----------:|---------:|----------:|-----------:|------------:|:---------:|
| 1,000      | 10.7783  | 0.4913    | 9.8152     | 11.7413     | yes       |
| 10,000     | 10.6298  | 0.1482    | 10.3394    | 10.9202     | yes       |
| 100,000    | 10.5120  | 0.0469    | 10.4202    | 10.6039     | yes       |
| 1,000,000  | 10.4723  | 0.0147    | 10.4434    | 10.5013     | yes       |
| 10,000,000 | 10.4501  | 0.0047    | 10.4410    | 10.4592     | yes       |

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
