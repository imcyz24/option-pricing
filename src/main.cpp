#include <iomanip>
#include <iostream>

#include "pricing.hpp"

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.2, T = 1.0;
    const double exact = bs_call(S, K, r, sigma, T);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Black-Scholes: " << exact << "\n\n";
    std::cout << "n          MC price   std err   95% CI low  95% CI high  BS in CI?\n";
    for (int n = 1000; n <= 10000000; n *= 10) {
        McResult mc = mc_call(S, K, r, sigma, T, n);
        double lo = mc.price - 1.96 * mc.std_err;
        double hi = mc.price + 1.96 * mc.std_err;
        bool inside = (exact >= lo && exact <= hi);
        std::cout << std::setw(10) << n << " " << mc.price << "   " << mc.std_err
                  << "   " << lo << "      " << hi << "      "
                  << (inside ? "yes" : "no") << "\n";
    }

    // Coverage test: how often does the 95% CI contain the true price?
    const int n_test = 10000;
    const unsigned trials = 100;
    int hits = 0;
    for (unsigned seed = 1; seed <= trials; ++seed) {
        McResult mc = mc_call(S, K, r, sigma, T, n_test, seed);
        double lo = mc.price - 1.96 * mc.std_err;
        double hi = mc.price + 1.96 * mc.std_err;
        if (exact >= lo && exact <= hi) ++hits;
    }
    std::cout << "\nCoverage test (n = " << n_test << ", " << trials << " seeds): "
              << hits << " / " << trials << " intervals contain the BS price\n";

    // Antithetic vs plain Monte Carlo, same number of payoff evaluations
    const int n_total = 100000;
    McResult plain = mc_call(S, K, r, sigma, T, n_total);
    McResult av = mc_call_av(S, K, r, sigma, T, n_total / 2);
    std::cout << "\nPlain MC      (" << n_total << " paths): price " << plain.price
              << ", std err " << plain.std_err << "\n";
    std::cout << "Antithetic MC (" << n_total / 2 << " pairs): price " << av.price
              << ", std err " << av.std_err << "\n";
    std::cout << "Std err ratio (antithetic / plain): " << av.std_err / plain.std_err << "\n";

    // Is the variance reduction stable across seeds?
    double ratio_sum = 0.0;
    const unsigned ratio_trials = 20;
    for (unsigned seed = 1; seed <= ratio_trials; ++seed) {
        McResult p = mc_call(S, K, r, sigma, T, n_total, seed);
        McResult a = mc_call_av(S, K, r, sigma, T, n_total / 2, seed);
        ratio_sum += a.std_err / p.std_err;
    }
    std::cout << "Average std err ratio over " << ratio_trials << " seeds: "
              << ratio_sum / ratio_trials << "\n";

    // Delta
    const double bs_d = bs_delta(S, K, r, sigma, T);
    std::cout << "\nBS Delta: " << bs_d << "\n";

    McResult d = mc_delta_fd(S, K, r, sigma, T, 100000, 1.0);
    double d_lo = d.price - 1.96 * d.std_err;
    double d_hi = d.price + 1.96 * d.std_err;
    std::cout << "MC Delta (finite diff, h = 1): " << d.price
              << ", std err " << d.std_err
              << ", 95% CI [" << d_lo << ", " << d_hi << "]"
              << ", BS Delta in CI? " << ((bs_d >= d_lo && bs_d <= d_hi) ? "yes" : "no") << "\n";

    McResult pw = mc_delta_pw(S, K, r, sigma, T, 100000);
    double pw_lo = pw.price - 1.96 * pw.std_err;
    double pw_hi = pw.price + 1.96 * pw.std_err;
    std::cout << "MC Delta (pathwise): " << pw.price
              << ", std err " << pw.std_err
              << ", 95% CI [" << pw_lo << ", " << pw_hi << "]"
              << ", BS Delta in CI? " << ((bs_d >= pw_lo && bs_d <= pw_hi) ? "yes" : "no") << "\n";
}