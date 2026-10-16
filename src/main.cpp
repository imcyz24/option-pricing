#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

struct McResult {
    double price;
    double std_err;
};

double norm_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

double bs_call(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}
double bs_delta(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return norm_cdf(d1);
}

McResult mc_call(double S, double K, double r, double sigma, double T, int n, unsigned seed = 42) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> N(0.0, 1.0);
    const double disc = std::exp(-r * T);
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        double ST = S * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * N(gen));
        double payoff = disc * std::max(ST - K, 0.0);
        sum += payoff;
        sum_sq += payoff * payoff;
    }
    double mean = sum / n;
    double var = (sum_sq - n * mean * mean) / (n - 1);
    return {mean, std::sqrt(var / n)};
}
McResult mc_call_av(double S, double K, double r, double sigma, double T, int n_pairs, unsigned seed = 42) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> N(0.0, 1.0);
    const double disc = std::exp(-r * T);
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < n_pairs; ++i) {
        double z = N(gen);
        double up = S * std::exp(drift + vol * z);
        double down = S * std::exp(drift - vol * z);
        double y = disc * 0.5 * (std::max(up - K, 0.0) + std::max(down - K, 0.0));
        sum += y;
        sum_sq += y * y;
    }
    double mean = sum / n_pairs;
    double var = (sum_sq - n_pairs * mean * mean) / (n_pairs - 1);
    return {mean, std::sqrt(var / n_pairs)};
}
McResult mc_delta_fd(double S, double K, double r, double sigma, double T,
                     int n, double h, unsigned seed = 42) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> N(0.0, 1.0);
    const double disc = std::exp(-r * T);
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        double growth = std::exp(drift + vol * N(gen));  // same Z for both bumps
        double up   = std::max((S + h) * growth - K, 0.0);
        double down = std::max((S - h) * growth - K, 0.0);
        double d = disc * (up - down) / (2.0 * h);
        sum += d;
        sum_sq += d * d;
    }
    double mean = sum / n;
    double var = (sum_sq - n * mean * mean) / (n - 1);
    return {mean, std::sqrt(var / n)};
}
McResult mc_delta_pw(double S, double K, double r, double sigma, double T,
                     int n, unsigned seed = 42) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> N(0.0, 1.0);
    const double disc = std::exp(-r * T);
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        double growth = std::exp(drift + vol * N(gen));   // S_T / S
        double d = (S * growth > K) ? disc * growth : 0.0;
        sum += d;
        sum_sq += d * d;
    }
    double mean = sum / n;
    double var = (sum_sq - n * mean * mean) / (n - 1);
    return {mean, std::sqrt(var / n)};
}

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

    const int n_total = 100000;
    McResult plain = mc_call(S, K, r, sigma, T, n_total);
    McResult av = mc_call_av(S, K, r, sigma, T, n_total / 2);
    std::cout << "\nPlain MC      (" << n_total << " paths): price " << plain.price
              << ", std err " << plain.std_err << "\n";
    std::cout << "Antithetic MC (" << n_total / 2 << " pairs): price " << av.price
              << ", std err " << av.std_err << "\n";
    std::cout << "Std err ratio (antithetic / plain): " << av.std_err / plain.std_err << "\n";

    double ratio_sum = 0.0;
    const unsigned ratio_trials = 20;
    for (unsigned seed = 1; seed <= ratio_trials; ++seed) {
        McResult p = mc_call(S, K, r, sigma, T, n_total, seed);
        McResult a = mc_call_av(S, K, r, sigma, T, n_total / 2, seed);
        ratio_sum += a.std_err / p.std_err;
    }
    std::cout << "Average std err ratio over " << ratio_trials << " seeds: "
              << ratio_sum / ratio_trials << "\n";
              // Is the variance reduction stable across seeds?
 std::cout << "\nBS Delta: " << bs_delta(S, K, r, sigma, T) << "\n";
 McResult d = mc_delta_fd(S, K, r, sigma, T, 100000, 1.0);
    double bs_d = bs_delta(S, K, r, sigma, T);
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