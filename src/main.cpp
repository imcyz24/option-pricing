#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>

double norm_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

double bs_call(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}

double mc_call(double S, double K, double r, double sigma, double T, int n, unsigned seed = 42) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> N(0.0, 1.0);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        double ST = S * std::exp((r - 0.5 * sigma * sigma) * T + sigma * std::sqrt(T) * N(gen));
        sum += std::max(ST - K, 0.0);
    }
    return std::exp(-r * T) * sum / n;
}

int main() {
    const double S = 100, K = 100, r = 0.05, sigma = 0.2, T = 1.0;
    const double exact = bs_call(S, K, r, sigma, T);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Black-Scholes: " << exact << "\n\n";
    std::cout << "n        MC price   abs error\n";
    for (int n = 1000; n <= 10000000; n *= 10) {
        double mc = mc_call(S, K, r, sigma, T, n);
        std::cout << std::setw(8) << n << " " << mc << "   " << std::abs(mc - exact) << "\n";
    }
}