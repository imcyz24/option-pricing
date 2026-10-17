#include "pricing.hpp"

#include <algorithm>
#include <cmath>
#include <random>

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

McResult mc_call(double S, double K, double r, double sigma, double T, int n, unsigned seed) {
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

McResult mc_call_av(double S, double K, double r, double sigma, double T, int n_pairs, unsigned seed) {
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
                     int n, double h, unsigned seed) {
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
                     int n, unsigned seed) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> N(0.0, 1.0);
    const double disc = std::exp(-r * T);
    const double drift = (r - 0.5 * sigma * sigma) * T;
    const double vol = sigma * std::sqrt(T);
    double sum = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < n; ++i) {
        double growth = std::exp(drift + vol * N(gen));  // S_T / S
        double d = (S * growth > K) ? disc * growth : 0.0;
        sum += d;
        sum_sq += d * d;
    }
    double mean = sum / n;
    double var = (sum_sq - n * mean * mean) / (n - 1);
    return {mean, std::sqrt(var / n)};
}
