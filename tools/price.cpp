#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "pricing.hpp"

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: price S K r sigma T\n"
                  << "Example: price 100 100 0.05 0.2 1\n";
        return 1;
    }
    const double S = std::atof(argv[1]);
    const double K = std::atof(argv[2]);
    const double r = std::atof(argv[3]);
    const double sigma = std::atof(argv[4]);
    const double T = std::atof(argv[5]);
    if (S <= 0 || K <= 0 || sigma <= 0 || T <= 0) {
        std::cerr << "S, K, sigma and T must be positive numbers\n";
        return 1;
    }

    McResult mc = mc_call_av(S, K, r, sigma, T, 200000);

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Call (Black-Scholes): " << bs_call(S, K, r, sigma, T) << "\n";
    std::cout << "Put  (Black-Scholes): " << bs_put(S, K, r, sigma, T) << "\n";
    std::cout << "Call Delta:           " << bs_delta(S, K, r, sigma, T) << "\n";
    std::cout << "Call (Monte Carlo):   " << mc.price << "  95% CI ["
              << mc.price - 1.96 * mc.std_err << ", "
              << mc.price + 1.96 * mc.std_err << "]\n";
}
