#include <gtest/gtest.h>

#include <cmath>

#include "pricing.hpp"

namespace {
const double S = 100, K = 100, r = 0.05, sigma = 0.2, T = 1.0;
}

// Independent check: known textbook values, not derived from our own code
TEST(BlackScholes, KnownValues) {
    EXPECT_NEAR(bs_call(S, K, r, sigma, T), 10.4506, 1e-4);
    EXPECT_NEAR(bs_put(S, K, r, sigma, T), 5.5735, 1e-4);
    EXPECT_NEAR(bs_delta(S, K, r, sigma, T), 0.6368, 1e-4);
}

// C - P = S - K e^{-rT}, for several parameter sets
TEST(BlackScholes, PutCallParity) {
    const double strikes[] = {80, 100, 120};
    const double vols[] = {0.1, 0.2, 0.5};
    for (double k : strikes) {
        for (double v : vols) {
            double lhs = bs_call(S, k, r, v, T) - bs_put(S, k, r, v, T);
            double rhs = S - k * std::exp(-r * T);
            EXPECT_NEAR(lhs, rhs, 1e-9) << "K=" << k << " sigma=" << v;
        }
    }
}

TEST(MonteCarlo, PriceWithinFourStdErrOfBS) {
    McResult mc = mc_call(S, K, r, sigma, T, 100000, 42);
    EXPECT_NEAR(mc.price, bs_call(S, K, r, sigma, T), 4.0 * mc.std_err);
}

TEST(MonteCarlo, AntitheticWithinFourStdErrOfBS) {
    McResult mc = mc_call_av(S, K, r, sigma, T, 50000, 42);
    EXPECT_NEAR(mc.price, bs_call(S, K, r, sigma, T), 4.0 * mc.std_err);
}

TEST(MonteCarlo, AntitheticReducesStdErr) {
    McResult plain = mc_call(S, K, r, sigma, T, 100000, 42);
    McResult av = mc_call_av(S, K, r, sigma, T, 50000, 42);
    EXPECT_LT(av.std_err, plain.std_err);
}

TEST(MonteCarlo, DeltaEstimatorsWithinFourStdErrOfBS) {
    const double exact = bs_delta(S, K, r, sigma, T);
    McResult fd = mc_delta_fd(S, K, r, sigma, T, 100000, 1.0, 42);
    McResult pw = mc_delta_pw(S, K, r, sigma, T, 100000, 42);
    EXPECT_NEAR(fd.price, exact, 4.0 * fd.std_err);
    EXPECT_NEAR(pw.price, exact, 4.0 * pw.std_err);
}