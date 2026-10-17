#pragma once

struct McResult {
    double price;
    double std_err;
};

double norm_cdf(double x);
double bs_call(double S, double K, double r, double sigma, double T);
double bs_delta(double S, double K, double r, double sigma, double T);

McResult mc_call(double S, double K, double r, double sigma, double T,
                 int n, unsigned seed = 42);
McResult mc_call_av(double S, double K, double r, double sigma, double T,
                    int n_pairs, unsigned seed = 42);
McResult mc_delta_fd(double S, double K, double r, double sigma, double T,
                     int n, double h, unsigned seed = 42);
McResult mc_delta_pw(double S, double K, double r, double sigma, double T,
                     int n, unsigned seed = 42);
                     