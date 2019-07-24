#pragma once

#include "common.h"

class Lognum {
private:
	double log_value;

	static const std::array<double, 16384> log_factorial_table;
	static double log_factorial(int n) {
		if(n >= 0 && n < 16384) {
			return log_factorial_table[n];
		} else {
			return (double)lgamma((double)n + 1.0);
		}
	}

	Lognum(double log_value) : log_value(log_value) {}

public:
	Lognum() : Lognum(-INFINITY) {}

	static Lognum from_double(double val) {
		return Lognum(std::log(val));
	}
	static Lognum from_log(double val) {
		return Lognum(val);
	}

	static Lognum zero() {
		return Lognum();
	}
	static Lognum one() {
		return Lognum(0.0);
	}

	// Sample from [0, 1]
	static Lognum uniform_rand() {
		return Lognum::from_log(-std::exponential_distribution<double>(1.0)(rng));
	}
	
	// Sample from [0, upper_bound]
	static Lognum uniform_rand(Lognum upper_bound) {
		return uniform_rand() * upper_bound;
	}
	
	static Lognum binomial(int n, int k) {
		return Lognum(log_factorial(n) - log_factorial(k) - log_factorial(n - k));
	}

	// Raise to integer power
	Lognum powi(int exponent) const {
		return Lognum((double)exponent * log_value);
	}

	Lognum operator*(Lognum other) const {
		return Lognum(this->log_value + other.log_value);
	}

	Lognum operator+(Lognum l2) const {
		double log1 = this->log_value;
		double log2 = l2.log_value;

		if(log1 == -INFINITY && log2 == -INFINITY)  {
			return Lognum();
		} else {
			if(log1 >= log2) {
				return Lognum(log1 + std::log(1.0 + std::exp(log2 - log1)));
			} else {
				return Lognum(log2 + std::log(1.0 + std::exp(log1 - log2)));
			}
		}
	}

	bool operator<(Lognum o) const {
		return log_value < o.log_value;
	}
	bool operator>(Lognum o) const {
		return log_value > o.log_value;
	}
	bool operator<=(Lognum o) const {
		return log_value <= o.log_value;
	}
	bool operator>=(Lognum o) const {
		return log_value >= o.log_value;
	}
};
