#include <iostream>
#include <math.h>
#include "fadiff.h"
#include "badiff.h"
#include <vector>
#include <vector>
#include <random>
#include <chrono>

using namespace std;
using  ns = chrono::nanoseconds;
using get_time = chrono::steady_clock;

//random numbers generator
random_device generator;

//normal distribution
normal_distribution<double> distribution(0, 1);

F<double> ExampleFunction(F<double> x, F<double> y){
	return 2 * x + exp(y);
}

F<double> Max(F<double> a, F<double> b){
	return a < b ? b : a;
}

B<double> Max(B<double> a, B<double> b){
	return a < b ? b : a;
}

double Max(double a, double b){
	return a < b ? b : a;
}

//Monte Carlo Example
F<double> MonteCarlo(F<double> S_0, F<double> mu, F<double> sigma, F<double> dt, F<double> T, F<double> K){
	F<double> RandNum;
	F<double> n = T / dt;
	F<double> dS;
	F<double> S_t = S_0;
	for (int i = 0; i < n; ++i){
		RandNum = distribution(generator);
		dS = S_t*mu*dt + S_t*sigma*sqrt(dt)*RandNum;
		S_t += dS;
	}
	F<double> PayOff = S_t - K;
	return Max(PayOff, 0);
}

B<double> MonteCarlo(B<double> S_0, B<double> mu, B<double> sigma, B<double> dt, B<double> T, B<double> K){
	B<double> RandNum;
	B<double> n = T / dt;
	B<double> dS;
	B<double> S_t = S_0;
	for (int i = 0; i < n; ++i){
		RandNum = distribution(generator);
		dS = S_t*mu*dt + S_t*sigma*sqrt(dt)*RandNum;
		S_t += dS;
	}
	B<double> PayOff = S_t - K;
	return Max(PayOff, 0);
}

double MonteCarlo(double S_0, double mu, double sigma, double dt, double T, double K){
	double RandNum;
	double n = T / dt;
	double dS;
	double S_t = S_0;
	for (int i = 0; i < n; ++i){
		RandNum = distribution(generator);
		dS = S_t*mu*dt + S_t*sigma*sqrt(dt)*RandNum;
		S_t += dS;
	}
	double PayOff = S_t - K;
	return Max(PayOff, 0);
}

int main(){
	int n_MC = 10000;

	/******************************************************************Forward AD******************************************************************/
	auto start_FAD = get_time::now();

	F<double> S_0_FAD = 100;
	F<double> mu_FAD = 0.05;
	F<double> sigma_FAD = 0.3;
	F<double> dt_FAD = 0.01;
	F<double> T_FAD = 1;
	F<double> K_FAD = 100;
	F<double> Payoff_FAD = 0;

	double Price_FAD = 0;
	double dPdS_FAD = 0;
	double dPdV_FAD = 0;

	for (int i = 0; i < n_MC; i++){  //loop on the Monte Carlo simulations

		Payoff_FAD = MonteCarlo(S_0_FAD, mu_FAD, sigma_FAD, dt_FAD, T_FAD, K_FAD); //computes the price associated to one path

		Price_FAD += Payoff_FAD.x() / n_MC;	//evaluates the value of the price and converts it into a double

		S_0_FAD.diff(0, 6);
		sigma_FAD.diff(2, 6);

		dPdS_FAD += Payoff_FAD.d(0) / n_MC;
		dPdV_FAD += Payoff_FAD.d(2) / n_MC;

	}

	std::cout << "FAD Delta: " << dPdS_FAD << endl;
	std::cout << "FAD Vega: " << dPdV_FAD << endl;

	auto end_FAD = get_time::now();
	auto diff_FAD = end_FAD - start_FAD;
	std::cout << "FAD time is :  " << chrono::duration_cast<chrono::milliseconds>(diff_FAD).count() << " ms " << endl;




	/******************************************************************Backward AD******************************************************************/
	auto start_BAD = get_time::now();
	B<double> S_0_BAD = 100;
	B<double> mu_BAD = 0.05;
	B<double> sigma_BAD = 0.3;
	B<double> dt_BAD = 0.01;
	B<double> T_BAD = 1;
	B<double> K_BAD = 100;
	B<double> Payoff_BAD = 0;

	double dPdS_BAD = 0;
	double dPdV_BAD = 0;
	double Price_BAD = 0;

	for (int i = 0; i < n_MC; i++){  //loop on the Monte Carlo simulations

		Payoff_BAD = MonteCarlo(S_0_BAD, mu_BAD, sigma_BAD, dt_BAD, T_BAD, K_BAD); //computes the price associated to one path

		Payoff_BAD.diff(0, 1);

		Price_BAD += Payoff_BAD.x() / n_MC;
		dPdS_BAD += S_0_BAD.d(0) / n_MC;
		dPdV_BAD += sigma_BAD.d(0) / n_MC;

	}

	std::cout << "BAD Delta: " << dPdS_BAD << endl;
	std::cout << "BAD Vega: " << dPdV_BAD << endl;

	auto end_BAD = get_time::now();
	auto diff_BAD = end_BAD - start_BAD;
	std::cout << "BAD time is :  " << chrono::duration_cast<chrono::milliseconds>(diff_BAD).count() << " ms " << endl;


	/******************************************************************Bump and Revalue******************************************************************/
	auto start_BNV = get_time::now();

	double S_0_double = 100;
	double mu_double = 0.05;
	double sigma_double = 0.3;
	double dt_double = 0.01;
	double T_double = 1;
	double K_double = 100;

	double F_payoff_double_S_Up = 0;
	double F_payoff_double_S_Down = 0;
	double F_payoff_double_V_Up = 0;
	double F_payoff_double_V_Down = 0;

	double Price_S_Up = 0;
	double Price_S_Down = 0;
	double Price_V_Up = 0;
	double Price_V_Down = 0;

	double dPdS_double = 0;
	double dPdV_double = 0;

	for (int i = 0; i < n_MC; i++){  //loop on the Monte Carlo simulations
		F_payoff_double_S_Up = MonteCarlo(S_0_double + 1, mu_double, sigma_double, dt_double, T_double, K_double); //computes the price associated to one path
		Price_S_Up += F_payoff_double_S_Up / n_MC;	//evaluates the value of the price and converts it into a double

		F_payoff_double_S_Down = MonteCarlo(S_0_double - 1, mu_double, sigma_double, dt_double, T_double, K_double); //computes the price associated to one path
		Price_S_Down += F_payoff_double_S_Down / n_MC;	//evaluates the value of the price and converts it into a double

		F_payoff_double_V_Up = MonteCarlo(S_0_double, mu_double, sigma_double + 0.01, dt_double, T_double, K_double); //computes the price associated to one path
		Price_V_Up += F_payoff_double_V_Up / n_MC;	//evaluates the value of the price and converts it into a double

		F_payoff_double_V_Down = MonteCarlo(S_0_double, mu_double, sigma_double - 0.01, dt_double, T_double, K_double); //computes the price associated to one path
		Price_V_Down += F_payoff_double_V_Down / n_MC;	//evaluates the value of the price and converts it into a double

	}

	dPdS_double += (Price_S_Up - Price_S_Down) / 2;
	dPdV_double += (Price_V_Up - Price_V_Down) / 0.02;

	std::cout << "BNV Delta: " << dPdS_double << endl;
	std::cout << "BNV Vega: " << dPdV_double << endl;

	auto end_BNV = get_time::now();
	auto diff_BNV = end_BNV - start_BNV;
	std::cout << "BNV time is :  " << chrono::duration_cast<chrono::milliseconds>(diff_BNV).count() << " ms " << endl;

	return 0;
}