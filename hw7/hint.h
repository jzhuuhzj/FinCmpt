/*
 *  alice_and_bob.h
 *  Loosing as little as possible
 *
 *  Created by Ramavarapu Sreenivas on 9/2/12.
 *  Copyright 2012 University of Illinois at Urbana-Champaign. All rights reserved.
 *
 */
#ifndef ALICE_AND_BOB
#define ALICE_AND_BOB

#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

class I_have_nothing_apropos_for_this_class
{
private:
	double alice_probability, bob_probability;
	int number_of_trials;
	
	// private member function: uniform RV generator
	double get_uniform()
	{
		// write the appropriate code here
		return (((double) random())/(pow(2.0,31.0)-1.0));
	}
	
	// private member function: nCi (i.e. n-take-i) 
	int take(int n, int i)
	{
		// write a **RECURSIVE** implementation of n-take-i. 
		// If you made it non-recurisive (i.e. n!/((n-i)!i!)) -- it 
		// will take too long for large sizes 
		
		// you cannot take a greater number of quantity from a smaller number of quantity
		// return 0
		if (i > n) {
			return 0;
		}
		// if you take 0 or i-many quantity from i-many quantity, you only have 1 way to do that
		if (i == n || i == 0 ){
			return 1;
		}
		// if you take 1 from n then you have n-many ways to do it
		if (i == 1){
			return n;
		}
		// C(n,i) = ((n-i+1)*(n-i+2)*....*(n))/i!
		return ((n - i + 1) * take(n, i - 1)) / i;

		// write a **RECURSIVE** implementation of n-take-i. 
		// If you made it non-recurisive (i.e. n!/((n-i)!i!)) -- it 
		// will take too long for large sizes 
	
	}
	
	// this routine implements the probability that Alice has more 
	// heads than Bob after n-many coin tosses
	double theoretical_value(double q, double p, int n)
	{
		// implement equation 1.1 of Addona-Wagon-Wilf paper
		double theo_value = 0.0;
		for (int r = 0; r < n; r++){
			double sum = 0.0;
			for (int s = r + 1; s <= n; s++){
				sum += take(n, s)*pow(q, s)*pow(1 - q, n - s);
			}
			theo_value += take(n, r)*pow(p, r)*pow(1 - p, n - r)*sum; 
		}
		return theo_value;
	}

public: 
	// public function: 
	static const int count = 30;

	void set_probability(double alice_p, double bob_p)
	{
		alice_probability = alice_p;
		bob_probability = bob_p;
	}

	void set_number_of_trials(int num_trials){
		number_of_trials = num_trials;
	}
	
	// probability of Alice winning the game.
	double simulated_value(int number_of_coin_tosses_in_each_game, int no_of_trials)
	{
		int no_of_wins_for_alice = 0;
		for (int i = 0; i < no_of_trials; i++) 
		{
			int number_of_heads_for_alice = 0;
			int number_of_heads_for_bob = 0;
			for (int j = 0; j < number_of_coin_tosses_in_each_game; j++) 
			{
				if (get_uniform() < alice_probability) 
					number_of_heads_for_alice++;
				if (get_uniform() < bob_probability)
					number_of_heads_for_bob++;
			}
			if (number_of_heads_for_alice > number_of_heads_for_bob)
				no_of_wins_for_alice++;
		}
		return (((double) no_of_wins_for_alice)/((double) no_of_trials));
	}
		


	int search_result()
	{
		// implememt a discrete-search procedure for the optimal n-value. 
		// start with n = 1 and find the discrete-value of n that has 
		// the largest probability for Alice winning.  Why would this work?
		// See Theorem 2.2 of the paper for the reason!
		
		int result = 0;
		for (int n = 1; n <= count; n++){
			// if f(n) >= f(n-1)
			if ((theoretical_value(alice_probability, bob_probability, n) 
					>= theoretical_value(alice_probability, bob_probability, n - 1)) &&
					// f(n+1) <= f(n)
					(theoretical_value(alice_probability, bob_probability, n + 1) 
						<= theoretical_value(alice_probability, bob_probability, n))){

				result = n;
				break;
			}
		}
		
		return result;
	}

	void create_output_file()
	{
		vector<double> theo_values;
		vector<double> simu_values;
		
		for (int i = 1; i <= count; i++)
		{
			// push back theoretical values
			theo_values.push_back(theoretical_value(alice_probability, bob_probability, i));
			// push back simulated values
			simu_values.push_back(simulated_value(i, number_of_trials));	
		}

		ofstream theoretical_values_file;
		ofstream simulated_values_file;

		// open files
		theoretical_values_file.open("theoretical_values");
		simulated_values_file.open("simulated_values");

		// write data into files
		for (int i = 0; i < count; i++) {
			theoretical_values_file << theo_values[i] << endl;
			simulated_values_file << simu_values[i] << endl;

		}
	}

};
#endif








