// IE523: Financial Computation
// "How to lose as little as possible" by Addona, Wagon and Wilf
// Written by Prof. Sreenivas
// 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "hint.h"
#include <time.h>

using namespace std;
using namespace std::chrono;
	
int main(int argc, char* argv[])
{

	clock_t t1,t2;
    t1=clock();

	I_have_nothing_apropos_for_this_class x;
	double alice_success_prob, bob_success_prob;
	int number_of_trials;

	sscanf (argv[1], "%lf", &alice_success_prob);
	sscanf (argv[2], "%lf", &bob_success_prob);
	sscanf (argv[3], "%d", &number_of_trials);

	
	cout << "Probability of success for Alice = " << alice_success_prob << endl;
	cout << "Probability of success for Bob = " << bob_success_prob << endl;
	cout << "Number of trials = " << number_of_trials << endl;
	
	x.set_probability(alice_success_prob, bob_success_prob);
	x.set_number_of_trials(number_of_trials);
	
	int optimal = x.search_result();
	
	if (optimal > 0)
		cout << "The optimal number of coin tosses in each game is " << optimal << endl;
	else {
		cout << "The optimal number of coin tosses in each game exceeds 100... Quitting" << endl;
	}

	x.create_output_file();

	t2=clock();
    float seconds = ((float)t2-(float)t1)/CLOCKS_PER_SEC;; 
    cout << "This program took me " << seconds << " seconds to run" << endl;

    return 0;

}
	