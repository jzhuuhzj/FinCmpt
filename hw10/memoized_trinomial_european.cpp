#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "normdist.h"

using namespace std;

//vector of vector 
vector<vector<double> > european_call;
vector<vector<double> > european_put;

double up_factor, uptick_prob, downtick_prob, notick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

//initialize a matrix used to store recursion result
void initialize()
{
    for (int i = 0; i < no_of_divisions + 1; i++)
    {
        vector<double> temp;
        european_call.push_back(temp);
        european_put.push_back(temp);
    }

    for (int i = 0; i < no_of_divisions + 1; i++)
    {
        for (int j = 0; j < 2 * no_of_divisions + 1; j++)
        {
            european_call[i].push_back(-1);
            european_put[i].push_back(-1);
        }
    }
}


//get max
double max(double a, double b) {
    return (b < a )? a:b;
}

//for european call
double european_call_option(int k, int i, double current_stock_price) {
	//max(S - K)
    if (european_call[k][no_of_divisions + i] == -1){
        if (k == no_of_divisions){
            european_call[k][no_of_divisions + i] = max(0.0, (current_stock_price - strike_price));
        } else{
             european_call[k][no_of_divisions + i] = ((uptick_prob*european_call_option(k+1, i+1, current_stock_price*up_factor)+
                    (1-uptick_prob-downtick_prob)*european_call_option(k+1, i, current_stock_price)+
               (downtick_prob*european_call_option(k+1, i-1, current_stock_price/up_factor)))/R);
        }
    }
    
    return european_call[k][no_of_divisions + i];
}

//for european put
double european_put_option(int k, int i, double current_stock_price) {
    //max(S - K)
    if (european_put[k][no_of_divisions + i] == -1){
        if (k == no_of_divisions){
            european_put[k][no_of_divisions + i] = max(0.0, (-current_stock_price + strike_price));
        } else{
             european_put[k][no_of_divisions + i] = ((uptick_prob*european_put_option(k+1, i+1, current_stock_price*up_factor)+
                    (1-uptick_prob-downtick_prob)*european_put_option(k+1, i, current_stock_price)+
               (downtick_prob*european_put_option(k+1, i-1, current_stock_price/up_factor)))/R);
        }
    }
    
    return european_put[k][no_of_divisions + i];
}

double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time){
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility
                                       const double& time) {  // time to maturity
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

//normal distribution
double N(const double& z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a=fabs(z);
    double t = 1.0/(1.0+a*p);
    double b = c2*exp((-z)*(z/2.0));
    double n = ((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
    n = 1.0-b*n;
    if ( z < 0.0 ) n = 1.0 - n;
    return n;
};


int main (int argc, char* argv[])
{
    
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%d", &no_of_divisions);
    sscanf (argv[3], "%lf", &risk_free_rate);
    sscanf (argv[4], "%lf", &volatility);
    sscanf (argv[5], "%lf", &initial_stock_price);
    sscanf (argv[6], "%lf", &strike_price);
    
    up_factor = exp(volatility*sqrt(2*expiration_time/((double) no_of_divisions)));
    R = exp(risk_free_rate*expiration_time/((double) no_of_divisions));
    uptick_prob = pow((sqrt(R) - (1/sqrt(up_factor)))/(sqrt(up_factor) - (1/sqrt(up_factor))), 2);
    downtick_prob = pow((sqrt(up_factor) - sqrt(R))/(sqrt(up_factor) - (1/sqrt(up_factor))), 2);
    notick_prob = 1 - uptick_prob - downtick_prob;
    
    //initialize the matrix
    initialize();

    //output
    cout << "(Memoized) Recursive Trinomial European Option Pricing" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "--------------------------------------" << endl;
    //cout << "R = " << R << endl;
    cout << "Up Factor = " << up_factor << endl;
    cout << "Uptick Probability = " << uptick_prob << endl;
    //cout << "Downtick Probability = " << downtick_prob << endl;
    //cout << "Notick Probability = " << notick_prob << endl;
    cout << "--------------------------------------" << endl;
    double european_call_price = european_call_option(0, 0, initial_stock_price);
    cout << "Trinomial Price of an European Call Option = " << european_call_price << endl;
    cout << "Call Price according to Black-Scholes = " << option_price_call_black_scholes(initial_stock_price, 
                                                                                            strike_price, 
                                                                                            risk_free_rate,
                                                                                            volatility, expiration_time) << endl;
    cout << "--------------------------------------" << endl;
    double european_put_price = european_put_option(0, 0, initial_stock_price);
    cout << "Trinomial Price of an European Put Option = " << european_put_price << endl;
    cout << "Put Price according to Black-Scholes = " << option_price_put_black_scholes(initial_stock_price, 
                                                                                            strike_price, 
                                                                                            risk_free_rate,
                                                                                            volatility, expiration_time) << endl;

    cout << "--------------------------------------" << endl;
    cout << "Verifying the Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
    cout <<  initial_stock_price << " + " << european_put_price << " - " << european_call_price;
    cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
    cout << initial_stock_price + european_put_price - european_call_price << " = " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
    cout << "--------------------------------------" << endl;

    
}
