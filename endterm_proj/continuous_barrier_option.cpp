// Adjusted from the code examples of Prof. Sreenivas
// Jingxia Zhu

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include "normdist.h"
using namespace std;

double risk_free_rate, strike_price, barrier_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions, no_of_trials;

double put_option_price_explicit = 0.0;
double call_option_price_explicit = 0.0;
double put_option_price_adjusted = 0.0;
double call_option_price_adjusted = 0.0;
double call_option_price_theo = 0.0;
double put_option_price_theo = 0.0;


unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();

std::default_random_engine generator(seed);

double get_uniform() 
{
    std::uniform_real_distribution <double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return (number);
}

double max(double a, double b) {
    return (b < a )? a:b;
}

double option_price_put_black_scholes(const double& S,      // spot price
                                      const double& K,      // Strike (exercise) price,
                                      const double& r,      // interest rate
                                      const double& sigma,  // volatility
                                      const double& time)
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt) + 0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
                                       const double& K,       // strike (exercise) price,
                                       const double& r,       // interest rate
                                       const double& sigma,   // volatility 
                                       const double& time)    // time to maturity 
{  
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+r*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt; 
    double d2 = d1-(sigma*time_sqrt);
    return S*N(d1) - K*exp(-r*time)*N(d2);
};

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

double closed_form_down_and_out_european_call_option()
{
    // I took this formula from Wilmott, Howison and Dewynne, "The Mathematics of Financial Derivatives"
    double K = (2*risk_free_rate)/(volatility*volatility);
    double A = option_price_call_black_scholes(initial_stock_price, strike_price, 
                                              risk_free_rate, volatility, expiration_time);
    double B = (barrier_price*barrier_price)/initial_stock_price;
    double C = pow(initial_stock_price/barrier_price, -(K-1));
    double D = option_price_call_black_scholes(B, strike_price, risk_free_rate, volatility, expiration_time);
    return (A - D*C);
}

double closed_form_down_and_in_european_put_option() 
{
    // just making it easier by renaming the global variables locally
    double S = initial_stock_price;
    double r = risk_free_rate;
    double T = expiration_time;
    double sigma = volatility;
    double H = barrier_price;
    double X = strike_price;
    
    // Took these formulae from some online reference
    double lambda = (r+((sigma*sigma)/2))/(sigma*sigma);
    double temp = 2*lambda - 2.0;
    double x1 = (log(S/H)/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    double y = (log(H*H/(S*X))/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    double y1 = (log(H/S)/(sigma*sqrt(T))) + (lambda*sigma*sqrt(T));
    return (-S*N(-x1) + X*exp(-r*T)*N(-x1 + sigma*sqrt(T)) + 
            S*pow(H/S, 2*lambda)*(N(y)-N(y1)) -
            X*exp(-r*T)*pow(H/S, temp)*(N(y-sigma*sqrt(T))-N(y1-sigma*sqrt(T))));
}

double closed_form_down_and_out_european_put_option()
{
    double vanilla_put = option_price_put_black_scholes(initial_stock_price, strike_price, 
                                                       risk_free_rate, volatility, expiration_time);
    double put_down_in = closed_form_down_and_in_european_put_option();
    return (vanilla_put - put_down_in);
}


// return S
double barrier_option(double S, double barrier_price){
    if (S <= barrier_price){
        return S = 0;
    } else{
        return S;
    }
}

// return S or p_c
double prob(double S, double barrier_price){
    if (S <= barrier_price){
        S = 0;
        return S;
    } else{
        return (double)exp(-(2 * log(initial_stock_price / barrier_price)*log(S / barrier_price)) / (expiration_time*pow(volatility, 2)));
    }
}


//get-four-paths-for-the-price-of-one-path
//The greeks by automatic differentiation of the simulation code
void explicit_simulation(){
    double S1 = 0.0;
    double S2 = 0.0;
    double S3 = 0.0;
    double S4 = 0.0;

    double S1_call_val;
    double S2_call_val;
    double S3_call_val;
    double S4_call_val;


    double R = (risk_free_rate - 0.5*pow(volatility,2))*(expiration_time/(double)no_of_divisions);
    double SD = volatility*sqrt(expiration_time/(double)no_of_divisions); 

    for (int i = 0; i < no_of_trials; i++) {
            S1 = initial_stock_price;
            S2 = initial_stock_price;
            S3 = initial_stock_price;
            S4 = initial_stock_price;
            

        for (int j = 0; j < no_of_divisions; j++){

            // generate unit-normals using Box-Muller Transform
            double x = get_uniform();
            double y = get_uniform();
            double a =  sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
            double b =  sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

            // check at time t (0<= t < T)
            barrier_option(S1, barrier_price);
            barrier_option(S2, barrier_price);
            barrier_option(S3, barrier_price);
            barrier_option(S4, barrier_price);

            S1 = S1 * exp(R + SD*a);
            S2 = S2 * exp(R - SD*a);
            S3 = S3 * exp(R + SD*b); 
            S4 = S4 * exp(R - SD*b);

        }

        // check at time T
        barrier_option(S1, barrier_price);
        barrier_option(S2, barrier_price);
        barrier_option(S3, barrier_price);
        barrier_option(S4, barrier_price);


        call_option_price_explicit += (((S1>0) ? max(0.0, S1 - strike_price):0) + 
                                  ((S2>0) ? max(0.0, S2 - strike_price) : 0) + 
                                  ((S3>0) ? max(0.0, S3 - strike_price) : 0) + 
                                  ((S4>0) ? max(0.0, S4 - strike_price) : 0))/4.0;

        put_option_price_explicit += (((S1>0) ? max(0.0, strike_price - S1):0) + 
                                 ((S2>0) ? max(0.0, strike_price - S2) : 0) + 
                                 ((S3>0) ? max(0.0, strike_price - S3) : 0) + 
                                 ((S4>0) ? max(0.0, strike_price - S4) : 0))/4.0;
    }
    call_option_price_explicit = exp(-risk_free_rate*expiration_time)*(call_option_price_explicit/((double) no_of_trials));
    put_option_price_explicit = exp(-risk_free_rate*expiration_time)*(put_option_price_explicit/((double) no_of_trials));
}

void prob_correction(){
    double S1 = 0.0;
    double S2 = 0.0;
    double S3 = 0.0;
    double S4 = 0.0;
    
    double R = (risk_free_rate - 0.5 * pow(volatility, 2)) * expiration_time;
    double SD = volatility * sqrt(expiration_time);

    for (int i = 0; i < no_of_trials; i++)
    {
        double S1 = initial_stock_price;
        double S2 = initial_stock_price;
        double S3 = initial_stock_price;
        double S4 = initial_stock_price;

        // check at time t ( 0 <= t < T)
        barrier_option(S1, barrier_price);
        barrier_option(S2, barrier_price);
        barrier_option(S3, barrier_price);
        barrier_option(S4, barrier_price);

        double x = get_uniform();
        double y = get_uniform();
        double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
        double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);
        
        S1 = S1*exp(R + SD*a);
        S2 = S2*exp(R - SD*a);
        S3 = S3*exp(R + SD*b);
        S4 = S4*exp(R - SD*b);

        // check at time T
        barrier_option(S1, barrier_price);
        barrier_option(S2, barrier_price);
        barrier_option(S3, barrier_price);
        barrier_option(S4, barrier_price);
       
        // calculate call option value
        call_option_price_adjusted += (((1-prob(S1, barrier_price)) * ((S1>0) ? max(0.0, S1 - strike_price) : 0)) + 
                                  ((1-prob(S2, barrier_price)) * ((S2>0) ? max(0.0, S1 - strike_price) : 0)) + 
                                  ((1-prob(S3, barrier_price)) * ((S3>0) ? max(0.0, S1 - strike_price) : 0)) + 
                                  ((1-prob(S4, barrier_price)) * ((S4>0) ? max(0.0, S1 - strike_price) : 0)))/4.0;

        // calculate put option value
        put_option_price_adjusted += (((1-prob(S1, barrier_price)) * ((S1>0) ? max(0.0, strike_price - S1) : 0)) + 
                                  ((1-prob(S2, barrier_price)) * ((S2>0) ? max(0.0, strike_price - S2) : 0)) + 
                                  ((1-prob(S3, barrier_price)) * ((S3>0) ? max(0.0, strike_price - S3) : 0)) + 
                                  ((1-prob(S4, barrier_price)) * ((S4>0) ? max(0.0, strike_price - S4) : 0)))/4.0; 
    }    
    // give the call and put option price after no_of trials          
    call_option_price_adjusted = exp(-risk_free_rate * expiration_time) * (call_option_price_adjusted / ((double)no_of_trials));
    put_option_price_adjusted = exp(-risk_free_rate * expiration_time) * (put_option_price_adjusted / ((double)no_of_trials));

}



int main (int argc, char* argv[])
{
    
    sscanf (argv[1], "%lf", &expiration_time);
    sscanf (argv[2], "%lf", &risk_free_rate);
    sscanf (argv[3], "%lf", &volatility);
    sscanf (argv[4], "%lf", &initial_stock_price);
    sscanf (argv[5], "%lf", &strike_price);
    sscanf (argv[6], "%d", &no_of_trials);
    sscanf (argv[7], "%d", &no_of_divisions);
    sscanf (argv[8], "%lf", &barrier_price);

    explicit_simulation();
    prob_correction();
    
    cout << "--------------------------------------" << endl;
    cout << "European Down-and-out Continuous Barrier Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike price = " << strike_price << endl;
    cout << "Barrier price = " << barrier_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "--------------------------------------" << endl;
    cout << "--------------------------------------" << endl;
    cout << "The average Call Price via explicit simulation = " << call_option_price_explicit << endl;
    cout << "The call price using (1-p)-adjusted term       = " << call_option_price_adjusted << endl;
    cout << "The theoratical Call price                     = " << closed_form_down_and_out_european_call_option() << endl;
    cout << "--------------------------------------" << endl;
    cout << "The average Put Price via explicit simulation  = " << put_option_price_explicit << endl;
    cout << "The Put price using (1-p)-adjusted term        = " << put_option_price_adjusted << endl;
    cout << "The Theoratical Put price                      = " << closed_form_down_and_out_european_put_option() << endl;
    cout << "--------------------------------------" << endl;
    
}
