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
double put_option_price_brownian_bridge = 0.0;
double call_option_price_brownian_bridge = 0.0;

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


double barrier_option(double S, double barrier_price){
    if (S <= barrier_price){
        S = 0;
        return S;
    } else{
        return S;
    }
}

double prob(double S, double barrier_price, double i){
    if (S <= barrier_price){  
        S = 0;
        return S;
    } else{
        double mu = initial_stock_price + ((double)i / (double)no_of_divisions)*(S - initial_stock_price);
        double sigma = ((double)i / (double)no_of_divisions)*((double)no_of_divisions - (double)i);
        return (1 - N((barrier_price - mu) / sqrt(sigma)));
    }
}


//get-four-paths-for-the-price-of-one-path
//The greeks by automatic differentiation of the simulation code
void explicit_simulation(){
    double S1 = 0.0;
    double S2 = 0.0;
    double S3 = 0.0;
    double S4 = 0.0;

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

            barrier_option(S1, barrier_price);
            barrier_option(S2, barrier_price);
            barrier_option(S3, barrier_price);
            barrier_option(S4, barrier_price);

            S1 = S1 * exp(R + SD*a);
            S2 = S2 * exp(R - SD*a);
            S3 = S3 * exp(R + SD*b); 
            S4 = S4 * exp(R - SD*b);

        }
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

void brownian_bridge(){
    double S1 = 0.0;
    double S2 = 0.0;
    double S3 = 0.0;
    double S4 = 0.0;

    double mu1 = 0.0;
    double mu2 = 0.0;
    double mu3 = 0.0;
    double mu4 = 0.0;

    double sigma1 = 0.0;
    double sigma2 = 0.0;
    double sigma3 = 0.0;
    double sigma4 = 0.0;
    
    double R = (risk_free_rate - 0.5 * pow(volatility, 2)) * expiration_time;
    double SD = volatility * sqrt(expiration_time);

    for (int i = 1; i <= no_of_trials; i++)
    {
        double S1 = initial_stock_price;
        double S2 = initial_stock_price;
        double S3 = initial_stock_price;
        double S4 = initial_stock_price;

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

        barrier_option(S1, barrier_price);
        barrier_option(S2, barrier_price);
        barrier_option(S3, barrier_price);
        barrier_option(S4, barrier_price);

        double p_d_1 = 1.0;
        double p_d_2 = 1.0;
        double p_d_3 = 1.0;
        double p_d_4 = 1.0;
        
        for (int i = 0; i < no_of_divisions; i++)
        {

            p_d_1 *= prob(S1, barrier_price, i);
            p_d_2 *= prob(S2, barrier_price, i);
            p_d_3 *= prob(S3, barrier_price, i);
            p_d_4 *= prob(S4, barrier_price, i);

        }

        call_option_price_brownian_bridge += (p_d_1 * ((S1>0) ? max(0.0, S1 - strike_price) : 0) + 
                                  p_d_2 * ((S2>0) ? max(0.0, S2 - strike_price) : 0) + 
                                  p_d_3 * ((S3>0) ? max(0.0, S3- strike_price) : 0) + 
                                  p_d_4 * ((S4>0) ? max(0.0, S4 - strike_price) : 0))/4.0;
        
        put_option_price_brownian_bridge += (p_d_1 * ((S1>0) ? max(0.0, strike_price - S1) : 0) + 
                                  p_d_2 * ((S2>0) ? max(0.0, strike_price - S2) : 0) + 
                                  p_d_3 * ((S3>0) ? max(0.0, strike_price - S3) : 0) + 
                                  p_d_4 * ((S4>0) ? max(0.0, strike_price - S4) : 0))/4.0;   
      
    }          
    call_option_price_brownian_bridge = exp(-risk_free_rate * expiration_time) * (call_option_price_brownian_bridge / ((double)no_of_trials));
    put_option_price_brownian_bridge = exp(-risk_free_rate * expiration_time) * (put_option_price_brownian_bridge / ((double)no_of_trials));

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
    brownian_bridge();
    
    cout << "--------------------------------------" << endl;
    cout << "European Down-and-out Discrete Barrier Options Pricing via Monte Carlo Simulation" << endl;
    cout << "Expiration time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike price = " << strike_price << endl;
    cout << "Barrier price = " << barrier_price << endl;
    cout << "Number of Trials = " << no_of_trials << endl;
    cout << "Number of Discrete Barriers = " << no_of_divisions << endl;
    cout << "--------------------------------------" << endl;
    cout << "The average Call Price via explicit simulation of price paths             = " << call_option_price_explicit << endl;
    cout << "The average Call Price with Brownian-Bridge correction on the final price = " << call_option_price_brownian_bridge << endl;
    cout << "The average Put Price via explicit simulation of price paths              = " << put_option_price_explicit << endl;
    cout << "The average Put Price with Brownian-Bridge correction on the final price  = " << put_option_price_brownian_bridge << endl;
    cout << "--------------------------------------" << endl;
    
}
