#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <vector>


using namespace std;

//vector of vector 
vector<vector<double> > american_call;
vector<vector<double> > american_put;

double up_factor, uptick_prob, downtick_prob, notick_prob, risk_free_rate, strike_price;
double initial_stock_price, expiration_time, volatility, R;
int no_of_divisions;

//initialize a matrix used to store recursion result
void initialize()
{
    for (int i = 0; i < no_of_divisions + 1; i++)
    {
        vector<double> temp;
        american_call.push_back(temp);
        american_put.push_back(temp);
    }

    for (int i = 0; i < no_of_divisions + 1; i++)
    {
        for (int j = 0; j < 2 * no_of_divisions + 1; j++)
        {
            american_call[i].push_back(-1);
            american_put[i].push_back(-1);
        }
    }
}


double max(double a, double b) {
    return (b < a )? a:b;
}

//for american call
double american_call_option(int k, int i, double current_stock_price){
    //max(S - K)
    if (american_call[k][no_of_divisions + i] == -1){
        if (k == no_of_divisions){
            american_call[k][no_of_divisions + i] = max(0.0, (current_stock_price - strike_price));
        } else{
             american_call[k][no_of_divisions + i] = max((current_stock_price - strike_price),
                   (uptick_prob*american_call_option(k+1, i+1, current_stock_price*up_factor)+
                    (1-uptick_prob-downtick_prob)*american_call_option(k+1, i, current_stock_price)+
               (downtick_prob*american_call_option(k+1, i-1, current_stock_price/up_factor)))/R);
        }
    }
        return american_call[k][no_of_divisions + i];
}

//for american put
double american_put_option(int k, int i, double current_stock_price){
    //max(-S + K)
    if (american_put[k][no_of_divisions + i] == -1){
        if (k == no_of_divisions){
            american_put[k][no_of_divisions + i] = max(0.0, (-current_stock_price + strike_price));
        } else{
            american_put[k][no_of_divisions + i] = max((-current_stock_price + strike_price),
                       (uptick_prob*american_put_option(k+1, i+1, current_stock_price*up_factor)+
                        (1-uptick_prob-downtick_prob)*american_put_option(k+1, i, current_stock_price)+
                        (downtick_prob*american_put_option(k+1, i-1, current_stock_price/up_factor)))/R);
        }
    }
    return american_put[k][no_of_divisions + i];
}


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
    cout << "(Memoized) Recursive Trinomial American Option Pricing" << endl;
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
    cout << "Notick Probability = " << notick_prob << endl;
    cout << "Downtick Probability = " << downtick_prob << endl;
    cout << "--------------------------------------" << endl;
    double american_call_price = american_call_option(0, 0,initial_stock_price);
    cout << "Trinomial Price of an American Call Option = " << american_call_price << endl;
    double american_put_price = american_put_option(0, 0, initial_stock_price);
    cout << "Trinomial Price of an American Put Option = " << american_put_price << endl;
    cout << "--------------------------------------" << endl;

    
}
