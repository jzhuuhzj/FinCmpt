// Written by Prof. Sreenivas for IE523: Financial Computing

// Modified from the work of Prof. Sreenivas

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include "lp_lib.h"

using namespace std;

const double ERROR = 1e-10;
int number_of_cash_flows;
vector <double> price_list;
vector <int> maturity_list;
vector <double> yield_to_maturity;
vector <double> duration;
vector <double> lp_duration;
vector <double> convexity;
vector <double> lp_convexity;
double debt_obligation_amount;
double time_when_debt_is_due;
vector <double> percentage_of_cash_flow_to_meet_debt_obligation;

// From the present value of the input cash flow's equation: The present value of bond = sum of discounted cash flows.
// We then can rewrite it to a function f(r), which can prepare us to solve for r(YTM) using Newton-Raphson Method
double function_r(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that computes f(r) in page 2 of lesson 3 of my notes
    
    // initialize sum
    double sum = 0.0;
    for (int t = 0; t < maturity; t++){
        sum += (cash_flow[t])*pow((1+rate),(maturity-t-1));
    }
    double f = price*pow((1+rate), maturity) - sum;
    return f;
    
}

// Find the derivative of f(r) from above, which we will be using in Newton-Raphson Method
double derivative_function(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that computes f'(r) in the bottom of page 2 of lesson 3
    // of my notes
    
    // initialize sum
    double sum = 0.0;
    for (int t = 1; t <= maturity; t++){
        sum += (cash_flow[t-1])*(maturity - t)*pow((1+rate),(maturity-t-1));
    }
    double f_derivative = maturity*price*pow((1+rate),(maturity-1)) - sum;
    return f_derivative;
}

// Here, we use the Newton-Raphson Method to find out our only positive root.
// We continue the root finding process until we get an answer that is within the desired error range.
double Newton_Raphson(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that finds the (only) +ve root of f(r) of page 2 of
    // lesson 3 using Newton-Raphson method
    while(abs(function_r(cash_flow, price, maturity, rate)) > ERROR){
        rate -= function_r(cash_flow, price, maturity, rate)/derivative_function(cash_flow, price, maturity, rate);
    }
    return rate;
}

// Find the duration of the input cash flow
double get_duration(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that computes the duration of a cash flow
    
    // initialize sum
    double sum = 0.0;
    for (int t = 1; t <= maturity; t++){
         sum += t*cash_flow[t-1]/pow((1+rate),t);
    }
    double duration = sum/price;
    return duration;
}

// Find the convexity of the input cash flow
double get_convexity(vector <double> cash_flow, double price, int maturity, double rate)
{
    // write a function that computes the convexity of a cash flow
    
    // initialize sum
    double sum = 0.0;
    for (int t = 1; t <= maturity; t++){
        sum += t*(t+1)*cash_flow[t-1]/pow((1+rate),(t+2));
    }
    double convexity = sum/price;
    return convexity;
}

// Helper method to find the average of yield to maturities
double mean()
{
    // initialize average and sum
    double average = 0.0;
    double sum = 0.0;
    
    for (int x = 0; x < yield_to_maturity.size(); x++)
    {
        sum += yield_to_maturity[x];
        average = sum / number_of_cash_flows;
    }
    return average;
}

// Using the average of YTMs computed from the mean() function, find out the present value of the future debt obligation
double present_value_of_debt()
{
    // compute PV of future debt obligation
    // using the average-value-of-the-YTMs
    
    //initialize presentValue
    double presentValue = 0.0;
    presentValue = (debt_obligation_amount)/pow((1 + mean()), time_when_debt_is_due);
    return presentValue;
    
}

// print out output
void print_data(char *filename)
{
    // print out the initial setting
    cout << "Input File: " << filename << endl;
    cout << "We owe " << debt_obligation_amount << " in " << time_when_debt_is_due << " years" << endl;
    cout << "Number of Cash Flows: " << number_of_cash_flows << endl;
    
    // print out result for every input cash flow
    for (int i = 0; i < number_of_cash_flows; i++)
    {
        cout << "------------------------------------------------------" << endl;

        cout << "Cash Flow #" << i+1 << endl;
        cout << "Price = " << price_list[i] << endl;
        cout << "Maturity = " << maturity_list[i] << endl;
        cout << "Percentage of Face Value that would meet the obligation = " <<
            percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
        cout << "Yield to Maturity = " << yield_to_maturity[i] << endl;
        cout << "Duration = " << duration[i] << endl;
        cout << "Duration (to be used in LP-formulation below) = " << lp_duration[i] << endl;
        cout << "(Note) " << duration[i] << " = " << lp_duration[i] << " x " << percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
        cout << "Convexity = " << convexity[i] << endl;
        cout << "Convexity (to be used in LP-formulation below) = " << lp_convexity[i] << endl;
        cout << "(Note) " << convexity[i] << " = " << lp_convexity[i] << " x " << percentage_of_cash_flow_to_meet_debt_obligation[i] << endl;
        
    }
    
    cout << "******************************************************" << endl;
    cout << "Average YTM (which I use to compute PV of Debt = " << mean() << endl;
    cout << "Present value of debt = " << present_value_of_debt() << endl;
    cout << "******************************************************" << endl;
}

// get data from input file
void get_data(char* argv[])
{
    // write the code that reads the data from the file identified
    // on the command-line.
    ifstream input_file(argv[1]);
    cout << std::endl << "Input File Name:" << argv[1]<< endl;
    
    // initialize
    double price = 0.0;
    int maturity = 0;
    double cash_flow = 0.0;
    double ytm = 0.0;
    double rate = 0.0;
    double dur = 0.0;
    double cov = 0.0;
    
    //if the input file exists
    if (input_file.is_open()){
        // get number of cash flows
        input_file >> number_of_cash_flows;
        for (int i = 0; i < number_of_cash_flows; i++){
            vector <double> cash_flow_list;
            
            // get next input：price
            input_file >> price;
            // push back to price list
            price_list.push_back(price);
            
            // get next input: maturity
            input_file >> maturity;
            // push back to maturity list
            maturity_list.push_back(maturity);
            
            // for every maturity, get the corresponding cash flows
            for (int j = 0; j < maturity; j++){
                input_file >> cash_flow;
                // push back to the cash flow list
                cash_flow_list.push_back(cash_flow);
            }
            
            // get ytm computed from Newton_Raphson Method
            ytm = Newton_Raphson(cash_flow_list, price, maturity, rate);
            // then push back to the list of yield_to_maturities
            yield_to_maturity.push_back(ytm);
            // get duration
            dur = get_duration(cash_flow_list, price, maturity, ytm);
            // then push back to the list of durations
            duration.push_back(dur);
            // get convexity
            cov = get_convexity(cash_flow_list, price, maturity, ytm);
            // then push back to the lit of convexities
            convexity.push_back(cov);
            
        }
        // get the debt obligation amount from the input file
        input_file >> debt_obligation_amount;
        // get timewhen debt is due from the input file
        input_file >> time_when_debt_is_due;
        
        for (int i = 0; i < number_of_cash_flows; i++){
            // compute the percentage of cash flow to meet the debt obligation amount
            double percent = present_value_of_debt()/price_list[i];
            // push back to the list of percentages
            percentage_of_cash_flow_to_meet_debt_obligation.push_back(percent);
            // compute the duration that will be used in lp_solve and then push back
            lp_duration.push_back(duration[i]/percentage_of_cash_flow_to_meet_debt_obligation[i]);
            // compute the convexity that will be used in lp_solve and then puah back
            lp_convexity.push_back(convexity[i]/percentage_of_cash_flow_to_meet_debt_obligation[i]);
        }
        
    } else{
        //if the input file does not exist, then exit the program
        cout << "Error: Input file does not exist in PWD" << endl;
        exit(0);
    }
   
}

void get_optimal_portfolio()
{
    lprec *lp;
    // write the lp_solve specific material that
    // computes the optimal_portfolio
    
    // initialize a solution array
    REAL solution[number_of_cash_flows];
    
    lp = make_lp(0, number_of_cash_flows);
    
    // keep the message reporting of lp_solve to a minimum
    set_verbose(lp, 3);
    
    // add constraint 1: all discounted cash flows will be equal to the present value of debt
    {
        double row[number_of_cash_flows+1];
        for (int i = 0; i < number_of_cash_flows+1; i++){
            row[i+1] = price_list[i];
        }
        row[0] = 0;
        add_constraint(lp, row, EQ, present_value_of_debt());
    }
    
    // add constraint 2: sum of (lambda * D) = N
    {
        double row[number_of_cash_flows+1];
        for (int i = 0; i < number_of_cash_flows+1; i++){
            row[i+1] = lp_duration[i];
        }
        row[0] = 0;
        add_constraint(lp, row, EQ, time_when_debt_is_due);
    }
    
    // set objective function
    {
        double row[number_of_cash_flows+1];
        for (int i = 0; i < number_of_cash_flows;i++){
            row[i+1]= -lp_convexity[i];
        }
        set_obj_fn(lp, row);
    }
    
    // print lp
    print_lp(lp);
    
    // solve the lp
    double return_value = solve(lp);
    
    // get the optimizing values of the variables
    get_variables(lp, solution);
    
    // if the constraints are met
    if (return_value == 0){
        // print optimal value
        cout << "Largest convexiry we can get is: " << -get_objective(lp) << endl;
        
        // print out the optimal portfolio results
        cout << "Optimal portfolio:" << endl;
        for (int i = 0; i < number_of_cash_flows; i++){
            cout << "%Cash Flow:" <<  i+1 << "  " << solution[i] << " ";
            cout << endl;
        }
        
        cout << endl << "That is, buy" << endl;
        
        // print out the actual dollar amount to buy for the optimal portfolio
        for (int i = 0; i < number_of_cash_flows; i++){
            if(solution[i] != 0){
                cout << "$" << price_list[i] * solution[i] << " of Cash Flow#" << i+1 << endl;
            }
        }
    } else{
        // if constraint not met
        cout << "There is no portfolio that meets the duration constraint of " <<
            time_when_debt_is_due << " years" << endl;
    }
    
    // delete the lp to release memory
    delete_lp(lp);
}

int main (int argc, char* argv[])
{
    if (argc == 1) {
        cout << "Input filename missing" << endl;
    }
    else
    {
        get_data(argv);
        
        print_data(argv[1]);
        
        get_optimal_portfolio();
    }
    return (0);

}
