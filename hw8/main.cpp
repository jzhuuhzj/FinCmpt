//
//  main.cpp
//  hw8
//
//  Created by zhujx on 11/6/18.
//  Copyright © 2018 zhujx. All rights reserved.
//

#include <iostream>
#include <unistd.h>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include "newmatap.h"
#include "newmat.h"
#include "newmatio.h"
#include <time.h>
#include <ctime>
#include <random>
#include <vector>



using namespace std;

float get_random(){
    // generate random number between -5 to 5
    float random = rand() % 10 - 5;
    // cout << random << endl;
    return random;
}


Matrix initialize(int no_rows){
    Matrix A(no_rows, no_rows);
    for (int i = 1; i <= no_rows; i++)
    {
        for (int j = 1; j <= no_rows; j++){
            // fill the matrix with random numbers from -5 to 5
            A(i,j) = get_random();
        }
    }
    return A;
}

Matrix repeated_squaring (Matrix A, int exponent, int no_rows){
    // if exponent is 0
    if (exponent == 0){
        IdentityMatrix I(no_rows);
        return I;
    }
    // if exponent is odd
    if ((exponent % 2) == 1)
        return (A * repeated_squaring (A * A, (exponent - 1) / 2, no_rows));
    // if exponent is even
    else
        return (repeated_squaring (A*A, exponent / 2, no_rows));
    
}

Matrix direct_multiplication (Matrix A, int exponent, int no_rows){
    Matrix C(no_rows,no_rows);
    // direct multiplication to calculate matrix's exponent result
    C = A;
    for (int i = 0; i < exponent - 1; i++)
        C = A*C;
    return C;
}


int main(int argc, const char * argv[]){
    
    // initialize random seed
    srand( (unsigned int) time(NULL) );
    
    // initialize
    int exponent = 0;
    int no_rows = 0;
    
    // input 2 arguments
    sscanf(argv[1], "%d", &exponent);
    sscanf(argv[2], "%d", &no_rows);
    
    Matrix A = initialize(no_rows);
    
    cout << "The number of rows/columns in the square matrix is: " << no_rows << endl;
    cout << "The exponent is: " << exponent << endl;
    cout << "Repeated Squaring Result: " << endl;
    
    // calculate repeated squaring result
    float time_before_repeated = clock();
    repeated_squaring(A, exponent, no_rows);
    float time_after_repeated = clock();
    float diff_time_repeated = (float(time_after_repeated)-float(time_before_repeated));
    cout << "It took " << diff_time_repeated/CLOCKS_PER_SEC << " seconds to complete" << endl;
    
    // calculate direct multiplication result
    cout << "Direct Multiplication Result: " << endl;
    float time_before_direct = clock();
    direct_multiplication (A, exponent, no_rows);
    float time_after_direct = clock();
    float diff_time_direct = (float(time_after_direct)-float(time_before_direct));
    cout << "It took " << diff_time_direct/CLOCKS_PER_SEC << " seconds to complete" << endl;
    
    
    // try exponent from 1 to 300
    int count = 300;
    // use vector to save time spent
    vector <float> repeated_squaring_time;
    vector <float> direct_multiplication_time;
    
    ofstream repeated_squaring_file;
    ofstream direct_multiplication_file;
    
    // open files
    repeated_squaring_file.open("repeated_squaring_time");
    direct_multiplication_file.open("direct_multiplication_time");
    
    // calculate repeated squaring result
    for (int i = 1; i <= count; i++) {
        float t1_repeated = clock();
        repeated_squaring(A, i, no_rows);
        float t2_repeated = clock();
        float time_spent_repeated = (float(t2_repeated) - float(t1_repeated))/CLOCKS_PER_SEC;
        repeated_squaring_time.push_back(time_spent_repeated);
    }

    // calculate direct multiplication result
    for (int j = 1; j <= count; j++) {
        float t1_direct = clock();
        direct_multiplication(A, j, no_rows);
        float t2_direct = clock();
        float time_spent_direct= (float(t2_direct) - float(t1_direct))/CLOCKS_PER_SEC;
        direct_multiplication_time.push_back(time_spent_direct);
    }
    
    // write into files
    for (int k = 0; k < count; k++) {
        repeated_squaring_file << repeated_squaring_time[k] << endl;
        direct_multiplication_file << direct_multiplication_time[k] << endl;
    }
    
    
}







