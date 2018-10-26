// IE523: Financial Computation
// Processing tick-data using the Fast Fourier Transform (FFT)
// Written by Prof. R.S. Sreenivas
// Modified from the work of Professor Sreenivas

#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <complex>
#include "newmat.h"
#include "newmatap.h"  // Need this for the FFT algorithm

using namespace std;

// function that is used for sorting a 2-col array based on the value of the entries in the 2nd column.
// taken from http://stackoverflow.com/questions/3041897/sorting-a-2-dimensional-array-on-multiple-columns
bool compareTwoRows2(double* rowA, double* rowB){
    return ( (rowA[1]>rowB[1]) || ((rowA[1]==rowB[1])&&(rowA[0]>rowB[0])) );
}

class Filtering_Instance
{
    // Private
    int no_of_terms, no_of_data_points;
    
    // Private member function that computes the mean of a data array
    double compute_mean(ColumnVector data)
    {
        // write the code to compute the mean of "data"
        // initialize average and sum
        double average = 0.0;
        double sum = 0.0;
        
        //calculate mean
        for (int x = 1; x <= data.nrows(); x++)
        {
            sum += data(x);
//            average = sum / data.nrows();
        }
        
        average = sum/data.nrows();
        return average;
        
    }
    
    // Private member function that computes the magnitude of (an array of) complex #s
    void compute_magnitude(ColumnVector &magnitude, ColumnVector real_part, ColumnVector imag_part)
    {
        // write the code to compute sqrt(real_part(i)^2 + imag_part(i)^2)
//        double sum = 0.0;
        for (int i = 1; i <= real_part.nrows(); i++){
            magnitude(i) = sqrt(pow(real_part(i), 2) + pow(imag_part(i), 2));
        }
        
    }
    
    // Private member function that reads the data from the input file
    // and stores it in "data"
    void get_data(char* file_name, ColumnVector &data)
    {
        // write code that reads the ticker-data from the input file
        // and stores it in "data"
        ifstream input_file(file_name);
        if (input_file.is_open()) // If the input file exists in the local directory
        {
            for (int i = 1; i <= no_of_data_points; i++){
                input_file >> data(i);
            }
        } else {
            cout << "Input file: " << file_name << " does not exist in present folder"  << endl;
        }
    }
    
    // private member function that writes the data file into a file
    void write_data(char* file_name, ColumnVector &data)
    {
        // write code that writes "data" to file_name.
        ofstream output_file(file_name);
        if (output_file.is_open()) // If the output file exists in the local directory
        {
            for (int i = 1; i <= data.nrows(); i++){
                output_file << data(i) << endl;
            }
        } else {
            cout << "Output file: " << file_name << " does not exist in present folder"  << endl;
        }
    }
    
    // private member function that filters data using the FFT
    // The filtered data is computed using the top "no_of_terms"-many
    // magnitude-components of the orginal data
    void filter_the_data(ColumnVector &data, ColumnVector &filtered_data, int no_of_terms)
    {
        ColumnVector fft_real_part(data.nrows()), fft_imag_part(data.nrows());
        ColumnVector mean_adjusted_data(data.nrows()), magnitude(data.nrows());
        
        double mean = compute_mean(data);
        for (int i = 1; i <= data.nrows(); i++)
            mean_adjusted_data(i) = data(i) - mean;
        
        RealFFT(mean_adjusted_data, fft_real_part, fft_imag_part);
        compute_magnitude(magnitude, fft_real_part, fft_imag_part);
        
        // creating a two dimensional array: first col is the index; second col is the
        // magnitude.  The plan is to have this 2-D array sorted using the 2nd col (ie.
        // sorted based on magnitude). Then we pick the top "no_of_terms"-many of these
        // components to reconstitute/reconstruct the signal back.
        
        double** two_dimensional_array = new double*[fft_imag_part.nrows()];
        for (int i = 0; i < fft_imag_part.nrows(); i++)
            two_dimensional_array[i] = new double[2];
        
        for (int i = 0; i < fft_imag_part.nrows(); i++)
        {
            two_dimensional_array[i][0] = i;
            two_dimensional_array[i][1] = magnitude(i+1);
        }
        std::sort(two_dimensional_array,two_dimensional_array+fft_imag_part.nrows(),&compareTwoRows2);
        
        // if do_we_pick_this(i) == 1, then we keep that component for reconstruction
        // of the filtered signal.  The rest of the array-names should be self-explanatory
        ColumnVector do_we_pick_this(fft_imag_part.nrows());
        ColumnVector filtered_fft_real_part(fft_imag_part.nrows());
        ColumnVector filtered_fft_imag_part(fft_imag_part.nrows());
        
        // write the code for picking the top "no_of_terms" many magnitudes
        // put the real-part in "filtered_fft_real_part" and imaginary-part in
        // "filtered_fft_imag_part" -- and reconstruct the filtered signal as
        // shown below.
        
        for (int i = 0; i < no_of_terms; i++){
            if (i < no_of_terms){
                //pick the top no_of_terms index (in descending order of magnitude and marked as 1
                do_we_pick_this(two_dimensional_array[i][0]+1) = 1;
            } else{
                //otherwise
                do_we_pick_this(two_dimensional_array[i][0]+1) = 0;
            }
            
        }
        
        for (int i = 1; i<=fft_imag_part.nrows(); i++){
            if(do_we_pick_this(i)==1){
                // keep the original value for the filtered real part
                filtered_fft_real_part(i) = fft_real_part(i);
                // keep the original value for the filtered imaginary part
                filtered_fft_imag_part(i) = fft_imag_part(i);
            } else{
                //zero out all remaining real part
                filtered_fft_real_part(i) = 0;
                //zero out all remaining imaginary part
                filtered_fft_imag_part(i) = 0;
            }
        }
        // reconstructed signal using just the "no_of_terms"-many top-magnitude
        // components.
        RealFFTI(filtered_fft_real_part, filtered_fft_imag_part, filtered_data);
        
        // write code to add the mean-back to the reconstructed, filtered-signal
        
        for(int i = 1; i <= data.nrows(); i++){
            filtered_data(i) += mean;
        }
    }
    
public:
    // Public member function that reads the incomplete puzzle
    // we are not doing any checks on the input puzzle -- that is,
    // we are assuming they are indeed valid
    void run_the_filter(int argc, char * const argv[])
    {
        sscanf (argv[1], "%d", &no_of_terms);
        sscanf (argv[2], "%d", &no_of_data_points);
        
        std::cout << "Input File Name: " << argv[3] << std::endl;
        std::cout << "Number of data points in the input file = " << no_of_data_points << endl;
        std::cout << "Number of dominant terms in the FFT = " << no_of_terms << endl;
        std::cout << "Output File Name: " << argv[4] << std::endl;
        
        ColumnVector data(no_of_data_points), filtered_data(no_of_data_points);
        
        // get ticker data
        get_data(argv[3], data);
        
        // filter the ticker data
        filter_the_data(data, filtered_data, no_of_terms);
        
        // write the filtered data
        write_data(argv[4], filtered_data);
    }
};


int main (int argc, char* argv[])
{
    Filtering_Instance x;
    x.run_the_filter(argc,argv);
}
