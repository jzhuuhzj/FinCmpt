/*
 *  sudoku.h
 *  Sudoku
 *  Created by Prof. Ramavarapu Sreenivas 
 *  Inspired by: http://web.eecs.utk.edu/courses/spring2012/cs140/Notes/Sudoku/index.html
 *  Modified from the work of Prof. Sreenivas
*/
#ifndef sudoku
#define sudoku

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using std::vector;
using namespace std;

class Sudoku 
{ 
	// Private
	int puzzle[9][9];
	int count;
	
	// Private member function that checks if the named row is valid
	bool row_valid(int row)
	{
		// write code that checks if "row" is valid

		//create a new array called row_count[] with size 10, and initialize to all 0s
        int row_count[10] = {0,0,0,0,0,0,0,0,0,0};

		//interate through all 9 numbers of a row 
        for (int i = 0; i < 9; i++){
        	//for any number that is not 0 but within [1,9]
            if (puzzle[row][i] >= 1 && puzzle[row][i] <= 9){
            	//find the corresponding index in row_count[] and increment the value there by 1
                row_count[puzzle[row][i]]++;
                //if a particular number appears more than once, return false
                if (row_count[puzzle[row][i]] > 1) {
                    return false;
                }
            }
        }
        return true;
	}
	
	// Private member function that checks if the named column is valid
	bool col_valid(int col)
	{
		// check validity of "col" 

		//create a new array called col_count[] with size 10, and initialize to all 0s
        int col_count[10] = {0,0,0,0,0,0,0,0,0,0};

        //interate through all 9 numbers of a column
        for (int i = 0; i < 9; i++) {
        	//for any number that is not 0 but within [1,9]
            if (puzzle[i][col] >= 1 && puzzle[i][col] <= 9){ 
            	//find the corresponding index in col_count[] and increment the value there by 1
                col_count[puzzle[i][col]]++;
                //if a particular number appears more than once, return false
                if (col_count[puzzle[i][col]] > 1) {
                    return false;
                }
            }
        }
        return true;
	}
	
	// Private member function that checks if the named 3x3 block is valid
	bool block_valid(int row, int col)
	{
		// check 3 x 3 block validity

		//create a new array called col_count[] with size 10, and initialize to all 0s
		int block_count[10] = {0,0,0,0,0,0,0,0,0,0};

		//interate through all 9 numbers of a 3*3 block
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				//for any number that is not 0 but within [1,9]
				if (puzzle[3*(row/3)+i][3*(col/3)+j] >= 1 && puzzle[3*(row/3)+i][3*(col/3)+j] <= 9){ 
					//find the corresponding index in col_count[] and increment the value there by 1
                	block_count[puzzle[3*(row/3)+i][3*(col/3)+j]]++;
                	//if a particular number appears more than once, return false
                	if (block_count[puzzle[3*(row/3)+i][3*(col/3)+j]] > 1) {
                    	return false;
                    }
                }
            }
		}
		return true;
	}

	//check if any blank space left
	bool check_blank(int row, int col)
    {
        for (int i = 0; i <9; i++)
        {
            for (int j = 0; j < 9; j++)
            {
             	if(puzzle[i][j] == 0){
                	return true;
             	}
            }
        }
        return false;
    }


	
public:
	// Public member function that reads the incomplete puzzle
	// we are not doing any checks on the input puzzle -- that is,
	// we are assuming they are indeed valid
	void read_puzzle(int argc, char * const argv[])
	{
		// write code that reads the input puzzle using the 
		// guidelines of figure 23 of the bootcamp material

		ifstream input_file(argv[1]);
		std::cout << std::endl << "Input File Name:" << argv[1]<< std::endl;

		//if the input file exists
		if (input_file.is_open()){ 
			//read input file and write it into our puzzle
			for (int i = 0; i < 9; i++){
				for (int j = 0; j < 9; j++){
						input_file >> puzzle[i][j];
				}
			}
		} else{
			//if the input file does not exist
			cout << "Error: Input file does not exist in PWD" << endl;
			exit(0);
		}

	}

	
	// Public member function that prints the puzzle when called
	void print_puzzle()
	{



		if (count >= 1){
            std::cout << std::endl << "Solution #" << count;
            std::cout << std::endl << "Board Position" << std::endl;
        }
		
		for (int i = 0; i < 9; i++)
		{
			for (int j = 0; j < 9; j++)
			{
				// check if we have a legitimate integer between 1 and 9
				if ((puzzle[i][j] >= 1) && (puzzle[i][j] <= 9))
				{
					// printing initial value of the puzzle with some formatting
					std::cout << puzzle[i][j] << " ";
				}
				else {
					// printing initial value of the puzzle with some formatting
					std::cout << "X ";
				}
			}
			std::cout << std::endl;
		}

		count++;
		
	}
	
	// Public member function that (recursively) implements the brute-force 
	// search for possible solutions to the incomplete Sudoku puzzle
	bool Solve(int row, int col)
	{
		// this part of the code identifies the row and col number of the 
		// first incomplete (i.e. 0) entry in the puzzle.  If the puzzle has
		// no zeros, the variable row will be 9 => the puzzle is done, as 
		// each entry is row-, col- and block-valid...
		
		// use the pseudo code of figure 3 of the description

		if (!check_blank(row, col)){
            print_puzzle();
            return false;
        } else{

			for (int i = row; i < 9; i++){
	        	for (int j = 0; j < 9; j++){
	        		//find a blank space in the puzzle to be solved
	                if (puzzle[i][j] == 0){
	                	//iterate through numbers from 1 to 9
	                	for (int k = 1; k < 10; k++){
	                		//fill one nuber in
	        				puzzle[i][j] = k;
	        				//if all theree row, column, black-assignment checks are satisfied, return true
	        				if(row_valid(i) && col_valid(j) && block_valid(i, j) && Solve(i, j)){
	        					return true;
	        				} 
	        				puzzle[i][j] = 0;
	        			}
	        			return false;
	        			std::cout << std::endl << "The Sukoku has no solution." << std::endl;
					}
	            }
	        }
	        //sudoku solved
	        return true;
    	}
	}

};

#endif
