// k_hanoi.h
// Written by Prof. Sreenivas for IE523
// Modified from the work of Prof. Sreenivas in hanoi.h which is uploaded on Compass

#ifndef K_HANOI_H
#define K_HANOI_H

#include <cmath>
#include <deque>

using std::deque;

class Towers 
{
	// private data member: #steps 
	int number_of_steps;
	
	// private data member: #pegs
	int number_of_pegs;
	
	// private data member: #disks on 1st Peg
	int number_of_disks;
	
	// private data member: k-Deques-of-Deques
	// towers[0] is deque for the 1st Peg, towers[1] is the deque 
	// for the 2nd Peg,..., towers[k-1] is the deque for the kth Peg
	deque < deque <int> > towers; 
	
	// private member function:  returns 'false' if a larger disk
	// is moved on top of a smaller disk during a move
	bool is_everything_legal() 
	{ 
		// write the code for this part
		bool result = true;
		for (int i = 0; i < number_of_pegs; i++)
			for (int j = 0; j < towers[i].size(); j++)
				for (int k = j; k < towers[i].size(); k++)
					 //if the larger disk is put on the top of a smaller one, the step is illegal
					if (towers[i][k] < towers[i][j]) {
						result = false;
					}
		return (result);
	}
	
	// private member function: initializes the pegs and #disks on Peg 1
	void initialize(int m, int n)
	{
		// write the code for this part

		// initialize the number of pegs, steps, ans disks
		number_of_pegs = m; 
		number_of_steps = 0;
		number_of_disks = n;

		//m pegs
		for (int i = 0; i < m; i++)
		{
			deque <int> peg;
			towers.push_back(peg);
		}

		//n disks on peg 1 at first
		for (int i = 0; i < n; i++) 
			towers[0].push_back(i+1);
	}
	
	// private member function:  picks the #disks to be moved based on
	// the #free pegs that are available at a given instant
	int pick_the_right_number_to_move(int number_of_disks, int number_of_free_pegs)
	{
		if (1 == number_of_free_pegs)
			return (number_of_disks-1);
		else 
			return (number_of_disks/2);
	}
	
	// private member function: primitive/basic move of the top disk
	// from the source peg to the destination peg -- and legality-check
	void move_top_disk(int source, int dest)
	{
		// write the code for this part
		// modified from "hanoi.h"
		std::cout << "Move disk " << towers[source].front() << " from Peg ";
		std::cout << source+1 << " to Peg " << dest+1;
		towers[dest].push_front(towers[source].front());
		towers[source].pop_front();
		if (true == is_everything_legal())
			std::cout << " (Legal)" << std::endl;
		else 
			std::cout << " (Illegal)" << std::endl;
	}
	
	// private member function: recursive solution to the 3 Peg Tower of Hanoi
	void Move(int number_of_disks, int source_peg_number, int destination_peg_number, deque <int> free_peg_numbers)
	{
		if (free_peg_numbers.size() != 0)
		{
			if (number_of_disks > 1) 
			{
				int m = pick_the_right_number_to_move(number_of_disks, free_peg_numbers.size());
				// write the code for this part
				
				//keep track of the intermediate peg, the last in the free peg deque but the second last among the total pegs
				int intermediate_peg_number = free_peg_numbers.back();

				//intermediate peg gets popped from the deque
				free_peg_numbers.pop_back();

				//push the destination peg onto the free_pegs deque
				free_peg_numbers.push_back(destination_peg_number);

				//move m disks from source_peg to intermediate_peg
				Move(m, source_peg_number, intermediate_peg_number, free_peg_numbers); 
			
				//the destination peg gets popped from the free_pegs deque
				free_peg_numbers.pop_back();

				//move (number of disks - m) disks from source peg to the destination peg
				Move(number_of_disks - m, source_peg_number, destination_peg_number, free_peg_numbers);
				
				//the sourse_peg is now free so push that back to the free_pegs deque
				free_peg_numbers.push_back(source_peg_number);

				//move the rest m disks from the intermediate_peg to the destination_peg
				Move(m, intermediate_peg_number, destination_peg_number, free_peg_numbers);
			}
			else {
				// number_of_disks == 1
				number_of_steps++;
				move_top_disk(source_peg_number, destination_peg_number);
			}
		}
		else {
			// there are no free pegs
			number_of_steps++;
			move_top_disk(source_peg_number, destination_peg_number);
		}
	}

	
	// private member function: prints the state of Peg m (m = 0, 1, 2)
	void print_peg_state(int m)
	{
		// write the code for this part
		// taken from "hanoi.h"
		std::cout << "-----------------------------" << std::endl;
		std::cout << "State of Peg " << m+1 << " (Top to Bottom): ";
		for (int i = 0; i < towers[m].size(); i++) 
			std::cout << towers[m][i] << " ";
		std::cout << std::endl;
		std::cout << "Number of Steps = " << number_of_steps << std::endl;
		std::cout << "-----------------------------" << std::endl;
	}
public:
	// Solves the K Peg Tower of Hanoi
	void solve(int number_of_pegs, int number_of_disks)
	{
		initialize(number_of_pegs, number_of_disks);
		print_peg_state(0);
		
		deque <int> free_peg_numbers;
		// source peg # is 0; destination peg # is no_of_pegs-1; 
		// free pegs {1,.., no_of_pegs-2}
		for (int i = 1; i <= number_of_pegs-2; i++)
			free_peg_numbers.push_back(i);
		Move(number_of_disks, 0, number_of_pegs-1, free_peg_numbers);
		print_peg_state(number_of_pegs-1);
	}
};
#endif