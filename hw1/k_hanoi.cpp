// K Peg Tower of Hanoi Solution
// Written by Prof. Sreenivas for IE523: Financial Computing

#include <iostream>
#include "k_hanoi.h"

int main (int argc, char * const argv[]) 
{
	int no_of_disks, no_of_pegs;
	
	// See 'k_hanoi.h' for the Towers class definition
    Towers x;
	
	std::cout << "Number of pegs? ";
	std::cin >> no_of_pegs;
	std::cout << "Number of disks? ";
	std::cin >> no_of_disks;

	
	x.solve(no_of_pegs, no_of_disks);

}
