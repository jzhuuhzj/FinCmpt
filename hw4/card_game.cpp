#include <iostream>
#include "card_game.h"
#include <ctime>
#include <cmath>

using namespace std;
	
int main (int argc, char * const argv[]) 
{
	int total_number_of_cards;	
	sscanf (argv[1], "%d", &total_number_of_cards); 

	cout << "Total Number of Cards = " << total_number_of_cards << endl;
	cout << "Value of the game = " << value(total_number_of_cards/2,total_number_of_cards/2) << endl;

    return 0;
}