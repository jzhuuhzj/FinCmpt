/*
 *  card_game.h
 *  Card Game
 *
 *  Created by Ramavarapu Sreenivas  
*/

#ifndef	CARDGAME_H
#define CARDGAME_H
#include <algorithm>
#include <iterator>
#include <map>

double value(int r, int b){	
	//key:number of red cards and blue cards; value: the value to you
	static std::map<std::pair<int, int>, double> cardMap;

	//if no red card left
	if (0 == r)
        return ((double) b);

    //if no blue card left
    if (0 == b)
        return (0);

    //if the the key is already in the map, just output the corresponding value
    if (cardMap.find(std::make_pair(r, b)) != cardMap.end()){
    	return (double)cardMap[std::make_pair(r, b)];
    
    //otherwise, calculate the value for the new pair
    } else{
    	double term1 = ((double) r/(r+b)) * value(r-1, b);
	
		double term2 = ((double) b/(r+b)) * value(r, b-1);
	
		double val = std::max((term1 + term2), (double) (b - r));

	    cardMap[std::make_pair(r, b)] = val;

		return val;
    }
}

#endif