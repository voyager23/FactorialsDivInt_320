/*
 * inverse_legendre_factorial.cxx
 * 
 * Copyright 2024 mike <mike@Fedora37>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <iostream>
#include <vector>
#include <utility>
#include "../inc/toolbox.hxx"

using namespace std;

// Uses bisection method
ul inverse_legendre_factorial(ul p, ul e);	//Finds smallest number, n, such that p^e divides n!
ul modified_legendre_factorial(ul n, ul p);	// find number of factors of p in n!

ul inverse_legendre_factorial(ul p, ul e){
	//Finds smallest number, n, such that p^e divides n!
	ul low = p;
	ul high = p*e;
	ul mid = 0;
	while((high-low) > 1){
		mid = (high+low)/2;
		if (modified_legendre_factorial(mid,p) > e)
			high = mid;
		else
			low = mid;
	}
	if (modified_legendre_factorial(low,p) >= e)
		return low;
	else
		return high;
}

ul modified_legendre_factorial(ul n, ul p){
	// find number of factors of p in n!
	ul r, sum = 0;
	if(p <= n){
		do{
			ul r = floor(n/p);
			p *= p;
			sum += r;
		}while(r != 0);
	}
	return sum;
}

int main(int argc, char **argv)
{
	
	return 0;
}

