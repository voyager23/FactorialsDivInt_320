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
#include <cmath>
#include <tuple>
#include <algorithm>
#include "../inc/toolbox.hxx"

using namespace std;

// Uses bisection method
ul inverse_legendre_factorial(ul p, ul e);	//Finds smallest number, n, such that p^e divides n!
ul modified_legendre_factorial(const ul n, ul p);	// find number of factors of p in n!
bool comp(tuple<ul,ul,ul> a,tuple<ul,ul,ul>);
 
ul inverse_legendre_factorial(ul p, ul e){
	//Finds smallest number, n, such that p^e divides n!
	ul low = p;
	ul high = p*e;
	ul mid = 0;
	while((high-low) > 1){
		mid = (high+low)/2;
		if (modified_legendre_factorial(mid,p) >= e)
			high = mid;
		else
			low = mid;
	}
	
	if (modified_legendre_factorial(low,p) >= e)
		return low;
	else
		return high;
}

ul modified_legendre_factorial(const ul n, ul p){
	// find number of factors of p in n!
	ul r, sum = 0;
	ul divisor = p;
	if(p <= n){
		do{
			r = (ul)(n/divisor);	// integer arithmetic
			divisor *= p;
			sum += r;
			//cout << n << " " << divisor << " " << r << endl;
		}while(r != 0);
	}
	return sum;
}

bool comp(tuple<ul,ul,ul> a,tuple<ul,ul,ul> b){
	//return true iff a < b
	return (get<2>(a) < get<2>(b));
	}

int main(int argc, char **argv)
{
	tuple<ul,ul,ul> temp;
	vector<tuple<ul,ul,ul>> current_nfact;
	const ul power = 1234567890;
	
	cout << inverse_legendre_factorial(2, power*7) << endl;
	temp = make_tuple(2,7,inverse_legendre_factorial(2, power*7));
	current_nfact.push_back(temp);
	
	cout << inverse_legendre_factorial(3, power*4) << endl;
	temp = make_tuple(3,4,inverse_legendre_factorial(3, power*4));
	current_nfact.push_back(temp);

	cout << inverse_legendre_factorial(5, power*1) << endl;
	temp = make_tuple(5,1,inverse_legendre_factorial(5, power*1));
	current_nfact.push_back(temp);

	cout << inverse_legendre_factorial(7, power*1) << endl;
	temp = make_tuple(7,1,inverse_legendre_factorial(7, power*1));
	current_nfact.push_back(temp);
	
	for(auto a : current_nfact)
		cout << get<0>(a) << " " << get<1>(a) << " " << get<2>(a) << endl;
	vector<tuple<ul,ul,ul>>::iterator result;
	result = max_element(current_nfact.begin(), current_nfact.end(), comp);
	cout << "maximum: " << get<2>(*result) << endl;
			

	return 0;
}

