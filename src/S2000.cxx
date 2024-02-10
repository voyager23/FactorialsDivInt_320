/*
 * S2000.cxx
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
 * S(100000)   = 6172360737804266019 	19 digits
 * 
 * S(1000000)  = 617 278157919195482643
 * 
 * S(10000000) = 61728327266010338190906
 * 
 */


#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <cstdint>
#include <array>
#include <unordered_map>
#include "../inc/toolbox.hxx"
using namespace std;

// Global constants
const ul power = 1234567890;
const ul modulus = 1000000000000000000;

// map key prime, value array<2, prime_exponent, minimum_n>
typedef unordered_map<ul, array<ul, 2>> mymap;

ul inverse_legendre_factorial(ul p, ul e);	//Finds smallest number, n, such that p^e divides n!
ul modified_legendre_factorial(const ul n, ul p);	// find number of factors of p in n!
ul advance_current_factorial(vector<ul> &primes,  mymap foobar, ul &current_factorial);

// TODO - This uses a vector of prime,exponent and n_min as a database.

ul advance_current_factorial(vector<ul> &primes, vector<pepf> &current_nfact, ul &current_factorial){
	PfactOfN pfofn;	// vector of PrimePowers
	PfactOfN::iterator pp;
	vector<pepf>::iterator p;
	
	current_factorial += 1;	// next factorial
	generate_descriptors(primes, current_factorial, pfofn);
	
	for(auto d = pfofn.begin(); d != pfofn.end(); ++d){ // vector of prime/powers

		for(p = current_nfact.begin(); p != current_nfact.end(); ++p) if (get<0>(*d) == get<0>(*p)) break;
		
		if (p == current_nfact.end()){ // new prime
			// calc the new power factor using prime and power				
			current_nfact.push_back(make_tuple(d->first, d->second, inverse_legendre_factorial(d->first, d->second*power)));
		} else {
			get<1>(*p) += get<1>(*d); //add exponent
			// calc new power factor using prime and adjusted exponent 
			get<2>(*p) = inverse_legendre_factorial(get<0>(*p), get<1>(*p)*power);
		}
		
	}
	vector<pepf>::iterator result;
	result = max_element(current_nfact.begin(), current_nfact.end(), comp);
	return get<2>(*result); 
}


int main(int argc, char **argv){
	
	// Checking a range of values
	vector<ul>primes;	// Required by Sieve
	SieveOfEratosthenes(primes, 1000001);	
	
	return 0;
}
