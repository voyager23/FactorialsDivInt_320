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
typedef std::tuple<ul,ul,ul> pepf; // prime, exponent, power_factor

// Uses bisection method
ul inverse_legendre_factorial(ul p, ul e);	//Finds smallest number, n, such that p^e divides n!
ul modified_legendre_factorial(const ul n, ul p);	// find number of factors of p in n!
ul advance_current_factorial(vector<ul> &primes, vector<pepf> &current_nfact, ul &current_factorial);
bool comp(pepf a,pepf);
bool find_prime(pepf a, pepf b);

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

bool comp(pepf a,pepf b){
	//return true iff a < b
	return (get<2>(a) < get<2>(b));
	}
	
bool find_prime(pepf a, pepf b){
	// return true iff prime(a) == prime(b)
	return (get<0>(a) == get<0>(b));
	}
	
ul advance_current_factorial(vector<ul> &primes, vector<pepf> &current_nfact, ul &current_factorial){
	PfactOfN pfofn;	// vector of PrimePowers
	PfactOfN::iterator pp;
	current_factorial += 1;	// next factorial
	generate_descriptors(primes, current_factorial, pfofn);
	// for each descriptor:
	//		if prime not found in current_nfact:  vector prime, exponent, power_factor
	//			current_nfact.push_back(tuple(prime, exponent, 0))
	//		else
	//			add exponents
	//			set power factor = 0	{changed}
	//
	// finally return maximum power factor
	for(auto d = pfofn.begin(); d != pfofn.end(); ++d){ // vector of prime/powers
		vector<pepf>::iterator p;
		for(p = current_nfact.begin(); p != current_nfact.end(); ++p) if (get<0>(*d) == get<0>(*p)) break;
		
		if (p == current_nfact.end()){ // new prime
			// calc the new power factor using prime and power				
			current_nfact.push_back(make_tuple(d->first, d->second, inverse_legendre_factorial(d->first, d->second*1234567890)));
		} else {
			get<1>(*p) += get<1>(*d); //add exponent
			// calc new power factor using prime and adjusted exponent 
			get<2>(*p) = inverse_legendre_factorial(get<0>(*p), get<1>(*p)*1234567890);
		}
		
	}
	vector<pepf>::iterator result;
	result = max_element(current_nfact.begin(), current_nfact.end(), comp);
	return get<2>(*result); 
}

int main(int argc, char **argv)
{
	pepf temp;
	vector<pepf> current_nfact;
	ul current_factorial = 9;
	const ul power = 1234567890;
	
	// Checking a range of values
	vector<ul>primes;	// Required by Sieve
	SieveOfEratosthenes(primes, 1000001);
	
	// Setup currrent_nfact for 9! 2^7, 3^4, 5^1, 7^1
	//cout << inverse_legendre_factorial(2, power*7) << endl;
	temp = make_tuple(2,7,inverse_legendre_factorial(2, power*7));
	current_nfact.push_back(temp);
	//cout << inverse_legendre_factorial(3, power*4) << endl;
	temp = make_tuple(3,4,inverse_legendre_factorial(3, power*4));
	current_nfact.push_back(temp);
	//cout << inverse_legendre_factorial(5, power*1) << endl;
	temp = make_tuple(5,1,inverse_legendre_factorial(5, power*1));
	current_nfact.push_back(temp);
	//cout << inverse_legendre_factorial(7, power*1) << endl;
	temp = make_tuple(7,1,inverse_legendre_factorial(7, power*1));
	current_nfact.push_back(temp);
	// 9! - not included in sum
	vector<pepf>::iterator result;
	result = max_element(current_nfact.begin(), current_nfact.end(), comp);
	cout << current_factorial << "! has maximum_power_factor: " << get<2>(*result) << endl;
	
	ul S = 0;
	for(ul i = 10; i != 1001; ++i){
		advance_current_factorial(primes, current_nfact, current_factorial);
		result = max_element(current_nfact.begin(), current_nfact.end(), comp);
		cout << current_factorial << "! has maximum_power_factor: " << get<2>(*result) << endl;
		S += get<2>(*result);
	}
	
	cout << "Final S() = " << S << endl;

	
	return 0;
}

