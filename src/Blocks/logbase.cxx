/*
 * logbase.cxx
 * 
 * Copyright 2024 Mike <mike@fedora38-2.home>
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
 #include <cmath>
#include "../../ToolBox/toolbox.hxx"
#include <unordered_map>

using namespace std;
using ul=uint64_t;
 
 double loga_baseb(double a, double b){
	 return log2(a)/log2(b);
}

#define MYLOG(a,b) log2(a)/log2(b)
/* 
 * examine the issues around creating a std::map<ul,ul> i_fact(999000)
 * consider 11! => 2^8 * 3^4 * 5^2 * 7^1 * 11^1
 * primes > floor(11/2) have exponent 1
 * use legrendre for primes <= floor(i/2)
 * for 999000! use primes < 499500
 */
 
 // Globals
const ul power = 1234567890;
const ul modulus = 1e18;

 
// Uses bisection method
ul inverse_legendre_factorial(ul p, ul e);	//Finds smallest number, n, such that p^e divides n!
ul modified_legendre_factorial(const ul n, ul p);	// find number of factors of p in n!

//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------ 
 void legendre(vector<ul> &primes, ul n, vector<pair<ul,ul>>&pp_out){
	// Returns a vector of PrimePowers, each of which divides n!
	// Corrected version 27/01/24
	ul sum, divisor, r;
	pp_out.clear();
	for(auto pp = primes.begin(); *pp <= n; ++pp){
		sum = 0;
		divisor = *pp;
		if (divisor > n) break;
		do {
			r = floor(n/divisor);	//floor imported via toolbox.hxx
			divisor *= (*pp);
			sum += r;			
		} while (r != 0);
		pp_out.push_back(make_pair(*pp, sum));
	}
}

int main(int argc, char **argv)
{
	// Basic list of primes
	vector<ul>primes;	// Required by Sieve
	SieveOfEratosthenes(primes, 1000001);
	
	cout << primes.back() << endl;
	
	vector<pair<ul,ul>> descriptors;	// pair< prime, exponent >
	vector<pair<ul,ul>> i_fact;
	
	legendre(primes, 999000, i_fact);
	
	// Using i_fact determine the minimum value for n! such that it is divisible by i_fact^1234567890
	ul n_min = 0;	// maintains the current minimum value of n such that (i!^power) | n!
	ul local_n_min;
	
	for(auto p = i_fact.begin(); p != i_fact.end(); ++p){
		local_n_min = inverse_legendre_factorial(p->first, p->second * 1234567890);
		if ( local_n_min > n_min) n_min = local_n_min;
	}
	
	cout << n_min << endl;
	
	return 0;
}

