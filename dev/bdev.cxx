/*
 * bdev.cxx
 * 
 * Copyright 2024 Mike <mike@Aorus39>
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

void legendre(vector<PrimePower>&pp_in, ul n, vector<PrimePower>&pp_out);
bool num_isdivisible_den(vector<PrimePower> num, vector<PrimePower> den);

void legendre(vector<PrimePower>&pp_in, ul n, vector<PrimePower>&pp_out){
	ul sum, divisor, r;
	pp_out.clear();
	for(auto pp : pp_in){
		sum = 0;
		divisor = pp.first;
		do {
			r = floor(n/divisor);	//floor imported via toolbox.hxx
			divisor *= pp.first;
			sum += r;			
		} while (r != 0);
		pp_out.push_back(make_pair(pp.first, sum));
	}
}

bool is_divisible(vector<PrimePower> num, vector<PrimePower> den){
	vector<PrimePower>::iterator inum, iden;
	iden = den.begin();
	while(iden != den.end()){ // Outer loop
		inum = num.begin();
		while((inum != num.end())and(inum->first != iden->first)) ++inum;
		if (inum == num.end()) return false;
		if (inum->second < iden->second) return false;
		// check next prime in denominator
		++iden;
	}
	return true;
}

int main(int argc, char **argv)
{
	const ul modulus = 1e18;
	vector<ul>primes;	// Required by Sieve
	SieveOfEratosthenes(primes, 1001);
	
	vector<PrimePower> pf_n, pp_out;
	vector<PrimePower> num, den;

	ul n=138;
	ul power=1234567890;
	generate_descriptors(primes, n, pf_n);	
	legendre(pf_n, n, pp_out);
	den = pp_out;
	// Raise demoninator powers
	for(auto pp = den.begin(); pp != den.end(); ++pp){
		pp->second = pp->second * power;
	}
	// Find a value for where n! is divisible by den.
	for(n = (138*power)-1; n != (138*power)+2; ++n){
		//cout << n << " ";
		generate_descriptors(primes, n, pf_n);
		legendre(pf_n, n, num);
		if(is_divisible(num, den)) {
			cout << endl << n << "! is divisible." << endl;
			break;
		}
	}		
	
	return 0;
}

