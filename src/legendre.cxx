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
typedef vector<PrimePower> vpp;
void legendre(vector<PrimePower>&pp_in, ul n, vector<PrimePower>&pp_out);
bool num_isdivisible_den(vector<PrimePower> num, vector<PrimePower> den);

// debug functions
void prt_vpp(vpp z);

// -----------------------------------------------------------------------

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
	// test if the numerator is divisble by the denominator
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
void prt_vpp(vpp z){
	for(auto pp : z) cout << pp.first << "^" << pp.second << "  ";
	cout << endl;
}

//========================================================================

int main(int argc, char **argv)
{
	// Checking a range of values
	vector<ul>primes;	// Required by Sieve
	SieveOfEratosthenes(primes, 1001);
	vpp pp_in, pp_out, denominator;
	
	// inputs are i (62) and power
	ul i = 62;
	ul power = 60;
	
	generate_descriptors(primes, i, pp_in);
	legendre(pp_in, i, pp_out);
	
	cout << "Prime factors of " << i << endl;
	prt_vpp(pp_in);
	cout << i << "! ^ " << power << endl;
	prt_vpp(pp_out);
	// save denominator and raise exponents by power
	denominator = pp_out;
	
	for(vpp::iterator pp = denominator.begin(); pp != denominator.end(); ++pp) pp->second = pp->second*power;
	
	ul hi_prime = pp_out.back().first;
	cout << "hi prime " << hi_prime << endl;
	
	// initial value for n will be i * power = 3720 - hi_prime
	// search for value < n | pp_out
	
	for(ul n = ((i*power) - hi_prime); n >= 3500; n -= hi_prime){
		generate_descriptors(primes, n, pp_in);
		legendre(pp_in, n, pp_out);	
		cout << "Prime factors of " << n << endl;
		prt_vpp(pp_out);
		cout << "denominator" << endl;
		prt_vpp(denominator);		
		cout << endl;
	}
		
	return 0;	// as required by standard
}









