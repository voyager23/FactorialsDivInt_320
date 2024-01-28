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
void legendre(vector<ul> &primes, ul n, vector<PrimePower>&pp_out);
bool isdivisible(vector<PrimePower> num, vector<PrimePower> den);
int test_num_den(vpp num, vpp den);

// debug functions
void prt_vpp(vpp z);

// ----------------------------------------------------------------------

void legendre(vector<ul> &primes, ul n, vector<PrimePower>&pp_out){
	// Corrected version 27/01/24
	ul sum, divisor, r;
	pp_out.clear();
	for(auto pp : primes){
		if (pp > n) break;
		sum = 0;
		divisor = pp;
		do {
			r = floor(n/divisor);	//floor imported via toolbox.hxx
			divisor *= pp;
			sum += r;			
		} while (r != 0);
		pp_out.push_back(make_pair(pp, sum));
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

int test_num_den(vpp num, vpp den, ul hi_prime){

	// if hi_prime_exp_num < hi_prime_exp_den return +1 hard fail
	vpp::iterator hip;
	for(hip = num.begin(); hip != num.end(); ++hip) if (hip->first == hi_prime) break;
	if((hip != num.end())and(hip->second < den.back().second )) return +1; // hard fail
	if(is_divisible(num,den)) return 0; //update solution
	return -1; // soft fail - continue to search
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
	vpp pp_in, pp_out, const_denom;
	
	// inputs are i and power
	ul i = 16;
	ul power = 2;
	legendre(primes, i, const_denom);
	//multiply exponents by power -> const_denom => (i!)^power
	for(auto pp = const_denom.begin(); pp != const_denom.end(); ++pp) pp->second *= power;
	
	for(ul j = 48; j > 25; --j){
		legendre(primes, j, pp_out);
		cout << j << "!" << endl;
		prt_vpp(pp_out);
		cout << "------denominator------"<<endl;
		prt_vpp(const_denom);
		cout << endl;
	}
	

	return 0;	// as required by standard
}

