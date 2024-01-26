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
bool isdivisible(vector<PrimePower> num, vector<PrimePower> den);
int test_num_den(vpp num, vpp den);

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
	
	// inputs are i (62) and power
	ul i = 231;
	ul power = 60;
	
	generate_descriptors(primes, i, pp_in);
	legendre(pp_in, i, pp_out);
	
	const_denom = pp_out;
	
	for(vpp::iterator pp = const_denom.begin(); pp != const_denom.end(); ++pp) pp->second = pp->second*power;
	
	ul hi_prime = pp_out.back().first;
	cout << "hi prime " << hi_prime << endl << endl;
	
	ul n = (i*power);
	ul minimum_soln = n;
	do{
		generate_descriptors(primes,n,pp_in);
		legendre(pp_in,n,pp_out);
		
		cout << "Prime factors of " << n << endl;
		prt_vpp(pp_out);
		cout << "const_denom" << endl;
		prt_vpp(const_denom);		
		cout << endl;
		vpp::iterator hip;
		for(hip = pp_out.begin(); hip != pp_out.end(); ++hip) if (hip->first == hi_prime) break;
		
		if(hip != pp_out.end()){ // found hi_prime
			if(is_divisible(pp_out, const_denom)and(hip->second >= const_denom.back().second)){				
				minimum_soln = n;
			} else {
				if(hip->second < const_denom.back().second){
					cout << "Search limit. Min n = " << minimum_soln << endl;
					break;
				}
			}
		}
		n -= hi_prime;
	}while(1);

	return 0;	// as required by standard
}

