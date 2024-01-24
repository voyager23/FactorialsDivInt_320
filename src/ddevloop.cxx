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
// debug function
void prt_vpp(vpp z);

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

int main(int argc, char **argv)
{
	// Checking a range of values
	vector<ul>primes;	// Required by Sieve
	SieveOfEratosthenes(primes, 1001);

	ul d = 17;
	while(d < 100){
		ul power = 60;
		vpp pf_n, numerator, denom, denominator;
		//vpp denom = {make_pair(2,16),make_pair(3,8),make_pair(5,3),make_pair(7,2),make_pair(11,1),make_pair(13,1),make_pair(17,1)}; // 18!
		//for(auto pp : denom) pp.second *= power;	// 18! ^ power
		generate_descriptors(primes, d, denom);	// get the prime powers
		legendre(denom, d, denominator);	// get prime powers of factorial
		for(vpp::iterator pp = denominator.begin(); pp != denominator.end(); ++pp) pp->second = pp->second*power;	// raise factorial to power
		
		ul n = d+1;
			n = 34;
		do{
			//cout << "checking n = " << n << "  ";
			generate_descriptors(primes, n, pf_n);  // generate the prime factors of n
			legendre(pf_n, n, numerator);			// generate the prime factors of n!
			if(is_divisible(numerator,denominator)) {		// check if denominator | numerator
				cout << endl << n << "! is divisible by " << d <<"! ^ " << power << endl;
				prt_vpp(numerator);
				prt_vpp(denominator);
				if(n != d*power){
					cout <<"FAIL expected " << d*power << "!  delta:" << (d*power - n)/d << endl << endl;
				}
				break;
			} else {
				//cout << endl;
				n += 1;
			}
		} while(1);
		d += 1;
	}
	return 0;	// as required by standard
}









