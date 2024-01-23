/*
 * adev.cxx
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
#include <vector>
#include <utility>
#include "../inc/toolbox.hxx"
using namespace std;

void legendre(vector<PrimePower>&pp_in, ul n, vector<PrimePower>&pp_out);

void legendre(vector<PrimePower>&pp_in, ul n, vector<PrimePower>&pp_out){
	//~ Outline:
	//~ declare ul sum
	//~ for each prime in pp_in
		//~ clear sum
		//~ Sum Floor n/p^1,2,3 etc
		//~ pp.push_back(make_pair(prime,sum))
	//~ stop
	ul sum, divisor, r;
	pp_out.clear();
	for(auto pp : pp_in){
		sum = 0;
		divisor = pp.first;
		do {
			r = floor(n/divisor);	//floor imported from toolbox.hxx
			divisor *= pp.first;
			sum += r;			
		} while (r != 0);
		pp_out.push_back(make_pair(pp.first, sum));
	}
}

int main(int argc, char **argv)
{
	const ul modulus = 1e18;
	vector<ul>primes;
	SieveOfEratosthenes(primes, 1001);
	PfactOfN pf_n;
	
	//~ generate_descriptors(primes, 138, pf_n);
	//~ for(auto pp : pf_n) cout << pp.first << "^" << pp.second << "  ";
	//~ cout << endl;
	vector<PrimePower> pp_out;
	//~ legendre(pf_n, 138, pp_out);
	//~ cout << "138 factorial prime/powers" << endl;
	//~ for(auto pp : pp_out) cout << pp.first << "^" << pp.second << "  ";
	//~ cout << endl;
	
	const ul n=276;
	generate_descriptors(primes, n, pf_n);
	for(auto pp : pf_n) cout << pp.first << "^" << pp.second << "  ";
	cout << endl;
	legendre(pf_n, n, pp_out);
	cout << n << " factorial prime/powers" << endl;
	for(auto pp : pp_out) cout << pp.first << "^" << pp.second << "  ";
	cout << endl;
	
	return 0;
}

