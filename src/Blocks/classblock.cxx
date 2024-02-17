/*
 * classblock.cxx
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
#include <cstdint>
#include <thread>
#include <array>
#include <vector>
#include <cmath>
#include "../../ToolBox/toolbox.hxx"

using namespace std;
using ul = uint64_t;

// Basic list of primes

const ul modulus = 1e18;
const ul power = 1234567890;
vector<ul> vprime;	// Required by Sieve
		
typedef struct
{
	unsigned idx;
	ul start, sentinel;
	vector<ul> &primes = vprime;
}Tdata;

class Legendre
{
	public:
	
		Legendre() {};
		void run(Tdata &x) {
				if(x.idx > 0) std::this_thread::sleep_for(x.idx * 1000ms);	//debug
				
				cout << "class legendre " << x.start << " " << x.sentinel << " " << endl;
				// initialise to 100! - Assume that smaller primes will be overwritten
				generate_descriptors(x.primes, x.start, i_fact);
				//DEBUG
				for(auto pp : i_fact) cout << pp.first << "^" << pp.second << endl;
								
				//~ for(auto d = descriptors.begin(); d != descriptors.end(); ++d){
					//~ i_fact.push_back(make_pair(d->first, d->second));
					//~ local_n_min = inverse_legendre_factorial(d->first, d->second*1234567890);
					//~ if ( local_n_min > n_min) n_min = local_n_min;
				//~ } // Next descriptor prime;
		}
		
	private:
	
		vector<pair<ul,ul>> descriptors;	// pair< prime, exponent >
		vector<pair<ul,ul>> i_fact;			
		ul local_n_min, n_min;
			
		void legendre(vector<ul> &primes, ul n, vector<pair<ul,ul>> pp_out){
			// Returns a vector of PrimePowers, each of which divides n!
			// Corrected version 27/01/24
			ul sum, divisor, r;
			pp_out.clear();
			for(auto pp = primes.begin(); *pp <= n; ++pp){
				sum = 0;
				divisor = *pp;
				if (divisor > n) break;
				do {
					r = floor(n/divisor);	//floor imported via cmath
					divisor *= (*pp);
					sum += r;			
				} while (r != 0);
				pp_out.push_back(make_pair(*pp, sum));
			}
		}

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
};


int main(int argc, char **argv)
{
	SieveOfEratosthenes(vprime, 1000001);

    const ul u_lo = 10;
    const ul u_hi = 1000;
    const ul n_threads = 10;
    
    const ul elements = u_hi - u_lo + 1;
    const ul width = elements / n_threads;
    ul residue = elements % n_threads;
    ul begin = u_lo;
    ul end; // Sentinel one-past-end value
    
    // Build the vector of thread data
    // td expects {idx, start, sentinel}
    // vector of Thread data
    
    vector<Tdata> td;
    td.reserve(n_threads);
    for(unsigned t = 0; t != n_threads; t++){
		end = begin + width;
		if(residue > 0) {end += 1; residue -= 1;}
		cout << begin << "\t -> \t" << (end -1) << endl;
		td.push_back(Tdata{t, begin, end});
		begin = end;
	}
	
    // vector of class
    vector<Legendre> vlg;
 
    // vector of threads
    vector<thread> vt;
    for(unsigned t = 0; t != n_threads; t++){
		vlg.push_back(Legendre());
		vt.emplace_back(std::thread(&Legendre::run, &(vlg.back()), std::ref(td[t])));
	}
	
	for(auto x = vt.begin(); x != vt.end(); ++x) (*x).join();
    
	cout << "Complete\n";
	return 0;
}

