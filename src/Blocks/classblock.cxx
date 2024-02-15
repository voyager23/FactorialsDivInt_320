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
vector<ul> vprime;	// Required by Sieve
const ul modulus = 1e18;
const ul power = 1234567890;
		
typedef struct
{
	ul n,i,S;
	vector<ul> &primes = vprime;
}Tdata;

class Legendre
{
	public:
	
		Legendre() {};
		void run(Tdata &x) {
				cout << "class legendre " << x.n << " " << x.i << " " << endl;
				// initialise to 100! - Assume that smaller primes will be overwritten
				generate_descriptors(x.primes, x.i, descriptors);
				for(auto d = descriptors.begin(); d != descriptors.end(); ++d){
					i_fact.push_back(make_pair(d->first, d->second));
					local_n_min = inverse_legendre_factorial(d->first, d->second*1234567890);
					if ( local_n_min > n_min) n_min = local_n_min;
				} // Next descriptor prime;
				
				//DEBUG
				for(auto pp : i_fact) cout << pp.first << "^" << pp.second << endl;	
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
	
    Legendre ct;
    Tdata td1;
    Tdata td2;
    Tdata td3 = {1, 100, 0};
    std::thread t3(&Legendre::run, &ct, std::ref(td3));
    
    t3.join();
	cout << "\ncomplete\n";
	return 0;
}

