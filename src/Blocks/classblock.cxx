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
#include <cstdint>
#include <unordered_map>
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
	const ul modulus = 1e18;
	ul Sum; // Computed thread sum mod 1e18
	vector<ul> &primes = vprime;
}Tdata;

typedef vector<pair<ul,ul>> pair_vect;

//---------------------------------------------------------------------------------------------
class Legendre
{
	public:
		Legendre() {};
		void run(Tdata &x) {
			/* outline
			 * 
			 * On Entry -  populate a map of prime/exponent pairs for factorial x.start.
			 * Scan complete map to find minimum n for which (i!)^1234567890 | n!. This is n_min
			 * 
			 * for(i = x.start+1; i != sentinel; ++i) {
			 * 		local_n_min = 0
			 * 		find prime/exponent pairs for i
			 * 		for each pair {
			 * 			update map
			 * 			update local_n_min inverse_legendre_factorial
			 * 		}
			 * 		update n_min using local_n_min
			 * 		update Sum using n_min
			 * 	} // next value of i
			 * 	finally - record (Sum % 1e18) in Tdata block
			 * 	Stop
			 */
			 
			// output formatting pause
			// if(x.idx > 0) std::this_thread::sleep_for(x.idx * 100ms);
				
				
			// Entry	
			i_fact.clear();
			legendre(std::ref(x.primes), x.start, std::ref(i_fact));

			n_min = 0;
			for(auto d = i_fact.begin(); d != i_fact.end(); ++d){
				local_n_min = inverse_legendre_factorial(d->first, d->second*1234567890);
				if ( local_n_min > n_min) n_min = local_n_min;
			} // Next descriptor prime;
			// n_min is now correct for the first value of n - set value of Sum
			x.Sum = n_min % x.modulus;
			
			// evaluate the remaining values of i in this block
			for(ul i = x.start+1; i != x.sentinel; ++i){
				// prime/exponent pairs for i
				gen_descript(std::ref(x.primes), i, descriptor);
				// for each pair - update the map and adjust the value of local_n_min 
				for(auto d = descriptor.begin(); d != descriptor.end(); ++d){
					auto p = i_fact.find(d->first);
					if(p != i_fact.end()){ // found existing prime update exponent
						p->second += d->second;
					} else { // new prime
						auto rc = i_fact.emplace(make_pair(d->first, d->second));
						// update p to point to new entry
						if(rc.second){
							p = rc.first;
						} else {
							cout << "Emplace new value failed" << endl;
						}
					}
					// find local_n_min for this prime/exponent pair
					local_n_min = inverse_legendre_factorial(p->first, p->second * 1234567890);
					// Adjust the n_min value as required
					if ( local_n_min > n_min) n_min = local_n_min; 
				} // next descriptor prime
				x.Sum = n_min % x.modulus;							
			} // for i...
		} // run
		
	private:
		unordered_map<ul,ul> i_fact;			
		pair_vect descriptor;			
		ul local_n_min, n_min;
		
		void gen_descript(vector<ul> &primes, ul n, pair_vect & ppv){
			// pair_vect ppv;
			ppv.clear();
			// Generates a vector of descriptors for integer n. Each element is a pair<prime, power>
			std::pair<ul,ul> temp;
			for(auto i = primes.begin(); i != primes.end(); ++i){
				ul p = *i;
				if(p > n) break;
				temp = {p,0};
				while((n % p)==0){
					temp.second += 1;
					n /= p;
				}
				// save <prime,power> to vector
				if(temp.second > 0) ppv.push_back(temp);
			}
			//return ppv;
		}
			
		void legendre(vector<ul> &primes, ul n, unordered_map<ul,ul> &pp_out){
			// Returns a vector of PrimePowers, each of which divides n!
			// Corrected version 27/01/24
			// pair_vect pp_out;
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
				auto rc = pp_out.emplace(*pp, sum);
				if(rc.second == false){
					cout << "Emplace failed. Prime:" << *pp << " sum:" << sum << endl;
					exit(1);
				}
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
//------------------------------------------------------------------------------------------

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
	vector<Tdata> td;
	td.reserve(n_threads);
	for(unsigned t = 0; t != n_threads; t++){
		end = begin + width;
		if(residue > 0) {end += 1; residue -= 1;}
		//cout << begin << "\t -> \t" << (end -1) << endl;
		td.push_back(Tdata{t, begin, end});
		begin = end;
	}
	
    // vector of classes
    std::vector<Legendre> lg;
    for(ul n = 0; n != n_threads; ++n) {
		lg.push_back(Legendre());
	}
	
    // vector of threads
    std::vector<std::thread> workers;
    for(ul n = 0; n != n_threads; ++n){
		auto v = std::thread(&Legendre::run, std::ref(lg[n]), std::ref(td[n]));
		workers.push_back(move(v));
	}
	
	ul S = 0;
    for(ul n = 0; n != n_threads; ++n) {
		workers[n].join();
		// output index and sum for each thread
		cout << "idx:" << td[n].idx << "  thread_Sum:" << td[n].Sum << endl;
		S += td[n].Sum;
	}
    
	cout << "Complete Sum:" << S << endl;
	return 0;
}
	
