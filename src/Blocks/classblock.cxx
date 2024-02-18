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

typedef vector<pair<ul,ul>> pair_vect;

class Legendre
{
	public:
		Legendre() {};
		void run(Tdata &x) {
				
				i_fact.clear();
				i_fact = legendre(x.primes, x.start);
				
				if(x.idx > 0) std::this_thread::sleep_for(x.idx * 1000ms);	//debug
				cout << "class legendre " << x.idx << " " << x.start << " " << x.sentinel << " " << i_fact.size() << endl;				
		}
		
	private:
		pair_vect i_fact;			
		ul local_n_min, n_min;
		
		pair_vect gen_descript(vector<ul> &primes, ul n){
			pair_vect ppv;
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
			return ppv;
		}
			
		pair_vect legendre(vector<ul> &primes, ul n){
			// Returns a vector of PrimePowers, each of which divides n!
			// Corrected version 27/01/24
			pair_vect pp_out;
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
			return pp_out;
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
    const ul n_threads = 4;
    
    const ul elements = u_hi - u_lo + 1;
    const ul width = elements / n_threads;
    ul residue = elements % n_threads;
    ul begin = u_lo;
    ul end; // Sentinel one-past-end value
    
    // Build the vector of thread data
    // td expects {idx, start, sentinel}
    
    // Thread data
    Tdata td1 = {1,600000,301};
    Tdata td2 = {2,700000,701};
    Tdata td3 = {3,800000,8001};
    Tdata td4 = {4,900000,9001};

	
    // classes
    Legendre lg1 = Legendre();
    Legendre lg2 = Legendre();
    Legendre lg3 = Legendre();
    Legendre lg4 = Legendre();
 
    // threads

	auto w = std::thread(&Legendre::run, &lg1, std::ref(td1));
	auto x = std::thread(&Legendre::run, &lg2, std::ref(td2));
	auto y = std::thread(&Legendre::run, &lg3, std::ref(td3));
	auto z = std::thread(&Legendre::run, &lg4, std::ref(td4));
	
	w.join();
	x.join();
	y.join();
	z.join();
    
	cout << "Complete\n";
	return 0;
}

