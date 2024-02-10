#include <iostream>
#include <unordered_map>
#include <array>
#include <vector>
#include <cstdint>
#include "../inc/toolbox.hxx"

using namespace std;
using ul=uint64_t;
typedef	std::unordered_map<ul,ul> MapPrimeExponent;	// simple map of key prime, value basic exponent
 
// Globals
const ul power = 1234567890;
ul n_min;	// maintains the current minimum value of n such that (i!^power) | n!
 
// Uses bisection method
ul inverse_legendre_factorial(ul p, ul e);	//Finds smallest number, n, such that p^e divides n!
ul modified_legendre_factorial(const ul n, ul p);	// find number of factors of p in n!
ul advance_current_factorial(vector<ul> &primes, MapPrimeExponent &current_nfact, ul &current_factorial);

//------------------------------------------------------------------------------
ul advance_current_factorial(vector<ul> &primes, MapPrimeExponent &current_nfact, ul &current_factorial){
	PfactOfN pfofn;			// vector of PrimePowers
	PfactOfN::iterator pp;	// iterator
	//vector<pepf>::iterator p;
	
	current_factorial += 1;	// next factorial
	generate_descriptors(primes, current_factorial, pfofn);
	
	for(auto d = pfofn.begin(); d != pfofn.end(); ++d){ // vector of prime/powers

		for(auto pp = current_nfact.begin(); pp != current_nfact.end(); ++pp) if (get<0>(*d) == get<0>(*pp)) break;
		
		if (pp == current_nfact.end()){ // new prime
			// calc the new power factor using prime and power				
			current_nfact.push_back(make_tuple(d->first, d->second, inverse_legendre_factorial(d->first, d->second*power)));
		} else {
			get<1>(*pp) += get<1>(*d); //add exponent
			// calc new power factor using prime and adjusted exponent 
			get<2>(*pp) = inverse_legendre_factorial(get<0>(*p), get<1>(*pp)*power);
		}
		
	}
	vector<pepf>::iterator result;
	result = max_element(current_nfact.begin(), current_nfact.end(), comp);
	return get<2>(*result); 
}

//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
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


//------------------------------------------------------------------------------
int main(int argc, char **argv){

	//~ std::vector<std::pair<ul,ul>> tenfact { make_pair(2,1), make_pair(5,1) };
	//~ MapPrimeArray mymap;
	
	//~ mymap.emplace(make_pair(2, array<ul,2> {7, 0}));
	//~ mymap.emplace(make_pair(3, array<ul,2> {4, 0}));
	//~ mymap.emplace(make_pair(5, array<ul,2> {1, 0}));
	//~ mymap.emplace(make_pair(7, array<ul,2> {1, 0}));
	
	//~ for(auto p : tenfact){
		//~ auto rc = mymap.find(p.first);
		//~ if(rc != mymap.end()){
			//~ cout << "found" << endl;
			//~ // add exponent to value in array
			//~ rc->second[0] += p.second;
			
			//~ // adjust min_n_fact in array
			//~ rc->second[1] = rc->second[0] + p.first;			
			
		//~ } else {
			//~ cout << "not found" << endl;			
		//~ }
	//~ }
		
	return 0;
}
	
