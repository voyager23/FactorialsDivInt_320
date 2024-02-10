
#include "../ToolBox/toolbox.hxx"
#include <unordered_map>

using namespace std;
using ul=uint64_t;
typedef	std::unordered_map<ul,ul> mapPrimeExp;	// simple map of key prime, value basic exponent
 
// Globals
const ul power = 1234567890;
const ul modulus = 1e18;

 
// Uses bisection method
ul inverse_legendre_factorial(ul p, ul e);	//Finds smallest number, n, such that p^e divides n!
ul modified_legendre_factorial(const ul n, ul p);	// find number of factors of p in n!

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

	const ul modulus = 1e18;
	const ul limit = 1000000; // loop limit
	ul n_min = 0;	// maintains the current minimum value of n such that (i!^power) | n!
	ul local_n_min;
	ul S = 0;	//Sum
	ul i = 0;	//loop variable
	
	
	mapPrimeExp i_fact;
	// initialise to 9!
	i_fact.emplace(make_pair(2, 7));
	i_fact.emplace(make_pair(3, 4));
	i_fact.emplace(make_pair(5, 1));
	i_fact.emplace(make_pair(7, 1));
	
	// Basic list of primes
	vector<ul>primes;	// Required by Sieve
	SieveOfEratosthenes(primes, 1000001);
	vector<PrimePower> descriptors;
	
	// main loop
	for(i = 10; i <= limit; ++i){
		generate_descriptors(primes, i, descriptors);
		for(auto d = descriptors.begin(); d != descriptors.end(); ++d){
			auto p = i_fact.find(d->first);
			if(p != i_fact.end()){ // found existing prime
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
			local_n_min = inverse_legendre_factorial(p->first, p->second * 1234567890);
			if ( local_n_min > n_min) n_min = local_n_min; 
		} // next descriptor prime
		S = (S + n_min) % modulus;
	} // next i value
		
	cout << "S(" << limit << ") = " << S << endl;
	cout << "S(1000000) = 278157919195482643" << endl;
	return 0;
}
	
