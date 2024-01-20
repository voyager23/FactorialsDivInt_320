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


int main(int argc, char **argv)
{
	const ul modulus = 1e18;
	vector<ul>primes;
	SieveOfEratosthenes(primes, 1001);
	PfactOfN pf_n;
	generate_descriptors(primes, 138, pf_n);
	for(auto pp : pf_n) cout << pp.first << "^" << pp.second << "  ";
	cout << endl;
	generate_descriptors(primes, 171, pf_n);
	for(auto pp : pf_n) cout << pp.first << "^" << pp.second << "  ";
	cout << endl;
	return 0;
}

