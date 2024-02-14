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

using namespace std;
using ul = uint64_t;

typedef struct
{
	ul n,i,S;
	const ul modulus = 1e18;
	const ul power = 1234567890;
}Tdata;

void threadCallback(Tdata & x)
{
	cout << x.n << " " << x.i << " " << x.S << endl;
}


int main(int argc, char **argv)
{
	Tdata td1 = {1,92,314};
	Tdata td2 = {121,22,33};
	vector<Tdata> v;
	v.push_back(td1);
	v.push_back(td2);
	
    std::thread t1(threadCallback,std::ref(v[0]));
    std::thread t2(threadCallback,std::ref(v[1]));
    
    t1.join();
    t2.join();
	cout << "\ncomplete\n";
	return 0;
}

