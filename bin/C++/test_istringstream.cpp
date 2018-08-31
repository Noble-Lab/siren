#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main(){

	string s("0 1 \n");
	istringstream iss( s );
	int moo;
	char foo;

	iss >> moo;
	printf("%i -%i-\n",0,moo);
	iss >> moo;
	printf("%i -%i-\n",1,moo);
	iss >> foo;
	printf("%i -%c-\n",2,foo);


	return 0;
}

