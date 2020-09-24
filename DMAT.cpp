//Primary cpp file for DMAt code including rk45 algorithm, correct differential form of boltzmann equation and main driver function

using namespace std;
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "DMAT.hpp"

double thermalXv, pi, m, Mp, gstar, Yeq;

//Define boltzmann equation in differential form 
double dydx(double x, double y){
thermalXv = 1;
pi = 3.14;
m=1;
Mp=1;
gstar=1;
Yeq=1;
		return sqrt(pi/45)*(m*Mp/(x*x))*gstar*thermalXv*((Yeq*Yeq)-(y*y));
		
}

//Use rk45 algorithm to solve ODE
int main() {
float x, y, z, h, k1, k2, k3, k4, k5, k6, yn;
int i, N;
x=3;
y=6;
h=0.1;
N=5;
	for(int i=0; i<=N; i++){
		k1 = h*dydx(x,y);
		k2 = h*dydx(x+(h/4),y+(k1/4));
		k3 = h*dydx(x+(3*h/8),y+(3*k1/32)+(9*k2/32));
		k4 = h*dydx(x+(12*h/13),y+(1932*k1/2197)-(7200*k2/2197)+(7296*k3/2197));
		k5 = h*dydx(x+h,y+(439*k1/216)-(8*k2)+(3680*k3/513)-(845*k4/4104));
		k6 = h*dydx(x+(h/2),y-(8*k1/27)+(2*k2)-(3544*k3/2565)+(1859*k4/4104)-(11*k5/40));
		yn = y+(25*k1/216)+(1408*k3/2565)+(2197*k4/4101)-(k5/5);
		z = y+(16*k1/135)+(6656*k3/12825)+(28561*k4/56430)-(9*k5/50)+(2*k6/55);
		//h = h*pow((tol*h/2*abs(z-y) ,1/4);
		cout<<x<<"\t"<<y<<"\t"<<yn<< endl;
		x = x+h;
		y = yn;
	}
}



