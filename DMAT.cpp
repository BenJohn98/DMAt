/*
C++ code for computing thermal abundance of dark matter candidates.
Initially written for assessment in PHYS4080 at The University of Queensland, semester 2, 2020.
Written by: Benjamin Oxley

*/

//include all required libraries
using namespace std;
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>


//define pi as a global variable
double pi = 3.14;	

//define prototype functions for main
int boltz(double x, const double y[], double f[], void *params_bolt);
double Yeq (double x);
double fthermalXv(double x, double s, void *params);

//main function uses GSL to solve boltzmann ODE. Outputs the values of the function thermal abundance at different values of x. 
//Stepping values as well as limits on x are defined in main.
int main() {


int dimension = 1; //dimension of ODE system
double eps_abs = 1.0E-8; //set absolute error
double eps_rel = 1.0E-10; //set relative error

const gsl_odeiv_step_type *type_bolt = gsl_odeiv_step_rkf45; //define method of ODE solver (in this case Runge-Kutta-Fehlberg or rk45)


gsl_odeiv_step *step_bolt = gsl_odeiv_step_alloc (type_bolt, dimension);
gsl_odeiv_control *control_bolt = gsl_odeiv_control_y_new (eps_abs, eps_rel);
gsl_odeiv_evolve *evolve_bolt = gsl_odeiv_evolve_alloc (dimension);

gsl_odeiv_system sys; //name system of ode's


double m = 1.6E-19; //set mass of particle
double Mp = 1.2E-19; //Plank mass: constant
double gstar = 1.0; //constant
double thermalXv = 0.001; //integer thermal cross section times velocity, becomes a function for scalar singlet dark matter plots (comment out this line)
double y[1]; //thermal abundance


double x, x_stepped; //define stepping values
double xmin, xmax, dx; //define bounds

double h = 1e-6; //set timestep for solver

//give function, dimension, and parameters to GSL
sys.function = boltz; 
sys.dimension = dimension;
sys.params = &m, &Mp, &gstar, &thermalXv;


//set bounds on ODE solver
xmin = 1;
xmax = 100;
dx = 1;

//set initial value
y[1] = 3.6E-9;

//set starting x
x=xmin;

//open data file
ofstream datafile;
datafile.open ("solution.txt");

//for loop to solve ODE. At each time step, the appropriate GSL function is called to solve the ODE using rk45.
	for (x_stepped = xmin + dx; x_stepped <= xmax; x_stepped += dx){
		while (x < x_stepped){
			gsl_odeiv_evolve_apply(evolve_bolt, control_bolt, step_bolt, &sys, &x, x_stepped, &h, y);
		}
		printf ("%.5e %.5e \n", x, log(y[0])/log(y[1])); //print normalised dimensionless thermal abundance to terminal(uncomment/comment out this line to toggle)
		//printf ("%.5e %.5e \n", m, y[0]/(3.6E-9)); //print actual thermal abundance to terminal (uncomment/comment out this line to toggle)
		datafile << ("%.5e %.5e \r\n", x, log(y[0])/log(y[1])); //write above data to file
		datafile << ("%.5e %.5e \n", m, y[0]/(3.6E-9)); // write data to file

	}

//close data file
datafile.close();	

//free up ode solver memory
gsl_odeiv_evolve_free (evolve_bolt);
gsl_odeiv_control_free (control_bolt);
gsl_odeiv_step_free (step_bolt);



return 0;

}

//define equilibrium abundance equation
double Yeq (double x) {
	double heff, gfactor, yeq; //define variables
	heff=86.25; //set heff constant
	gfactor=2; //set gfactor constant
	yeq = (45/(4*pi*pi*pi*pi))*((x*x)/heff)*gfactor*gsl_sf_bessel_Kn_scaled(2,x); //compute equilibrium abundance for each x
	return yeq;
	
}

/* ******************UNCOMMENT FOR SCALAR SINGLET DARK MATTER******************
float aDhs2(float s) {
	double mh =125;
	double s = ;
	return 1/((s-mh*mh)(s-(mh*mh)+mh*mh*dwidth(mh)*dwidth(mh))
	
}

float sigv_cms(float s) {
	
	return (2*lhs*lhs*v0*v0/sqrt(s))*aDhs2*dwidth(sqrt(s))
	
}

double fthermalXv(double x, double s){
	double m = *(double *) params;
	double f = (x*s*sqrt(s-4*m*m)*gsl_sf_bessel_i1_scaled(x*sqrt(s)/m)*sigv_cms(s))/(16*(m*m*m*m*m)*((gsl_sf_bessel_i2_scaled(x))*(gsl_sf_bessel_i2_scaled(x))));
	return f;
}

int interpolate (){
		int N=18; //length of data arrays
		double x[18] = {80.0,90.0,100.0,110.0,120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0,220.0,240.0,260.0,280.0,300.0}; //define x data in array
		double y[18] = {1.99,2.22,2.48,2.85,3.51,4.91,8.17,17.3,83.1,380,631,1040,1430,2310,3400,4760,6430,8430}; //define y data in array
		
		gsl_interp_accel *acc = gsl_interp_accel_alloc (); //gsl interpolator setup
		const gsl_interp_type *t = gsl_interp_cspline_periodic;
		gsl_spline *spline = gsl_spline_alloc (t, N);
		
		int i; //define values for interpolator
		double x, y;
		
		gsl_spline_init (spline, x, y, N);
		
		//for loop to interpolate over range
		for (i=0; i<=100; i++){
			x = (1 - i / 100.0) * x[0] + (i / 100.0) * x[N-1];
			y = gsl_spline_eval (spline, xi, acc); 
		}
	//free up interpolator things
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
	
	
  return 0;
}
*/

//Define boltzmann equation in differential form 
int boltz(double x, const double y[], double f[], void *params_bolt){
	double m = *(double *) params_bolt; //point to params_bolt to fetch required values
	double Mp = *(double *) params_bolt;
	double gstar = *(double *) params_bolt;
	double thermalXv = *(double *) params_bolt;

	f[0] = sqrt(pi/45)*(m*Mp/(x*x))*gstar*thermalXv*((Yeq(x)*Yeq(x))-(y[0]*y[0]));	//Boltzmann ODE to be solved (see documentation for full explanation)
	
	return GSL_SUCCESS;
	
}






