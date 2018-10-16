#include "BH_cross_sections.h"
#include "stdio.h"
#include <stdlib.h>
#include <cmath>
#include <time.h>
#include <unistd.h>

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>

#define MAXSTRLEN 200
using namespace std;

// define constants:
const double pi = 3.1415626;
const double m_e = 0.511;
const double e_gamma = 60;
const double Z = 92;
const double A = 238;
const double m_r = A * 931.5;

// Below integrations based on observation that:
// Cross section value sum over (phi_p,phi_q,E)
// for fixed (Theta_p, Theta_q) times spherical angle
// has contribution mostly coming from region which at
// least one of Theta_p, Theta_q is very small.
// So sum over (Theta_p,Theta_q) over 2 strips with one
// of Theta_p,Theta_q very small should be a good approximation

// min_polar: the lower bound of (theta_p,theta_q) strip region
const double min_polar = 90;
// max_polar: the upper bound of (theta_p,theta_q) strip region
const double max_polar = 90;
//const double min_polar = 4;
//const double max_polar = 90;

// max_polar_sqr: in addition to strip region, 
// theta_p*theta_q < max_polar_sqr region is also included for xsec computation
const double max_polar_sqr = 64;

// max_phi_diff:  phi in range 180-|phi_p-phi_q| < max_phi_diff
const double max_phi_diff = 45;
// min_E: Minimum energy for electron/positron, thus energy range (min_E , E_gamma-min_E)
const double min_E = 0.511;

// number of divisions .etc
// Note: For computation time:
// if div_phi increased n times, computation time will be n^2 times
// if div_E increased n times, compuatation time will be n times
// For div_phi == 90, div_E == 20, approx. 1 min/event

// del_th: intervals of theta
const double del_th = 1;
const int div_phi = 90;
const int div_E = 20;
// div_phi: number of phi divisions over 360 degrees
// div_E: number of Energy divisions over 60 MeV

//For given Theta_p, Theta_q, compute sum over Phi_p,Phi_q,E_p,E_q
double Sum_Phi_E(double Theta_p, double Theta_q, bool cut = true);

//Parameters relates to cut
// c_min_polar: theta_p/theta_q smaller than c_min_polar will be ignored
const double c_min_polar = 6;
// c_max_t_th: t_theta_p/q (horizonal projection angle) larger than c_max_th will be ignored
const double c_max_t_th = 15;
// c_max_t_phi: t_phi_p/q (vertical projection angle) larger than c_max_phi will be ignored
const double c_max_t_phi = 2.2;
// Apply cuts for Sum_Phi_E()
// Note: 1. Input takes rad unit angles
// 	 2. Caution when pick phi division, 
//          such that there is no big chuck 
//          of solid angle block ignored(picked)
bool Cut(double th[2], double phi[2]);


//Solve positron energy
inline double solve_e_positron(double);

// Deg Rad conversion
inline double DToR(double deg){
	return deg*pi/180.;
}

inline double RToD(double rad){
	return rad*180/pi;
}

// main() function, follow procedures:
// 1. Initialize BH cross section codes

// 2. For fixed (Theta_p,Theta_q), compute differetial cross 
//    sections sum over Phi_p,Phi_q,E through function Sum_Phi_E() by using BH code. 

// 3. Sum result of 2. by sin(theta_p)*dtheta_p*sin(theta_q)*dtheta_q, that results in
//    total cross sections.

// 4. Store data and total cross section in step 2. 3 in data file.

int main(int argc, char *argv[]){

	/*
	// Test Sum_Phi_E function, ignored from main procedure
	double th_p;
	double th_q;
	// use th_p and th_q both 5 deg
	th_p = DToR(5.);
	th_q = DToR(5.);
	//cerr << "Sum for th_p = " << RToD(th_p) <<", th_q = " << RToD(th_q) << " is " << Sum_Phi_E(th_p, th_q) << endl; 
	*/

	// Test Code for function Cut()
	/*
	double c_th[2] = {5., 15.5};	
	double c_phi[2] = {5., (180-5.)};	
	cerr << "c_th_p = " << c_th[0] << " , c_phi_p = " << c_phi[0] << endl
	     << "c_th_q = " << c_th[1] << " , c_phi_q = " << c_phi[1] << endl;
	for(int i = 0; i < 2; i++){
		c_th[i] = DToR(c_th[i]);
		c_phi[i] = DToR(c_phi[i]);
	}
	cerr << "Cut = " <<  Cut(c_th, c_phi) << endl;
 	*/

	// FIXME:Why result time is not(larger) in unit s(econd) ?
	clock_t t_s = clock();

	// open file to store data
	string tag;
	cerr << "Please enter tag:" << endl;
	cin >> tag;
	string filename_p = "para."+tag+".dat";
	string filename_d = "data."+tag+".dat";

	// div_th: number of theta_p/q division (equals to max_polar/del_th)
	int div_th = (int)(max_polar/del_th);
	//cerr << "div_th = " << div_th << endl;
	double th[div_th];
	// file storing run parameters
	ofstream of_p;
	of_p.open(filename_p.c_str());
	of_p << setw(10) << "min_polar" << setw(10) << "max_polar" << setw(10) << "div_th" << setw(20) << "theta_div(degree)" << endl
	     << setw(10) << min_polar << setw(10) << max_polar << setw(10) << div_th << setw(20) << 
		(max_polar)/(double)div_th << endl << endl
	     << setw(10) << "div_phi" << setw(20) << "phi_div(degree)" <<endl
	     << setw(10) << "div_E" << setw(20) << "energy_div(MeV)" << endl
	     << setw(10) << "Cut used:" << endl;
	of_p.close();

	// open file storing data
	// format for each line: theta_p(degree) , theta_q(degree) , Sum_Phi_E(theta_p, theta_q)
	ofstream of_d;
	of_d.open(filename_d.c_str());
	

	// Note: th[i] in unit rad
	for(int i = 0; (double)i*del_th < max_polar; i++){
		th[i] = ((double)i+.5)*DToR(del_th);
		cerr << "th[i] = " << RToD(th[i]) << endl;
	}
	// Total cross section
	double xsec=0;

	//double sum0 = Sum_Phi_E(th[0],th[0])*sin(th[0])*DToR(del_th)*sin(th[0])*DToR(del_th);
	double sum0 = 1*sin(th[0])*DToR(del_th)*sin(th[0])*DToR(del_th);
	for(int i = 0; i < div_th; i++){
		cerr << endl << endl;
		for(int j = 0; j < div_th; j++){
			// Below: integrated between 2 strips.
			if( (!( th[i] >= DToR(min_polar) && th[j] >= DToR(min_polar) )) || 
				th[i]*th[j] <= DToR(DToR(max_polar_sqr)) ){
				//double sum = Sum_Phi_E(th[i], th[j])*sin(th[i])*DToR(del_th)*sin(th[j])*DToR(del_th);
				double sum = Sum_Phi_E(th[i], th[j],true)*sin(th[i])*DToR(del_th)*sin(th[j])*DToR(del_th);
				//double sum = 1*sin(th[i])*DToR(del_th)*sin(th[j])*DToR(del_th);
				//double sum = 1;
				if(sum != 0){
					cerr << "(theta_p, theta_q) = (" << RToD(th[i]) << "," << RToD(th[j]) << ")" << " ,Sum = " << sum 
				     	<< " ,Relative contribution "<< sum*100./sum0 << endl;
					of_d << setw(10) << RToD(th[i]) << setw(10) << RToD(th[j]) <<  setprecision(5) << setw(20) << sum 
				     	<< setw(20) << sum*100./sum0 << endl;
					xsec += sum;
				}
			}
		}
	}

	of_d << endl << endl << setw(10) << xsec << endl;
	of_d << "Run time = %.2f (s)" << (double)(clock()-t_s)/CLOCKS_PER_SEC << endl;
	of_d.close();
	cerr << "xsec = " << xsec << endl;






	/*
	// Code below: sum over theta, using uniform sphericla angle (or uniform in cos(theta))
	// Sum over del_th_p and del_th_q as other loop
        double del_x = ( cos(DToR(max_polar)) - cos(DToR(min_polar)) )/(double)div_th;
        //cerr << "cos max polar is " << cos(DToR(max_polar)) << " , cos min polar is " << cos(DToR(min_polar)) << endl;
	//cerr << "max polar is " << max_polar << endl;

	double th[div_th];
	double th_div[div_th+1];
	th_div[0] = DToR(min_polar);
	for(int i = 0; i < div_th; i++){
		th[i] = acos( cos(DToR(min_polar)) + ((double)i+0.5)*del_x );
		th_div[i+1] = acos( cos(DToR(min_polar)) + (i+1)*del_x );
	}
	for(int i = 0; i < div_th; i++)
		cerr << "i th theta = " << RToD(th[i]) << " degrees." <<  "i th theta division = " << RToD(th_div[i+1]) << " degrees." << endl;

	// Now found all theta value and division values, sum over cos(theta_p) and cos(theta_q)
	double sum;
	cerr << "Number of division in theta is "  << div_th << endl;
	for(int i = 0; i < div_th; i++)
		for(int j = 0; j < div_th; j++){
			sum += 1*sin(th[i])*(th_div[i+1]-th_div[i])*sin(th[j])*(th_div[j+1]-th_div[j]);
			//double sum_phi_e = Sum_Phi_E(th[i], th[j]);
			//sum += sum_phi_e; 
			//cerr << "For division (" << i << "," << j << "), (th[i],th[j]) = (" << th[i] << "," << th[j] << ")"
			//	<< " ,sum over phi and E result is " << sum_phi_e << endl;
	}
	cerr << "Result is " << sum << endl;
	*/
}




// Simple approximation assumes that recoil has 0 energy
inline double solve_e_positron(double E_p){
	return e_gamma-E_p;
}




// Note: input th_p and th_q is in unit rad
double Sum_Phi_E(double th_p, double th_q, bool cut){
	// i: phi_p index
	// j: phi_q index
	// k: E_p/E_q index
	double phi_p = 0;
	double phi_q = 0;
	double phi = 0;
	double E_p = 0;
	double E_q = 0;

	// Input angles for BH_cross_section use rad as unit
	// phi ranges from 0 to 360 deg or 0 to 2pi
	double del_phi = DToR(360./div_phi);
	double del_E = 60./div_E;

	//initialize phi_p, phi_q, E_p, E_q that they're center value of bin under summation

	// sum is the result of summation over phi_p,phi_q, E_p, E_q
	double sum = 0;
	BH_cross_sections *xs = new BH_cross_sections(Z, e_gamma);

	for(int i = 0; i < div_phi; i++){
		for(int j = 0; j < div_phi; j++){
			phi_p = i*del_phi + del_phi/2.;
			phi_q = j*del_phi + del_phi/2.;
			phi = phi_p - phi_q;

			//cerr << "th_p = " << RToD(th_p) << ",th_q = " << RToD(th_q) << endl;
			//cerr << "phi_p = " << RToD(phi_p) << ",phi_q = " << RToD(phi_q) << ", phi = " << RToD(phi) << endl;
			//cerr << "Relative Phi diff is " <<  abs(180-abs(RToD(phi_p)-RToD(phi_q))) << endl;

			if( abs(180-abs(RToD(phi_p)-RToD(phi_q))) < max_phi_diff){
				for(int k = 0; k < div_E; k++){
					E_p = k*del_E + del_E/2. ;
					E_q = solve_e_positron(E_p);

					if(E_p > min_E && E_p < e_gamma-min_E){
					//cerr << "E_p = " << E_p << ",E_q = " << E_q << endl;
					//cerr << "xsec = " << xs->xsec_full(E_p, E_q, th_p, th_q, phi) << endl << endl;
					sum += xs->xsec_full(E_p, E_q, th_p, th_q, phi);
				}
	}}}}
	// Then multiply by unit volume
	cerr << "del_E = " << del_E << ", del_phi = " << del_phi << endl; 
	sum *= pow(del_phi,2)*del_E;
	//cerr << "del_phi = " << RToD(del_phi) << " , del_E = " << del_E << "sum is " << sum << endl;
	return sum;
}

// Apply acceptance cut to range of theta, phi
bool Cut(double th[2], double phi[2]){
	bool c_polar = (th[0] > DToR(c_min_polar)) && (th[1] > DToR(c_min_polar));

	double t_x[2];
	double t_y[2];
	double max_t_x = tan(DToR(c_max_t_th));
	double max_t_y = tan(DToR(c_max_t_phi));
	for(int i = 0; i < 2; i++){
		t_x[i] = tan(th[i])*cos(phi[i]);
		t_y[i] = tan(th[i])*sin(phi[i]);
	}

	bool c_xy[2];
	for(int i = 0; i < 2; i++){
		// for x (horizontal), electron/positron only goes to right/left.
		// for y (vertical), any value in [-max_t_y, max_t_y]
	 	c_xy[i]= ( pow(-1,(double)i)*t_x[i] > 0)&&( pow(-1,double(i))*t_x[i] < max_t_x)&&(t_y[i] > -max_t_y)&&(t_y[i] < max_t_y);
	}

	bool cut =  c_polar && c_xy[0] && c_xy[1];
	cerr << "In function Cut(): " << endl 
	     << "c_min_polar = " <<  c_min_polar << endl
	     << "th[0] = " << RToD(th[0]) << ", phi[0] = " << RToD(phi[0]) << endl
	     << "th[1] = " << RToD(th[1]) << ", phi[1] = " << RToD(phi[1]) << endl << endl
	     << "c_max_t_x = " << max_t_x << ", c_max_t_y = " << max_t_y << endl
	     << "t_x[0] = " << t_x[0] << ", t_y[0] = " << t_y[0] << endl
	     << "t_x[1] = " << t_x[1] << ", t_y[1] = " << t_y[1] << endl
	     << "Thus cut = " << cut;
	return cut;
}
