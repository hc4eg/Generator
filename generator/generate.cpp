#include "BH_cross_sections.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include <sstream>
using namespace std;
//solve positron energy
inline double solve_e_positron(double e_gamma, double m_e, double m_r,
	double e_p, double th_p, double th_q, double phi_q);
//From all the parameter values, compute recoil nucleus kinetic energy
inline double ke_recoil(double e_gamma, double e_p, double e_q, double m_e, double m_r,
	double th_p, double th_q, double phi_q);
//Generate random number in [0,1]
inline double randfloat(void);

int main(int argc, char *argv[])
{
	double pi=3.14159265358979323846264338;
	//If value is X(degrees), the X*deg is in unit rad, thus can be used for trigonometric functions
	double deg = pi/180.;

	double m_e = 0.511; //MeV
	double e_gamma = 60.; // MeV
	double Z = 92; // Z uranium
	double A = 238; // A u 
	double m_r = A*931.5; // recoil mass
	
	//th_p,th_q: polar angles of electron and positron
	//e_p,e_q: total energy of electron and positron
	//ke_r: recoil nulcleus kinetic energy
	//phi_q: azimuthal angle difference of electron and positron
	double th_p, th_q, delta, e_p, e_q, ke_r, phi_q;

	char filename[100];

	BH_cross_sections * xs = new BH_cross_sections(Z, e_gamma);

	srandom((unsigned int) time(NULL));
	long events_to_do = 1000000;
	if(argc > 1) events_to_do = atoi(argv[1]);
	long total_count = 0;
	long print_frequency = 1000;

	// experiment parameters
	double central_energy = 30.;
	double max_percent = 40.;
	double min_percent = -40.;
	double max_eng = central_energy * (1.+max_percent/100.);
	double min_eng = central_energy * (1.+min_percent/100.);
	// t_theta, t_phi: x and y projection angle
	// So only pairs wihtin a rectangular projection specified by max_t_theta and max_t_phi will be computed
	const double max_t_theta = 15*deg;
	//const double max_t_phi = 2.2*deg;
	const double max_t_phi = 4.*deg;
	
	//Random nubmer generation will generate pairs with polar angle
	//In range [min_polar, max_polar]
	//And posible max_polar will be direction on the corners of rectangular specified by max_t_theta and max_t_phi
	double min_polar = 4.*deg;
	double max_polar = atan( sqrt(tan(max_t_theta)*tan(max_t_theta)
				+ tan(max_t_phi)*tan(max_t_phi)) );
	double cos_min_polar = cos(min_polar);
	double rand_norm = cos_min_polar - cos(max_polar);
	double maximum_xsec = 5.e-6; // This needs to be determined
					// checked for min_polar = 4 deg

	ifstream runin;
	ofstream runout;
	int runno = 0;
	// if the file "events.runno" exists - use its contents to determine run number
	runin.open("events.runno");
	if(runin.good())
		{
		runin>>runno;
		runin.close();
		}
	runno++;
	runout.open("events.runno", ios::trunc);
	if(runout.good())
		{
		runout << runno << endl;
		runout.close();
		}
	else
		{
		cout << "Cannot open events.runno for output." << endl;
		}

	double t_theta_p, t_theta_q, t_phi_p, t_phi_q;
	double phi_p_actual, phi_q_actual;

	FILE *fp;
	sprintf(filename, "events.run%.3d.dat",runno);
	fp = fopen(filename, "w");

	FILE *stdfp;
	sprintf(filename, "output.run%.3d.dat",runno);
	stdfp = fopen(filename, "w");

	fprintf(stdfp, "Run %d: Events to do = %d\n", runno, events_to_do);
	fprintf(stdfp, "Run parameters:\n");
	fprintf(stdfp, "Central Energy = %.2d MeV\n", central_energy); 
	fprintf(stdfp, "Percent acceptance range = %.1f%% to %.1f%%\n", min_percent, max_percent); 
	fprintf(stdfp, "Maximum horizontal angle = %.2f degrees\n", max_t_theta/deg);
	fprintf(stdfp, "Maximum vertical angle = %.2f degrees\n", max_t_phi/deg);
	fprintf(stdfp, "Minimum polar angle = %.2f degrees\n", min_polar/deg);
	fprintf(stdfp, "Maximum polar angle = %.2f degrees\n", max_polar/deg);
	fprintf(stdfp, "Maximum Expected Cross Section = %.2g\n", maximum_xsec);
	fflush(stdfp);

	int found_in_last = 0;
	long loops = 0;
	long starttime = time(NULL);
	while(total_count < events_to_do)
		{
		loops++;
		// choose random energies and angles for electron and positron
		// choose electron angles within spectrometer acceptance
		do	{
			double cos_polar = cos_min_polar - randfloat()*rand_norm;
			th_p = acos(cos_polar);
			phi_p_actual = randfloat()*360.*deg;
			double xx = sin(th_p)*cos(phi_p_actual);
			double yy = sin(th_p)*sin(phi_p_actual);
			double zz = cos_polar;
			t_theta_p = atan(xx/zz);
			t_phi_p = atan(yy/zz);
			} while ( abs(t_theta_p) > max_t_theta
				|| abs(t_phi_p) > max_t_phi
				|| th_p < min_polar);
		// choose positron angles within spectrometer acceptance
		do	{
			double cos_polar = cos_min_polar - randfloat()*rand_norm;
			th_q = acos(cos_polar);
			phi_q_actual = randfloat()*360.*deg;
			double xx = sin(th_q)*cos(phi_q_actual);
			double yy = sin(th_q)*sin(phi_q_actual);
			double zz = cos_polar;
			t_theta_q = atan(xx/zz);
			t_phi_q = atan(yy/zz);
			} while ( abs(t_theta_q) > max_t_theta
				|| abs(t_phi_q) > max_t_phi
				|| th_q < min_polar);

		phi_q = phi_q_actual - phi_p_actual;

		// NOTE: phi_q ranges from 0 to 2PI
		if(phi_q < 0.) phi_q += 360.*deg;
		
		// choose electron energy within spectrometer acceptance
		double ke_p = min_eng + (max_eng-min_eng)*randfloat();
		e_p = ke_p + m_e;
		// solve for allowed positron energy
		double e_q = solve_e_positron(e_gamma, m_e, m_r, e_p, th_p, th_q, phi_q);
		
		// calculate the cross section
		double xsec = xs->xsec_full(e_p, e_q, th_p, th_q, phi_q);	

		double ratio = xsec/maximum_xsec;
		if(ratio > 1.)
			{
			fprintf(stdfp, "Cross Section = %.8g > maximum cross section = %.8g\n", xsec, maximum_xsec);
			maximum_xsec = xsec;
			fflush(fp);
			fflush(stdfp);
			}
		if( randfloat() > ratio) continue;

		// write the event to a file
		// all angles are printed in unit degrees
		double ke_q = e_q - m_e;
		fprintf(fp, "%.4f %.4f %.4f %.4f %.4f %.4f %.4g\n",
			ke_p, ke_q, t_theta_p/deg, t_phi_p/deg, t_theta_q/deg, t_phi_q/deg, xsec);
		fflush(fp);
		total_count++;

		if(total_count%print_frequency == 0)
			{
			double percent = double(total_count)/double(events_to_do)*100.;
			long now = time(NULL);
			long seconds = now - starttime;
			fprintf(stdfp, "Total Counts = %d (%.0f%%) Loops since last = %d Runtime = %d sec.\n", total_count, percent, loops, seconds);
			fflush(fp);
			fflush(stdfp);
			loops = 0;
			}
		}
	fclose(fp);
	
	fprintf(stdfp, "\nFinished: total_count = %d\n", total_count);
	long now = time(NULL);
	long seconds = now - starttime;
	fprintf(stdfp, "Run time = %d sec\n", seconds);
	fclose(stdfp);
}

double
solve_e_positron(double e_gamma, double m_e, double m_r,
	double e_p, double th_p, double th_q, double phi_q)
{
	// solve for allowed positron energy
	double ke_r = 0.1; // first guess
	double e_q = e_gamma - e_p - ke_r;
	// solve for actual ke_r and hence actual e_q
	double tol = 1.0e-7;
	double diff;
	double step = 0.01;
	// first guess
	double ke_r_cal = ke_recoil(e_gamma, e_p, e_q, m_e, m_r, th_p, th_q, phi_q);
	double old_diff = ke_r - ke_r_cal;
	//cout << "th_p_deg = " << th_p_deg << " delta = " << delta << endl;

	int N_itr = 0;
	while(1)
		{
		ke_r += step;
		e_q = e_gamma - e_p - ke_r;
		ke_r_cal = ke_recoil(e_gamma, e_p, e_q, m_e, m_r, th_p, th_q, phi_q);
		diff = ke_r - ke_r_cal;
		if(fabs(diff) < tol) break;
		if( diff*old_diff > 0.)
			{
			// same sign
			if(fabs(diff) > fabs(old_diff))
				{
				// we are going in the wrong direction
				step = -step;
				}
			}
		else
			{
			// we passed through the solution
			step = -step/2.;
			}
		old_diff = diff;
		N_itr++;
		}
	return(e_q);
}

double ke_recoil(double e_gamma, double e_p, double e_q, double m_e, double m_r,
	double th_p, double th_q, double phi_q)
	{
	//if(e_gamma - e_p - e_q < 0.) {
	//	printf("E_gamma = %.2f E_p = %.2f E_q = %.2f\n",e_gamma, e_p, e_q);
	//	return -1.;
	//	}
	//if(e_p < m_e) return -1.;
	//if(e_q < m_e) return -1.;
	double p_p = sqrt(e_p*e_p - m_e*m_e);
	double p_q = sqrt(e_q*e_q - m_e*m_e);
	double p_p_x = p_p * sin(th_p);
	double p_p_z = p_p * cos(th_p);
	double p_p_y = 0.;
	double p_q_x = p_q * sin(th_q) * cos(phi_q);
	double p_q_y = p_q * sin(th_q) * sin(phi_q);
	double p_q_z = p_q * cos(th_q);
	double p_r_z = e_gamma - p_p_z - p_q_z;
	double p_r_x = -(p_p_x + p_q_x);
	double p_r_y = -(p_p_y + p_q_y);
	double p_r_2 = p_r_x*p_r_x + p_r_y*p_r_y + p_r_z*p_r_z;
//	double p_r = sqrt(p_r_2);
//	printf("Electron Momentum = %f\n", p_p);
//	printf("  (x,y,z) = (%f,%f,%f)\n", p_p_x, p_p_y, p_p_z);
//	printf("Positron Momentum = %f\n", p_q);
//	printf("  (x,y,z) = (%f,%f,%f)\n", p_q_x, p_q_y, p_q_z);
//	printf("Recoil Momentum = %f\n", p_r);
//	printf("  (x,y,z) = (%f,%f,%f)\n", p_r_x, p_r_y, p_r_z);
	double ke_r = sqrt(p_r_2 + m_r*m_r) - m_r;	
	return(ke_r);
	}

double randfloat(void)
	{
	double r = random()/double(RAND_MAX);
	return r;
	}
