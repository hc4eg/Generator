#include "BH_cross_sections.h"
#include "stdio.h"
#include <iostream>
#include <fstream>
#include <iomanip>

#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <unistd.h>
#include <stdlib.h>

#define MAXSTRLEN 200
using namespace std;
//solve positron energy
inline double solve_e_positron(double e_gamma, double m_e, double m_r,
                        double e_p, double th_p, double th_q,
                        double phi_q);
//From all the parameter values, compute recoil nucleus kinetic energy
inline double ke_recoil(double e_gamma, double e_p, double e_q, double m_e,
                        double m_r, double th_p, double th_q, double phi_q);
//Generate random number in [0,1]
inline double randfloat(void);

// Angle cut
bool Cut(double th[2], double phi[2]);

const double pi = 3.14159265358979323846264338;
const double deg = pi / 180.;
// t_theta, t_phi: x and y projection angle
// So only pairs wihtin a rectangular projection specified by max_t_theta and max_t_phi will be computed
const double max_t_theta = 15 * deg;
const double max_t_phi = 2.2 * deg;

//const double min_polar[2] = { 4.*deg, (5.)*deg };
//const double max_polar[2] = { atan(sqrt(tan(max_t_theta) * tan(max_t_theta) + tan(max_t_phi) * tan(max_t_phi))), 
// atan(sqrt(tan(max_t_theta) * tan(max_t_theta) + tan(max_t_phi) * tan(max_t_phi)))  };

const double min_polar[2] = {  4*deg,  4*deg };
const double max_polar[2] = { 15*deg, 15*deg };
 
// experiment parameters
double max_eng = 42.;
double min_eng = 18.;


// argv[1]: Number of events to run
// argv[2]: tag of output data file

int main(int argc, char *argv[])
{
    //If value is X(degrees), the X*deg is in unit rad, thus can be used for trigonometric functions

    double m_e = 0.511;         //MeV
    double e_gamma = 60.;       // MeV
    double Z = 92;              // Z uranium
    double A = 238;             // A u 
    double m_r = A * 931.5;     // recoil mass

    //th_p,th_q: polar angles of electron and positron
    //e_p,e_q: total energy of electron and positron
    //ke_r: recoil nulcleus kinetic energy
    //phi_q: azimuthal angle difference of electron and positron
    double th_p, th_q, delta, e_p, e_q, ke_r, phi_q;

    char filename[MAXSTRLEN];
    char tag[MAXSTRLEN];

    BH_cross_sections *xs = new BH_cross_sections(Z, e_gamma);

    unsigned int seed;
    FILE *rf = fopen("/dev/urandom", "r");
    if(! rf ) {
        fprintf(stderr, "Can't open /dev/urandom.\n");
        exit(1);
    }
    if(fread(&seed, 1, sizeof(seed), rf) < 0) {
        fprintf(stderr, "Random seed not initialized correctly\n");
        exit(1);
    }
    fclose(rf);
    srandom(seed);
    strncpy(tag, "0", MAXSTRLEN);  // default tag

    // Number of pairs to generate: 1 Million
    long events_to_do = 1000000;
    if (argc > 1) {
        events_to_do = atoi(argv[1]);
    }
    if (argc > 2) {
        strncpy(tag, argv[2], MAXSTRLEN);
    }
    if (argc > 3) {
        fprintf(stderr, "Invalid argument list\n  %s <numEvents> [<unique_tag>]\n", argv[0]);
        exit(2);
    }
    long total_count = 0;
    long print_frequency = 1000;


    //Random nubmer generation will generate pairs with polar angle
    //In range [min_polar, max_polar]
    //And posible max_polar will be direction on the corners of rectangular specified by max_t_theta and max_t_phi
    // FIXME: The reason below min_polar, max_polar, cos_min_polar, rand_norm are in pairs is:
    // try to generate pairs in different range of theta_p and theta_q
    // rather than having theta_p, theta_q in same range
    double cos_min_polar[2] = { cos(min_polar[0]), cos(min_polar[1]) };
    double rand_norm[2] = { cos_min_polar[0] - cos(max_polar[0]), cos_min_polar[1] - cos(max_polar[1]) };
    //double maximum_xsec = 5.e-6;        // For xsec without factor, this needs to be determined
    double maximum_xsec = 1.3e3;        // For xsec with factor

    ifstream runin;
    ofstream runout;
    int runno = 0;
    double t_theta_p, t_theta_q, t_phi_p, t_phi_q;
    double phi_p_actual, phi_q_actual;

    //FILE *fp;
    //snprintf(filename, MAXSTRLEN, "events.run_%s.dat", tag);
    //fp = fopen(filename, "w");

    FILE *stdfp;
    snprintf(filename, MAXSTRLEN, "output.run_%s.dat", tag);
    stdfp = fopen(filename, "w");


    fprintf(stdfp, "Run %d: Events to do = %d\n", runno, events_to_do);
    fprintf(stdfp, "Run parameters:\n");
    //fprintf(stdfp, "Central Energy = %.2d MeV\n", central_energy);
    //fprintf(stdfp, "Percent acceptance range = %.1f%% to %.1f%%\n",
    //        min_percent, max_percent);
    fprintf(stdfp, "Maximum horizontal angle = %.2f degrees\n",
            max_t_theta / deg);
    fprintf(stdfp, "Maximum vertical angle = %.2f degrees\n",
            max_t_phi / deg);
    //fprintf(stdfp, "Minimum polar angle = %.2f degrees\n",
    //        min_polar / deg);
    //fprintf(stdfp, "Maximum polar angle = %.2f degrees\n",
    //        max_polar / deg);
    fprintf(stdfp, "Maximum Expected Cross Section = %.2g\n",
            maximum_xsec);
    fprintf(stdfp, "Random seed = %u\n", seed);
    fflush(stdfp);

    // Write table of max xsec in each region and its ralative total xsec
    bool b_table = true;
    if(b_table){
	// Ask for read pre.*.dat file to run simulation or start from scratch:
	cerr << "Read exsisting pre.*.data files? (Yy/Nn) ? " << endl;
	string s_read;
	cin >> s_read;
	bool b_read = ( ( s_read == "y" || s_read == "Y") ? true : false );
	cerr << "b_read = " << b_read << endl;
	

	cerr << "Enter tag for pre-run data file (either for read existing or open new file):" << endl;
	string otag = tag;
	cin >> otag;
	string ofname = "pre." + otag + ".dat";

	long starttime = time(NULL);

if(!b_read){
	ofstream* os = new ofstream;
	os->open(ofname.c_str());

	// b_preonly: true : Only generate prerun file ; false : generate pairs after prerun.
	bool b_pre_only = true;

	// Number of regions divided from min_polar to max_polar
	//int Nreg = 110;
	int Nreg = 22;
	double th_p_min[Nreg], th_p_max[Nreg], th_q_min[Nreg], th_q_max[Nreg];
	// max xsec values in each region, and their relative max xsec values
	double xsec_max[Nreg][Nreg];

	// Number of divisions in each p/q dimension when Area[][] in each region
	//double N_th_div = 100;
	//double N_phi_div = 180;
	double N_th_div = 25;
	double N_phi_div = 45;
	// In order to find area of cut region, which is hard to compute analytically, use numerical integration to compute Area
	double Area[Nreg][Nreg];
	double Area_max = 0;
	int N_xsec[Nreg][Nreg]; // number of points needed to compute total cross section, which is proportional to phase space area Area

	for(int i = 0; i < Nreg; i++)
		for(int j = 0; j < Nreg; j++){
			Area[i][j] = 0;
			N_xsec[i][j] = 0;
	}
	double N_xsec_max = 1e6; // N_xsec will scaled by N_xsec_max 

	// Compute th_p_min, th_p_max, th_q_min, th_q_max:
	for(int i = 0; i < Nreg; i++ ){
		th_p_min[i] = min_polar[0] + (double)i*(max_polar[0] - min_polar[0])/((double)Nreg);
		th_q_min[i] = min_polar[1] + (double)i*(max_polar[1] - min_polar[1])/((double)Nreg);
		th_p_max[i] = min_polar[0] + (double)(i+1)*(max_polar[0] - min_polar[0])/((double)Nreg);
		th_q_max[i] = min_polar[1] + (double)(i+1)*(max_polar[1] - min_polar[1])/((double)Nreg);
	}


	// Compute max and total cross section in each section:
	double xsec_tot[Nreg][Nreg];
	for(int i = 0; i < Nreg; i++)
		for(int j = 0; j < Nreg; j++){


			// N_loop: how many events to generate in order to find xsec_max and xsec_tot in each region
			// FIXME: the larger N_loop the better the precision of xsec_max and xsec_tot, needs larger values below
			//int N_loop = 1000;
			int N_loop = 1e4;
			// generate sphericla uniform angles for p,q
			double th[2];
			double ph[2];
			double phi;
			double cos_th[2];
			double energy[2];
			double xsec;
			xsec_tot[i][j] = 0;
			xsec_max[i][j] = 0;

			// FIXME: Since xsec_max always lies in th_p_min and th_q_min for the region(needs to be verified)
			// finding xsec_max on under this condition will make xsec_max finding much quicker
			for(int l = 0; l < N_loop; l++){

				do{
					th[0] = th_p_min[i];
					th[1] = th_q_min[j];
					for(int k = 0; k < 2; k++){ ph[k] = randfloat()*360.*deg; }
				}while( !Cut(th, ph) );
				// FIXME: for above line: since for finding xsec_max in certain region, th[0/1] using th_p/q_min instead of 
				// all possible theta for the region, it can be all being rejected by Cut(th,ph) in some cases.
					phi = ph[1] - ph[0];
					if(phi < 0.) phi+=360.*deg;
				//cerr << "In xsec_max finding: th_p = " << th[0] << ", th_q = " << th[1] << ", phi = " << phi << endl;

				energy[0] = min_eng + (max_eng - min_eng)*randfloat();
				//energy[1] = e_gamma - energy[0];
				energy[1] = solve_e_positron(e_gamma, m_e, m_r, energy[0], th[0], th[1], phi);

				xsec = 21.283*(pow(energy[0]*energy[1],2.)/pow(m_e,4.))*xs->xsec_full(energy[0], energy[1], th[0], th[1], phi);
				if( xsec_max[i][j] < xsec ) xsec_max[i][j] = xsec;
			}

			for(int l = 0; l < N_loop; l++){
				do{
					cos_th[0] = cos(th_p_min[i])-randfloat()*(cos(th_p_min[i]) - cos(th_p_max[i]));
					cos_th[1] = cos(th_q_min[j])-randfloat()*(cos(th_q_min[j]) - cos(th_q_max[j]));
					for(int k = 0; k < 2; k++){
						th[k] = acos( cos_th[k] );
						ph[k] = randfloat()*360.*deg;
					}
				}while( !Cut(th,ph) );

				phi = ph[1]-ph[0];
				if(phi < 0.) phi+=360.*deg;
				energy[0] = min_eng + (max_eng - min_eng)*randfloat();
				//energy[1] = e_gamma - energy[0];
				energy[1] = solve_e_positron(e_gamma, m_e, m_r, energy[0], th[0], th[1], phi);

				xsec = 21.283*(pow(energy[0]*energy[1],2.)/pow(m_e,4.))*xs->xsec_full(energy[0], energy[1], th[0], th[1], phi);
				xsec_tot[i][j] += xsec;
			}
		cerr << "Pre-run with index [i,j] = [" << i << "," << j << "] finished." << endl; 
		cerr << "Xsec_max[" << i << "," << j  << "] = " << xsec_max[i][j] << ", xsec_tot = " << xsec_tot[i][j] << endl;
	}


	// Compute Solid angle Area in each region:
	// FIXME: Important: When creating table, when finding xsec_max in each bin, we're generating same number of random pairs regardless of solid angle area of the portion is.
	(*os) << min_eng << " " << max_eng << endl;

	for(int i = 0; i < Nreg; i++ ){
		for(int j = 0; j < Nreg; j++){
			// divide both [ th_p_min[i], th_p_max[i] ] and [ th_q_min[j], th_q_max[j] ] into N_th_div==100 pieces and do numerical integration.
			double div_th_p = (th_p_max[i]-th_p_min[i])/((double)N_th_div);
			double div_th_q = (th_q_max[i]-th_q_min[i])/((double)N_th_div);

			double phi_range = 360*deg;
			double div_phi_p = phi_range/N_phi_div;
			double div_phi_q = div_phi_p;

			// In order to apply cut in th_p/q and phi_p/q, compute numerical integral by dividing both th_p/q and phi_p/q into small chunks.

			double temp_th_p, temp_th_q, temp_phi_p, temp_phi_q;
			for(int k = 0; k < N_th_div; k++){
				// th_p_min index i
				temp_th_p = th_p_min[i] + (((double)k) +0.5) * div_th_p;
				for(int l = 0; l < N_th_div; l++){
					// th_q_min index j
					temp_th_q = th_q_min[j] + (((double)k) +0.5) * div_th_q;
					for(int m = 0; m < N_phi_div; m++){
						// phi_p index m
						temp_phi_p = ((double)m +0.5) * (div_phi_p);
						for(int n = 0; n < N_phi_div; n++){
							// phi_q index n
							temp_phi_q = ( (double)n +0.5) * (div_phi_q);
							double th_a[2] = {temp_th_p, temp_th_q};
							double phi_a[2] = {temp_phi_p, temp_phi_q};
							// only angles statisfies cut will count in Area computation
							if( Cut(th_a, phi_a) )
								Area[i][j] += sin(temp_th_p) * div_th_p * sin(temp_th_q) * div_th_q * div_phi_p * div_phi_q;
						}
					}
				}
			}
			if(Area_max < Area[i][j]) Area_max = Area[i][j];
			cerr << "Area (in rad^2) unit of th_p = " << th_p_min[i]/deg << " to " << th_p_max[i]/deg << endl;
			cerr << "			 th_q = " << th_q_min[j]/deg << " to " << th_q_max[j]/deg << " is " << Area[i][j] << endl << endl; 
		}
	}
	cerr << "Maximum Area is " << Area_max << endl;
	



	// FIXME: Since the nature of rescaling computation of N_xsec below,
	// ofstream* os print data is delayed after all xsec_tot[][] is computed
	// How to circumvent this?
	// Print xsec_tot & Area in different column, and process data later when read them in pair generate process

	// For same number of Monte_Carlo sample points(UNIFORMLY DISTRIBUTED IN PHASE SPACE)
	// put into a portion when computing xsec_tot, it needs to scaled by Area to render N_xsec
	// (number of pairs needs to be generate in each portion)
	double weight  = 0;
	for(int i = 0; i < Nreg; i++)
		for(int j = 0; j < Nreg; j++){
			weight += xsec_tot[i][j]*Area[i][j];
	}
	int N_sum = 0;
	// Compute how many number of pairs need to be calculated for each portion.
	// FIXME: Since covertion from double to int will change result of N_sum, N_sum < N_xsec_max
	// and for small number of N_xsec_max, there could be discontinuiy of histogram result
	// from this effect?
	for(int i = 0; i < Nreg; i++)
		for(int j = 0; j < Nreg; j++){
			N_xsec[i][j] =(int)( N_xsec_max * (xsec_tot[i][j] * Area[i][j])/weight );
			cerr << "N_xsec[" << i << "][" << j << "] = " << N_xsec[i][j] << endl;
			N_sum+=N_xsec[i][j];
	}
	cerr << "Sum of N_xsec is " << N_sum <<endl;

	// pre-run data file stores in below format:
	// First line: min_eng max_eng
	// Following lines: th_p_min th_p_max th_q_min th_q_max xsec_max N_xsec
	for(int i = 0; i < Nreg; i++)
		for(int j = 0; j < Nreg; j++){
		(*os) << setw(10) << th_p_min[i]/deg << setw(10) << th_p_max[i]/deg 
		      << setw(10) << th_q_min[j]/deg << setw(10) << th_q_max[j]/deg
		      << setw(15) << xsec_max[i][j]  << setw(15) << N_xsec[i][j] <<
		      setw(15) << xsec_tot[i][j] <<  endl;
	}
	os->close();

if(!b_pre_only){ // Not only generate pre.*.dat file, but also use its result to generate pairs
	string pfilename = "pair." + otag + ".dat";
	os->open( pfilename.c_str() );
	
	int tot_count = 0;
	int print_freq = 1000;
	int events_to_do = 0;
	for(int i = 0; i < Nreg; i++)
		for(int j = 0; j < Nreg; j++)
			events_to_do += N_xsec[i][j];
	cerr << "Total number of pairs to be generated = " << events_to_do << endl << endl;
	long prevtime = time(NULL);

	// Now generate pairs based on information on xsec_max[][] and N_xsec[][] and print in pair.*.dat files
	for(int i  = 0; i < Nreg; i++)
		for(int j = 0; j < Nreg; j++){
			// Increase xsec_max, so that it's much harder to make ratio = xsec/xsec_max > 1 during pair generation
			xsec_max[i][j] *=2;
			double cos_th[2];
			double th[2];
			double ph[2];
			double energy[2];
			double phi;
			double xsec;
			double ratio;
			double count = 0;

			double t_th[2];
			double t_ph[2];
			int Nloop = 0;

			while(count < N_xsec[i][j]){
				Nloop++;
				// Add cut later when generate pairs in designated spherical angle region
				bool cut = true;
				// FIXME: This cut should be exactly the same as pre-run file
				do{
					cos_th[0] = cos(th_p_min[i])-randfloat()*(cos(th_p_min[i]) - cos(th_p_max[i]));
					cos_th[1] = cos(th_q_min[j])-randfloat()*(cos(th_q_min[j]) - cos(th_q_max[j]));
					for(int k = 0; k < 2; k++){
						th[k] = acos( cos_th[k] );
						ph[k] = randfloat()*360.*deg;
					}
				}while( !Cut(th,ph) );

				for(int k = 0; k < 2; k++){
					t_th[k] = atan( tan(th[k]) * cos(ph[k]) );
					t_ph[k] = atan( tan(th[k]) * sin(ph[k]) );
				}

				phi = ph[1]-ph[0];
				if(phi < 0.) phi+=360.*deg;

				energy[0] = min_eng + (max_eng - min_eng)*randfloat();
				//energy[1] = e_gamma - energy[0];
				energy[1] = solve_e_positron(e_gamma, m_e, m_r, energy[0], th[0], th[1], phi);

				xsec = 21.283*(pow(energy[0]*energy[1],2.)/pow(m_e,4.))*xs->xsec_full(energy[0], energy[1], th[0], th[1], phi);

				ratio = xsec/xsec_max[i][j];
				//cerr << "Current event: th_p = " << th[0]/deg << ", phi_p = " << ph[0]/deg << endl
				//     << "               th_q = " << th[1]/deg << ", phi_q = " << ph[1]/deg << endl
				//     << "               E_p  = " << energy[0] << ", E_q = " << energy[1] << endl;
				//cerr << " xsec = " << xsec << ", xsec_max = " << xsec_max[i][j] << ", ratio = " << ratio << endl;
				
				// Throwing error notes when ratio is over 1 (ie. current xsec value is over xsec_max simulated in pre-run)
				if(ratio > 1) cerr << "Exception, xsec exceed xsec_max in region th_p [" << th_p_min[i] << "," << th_p_max[i] <<
						"] , th_q [" << th_q_min[j] << "," << th_q_max[j] << "]." << endl;
				if( randfloat() > ratio ) continue;
				(*os) << energy[0]-m_e << " " << energy[1]-m_e << " " << t_th[0]/deg << " " << t_ph[0]/deg <<
				 " " << t_th[1]/deg << " " << t_ph[1]/deg << " " << xsec << endl;
				count++;
		 		cerr << "After " << Nloop << " trails, complete one event." << endl;
				Nloop = 0;
				tot_count++;

				if(tot_count % print_freq == 0){
					double percent = double(tot_count) / double(events_to_do) *100;
					long now  = time(NULL);
					long seconds = now - prevtime;
					prevtime = now;
					cerr << "Total counts = " << tot_count << "(" << percent << "%)" << ", these " << print_freq << " events run time = " << seconds << " secs" << endl;
				}
			}
			cerr << "Pair gen of [" << i << "," << j << "] portion has completed." << endl;
		}


	long now = time(NULL);
	long seconds = now -starttime;
	cerr << "Finished, total pairs = " << tot_count << " Total run time = "<< seconds << "sec" << endl;
}

}
else{
	ifstream* is = new ifstream;
	cerr << "Open file name " << ofname << endl;

	is->open(ofname.c_str());
	double e_min, e_max;
	int Nline = 0;

	// Step 1: find number of lines in the pre.*.dat file
	(*is) >> e_min >> e_max;
	is->ignore(500, '\n');
	while( !is->eof() ){
		is->ignore(500, '\n');
		Nline++;
	}
	Nline--;
	cerr << "Nline = "  << Nline << endl;

	// Step 2: read Nlines of th_p/q_min/max, xsec_max, N_xsec.
	double th_p_min[Nline], th_p_max[Nline], th_q_min[Nline], th_q_max[Nline], xsec_max[Nline], N_xsec[Nline];
	is->seekg(0);
	// skip 1st line of e_min, e_max
	is->ignore(500, '\n');
	for(int i = 0; i < Nline; i++){
		(*is) >> th_p_min[i] >> th_p_max[i] >> th_q_min[i] >> th_q_max[i] >> xsec_max[i] >> N_xsec[i];
		is->ignore(500, '\n');

		th_p_min[i] *= deg;
		th_q_min[i] *= deg;
		th_p_max[i] *= deg;
		th_q_max[i] *= deg;
	}
	is->close();




	// Step 3: open pair.*.dat file, generate random pairs and write into pair.*.dat file
	cerr << "Type run number of the current run:" << endl;
	cerr <<  "(it will add to tag of pair.*.dat file to distinguish exsisting pair.*.dat file) " << endl;
	string runNO;
	cin >> runNO;
	//char* runNO_char = itoa(runNO);

	ofstream* os = new ofstream;
	// Seems not C++ 11, so std::to_string() can't be used from int to string directly
	//string pfilename = "pair." + otag + "runNO" + string(runNO_char) + ".dat";

	string pfilename = "pair." + otag + "run" + runNO + ".dat";
	os->open(pfilename.c_str());

	int tot_count = 0;
	int print_freq = 1000;
	int events_to_do = 0;
	for(int i = 0; i < Nline; i++) events_to_do += N_xsec[i];
	cerr << "Total number of pairs to be generated = " << events_to_do << endl << endl;
	long prevtime = time(NULL);

	for(int i = 0; i < Nline; i++){
			// Generate random pairs:
			xsec_max[i] *=2.;
			double cos_th[2];
			double th[2];
			double ph[2];
			double energy[2];
			double phi;
			double xsec;
			double ratio;
			double count = 0;
			double t_th[2];
			double t_ph[2];
			double Nloop = 0;

			while(count < N_xsec[i]){
				Nloop ++;
				// FIXME: This cut should be exactly the same as pre-run file
				do{
					cos_th[0] = cos(th_p_min[i])-randfloat()*(cos(th_p_min[i]) - cos(th_p_max[i]));
					cos_th[1] = cos(th_q_min[i])-randfloat()*(cos(th_q_min[i]) - cos(th_q_max[i]));
					for(int k = 0; k < 2; k++){
						th[k] = acos( cos_th[k] );
						ph[k] = randfloat()*360.*deg;
					
					}
					//cerr << "Current event: th_p = " << th[0]/deg << ", phi_p = " << ph[0]/deg << endl
					//     << "               th_q = " << th[1]/deg << ", phi_q = " << ph[1]/deg << endl;
				}while( !Cut(th,ph) );
				for(int k = 0; k < 2; k++){
					t_th[k] = atan( tan(th[k]) * cos(ph[k]) );
					t_ph[k] = atan( tan(th[k]) * sin(ph[k]) );
				}
				phi = ph[1]-ph[0];
				if(phi < 0.) phi+=360.*deg;

				energy[0] = min_eng + (max_eng - min_eng)*randfloat();
				//energy[1] = e_gamma - energy[0];
				energy[1] = solve_e_positron(e_gamma, m_e, m_r, energy[0], th[0], th[1], phi);

				xsec = 21.283*(pow(energy[0]*energy[1],2.)/pow(m_e,4.))*xs->xsec_full(energy[0], energy[1], th[0], th[1], phi);
				ratio = xsec/xsec_max[i];
				//cerr << "Current event: th_p = " << th[0]/deg << ", phi_p = " << ph[0]/deg << endl
				//     << "               th_q = " << th[1]/deg << ", phi_q = " << ph[1]/deg << endl
				//     << "               E_p  = " << energy[0] << ", E_q = " << energy[1] << endl;
				//cerr << " xsec = " << xsec << ", xsec_max = " << xsec_max[i] << ", ratio = " << ratio << endl;
			
				// Throwing error notes when ratio is over 1 (ie. current xsec value is over xsec_max simulated in pre-run)
				if(ratio > 1) cerr << "Exception, xsec exceed xsec_max in region th_p [" << th_p_min[i] << "," << th_p_max[i] <<
						"] , th_q [" << th_q_min[i] << "," << th_q_max[i] << "]." << endl;
				if( randfloat() > ratio ) continue;

				// Print data into pair.*.dat file
				(*os) << energy[0]-m_e << " " << energy[1]-m_e << " " << t_th[0]/deg << " " << t_ph[0]/deg <<
				 " " << t_th[1]/deg << " " << t_ph[1]/deg << " " << xsec << endl;
				count++;
				//cerr << "After " << Nloop << " trails, complete one event." << endl;
				Nloop = 0;
				tot_count ++;

				if(tot_count % print_freq == 0){
					double percent = double(tot_count) / double(events_to_do) *100;
					long now  = time(NULL);
					long seconds = now - prevtime;
					prevtime = now;
					cerr << "Total counts = " << tot_count << "(" << percent << "%)" << ", these " << print_freq << " events run time = " << seconds << " sec " << endl;
				}
			}
	}
	long now = time(NULL);
	long seconds = now -starttime;
	cerr << "Finished, total pairs = " << tot_count << " Total run time = "<< seconds << "sec" << endl;
}


	}


}

double
solve_e_positron(double e_gamma, double m_e, double m_r,
                 double e_p, double th_p, double th_q, double phi_q)
{
    // solve for allowed positron energy
    double ke_r = 0.1;          // first guess
    double e_q = e_gamma - e_p - ke_r;
    // solve for actual ke_r and hence actual e_q
    double tol = 1.0e-7;
    double diff;
    double step = 0.01;
    // first guess
    double ke_r_cal =
        ke_recoil(e_gamma, e_p, e_q, m_e, m_r, th_p, th_q, phi_q);
    double old_diff = ke_r - ke_r_cal;
    //cout << "th_p_deg = " << th_p_deg << " delta = " << delta << endl;

    int N_itr = 0;
    while (1) {
        ke_r += step;
        e_q = e_gamma - e_p - ke_r;
        ke_r_cal =
            ke_recoil(e_gamma, e_p, e_q, m_e, m_r, th_p, th_q, phi_q);
        diff = ke_r - ke_r_cal;
        if (fabs(diff) < tol)
            break;
        if (diff * old_diff > 0.) {
            // same sign
            if (fabs(diff) > fabs(old_diff)) {
                // we are going in the wrong direction
                step = -step;
            }
        } else {
            // we passed through the solution
            step = -step / 2.;
        }
        old_diff = diff;
        N_itr++;
    }
    return (e_q);
}

double ke_recoil(double e_gamma, double e_p, double e_q, double m_e,
                 double m_r, double th_p, double th_q, double phi_q)
{
    //if(e_gamma - e_p - e_q < 0.) {
    //      printf("E_gamma = %.2f E_p = %.2f E_q = %.2f\n",e_gamma, e_p, e_q);
    //      return -1.;
    //      }
    //if(e_p < m_e) return -1.;
    //if(e_q < m_e) return -1.;
    double p_p = sqrt(e_p * e_p - m_e * m_e);
    double p_q = sqrt(e_q * e_q - m_e * m_e);
    double p_p_x = p_p * sin(th_p);
    double p_p_z = p_p * cos(th_p);
    double p_p_y = 0.;
    double p_q_x = p_q * sin(th_q) * cos(phi_q);
    double p_q_y = p_q * sin(th_q) * sin(phi_q);
    double p_q_z = p_q * cos(th_q);
    double p_r_z = e_gamma - p_p_z - p_q_z;
    double p_r_x = -(p_p_x + p_q_x);
    double p_r_y = -(p_p_y + p_q_y);
    double p_r_2 = p_r_x * p_r_x + p_r_y * p_r_y + p_r_z * p_r_z;
//      double p_r = sqrt(p_r_2);
//      printf("Electron Momentum = %f\n", p_p);
//      printf("  (x,y,z) = (%f,%f,%f)\n", p_p_x, p_p_y, p_p_z);
//      printf("Positron Momentum = %f\n", p_q);
//      printf("  (x,y,z) = (%f,%f,%f)\n", p_q_x, p_q_y, p_q_z);
//      printf("Recoil Momentum = %f\n", p_r);
//      printf("  (x,y,z) = (%f,%f,%f)\n", p_r_x, p_r_y, p_r_z);
    double ke_r = sqrt(p_r_2 + m_r * m_r) - m_r;
    return (ke_r);
}

double randfloat(void)
{
    double r = random() / double (RAND_MAX);
    return r;
}

// Note: input unit is rad
bool Cut(double th[2], double phi[2]){
	bool c_polar = (th[0] >= min_polar[0] && th[1] >= min_polar[1] );
	double t_x[2];
	double t_y[2];
	double max_t_x = tan( max_t_theta );	
	double max_t_y = tan( max_t_phi );

	for(int i = 0 ; i < 2; i++){
		t_x[i] = tan(th[i]) * cos(phi[i]);
		t_y[i] = tan(th[i]) * sin(phi[i]);
	}

	bool c_xy[2];
	for(int i = 0; i < 2; i++){
		c_xy[i] = ( pow(-1, (double)i) * t_x[i] > 0 ) && ( pow(-1,double(i))*t_x[i] < max_t_x  ) && (t_y[i] > -max_t_y) && (t_y[i] < max_t_y);
	}

	return c_polar && c_xy[0] && c_xy[1];
}
