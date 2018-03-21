#ifndef BH_cross_sections_h
#define BH_cross_sections_h

#include <cmath>
#include <complex>
#include "nr.h"
#include <ctime>
using namespace std;
using namespace NR;


class BH_cross_sections
{
	public:
	BH_cross_sections(); // constructor
	BH_cross_sections(double Zin, double photon_energy);
	~BH_cross_sections(); // destructor

	inline void set_Z(double val) { fZ = val; eta = fZ*falpha;}
	inline void set_photon_energy(double val) { fe_gamma = val; omega = fe_gamma/m_e;}

	double xsec_s(double e_p, double e_q, double th_p, double th_q, double phi_q);
	double xsec_a(double e_p, double e_q, double th_p, double th_q, double phi_q);
	double xsec_full(double e_p, double e_q, double th_p, double th_q, double phi_q);

	private:
	// private functions
	//DP MYqgaus(DP (BH_cross_sections::*func)(const DP,const DP,const DP,const DP,const DP),
	//	const DP a, const DP b,const DP p,const DP q,const DP phi,const DP xx);
	DP MYqgaus( const DP a, const DP b,const DP p,const DP q,const DP phi,const DP xx);
	complex<DP> gamma(complex<DP> xx);//gamma function on complex number xx
	void interp(void);    //make interpolation	
	DP  M(DP p,DP q,DP phi,DP x,DP lambda ); //calculate integrand 
	DP Mtot(DP p,DP q,DP phi,DP x,DP lambda);
	DP IntMtotGauss (DP p,DP q,DP phi,DP x); // calculates anit-symmetric part
	DP dsima_s(DP p,DP q,DP phi,DP x); // calculate symmetric part
	void deltag(DP p,DP q,DP phi,DP x,DP lambda,DP sign,complex<DP> *dg);
	void Mtot1(DP p,DP q,DP phi,DP x,DP lambda,complex<DP> *dgtot);
	void dgInt(const DP a, const DP b, DP p, DP q, DP phi, DP xx,  complex<DP>* res);
	DP dsigmas2 (DP p,DP q,DP phi,DP x);

	private:
	int NP;
	DP xmin,xmax;        //interpolation interval
	Vec_DP xa,yaF1,y2F1,yaF2,y2F2,
		yaG1Re,y2G1Re,yaG1Im,y2G1Im,
		yaG2Re,y2G2Re,yaG2Im,y2G2Im;

	
	DP fZ;		// Z
	DP fe_gamma;	// photon energy
	DP eta;		// Z*alpha
	DP m_e;		// electron mass
	DP omega;	// photon energy / m_e
	DP pi;
	DP falpha;	// fine structure constant
};

#endif
