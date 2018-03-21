#include "BH_cross_sections.h"

// constructor
BH_cross_sections::BH_cross_sections(double Zin, double photon_energy)
{
	// set up required values
	NP = 200;
	xmin=0;
	xmax=0.9992;        //interpolation interval
	falpha = 1./137.036;		// fine structure constant
	pi=3.14159265358979323846264338;
	m_e = 0.510999; // electron mass in MeV/c2
		//arrays for interpolation data
	xa = Vec_DP(NP);
	yaF1 = Vec_DP(NP);
	yaF2 = Vec_DP(NP);
	y2F1 = Vec_DP(NP);
	y2F2 = Vec_DP(NP);
	yaG1Re = Vec_DP(NP);
	yaG1Im = Vec_DP(NP);
	y2G1Re = Vec_DP(NP);
	y2G1Im = Vec_DP(NP);
	yaG2Re = Vec_DP(NP);
	yaG2Im = Vec_DP(NP);
	y2G2Re = Vec_DP(NP);
	y2G2Im = Vec_DP(NP);

	set_Z(Zin); // also sets eta
	set_photon_energy(photon_energy); // also sets omega
	interp(); // do the interpolation
}
BH_cross_sections::BH_cross_sections()
{
	// defaults:
	// uranium, e_gamma = 60 MeV
	BH_cross_sections(92., 60.);
}
BH_cross_sections::~BH_cross_sections()
	{}

double BH_cross_sections::xsec_s(double e_p, double e_q, double th_p, double th_q, double phi_q)
	{
	// symmetric part of cross section
	// inputs in MeV and radians
	DP p_p = sqrt(e_p*e_p - m_e*m_e);
	DP p_p_trans = p_p * sin(th_p)/m_e;
	DP p_q = sqrt(e_q*e_q - m_e*m_e);
	DP p_q_trans = p_q * sin(th_q)/m_e;
	DP phi = phi_q;
	DP x = e_p/fe_gamma;
	double xsec = dsima_s(p_p_trans,p_q_trans,phi,x);	
	xsec += dsigmas2(p_p_trans,p_q_trans,phi,x);	
	return xsec;
	}
double BH_cross_sections::xsec_a(double e_p, double e_q, double th_p, double th_q, double phi_q)
	{
	// anti-symmetric part of cross section
	// inputs in MeV and radians
	DP p_p = sqrt(e_p*e_p - m_e*m_e);
	DP p_p_trans = p_p * sin(th_p)/m_e;
	DP p_q = sqrt(e_q*e_q - m_e*m_e);
	DP p_q_trans = p_q * sin(th_q)/m_e;
	DP phi = phi_q;
	DP x = e_p/fe_gamma;
	double xsec = IntMtotGauss(p_p_trans,p_q_trans,phi,x);	
	return xsec;
	}
double BH_cross_sections::xsec_full(double e_p, double e_q, double th_p, double th_q, double phi_q)
	{
	// Total cross section
	// inputs in MeV and radians
	DP p_p = sqrt(e_p*e_p - m_e*m_e);
	DP p_p_trans = p_p * sin(th_p)/m_e;
	DP p_q = sqrt(e_q*e_q - m_e*m_e);
	DP p_q_trans = p_q * sin(th_q)/m_e;
	DP phi = phi_q;
	DP x = e_p/fe_gamma;
	double xsec = dsima_s(p_p_trans,p_q_trans,phi,x);	
	xsec += dsigmas2(p_p_trans,p_q_trans,phi,x);	
	xsec += IntMtotGauss(p_p_trans,p_q_trans,phi,x);	
	return xsec;
	}
	


///  below was part of MYfunc.h
//===========================================================
//DP BH_cross_sections::MYqgaus(DP (*func)(const DP,const DP,const DP,const DP,const DP),
//		const DP a, const DP b,const DP p,const DP q,const DP phi,const DP xx)
DP BH_cross_sections::MYqgaus( const DP a, const DP b,const DP p,const DP q,const DP phi,const DP xx)
{           //integrate a function by Gaussian quadratures
    static const DP x[]={0.1488743389816312,0.4333953941292472,
        0.6794095682990244,0.8650633666889845,0.9739065285171717};
    static const DP w[]={0.2955242247147529,0.2692667193099963,
        0.2190863625159821,0.1494513491505806,0.0666713443086881};
    int j;
    DP xr,xm,dx,s;

    xm=0.5*(b+a);
    xr=0.5*(b-a);
    s=0;
    for (j=0;j<5;j++) {
        dx=xr*x[j];
        //s += w[j]*(func(p,q,phi,xx,xm+dx)+func(p,q,phi,xx,xm-dx));
        s += w[j]*(Mtot(p,q,phi,xx,xm+dx)+Mtot(p,q,phi,xx,xm-dx));
    }
    return s *= xr;
}


complex<DP> BH_cross_sections::gamma(complex<DP> xx)//gamma function on complex number xx
{
    int j;
    if (xx.real()<=0)
        return 1;
    complex<DP> y,tmp,ser,x;
    complex<DP> cof[6]={76.18009172947146,-86.50532032941677,
         24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
         -0.5395239384953e-5};

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<6;j++)
    {
         y=y+complex<DP>(1,0);
         ser += cof[j]/y;
    }
    return exp(-tmp+log(2.5066282746310005*ser/x));
}

void BH_cross_sections::interp(void)    //make interpolation
{
    DP yF1p1,yF1pn,yF2p1,yF2pn,yG1Rep1,yG1Repn,yG1Imp1,yG1Impn,yG2Rep1,yG2Repn,yG2Imp1,yG2Impn;
    int i;
    for (i=0;i<=NP-1;i++) {
      xa[i]=(xmax-xmin)*(1-(NP-1.-i)*(NP-1.-i)/(NP-1)/(NP-1));
      yaF1[i]=NR::hypgeo(complex<DP>(0.,-eta),complex<DP>(0,eta),1,xa[i]).real();
      yaF2[i]=(1-xa[i])*NR::hypgeo(complex<DP>(1.,-eta),complex<DP>(1,eta),2,xa[i]).real();
      yaG1Re[i]=NR::hypgeo(complex<DP>(0.5,-eta),complex<DP>(0,eta),1,xa[i]).real();
      yaG1Im[i]=NR::hypgeo(complex<DP>(0.5,-eta),complex<DP>(0,eta),1,xa[i]).imag();
      yaG2Re[i]=NR::hypgeo(complex<DP>(1.5,-eta),complex<DP>(1,eta),2,xa[i]).real();
      yaG2Im[i]=NR::hypgeo(complex<DP>(1.5,-eta),complex<DP>(1,eta),2,xa[i]).imag();
    }

    yF1p1=eta*eta*NR::hypgeo(complex<DP>(1.,-eta),complex<DP>(1,eta),2,xa[0]).real();
    yF2p1=(1+eta*eta)/2*(1-xa[0])*NR::hypgeo(complex<DP>(2.,-eta),complex<DP>(2,eta),3,xa[0]).real()-NR::hypgeo(complex<DP>(1.,-eta),complex<DP>(1,eta),2,xa[0]).real();
    yG1Rep1=(complex<DP>(0.5,-eta)*complex<DP>(0,eta)*NR::hypgeo(complex<DP>(1.5,-eta),complex<DP>(1,eta),2,xa[0])).real();
    yG1Imp1=(complex<DP>(0.5,-eta)*complex<DP>(0,eta)*NR::hypgeo(complex<DP>(1.5,-eta),complex<DP>(1,eta),2,xa[0])).imag();
    yG2Rep1=(complex<DP>(1.5,-eta)*complex<DP>(1./2,eta/2)*NR::hypgeo(complex<DP>(2.5,-eta),complex<DP>(2,eta),3,xa[0])).real();
    yG2Imp1=(complex<DP>(1.5,-eta)*complex<DP>(1./2,eta/2)*NR::hypgeo(complex<DP>(2.5,-eta),complex<DP>(2,eta),3,xa[0])).imag();

    yF1pn=eta*eta*NR::hypgeo(complex<DP>(1.,-eta),complex<DP>(1,eta),2,xa[NP-1]).real();
    yF2pn=(1+eta*eta)/2*(1-xa[NP-1])*NR::hypgeo(complex<DP>(2.,-eta),complex<DP>(2,eta),3,xa[NP-1]).real()-NR::hypgeo(complex<DP>(1.,-eta),complex<DP>(1,eta),2,xa[NP-1]).real();
    yG1Repn=(complex<DP>(0.5,-eta)*complex<DP>(0,eta)*NR::hypgeo(complex<DP>(1.5,-eta),complex<DP>(1,eta),2,xa[NP-1])).real();
    yG1Impn=(complex<DP>(0.5,-eta)*complex<DP>(0,eta)*NR::hypgeo(complex<DP>(1.5,-eta),complex<DP>(1,eta),2,xa[NP-1])).imag();
    yG2Repn=(complex<DP>(1.5,-eta)*complex<DP>(1./2,eta/2)*NR::hypgeo(complex<DP>(2.5,-eta),complex<DP>(2,eta),3,xa[NP-1])).real();
    yG2Impn=(complex<DP>(1.5,-eta)*complex<DP>(1./2,eta/2)*NR::hypgeo(complex<DP>(2.5,-eta),complex<DP>(2,eta),3,xa[NP-1])).imag();

    // Call spline to get second derivatives
    NR::spline(xa,yaF1,yF1p1,yF1pn,y2F1);
    NR::spline(xa,yaF2,yF2p1,yF2pn,y2F2);
    NR::spline(xa,yaG1Re,yG1Rep1,yG1Repn,y2G1Re);
    NR::spline(xa,yaG1Im,yG1Imp1,yG1Impn,y2G1Im);
    NR::spline(xa,yaG2Re,yG2Rep1,yG2Repn,y2G2Re);
    NR::spline(xa,yaG2Im,yG2Imp1,yG2Impn,y2G2Im);

//    check interpolation
//      // Call splint for interpolations
//    double x,y,f;
//      cout << endl << setw(9) << "x" << setw(14) << "f(x)";
//      cout << setw(18) << "interpolation" << endl;
//      cout << fixed << setprecision(6);
//      for (i=0;i<10;i++) {
//        x=(-0.05+(i+1)/10.0)*1.0;
//        f=NR::hypgeo(complex<DP>(0.5,-eta),complex<DP>(0.,eta),1,x).imag();
//        NR::splint(xa,yaG1Im,y2G1Im,x,y);
//        cout << setw(12) << x << setw(13) << f;
//        cout << setw(13) << y << setw(13)  << (f-y)/f<<endl;
//      }
}


DP  BH_cross_sections::M(DP p,DP q,DP phi,DP x,DP lambda ) //calculate integrand (M2 see Methematica file Check_prog.m)
{                                                               //x=electron energy divided by omega
	DP m = 1;
    DP HF,HF1,HG1Re,HG1Im,HG2Re,HG2Im;
    DP chip=m*m/(p*p+m*m);
    DP chiq=m*m/(q*q+m*m);
    DP Q2=p*p+q*q+2*p*q*cos(phi);
    DP Q2par=pow((p*p/x+q*q/(1-x)+m*m/(x*(1-x)))/2/omega,2);
    DP u=1-Q2*chip*chiq;
    DP z=1-(1-u)*(1+lambda*lambda/sqrt(chip*chiq))/(1+chip*lambda*lambda/sqrt(chip*chiq))/(1+chiq*lambda*lambda/sqrt(chip*chiq));
    complex<DP> F,F1,G,G1,f1,f2,f3,M1,M2;
    NR::splint(xa,yaF1,y2F1,u,HF);
    NR::splint(xa,yaF2,y2F2,u,HF1);
    NR::splint(xa,yaG1Re,y2G1Re,z,HG1Re);
    NR::splint(xa,yaG1Im,y2G1Im,z,HG1Im);
    NR::splint(xa,yaG2Re,y2G2Re,z,HG2Re);
    NR::splint(xa,yaG2Im,y2G2Im,z,HG2Im);
    F=complex<DP>(HF,0);
    F1=complex<DP>(HF1*eta*eta,0);
    G=complex<DP>(HG1Re,HG1Im);
    G1=complex<DP>(HG2Re,HG2Im)*complex<DP>(0.5,-eta)*complex<DP>(0,eta);
    f1=(complex<DP>(0.5,-eta)*G-complex<DP>(1-z,0)*G1)*complex<DP>(1/(1+chip*lambda*lambda/sqrt(chip*chiq)),0);
    f2=(complex<DP>(0,eta)*G-complex<DP>(1-z,0)*G1)*complex<DP>(1/(1+chiq*lambda*lambda/sqrt(chip*chiq)),0);
    f3=complex<DP>((1-z)/(1+lambda*lambda/sqrt(chip*chiq)),0)*G1;
    M1=complex<DP>(-2*m*m*eta*eta*M_PI*eta/(sinh(M_PI*eta)*2*pow(M_PI*(Q2+Q2par),1.5))*pow(chip/(1+chip*lambda*lambda/sqrt(chip*chiq)),0.5)/(1-x),0)*gamma(complex<DP>(1,-eta))*gamma(complex<DP>(0.5,eta))*pow(complex<DP>((1+chip*lambda*lambda/sqrt(chip*chiq))/(1+chiq*lambda*lambda/sqrt(chip*chiq)),0),complex<DP>(0,eta));
    M2=(F*complex<DP>(0,(chip-chiq)*eta)+complex<DP>((1-chip-chiq),0)*F1)*(complex<DP>(4*x*(1-x)*chip+x*x+(1-x)*(1-x),0)*f1+complex<DP>(4*x*(1-x)*chiq+x*x+(1-x)*(1-x),0)*f2+complex<DP>(4*x*(1-x)+2*x*x+2*(1-x)*(1-x),0)*f3)+(f1-f2)*F*complex<DP>(0,eta*(1-u)*(x*x+(1-x)*(1-x)))-(f1+f2)*F1*complex<DP>(u*(x*x+(1-x)*(1-x)),0);
    return (M2*M1).imag()/pow(chip*chiq,0.25)/omega/omega;

}

DP BH_cross_sections::Mtot(const DP p,const DP q,const DP phi,const DP x,const DP lambda){
	//calculate integrand (Mtot see Methematica file Check_prog.m)
        return (M(p,q,phi,x,lambda/(1-lambda))-M(q,p,phi,1-x,lambda/(1-lambda)))/(1-lambda)/(1-lambda);
}


DP BH_cross_sections::IntMtotGauss (DP p,DP q,DP phi,DP x)//integrate Mtot by Gaussian quadratures
{
	static const DP a0 = 0., a25 = 0.25, a5 = 0.5, a75 = 0.75, a1 = 1.;
	DP sum = MYqgaus(a0,a25,p,q,phi,x);
	sum += MYqgaus(a25,a5,p,q,phi,x);
	sum += MYqgaus(a5,a75,p,q,phi,x);
	sum += MYqgaus(a75,a1,p,q,phi,x);
//    return MYqgaus(Mtot,0.,0.25,p,q,phi,x)+MYqgaus(Mtot,0.25,0.5,p,q,phi,x)+MYqgaus(Mtot,0.5,0.75,p,q,phi,x)+MYqgaus(Mtot,0.75,1.,p,q,phi,x);

	return sum;
}


DP BH_cross_sections::dsima_s(DP p,DP q,DP phi,DP x) //compute symmetric part of differential cross section
{
	DP m = 1;
    DP HF,HF1;
    DP chip=m*m/(p*p+m*m);
    DP chiq=m*m/(q*q+m*m);
    DP Q2=p*p+q*q+2*p*q*cos(phi);
    DP Q2par=pow((p*p/x+q*q/(1-x)+m*m/(x*(1-x)))/2/omega,2);
    DP u=1-Q2*chip*chiq;
    NR::splint(xa,yaF1,y2F1,u,HF);
    NR::splint(xa,yaF2,y2F2,u,HF1);
//    HF=NR::hypgeo(complex<DP>(0.,-eta),complex<DP>(0,eta),1,u).real();
//    HF1=(1-u)*NR::hypgeo(complex<DP>(1.,-eta),complex<DP>(1,eta),2,u).real();
    return 2*pow(eta/sinh(M_PI*eta),2)/(Q2+Q2par)/(Q2+Q2par)*(((1-u)*(x*x+(1-x)*(1-x))+2*x*(1-x)*(chip-chiq)*(chip-chiq))*eta*eta*HF*HF+(u*(x*x+(1-x)*(1-x))+2*x*(1-x)*(1-chip-chiq)*(1-chip-chiq))*HF1*HF1*eta*eta*eta*eta)/omega;
}

void  BH_cross_sections::deltag(DP p,DP q,DP phi,DP x,DP lambda,DP sign,complex<DP> *dg) //calculate delta g function (M2 see Methematica file Check_prog.m)
{                                                               //x=electron energy divided by omega
	DP m = 1;
    DP HF,HF1,HG1Re,HG1Im,HG2Re,HG2Im;
    DP chip=m*m/(p*p+m*m);
    DP chiq=m*m/(q*q+m*m);
    DP Q2=p*p+q*q+2*p*q*cos(phi);
    DP Q2par=pow((p*p/x+q*q/(1-x)+m*m/(x*(1-x)))/2/omega,2);
    DP u=1-Q2*chip*chiq;
    DP z=1-(1-u)*(1+lambda*lambda/sqrt(chip*chiq))/(1+chip*lambda*lambda/sqrt(chip*chiq))/(1+chiq*lambda*lambda/sqrt(chip*chiq));
    complex<DP> F,F1,G,G1,f1,f2,f3,factor;
    NR::splint(xa,yaF1,y2F1,u,HF);
    NR::splint(xa,yaF2,y2F2,u,HF1);
    NR::splint(xa,yaG1Re,y2G1Re,z,HG1Re);
    NR::splint(xa,yaG1Im,y2G1Im,z,HG1Im);
    NR::splint(xa,yaG2Re,y2G2Re,z,HG2Re);
    NR::splint(xa,yaG2Im,y2G2Im,z,HG2Im);
    F=complex<DP>(HF,0);
    F1=complex<DP>(HF1*eta*eta,0);
    G=complex<DP>(HG1Re,HG1Im*sign);
    G1=complex<DP>(HG2Re,HG2Im*sign)*complex<DP>(0.5,-eta*sign)*complex<DP>(0,eta*sign);
    f1=(complex<DP>(0.5,-eta*sign)*G-complex<DP>(1-z,0)*G1)*complex<DP>(1/(1+chip*lambda*lambda/sqrt(chip*chiq)),0);
    f2=(complex<DP>(0,eta*sign)*G-complex<DP>(1-z,0)*G1)*complex<DP>(1/(1+chiq*lambda*lambda/sqrt(chip*chiq)),0);
    f3=complex<DP>((1-z)/(1+lambda*lambda/sqrt(chip*chiq)),0)*G1;
    factor=complex<DP>(pow(pi,3./2)*eta*eta/(2*m*pow(Q2+Q2par,1./2))*2*pow(chip/(1+chip*lambda*lambda/sqrt(chip*chiq)),0.5)/(1-x)/pow(chip*chiq,0.25)/omega,0)*gamma(complex<DP>(1,-eta*sign))*gamma(complex<DP>(0.5,eta*sign))*pow(complex<DP>((1+chip*lambda*lambda/sqrt(chip*chiq))/(1+chiq*lambda*lambda/sqrt(chip*chiq)),0),complex<DP>(0,eta*sign));
    dg[0]=factor*(chip*f1+chiq*f2+f3);
    dg[1]=factor*(-chip*f1*cos(phi/2)*p+chiq*f2*cos(phi/2)*q);
    dg[2]=factor*(-chip*f1*sin(phi/2)*p-chiq*f2*sin(phi/2)*q);
//    cout<<"factor  "<<dg[2]<<endl;
}

void BH_cross_sections::Mtot1(DP p,DP q,DP phi,DP x,DP lambda,complex<DP> *dgtot)
{   //calculate integrand (Mtot1 see Methematica file Check_prog.m)
    complex<DP> a1[3],a2[3];
    deltag(p,q,phi,x,lambda/(1-lambda),1,a1);
    deltag(q,p,-phi,1-x,lambda/(1-lambda),-1,a2);

    dgtot[0]=(a1[0]+a2[0])/(1-lambda)/(1-lambda);
    dgtot[1]=(a1[1]-a2[1])/(1-lambda)/(1-lambda);
    dgtot[2]=(a1[2]-a2[2])/(1-lambda)/(1-lambda);
}


void BH_cross_sections::dgInt(const DP a, const DP b, DP p, DP q, DP phi, DP xx,  complex<DP>* res)
{                                                                       //integrate a function by Gaussian quadratures
    static const DP x[]={0.1488743389816312,0.4333953941292472,
        0.6794095682990244,0.8650633666889845,0.9739065285171717};
    static const DP w[]={0.2955242247147529,0.2692667193099963,
        0.2190863625159821,0.1494513491505806,0.0666713443086881};
    int j;
    DP xr,xm,dx;
    complex<DP> b1[3],b2[3];
    res[0]=complex<DP>(0,0);
    res[1]=complex<DP>(0,0);
    res[2]=complex<DP>(0,0);
    xm=0.5*(b+a);
    xr=0.5*(b-a);
    for (j=0;j<5;j++) {
        dx=xr*x[j];
        Mtot1(p, q, phi, xx, xm+dx,b1);
        Mtot1(p, q, phi, xx, xm-dx,b2);
        res[0]+= w[j]*(b1[0]+b2[0]);
        res[1]+= w[j]*(b1[1]+b2[1]);
        res[2]+= w[j]*(b1[2]+b2[2]);
//        s += w[j]*(func(p,q,phi,xx,xm+dx)+func(p,q,phi,xx,xm-dx));
    }
    res[0] *= xr;
    res[1] *= xr;
    res[2] *= xr;
}

DP BH_cross_sections::dsigmas2 (DP p,DP q,DP phi,DP x)//calculate modulo squared antisymmetric part of matrix element
{
	DP m = 1;
    complex<DP> c1[3],c2[3],c3[3],c4[3],gtot[3];
    dgInt(0,0.25,p,q,phi,x,c1);
    dgInt(0.25,0.5,p,q,phi,x,c2);
    dgInt(0.5,0.75,p,q,phi,x,c3);
    dgInt(0.75,1.,p,q,phi,x,c4);
    gtot[0]= c1[0]+c2[0]+c3[0]+c4[0];
    gtot[1]= c1[1]+c2[1]+c3[1]+c4[1];
    gtot[2]= c1[2]+c2[2]+c3[2]+c4[2];
    return pow(m/pi,4)/(2*omega)*((x*x+(1-x)*(1-x))*(gtot[2].imag()*gtot[2].imag()+gtot[2].real()*gtot[2].real()+gtot[1].imag()*gtot[1].imag()+gtot[1].real()*gtot[1].real())+gtot[0].imag()*gtot[0].imag()+gtot[0].real()*gtot[0].real());


}

