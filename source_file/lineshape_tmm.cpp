#include "qwdesign_tmm.h"
using namespace std;



// evaluation lineshape function using linear intepretation
double CLineShape::F_ev(double energy_in_ev)
{
	double PE=fabs(energy_in_ev);
	int k;
	double PE_k;
	if(PE<lsB1){
		k=int(PE/dStepPE1);
		//use linear intepretation 
		PE_k=k*dStepPE1;
		double h=1.0-(PE-PE_k)/dStepPE1;
		return (F1[k]-F1[k+1])*h+F1[k+1];
	}
	else if(PE<lsB2)
	{
		k=int((PE-lsB1)/dStepPE2);
		//use linear intepretation 
		PE_k=lsB1+k*dStepPE2;
		double h=1.0-(PE-PE_k)/dStepPE2;
		return (F2[k]-F2[k+1])*h+F2[k+1];
	}
	else if(PE<lsB3)
	{
		k=int((PE-lsB2)/dStepPE3);
		//use linear intepretation 
		PE_k=lsB2+k*dStepPE3;
		double h=1.0-(PE-PE_k)/dStepPE3;
		return (F3[k]-F3[k+1])*h+F3[k+1];
	}
	else
		return 0.0;
}

// evaluation lineshape function integral using linear intepretation
double CLineShape::intg_F_ev_to_inf(double dPE) // dPE=photon energy - energy difference in transition
{
	int k;
	double PE_k;
	if(dPE==0)
		return 0.5;

	double inte_F_from_0_plus=0.0;
	double PE=fabs(dPE); // PE is the magnitude of dPe!!!!!
	

	if(PE<lsB1){
		k=int(PE/dStepPE1);
		//use linear intepretation 
		PE_k=k*dStepPE1;
		double h=1.0-(PE-PE_k)/dStepPE1;
		inte_F_from_0_plus=(intgF1[k]-intgF1[k+1])*h+intgF1[k+1];
	}
	else if(PE<lsB2)
	{
		k=int((PE-lsB1)/dStepPE2);
		//use linear intepretation 
		PE_k=lsB1+k*dStepPE2;
		double h=1.0-(PE-PE_k)/dStepPE2;
		inte_F_from_0_plus=(intgF2[k]-intgF2[k+1])*h+intgF2[k+1];
	}
	else if(PE<lsB3)
	{
		k=int((PE-lsB2)/dStepPE3);
		//use linear intepretation 
		PE_k=lsB2+k*dStepPE3;
		double h=1.0-(PE-PE_k)/dStepPE3;
		inte_F_from_0_plus=(intgF3[k]-intgF3[k+1])*h+intgF3[k+1];
	}
	else
		inte_F_from_0_plus=0.5;
	
	if(dPE>0)
		return 0.5+inte_F_from_0_plus;
	else
		return 0.5-inte_F_from_0_plus;

}



void CLineShape::load_voigt_lineshape(double tL, double tG) // default is Voigt lineship, n=1 is lorentizian, and otherwise is Gaussian shape 
// tL is gammaL and tG is gammaG, they are in eV
{
	const double A=sqrt(log(2.0)); // A is a const sqrt(ln(2.0)), that is frequently used in voigt lineshape calculation
	const double B=2.0/sqrt(pi); // B is a const

	double mu=A*tL/tG;
	double v;
	
	double PE=0.0;
	
	F1[0]=exp(mu*mu)*erfc(mu);
	intgF1[0]=0;

	for(int i=1;i<=nW1;i++)
	{
		PE+=dStepPE1;
		v=A*PE/tG;

		// calculate the integration part of voigt shape
		double intg=0.0;
		double intg_step=min(1/(mu*200), v/100);
		for(double t=intg_step/2.0;t<v;t+=intg_step)
			intg+=intg_step*exp(t*t-v*v)*sin(2*mu*(v-t));
		F1[i]=F1[0]*exp(-v*v)*cos(2*mu*v)+B*intg; // calculate lineshape and F1[0]=exp(mu*mu)*erfc(mu)
		intgF1[i]=intgF1[i-1]+(F1[i]+F1[i-1])*0.5*dStepPE1; // calculate lineshap integral
	}

	PE=lsB1;
	F2[0]=F1[nW1];
	intgF2[0]=intgF1[nW1];

	for(int i=1;i<=nW2;i++)
	{		
		PE+=dStepPE2;
		v=A*PE/tG;
		// calculate the integration part of voigt shape
		double intg=0.0;
		double intg_step=min(1/(mu*200), v/100);
		for(double t=intg_step/2.0;t<v;t+=intg_step)
			intg+=intg_step*exp(t*t-v*v)*sin(2*mu*(v-t));
		F2[i]=F1[0]*exp(-v*v)*cos(2*mu*v)+B*intg;
		intgF2[i]=intgF2[i-1]+(F2[i]+F2[i-1])*0.5*dStepPE2; // calculate lineshap integral and and F1[0]=exp(mu*mu)*erfc(mu)

	}

	PE=lsB2;
	F3[0]=F2[nW2];
	intgF3[0]=intgF2[nW2];

	for(int i=1;i<=nW3;i++)
	{
		PE+=dStepPE3;
		v=A*PE/tG;
		// calculate the integration part of voigt shape
		double intg=0.0;
		double intg_step=min(1/(mu*200), v/100);
		for(double t=intg_step/2.0;t<v;t+=intg_step)
			intg+=intg_step*exp(t*t-v*v)*sin(2*mu*(v-t));
		F3[i]=F1[0]*exp(-v*v)*cos(2*mu*v)+B*intg;
		intgF3[i]=intgF3[i-1]+(F3[i]+F3[i-1])*0.5*dStepPE3; // calculate lineshap integral and and F1[0]=exp(mu*mu)*erfc(mu)
	}

	// Normalize the lineshape function
	// intgF3[nW3] should be equal to 0.5 in the normalized situation, thus
	double dNormalizeFactor=0.5/intgF3[nW3];
	
	for(int i=0;i<=nW1;i++)
	{
		F1[i]=F1[i]*dNormalizeFactor;
		intgF1[i]=intgF1[i]*dNormalizeFactor;
	}

	for(int i=0;i<=nW2;i++)
	{
		F2[i]=F2[i]*dNormalizeFactor;
		intgF2[i]=intgF2[i]*dNormalizeFactor;
	}
	for(int i=0;i<=nW3;i++)
	{
		F3[i]=F3[i]*dNormalizeFactor;
		intgF3[i]=intgF3[i]*dNormalizeFactor;
	}

}

void CLineShape::load_gaussian_lineshape(double tG) 
// tL is gammaL and tG is gammaG, they are in eV
{
	const double A=sqrt(log(2.0)/pi)/tG; // A and B are frequently used in gaussian lineshape calculation
	const double B=log(2.0)/(tG*tG);

	
	double PE=0.0;
	
	F1[0]=A;
	intgF1[0]=0;

	for(int i=1;i<=nW1;i++)
	{
		PE+=dStepPE1;
		F1[i]=A*exp(-B*PE*PE); 
		intgF1[i]=intgF1[i-1]+(F1[i]+F1[i-1])*0.5*dStepPE1; // calculate lineshap integral
	}

	PE=lsB1;
	F2[0]=F1[nW1];
	intgF2[0]=intgF1[nW1];

	for(int i=1;i<=nW2;i++)
	{		
		PE+=dStepPE2;
		F2[i]=A*exp(-B*PE*PE); 
		intgF2[i]=intgF2[i-1]+(F2[i]+F2[i-1])*0.5*dStepPE2; // calculate lineshap integral 

	}

	PE=lsB2;
	F3[0]=F2[nW2];
	intgF3[0]=intgF2[nW2];

	for(int i=1;i<=nW3;i++)
	{
		PE+=dStepPE3;
		F3[i]=A*exp(-B*PE*PE); 
		intgF3[i]=intgF3[i-1]+(F3[i]+F3[i-1])*0.5*dStepPE3; // calculate lineshap integral and and F1[0]=exp(mu*mu)*erfc(mu)
	}

	// Normalize the lineshape function
	// intgF3[nW3] should be equal to 0.5 in the normalized situation, thus
	double dNormalizeFactor=0.5/intgF3[nW3];
	
	for(int i=0;i<=nW1;i++)
	{
		F1[i]=F1[i]*dNormalizeFactor;
		intgF1[i]=intgF1[i]*dNormalizeFactor;
	}

	for(int i=0;i<=nW2;i++)
	{
		F2[i]=F2[i]*dNormalizeFactor;
		intgF2[i]=intgF2[i]*dNormalizeFactor;
	}
	for(int i=0;i<=nW3;i++)
	{
		F3[i]=F3[i]*dNormalizeFactor;
		intgF3[i]=intgF3[i]*dNormalizeFactor;
	}

}


void CLineShape::load_lorentzian_lineshape(double tL) 
// tL is gammaL and tG is gammaG, they are in eV
{
	const double A=tL/pi;

	double PE=0.0;
	
	F1[0]=1/(pi*tL);
	intgF1[0]=0;

	for(int i=1;i<=nW1;i++)
	{
		PE+=dStepPE1;
		F1[i]=A/(PE*PE+tL*tL); 
		intgF1[i]=intgF1[i-1]+(F1[i]+F1[i-1])*0.5*dStepPE1; // calculate lineshap integral
	}

	PE=lsB1;
	F2[0]=F1[nW1];
	intgF2[0]=intgF1[nW1];

	for(int i=1;i<=nW2;i++)
	{		
		PE+=dStepPE2;
		F2[i]=A/(PE*PE+tL*tL);; 
		intgF2[i]=intgF2[i-1]+(F2[i]+F2[i-1])*0.5*dStepPE2; // calculate lineshap integral 

	}

	PE=lsB2;
	F3[0]=F2[nW2];
	intgF3[0]=intgF2[nW2];

	for(int i=1;i<=nW3;i++)
	{
		PE+=dStepPE3;
		F3[i]=A/(PE*PE+tL*tL);; 
		intgF3[i]=intgF3[i-1]+(F3[i]+F3[i-1])*0.5*dStepPE3; // calculate lineshap integral and and F1[0]=exp(mu*mu)*erfc(mu)
	}

	// Normalize the lineshape function
	// intgF3[nW3] should be equal to 0.5 in the normalized situation, thus
	double dNormalizeFactor=0.5/intgF3[nW3];
	
	for(int i=0;i<=nW1;i++)
	{
		F1[i]=F1[i]*dNormalizeFactor;
		intgF1[i]=intgF1[i]*dNormalizeFactor;
	}

	for(int i=0;i<=nW2;i++)
	{
		F2[i]=F2[i]*dNormalizeFactor;
		intgF2[i]=intgF2[i]*dNormalizeFactor;
	}
	for(int i=0;i<=nW3;i++)
	{
		F3[i]=F3[i]*dNormalizeFactor;
		intgF3[i]=intgF3[i]*dNormalizeFactor;
	}

}
void CLineShape::load_modi_lorentzian_lineshape(double tL) 
// tL is gammaL and tG is gammaG, they are in eV
{
//	const double A=1.0/2.0/0.00339;
	double A=1.0/2.0/tL;
	double PE=0.0;
	
	F1[0]=A/2;
	intgF1[0]=0;

	for(int i=1;i<=nW1;i++)
	{
		PE+=dStepPE1;
		F1[i]=A/(1+exp(PE*A)); 
		intgF1[i]=intgF1[i-1]+(F1[i]+F1[i-1])*0.5*dStepPE1; // calculate lineshap integral
	}
	
	PE=lsB1;
	F2[0]=F1[nW1];
	intgF2[0]=intgF1[nW1];

	for(int i=1;i<=nW2;i++)
	{		
		PE+=dStepPE2;
		F2[i]=A/(1+exp(PE*A)); 
		intgF2[i]=intgF2[i-1]+(F2[i]+F2[i-1])*0.5*dStepPE2; // calculate lineshap integral 

	}

	PE=lsB2;
	F3[0]=F2[nW2];
	intgF3[0]=intgF2[nW2];

	for(int i=1;i<=nW3;i++)
	{
		PE+=dStepPE3;
		F3[i]=A/(1+exp(PE*A)); 
		intgF3[i]=intgF3[i-1]+(F3[i]+F3[i-1])*0.5*dStepPE3; // calculate lineshap integral and and F1[0]=exp(mu*mu)*erfc(mu)
	}

	// Normalize the lineshape function
	// intgF3[nW3] should be equal to 0.5 in the normalized situation, thus
	double dNormalizeFactor=0.5/intgF3[nW3];
	
	for(int i=0;i<=nW1;i++)
	{
		F1[i]=F1[i]*dNormalizeFactor;
		intgF1[i]=intgF1[i]*dNormalizeFactor;
	}

	for(int i=0;i<=nW2;i++)
	{
		F2[i]=F2[i]*dNormalizeFactor;
		intgF2[i]=intgF2[i]*dNormalizeFactor;
	}
	for(int i=0;i<=nW3;i++)
	{
		F3[i]=F3[i]*dNormalizeFactor;
		intgF3[i]=intgF3[i]*dNormalizeFactor;
	}

}

CLineShape::CLineShape(double LineWidth, double B1, double B2, double B3, int n1, int n2, int n3, int nType)
{

	nW1=n1;
	nW2=n2;
	nW3=n3;

	lsB1=B1;
	lsB2=B2;
	lsB3=B3;

	F1=new double [nW1+3];
	F2=new double [nW2+3];
	F3=new double [nW3+3];

	intgF1=new double [nW1+3];
	intgF2=new double [nW2+3];
	intgF3=new double [nW3+3];

	if(nType==0)
		load_lorentzian_lineshape(LineWidth);
	else if(nType==1)
		load_gaussian_lineshape(LineWidth);
	else if(nType==2)
		load_modi_lorentzian_lineshape(LineWidth);
	return;
}

CLineShape::CLineShape(double tL, double tG, double B1, double B2, double B3, int n1, int n2, int n3)
{

	nW1=n1;
	nW2=n2;
	nW3=n3;

	lsB1=B1;
	lsB2=B2;
	lsB3=B3;

	F1=new double [nW1+3];
	F2=new double [nW2+3];
	F3=new double [nW3+3];

	intgF1=new double [nW1+3];
	intgF2=new double [nW2+3];
	intgF3=new double [nW3+3];

	load_voigt_lineshape(tL,tG);

	return;
}

CLineShape::~CLineShape()
{
	if(F1!=NULL)
		delete F1;
	if(F2!=NULL)
		delete F2;
	if(F3!=NULL)
		delete F3;
	if(intgF1!=NULL)
		delete intgF1;
	if(intgF2!=NULL)
		delete intgF2;
	if(intgF3!=NULL)
		delete intgF3;
}



