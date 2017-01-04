#include "qwdesign_tmm.h"

using namespace std;


void CQuantumWell::f012(double alfa, double mredlp, double* cwave_mag_sq, double* hwave_mag_sq,  double& f0, double& f1, double& f2)
{
	double c=(m_dQWLength)*1E-10/(m_nPoints-1);	//c is the coefficient to make wavefunction from dimensionless to the dimension of meter

	int n=m_nPoints;

	double Q0=0.0;
	double Q1=0.0;
	double Q2=0.0;

	int nPoints=50;

	double rStep=16e-10;
//	double r_start=0.0;
	double r_End=rStep*nPoints;
	
	double r;

	for(r=rStep/2.0;r<r_End;r+=rStep)
	{
		double t2=0.0;

		for(int zh=0;zh<n;zh++)
		{
			double t1=0.0;
			
			for(int ze=0;ze<n;ze++)
			{
				t1+=cwave_mag_sq[ze]/sqrt(1+(ze-zh)*(ze-zh)*c*c/(r*r));
			}
			
			t2+=hwave_mag_sq[zh]*t1;

		}

		double temp=exp(-2.0*alfa*r)*rStep*t2;

		Q0+=temp;
		Q1+=temp*(-2.0*r);
		Q2+=temp*(4.0*r*r);
	}

	double Q0_p=0.0;
	double Q1_p=0.0;
	double Q2_p=0.0;

	double temp;

	for(r=rStep;r<=r_End;r+=rStep)
	{
		double t2=0.0;

		for(int zh=0;zh<n;zh++)
		{
			double t1=0.0;
			
			for(int ze=0;ze<n;ze++)
			{
				t1+=cwave_mag_sq[ze]/sqrt(1+(ze-zh)*(ze-zh)*c*c/(r*r));
			}
			
			t2+=hwave_mag_sq[zh]*t1;

		}

		temp=exp(-2.0*alfa*r)*rStep*t2;

		Q0_p+=temp;
		Q1_p+=temp*(-2.0*r);
		Q2_p+=temp*(4.0*r*r);
	}

	Q0_p-=temp/2.0;
	Q1_p-=temp*(-2.0*r)/2.0;
	Q2_p-=temp*(4.0*r*r)/2.0;

	Q0=(2.0*Q0+Q0_p)/3.0;
	Q1=(2.0*Q1+Q1_p)/3.0;
	Q2=(2.0*Q2+Q2_p)/3.0;

	f0=hbarj*hbarj*alfa*alfa/(2*mredlp)/q-q*alfa*alfa*Q0/(m_dielect*epsilon0*PI);
	f1=hbarj*hbarj*alfa/(mredlp)/q-q*alfa*(alfa*Q1+2*Q0)/(m_dielect*epsilon0*PI);
	f2=hbarj*hbarj/(mredlp)/q-q*(2*Q0+4*Q1*alfa+alfa*alfa*Q2)/(m_dielect*epsilon0*PI);
}

bool CQuantumWell::SeekExcitonBondEng_e_lh(double initGuess, CState* pEState, CState* pHState, double& maxAlfa, double& maxf ){
	
	double alfa=initGuess;
	double f0, f1, f2;
	
	double* cwave_mag_sq=new double [m_nPoints];
	double* hwave_mag_sq=new double [m_nPoints];

	CCOMPLEX* cwave=pEState->wave();
	CCOMPLEX* hwave=pHState->wave();

	for (int i=0;i<m_nPoints;i++)
	{
		cwave_mag_sq[i]=norm(cwave[i]);
		hwave_mag_sq[i]=norm(hwave[i]);
	}

	f012(alfa,m_dRedLp, cwave_mag_sq,hwave_mag_sq,f0,f1,f2);

	double alfa_next;

	if(abs(f2)<1e-60){
		cout<<"failed to find the max value of Q for E-H pair:\n";
		return 0;
	}

	alfa_next=alfa-f1/f2;

	int nIteration=1;
	const int MAXITERATION=40;

//	while(fabs(alfa_next-alfa)<maxError){
	for(;fabs(alfa_next-alfa)>maxError;){

		alfa=alfa_next;
		f012(alfa,m_dRedLp, cwave_mag_sq, hwave_mag_sq,f0,f1,f2);

		if(abs(f2)<1e-60 ||nIteration>MAXITERATION){
			cout<<"failed to find the max value of f for the E-H pair: \n";
			return 0;
		}

		alfa_next=alfa-f1/f2;
		nIteration++;
	}

	maxAlfa=alfa;
	maxf=f0;

	return 1;
}

bool CQuantumWell::SeekExcitonBondEng_e_hh(double initGuess, CState* pEState, CState* pHState, double& maxAlfa, double& maxf ){
	
	double alfa=initGuess;
	double f0, f1, f2;
	
	double* cwave_mag_sq=new double [m_nPoints];
	double* hwave_mag_sq=new double [m_nPoints];

	CCOMPLEX* cwave=pEState->wave();
	CCOMPLEX* hwave=pHState->wave();

	for (int i=0;i<m_nPoints;i++)
	{
		cwave_mag_sq[i]=norm(cwave[i]);
		hwave_mag_sq[i]=norm(hwave[i]);
	}

	f012(alfa,m_dRedHp, cwave_mag_sq,hwave_mag_sq,f0,f1,f2);

	double alfa_next;

	if(abs(f2)<1e-60){
		cout<<"failed to find the max value of Q for E-H pair:\n";
		return 0;
	}

	alfa_next=alfa-f1/f2;

	int nIteration=1;
	const int MAXITERATION=40;

//	while(fabs(alfa_next-alfa)<maxError){
	for(;fabs(alfa_next-alfa)>maxError;){

		alfa=alfa_next;
		f012(alfa,m_dRedHp, cwave_mag_sq, hwave_mag_sq,f0,f1,f2);

		if(abs(f2)<1e-60 ||nIteration>MAXITERATION){
			cout<<"failed to find the max value of f for the E-H pair: \n";
			return 0;
		}

		alfa_next=alfa-f1/f2;
		nIteration++;
	}

	maxAlfa=alfa;
	maxf=f0;

	return 1;
}

