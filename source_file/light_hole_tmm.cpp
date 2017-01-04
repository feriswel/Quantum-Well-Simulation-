#include "qwdesign_tmm.h"

using namespace std;

void CQuantumWell::FindLightHoleWaves(double E, double dEStep,  CTotalBoundedStates* pTotalStates, CEngLevels* pPrevLHEng ){
	
	CCOMPLEX* wave=new CCOMPLEX[m_nPoints];
	CCOMPLEX* dwave=new CCOMPLEX[m_nPoints];
	double* potential=new double[m_nPoints];

//	double dB1=m_dQWLength*E/2*1e-8;
//	double dB2=-m_dQWLength*E/2*1e-8;
	
	double dB1=-0.2; //m_dQWLength*E/2*1e-8;
	double dB2=-0.2;//-m_dQWLength*E/2*1e-8;

	double dVMax=E*m_dQWLength*1E-8; //E in unit of V/cm, length is in unit of A

	int i;
	for(i=0;i<m_nPoints;i++)
	{
		potential[i]=m_pVV[i]+dVMax/2-(dVMax/m_nPoints)*(i);
	}

	double dEngMax=m_dMaxHEnergy-E*m_dLHMax*1e-8;;
	double dEngMin=m_dMinHEnergy-E*m_dLHMin*1e-8;

	if(m_pLHEng==NULL)
	{
		FindLightHoleWavesInRange(dB1, dB2, potential, wave, dwave, dEngMax, dEngMin, dEStep, pTotalStates);
		m_pLHEng=new CEngLevels(pTotalStates);
	}
	else
	{
		//double* dRefEngLevels=pPrevLHEng->GetEngLevels();
		//for(int i=0;i<pPrevLHEng->GetNumberOfStates();i++)
		//{
		//	double dEngMaxLocal=min(dRefEngLevels[i]+0.001,dEngMax);
		//	double dEngMinLocal=max(dRefEngLevels[i]-0.005, dEngMin);
		//	if(dEngMaxLocal>dEngMinLocal){
		//		FindLightHoleWavesInRange(dB1, dB2, potential, wave, dwave, dEngMaxLocal, dEngMinLocal, dEStep/2, pTotalStates);
		//	}

		//}
		FindLightHoleWavesInRange(dB1, dB2, potential, wave, dwave, dEngMax, dEngMin, dEStep, pTotalStates);

		pPrevLHEng->refresh(pTotalStates);

	}

	if(wave!=NULL)
		delete wave;
	if(dwave!=NULL)
		delete dwave;
	if(potential!=NULL)
		delete potential;
}

void CQuantumWell::FindLightHoleWavesInRange(double dB1, double dB2, double* potential, CCOMPLEX* wave, CCOMPLEX* dwave, double dEngMax, double dEngMin, double dEStep, CTotalBoundedStates* pTotalStates )
{

	int nEngPoints=int((dEngMax-dEngMin)/dEStep);

	double* T=new double[nEngPoints+3];
	
	double incident_E=dEngMin;

	int i=0;

	for(i=0;i<3;incident_E+=dEStep)
	{
		T[i]=TMM_lh(potential, incident_E, wave, dwave, dB1, dB2);
		i++;
	}

	for(;incident_E<dEngMax;incident_E+=dEStep)
	{
		T[i]=TMM_lh(potential, incident_E, wave, dwave, dB1, dB2);
		if((T[i]>=T[i-1]) && (T[i-1]<=T[i-2]))
		{
			SeekMaxT_lh(potential,incident_E-dEStep-dEStep, incident_E-dEStep, incident_E, T[i-2], T[i-1], T[i], pTotalStates, dB1, dB2);
		}
		i++;
		
	}

	if(T!=NULL)
		delete T;

}

double CQuantumWell::TMM_lh(double *potential, double incident_E, CCOMPLEX *wave_function, CCOMPLEX *dwave_function, double dEBoundary1, double dEBoundary2)
{
	
	//wave_function		wave function of the partical
	//dwave_function    derivative of wave function
	//potential			pointer to potential distribution array in ev 
	//partical_mass		effective mass of the partical
	//incident_E		incident partical energy in ev
	//n					number of points
	//lenght			lenght of the structure

	int i;
	double kSq;

	double step_size=m_dStepSize*1e-10;
	double step_size_sq_half=step_size*step_size*0.5;
	double integ_wf;
	double normal_const;
	double T;
	int n=m_nPoints;


	if(incident_E<dEBoundary1) // half bounded state
	{
		double k_boundary1=sqrt(2*m_pMEFFlh[n-1]*(dEBoundary1-incident_E)*q)/hbarj;
		wave_function[0]=CCOMPLEX(1.0,0.0);	//initial condition
		dwave_function[0]=CCOMPLEX(k_boundary1,0);
	}
	else // unbounded state
	{
		double k_boundary1=sqrt(2*m_pMEFFlh[n-1]*(-dEBoundary1+incident_E)*q)/hbarj;
		wave_function[0]=CCOMPLEX(1.0,0.0);	//initial condition
		dwave_function[0]=CCOMPLEX(0,-k_boundary1);
	}
	

	for (i=1; i<n;i++)
	{
		kSq=4*m_pMEFFlh[i-1]*q*(incident_E-potential[i-1])*h_bar_sq_inv;
		wave_function[i]=(1-kSq*step_size_sq_half)*wave_function[i-1]+step_size*dwave_function[i-1];
		dwave_function[i]=(1-kSq*step_size_sq_half)*dwave_function[i-1]-kSq*step_size*wave_function[i-1];
		dwave_function[i]*=m_pMEFFlh[i]/(m_pMEFFlh[i-1]);
	}

	double k_boundary=sqrt(2*m_pMEFFlh[n-1]*(-dEBoundary2+incident_E)*q)/hbarj;
	CCOMPLEX temp=(wave_function[n-1]-dwave_function[n-1]/CCOMPLEX(0.0, k_boundary));

	T=4.0/norm(temp);
	
/****************************normalization******************************/


	integ_wf=0.0;					//integration of wavefunction
	
	for(i=1; i<n-1;i++)
		integ_wf+=step_size*real(wave_function[i]*conj(wave_function[i]));

	integ_wf=integ_wf+real(wave_function[0]*conj(wave_function[0]))*step_size*0.5+real(wave_function[n-1]*conj(wave_function[n-1]))*step_size*0.5;
	
	//T=integ_wf/real(wave_function[0]*conj(wave_function[0]));  //

	double gamma=fabs(4.0*integ_wf/norm(temp));
	double tao=gamma*m_pMEFFlh[n-1]/(hbarj*k_boundary);
	double dLineWidth=hbarj/(2.0*tao*q);


	normal_const=sqrt(step_size/integ_wf);

	for(i=0;i<n;i++)
	{
		wave_function[i]=wave_function[i]*normal_const;
	}

	return dLineWidth;
//	return -T;

}
double CQuantumWell::SeekMaxT_lh(double* potential, double a, double b, double c, double fa, double fb, double fc, CTotalBoundedStates* pTotalStates, double dB1, double dB2 ) {

	CCOMPLEX* wave=new CCOMPLEX[m_nPoints];
	CCOMPLEX* dwave=new CCOMPLEX[m_nPoints];

	double x;
	double b0;
	double fx;
	
	do
	{
		b0=b;

		double c1=(b-a);
		double c2=(b-c);
		double df1=(fb-fc);
		double df2=(fb-fa);

		double C=c1*df1;
		double D=c2*df2;
		double A=c1*C;
		double B=c2*D;

		x=b-0.5*(A-B)/(C-D);


		fx=TMM_lh(potential,x,wave,dwave, dB1, dB2);

		if(x<b)
		{
			if((fa>fx) && (fb>fx))
			{
				//a=a;
				c=b;
				b=x;

			//	fa=fa;
				fc=fb;
				fb=fx;
			}
			else
			{
				a=x;
			//	c=c;
			//	b=b;

				fa=fx;
			//	fc=fc;
			//	fb=fb;
			}
		}
		else
		{
			if((fx>fb) && (fa>fb))
			{
				c=x;
				fc=fx;
			}
			else
			{
				a=b;
				b=x;
			//	c=c;

				fa=fb;
				fb=fx;
			//	fc=fc;
			}
		}
	}while(fabs(x-b0)>ENGYTOLERANCE); //0.000001 is the tolerance

	if(fx<BOUNDEDLINEWIDTH)// check bounded state, we assume the bounded state has a linewidth of less than 0.1 meV
		pTotalStates->AddState(new CState(x,wave,m_nPoints,fx));
	
//	pTotalStates->AddState(new CState(x,wave,m_nPoints,0.0001));

	if(dwave!=NULL)
		delete dwave;

	return fx;
}


