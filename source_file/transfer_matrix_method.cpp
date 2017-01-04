#include "qwdesign_tmm.h"

using namespace std;

void LoadParticalMass(double* particalmass,double* layermass, double* width, double nLayers, int nPoints, double dQWLength)
{
	int j=0;
	int i=0;
	int k=0;

	double dStepSize=dQWLength/nPoints;

	for (j=0; (j<nPoints) && (i<nLayers); j++){

		particalmass[j]=layermass[i];

		if((j-k)>width[i]/dStepSize)
		{
			k=j;
			i++;
		}
	}
}




void LoadPotentialShape(double* cv, double* vv, double E, double* cptl, double* vptl, double* width, double nLayers, int nPoints, double dQWLength)
{
	double dVMax=E*dQWLength*1E-8; //E in unit of V/cm, length is in unit of A
	double dStepSize=dQWLength/nPoints;

	int i=0;
	int j=0;
	int k=0;

	for (j=0; (j<nPoints) && (i<nLayers); j++)
	{
		cv[j]=cptl[i]-dVMax/2+(dVMax/nPoints)*(j);			//always set the mid of the QW as the zero bias, with no change of potentials
		vv[j]=vptl[i]+dVMax/2-(dVMax/nPoints)*(j);			//cv and vv are for diagonal element for the matrix

		if ((j-k)>width[i]/dStepSize)
		{
			k=j;
			i++;
		}
	}

	if (j<nPoints)					//set the value of the rest few points, which are missed at the end boundary
		for (;j<nPoints;j++)
		{
			cv[j]=cv[j-1];
			vv[j]=vv[j-1];
		}
}




double CQuantumWell::TMM_e(double* potential, double incident_E, CCOMPLEX *wave_function, CCOMPLEX *dwave_function)
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
//	double eff_mass[n];
//	double r=;				//		band nonparabolicity
//	double E_np_inv=2*partical_mass*r*h_bar_sq_inv;
	double step_size=-m_dStepSize*1e-10;
	double step_size_sq_half=step_size*step_size*0.5;
	double integ_wf;
	double normal_const;
	double T;
	int n=m_nPoints;


	double dEBoundary1=1.0; //ev
	double k_boundary1=sqrt(-2*m_pMEFFe[n-1]*(-dEBoundary1+incident_E)*q)/hbarj;

	wave_function[n-1]=CCOMPLEX(1.0,0.0);	//initial condition
	dwave_function[n-1]=CCOMPLEX(-k_boundary1,0.0);////////////////////////
	

	for (i=n-2; i>=0;i--)
	{
		kSq=4*m_pMEFFe[i+1]*q*(incident_E-potential[i+1])*h_bar_sq_inv;
		wave_function[i]=(1-kSq*step_size_sq_half)*wave_function[i+1]+step_size*dwave_function[i+1];
		dwave_function[i]=(1-kSq*step_size_sq_half)*dwave_function[i+1]-kSq*step_size*wave_function[i+1];
		dwave_function[i]*=m_pMEFFe[i]/(m_pMEFFe[i+1]);
	}

	double dEBoundary=-0.2; //ev
	double k_boundary=sqrt(2*m_pMEFFe[n-1]*(-dEBoundary+incident_E)*q)/hbarj;

	CCOMPLEX temp=(wave_function[0]+dwave_function[0]/CCOMPLEX(0.0, k_boundary));
	double real_part=temp.real();
	double img_part=temp.imag();

	T=4.0/(real_part*real_part+img_part*img_part);
	
/****************************normalization******************************/


	integ_wf=0.0;					//integration of wavefunction
	
	for(i=1; i<n-1;i++)
		integ_wf+=step_size*real(wave_function[i]*conj(wave_function[i]));

	integ_wf=integ_wf+real(wave_function[0]*conj(wave_function[0]))*step_size*0.5+real(wave_function[n-1]*conj(wave_function[n-1]))*step_size*0.5;
	
	//T=integ_wf/real(wave_function[0]*conj(wave_function[0]));  //

	//CCOMPLEX show[1000];
	//for(int mm=0;mm<m_nPoints;mm++)
	//	show[mm]=wave_function[mm];

	double gamma=fabs(4.0*integ_wf/norm(temp));
	double tao=gamma*m_pMEFFe[n-1]/(hbarj*k_boundary);
	double dLineWidth=hbarj/(2.0*tao*q);


	normal_const=sqrt(step_size/integ_wf);

	for(i=0;i<n;i++)
	{
		wave_function[i]=wave_function[i]*normal_const;
	}

//	return -T;
	return dLineWidth;
}

double CQuantumWell::TMM_lh(double *potential, double incident_E, CCOMPLEX *wave_function, CCOMPLEX *dwave_function)
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


	double dEBoundary1=1.0; //ev
	double k_boundary1=sqrt(-2*m_pMEFFlh[n-1]*(-dEBoundary1+incident_E)*q)/hbarj;

	wave_function[0]=CCOMPLEX(1.0,0.0);	//initial condition
	dwave_function[0]=CCOMPLEX(-k_boundary1,0);
	

	for (i=1; i<n;i++)
	{
		kSq=4*m_pMEFFlh[i-1]*q*(incident_E-potential[i-1])*h_bar_sq_inv;
		wave_function[i]=(1-kSq*step_size_sq_half)*wave_function[i-1]+step_size*dwave_function[i-1];
		dwave_function[i]=(1-kSq*step_size_sq_half)*dwave_function[i-1]-kSq*step_size*wave_function[i-1];
		dwave_function[i]*=m_pMEFFlh[i]/(m_pMEFFlh[i-1]);
	}

	double dEBoundary=-0.05; //ev
	double k_boundary=sqrt(2*m_pMEFFlh[n-1]*(-dEBoundary+incident_E)*q)/hbarj;

	CCOMPLEX temp=(wave_function[n-1]+dwave_function[n-1]/CCOMPLEX(0.0, k_boundary));
	double real_part=temp.real();
	double img_part=temp.imag();

	T=4.0/(real_part*real_part+img_part*img_part);
	
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

//	return -T;
	return dLineWidth;

}
double CQuantumWell::TMM_hh(double *potential, double incident_E, CCOMPLEX *wave_function, CCOMPLEX *dwave_function)
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


	double dEBoundary1=1.0; //ev
	double k_boundary1=sqrt(-2*m_pMEFFhh[n-1]*(-dEBoundary1+incident_E)*q)/hbarj;

	wave_function[0]=CCOMPLEX(1.0,0.0);	//initial condition
	dwave_function[0]=CCOMPLEX(-k_boundary1,0);
	

	for (i=1; i<n;i++)
	{
		kSq=4*m_pMEFFhh[i-1]*q*(incident_E-potential[i-1])*h_bar_sq_inv;
		wave_function[i]=(1-kSq*step_size_sq_half)*wave_function[i-1]+step_size*dwave_function[i-1];
		dwave_function[i]=(1-kSq*step_size_sq_half)*dwave_function[i-1]-kSq*step_size*wave_function[i-1];
		dwave_function[i]*=m_pMEFFhh[i]/(m_pMEFFhh[i-1]);
	}

	double dEBoundary=-0.05; //ev
	double k_boundary=sqrt(2*m_pMEFFhh[n-1]*(-dEBoundary+incident_E)*q)/hbarj;

	CCOMPLEX temp=(wave_function[n-1]+dwave_function[n-1]/CCOMPLEX(0.0, k_boundary));
	double real_part=temp.real();
	double img_part=temp.imag();

	T=4.0/(real_part*real_part+img_part*img_part);
	
/****************************normalization******************************/


	integ_wf=0.0;					//integration of wavefunction
	
	for(i=1; i<n-1;i++)
		integ_wf+=step_size*real(wave_function[i]*conj(wave_function[i]));

	integ_wf=integ_wf+real(wave_function[0]*conj(wave_function[0]))*step_size*0.5+real(wave_function[n-1]*conj(wave_function[n-1]))*step_size*0.5;
	
	//T=integ_wf/real(wave_function[0]*conj(wave_function[0]));  //

	double gamma=fabs(4.0*integ_wf/norm(temp));
	double tao=gamma*m_pMEFFhh[n-1]/(hbarj*k_boundary);
	double dLineWidth=hbarj/(2.0*tao*q);


	normal_const=sqrt(step_size/integ_wf);

	for(i=0;i<n;i++)
	{
		wave_function[i]=wave_function[i]*normal_const;
	}

//	return -T;
	return dLineWidth;

}

void CQuantumWell::FindElectronWaves(double E, double dEStep, CTotalBoundedStates* pTotalStates ){
	
	CCOMPLEX* wave=new CCOMPLEX[m_nPoints];
	CCOMPLEX* dwave=new CCOMPLEX[m_nPoints];

	double* potential=new double[m_nPoints];

	double dVMax=E*m_dQWLength*1E-8; //E in unit of V/cm, length is in unit of A

	int i;

	int nB=int(m_dBLength/m_dStepSize);

	for(i=0;i<m_nPoints;i++)
	{
		potential[i]=m_pCV[i]-dVMax/2+(dVMax/m_nPoints)*(i);
	}
	
	double dEngMax=m_dMaxEEnergy;
	double dEngMin=m_dMinEEnergy;

	//dEngMax=potential[m_nPoints-nB+2];
	//dEngMin=potential[nB+2];
	//dEngMin+=0.001;
	
	int nEngPoints=int((dEngMax-dEngMin)/dEStep);

	double* T=new double[nEngPoints+3];
	
	double incident_E=dEngMin;

	i=0;
	
	for(i=0;i<3;incident_E+=dEStep)
	{
		T[i]=TMM_e(potential, incident_E, wave, dwave);
		i++;
	}

	for(;incident_E<dEngMax;incident_E+=dEStep)
	{
		T[i]=TMM_e(potential, incident_E, wave, dwave);
		if((T[i]>=T[i-1]) && (T[i-1]<=T[i-2]))
		{
			SeekMaxT_e(potential,incident_E-dEStep-dEStep, incident_E-dEStep, incident_E, T[i-2], T[i-1], T[i], pTotalStates);
		}
		i++;
		
	}

	if(wave!=NULL)
		delete wave;
	if(dwave!=NULL)
		delete dwave;
	if(potential!=NULL)
		delete potential;
	if(T!=NULL)
		delete T;
}

void CQuantumWell::FindLightHoleWaves(double E, double dEStep,  CTotalBoundedStates* pTotalStates ){
	
	CCOMPLEX* wave=new CCOMPLEX[m_nPoints];
	CCOMPLEX* dwave=new CCOMPLEX[m_nPoints];

	double* potential=new double[m_nPoints];

	double dVMax=E*m_dQWLength*1E-8; //E in unit of V/cm, length is in unit of A


	int i;
	for(i=0;i<m_nPoints;i++)
	{
		potential[i]=m_pVV[i]+dVMax/2-(dVMax/m_nPoints)*(i);

	}

	double dEngMax=m_dMaxHEnergy;
	double dEngMin=m_dMinHEnergy;

	//int nB=int(m_dBLength/m_dStepSize);
	//dEngMax=potential[nB-2];
	//dEngMin=potential[m_nPoints-nB-2];
	//dEngMin+=0.001;
	
	int nEngPoints=int((dEngMax-dEngMin)/dEStep);

	double* T=new double[nEngPoints+3];
	//double T[2000];
	
	double incident_E=dEngMin;

	i=0;
	
	for(i=0;i<3;incident_E+=dEStep)
	{
		T[i]=TMM_lh(potential, incident_E, wave, dwave);
		i++;
	}

	for(;incident_E<dEngMax;incident_E+=dEStep)
	{
		T[i]=TMM_lh(potential, incident_E, wave, dwave);
		if((T[i]>=T[i-1]) && (T[i-1]<=T[i-2]))
		{
			SeekMaxT_lh(potential,incident_E-dEStep-dEStep, incident_E-dEStep, incident_E, T[i-2], T[i-1], T[i], pTotalStates);
		}
		i++;
		
	}

	if(wave!=NULL)
		delete wave;
	if(dwave!=NULL)
		delete dwave;
	if(potential!=NULL)
		delete potential;
	if(T!=NULL)
		delete T;
}

double CQuantumWell::SeekMaxT_e(double* potential, double a, double b, double c, double fa, double fb, double fc, CTotalBoundedStates* pTotalStates){

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


		fx=TMM_e(potential,x,wave,dwave);

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

	if(fx<BOUNDEDLINEWIDTH) // check bounded state, we assume the bounded state has a linewidth of less than 0.1 meV
		pTotalStates->AddState(new CState(x,wave,m_nPoints,fx));

	if(dwave!=NULL)
		delete dwave;

	return fx;
}


double CQuantumWell::SeekMaxT_lh(double* potential, double a, double b, double c, double fa, double fb, double fc, CTotalBoundedStates* pTotalStates ) {

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


		fx=TMM_lh(potential,x,wave,dwave);

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

	if(dwave!=NULL)
		delete dwave;

	return fx;
}


