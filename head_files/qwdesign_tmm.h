#pragma once

#ifndef _WIN32_WINNT		// Allow use of features specific to Windows XP or later.                   
#define _WIN32_WINNT 0x0501	// Change this to the appropriate value to target other versions of Windows.
#endif						
#include<windows.h>
#include<process.h>
#include <stdio.h>
#include <tchar.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;


#include "string.h"


#ifndef QWDESIGN_HEAD_TMM 
#define QWDESIGN_HEAD_TMM

#include <complex>


using namespace std;
typedef complex <double> CCOMPLEX;
//typedef double[502][502] TWODARRAY ;
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

const double dMaxEEnergy=0.15; // max electron wave energy
const double dMaxHEnergy=0.11; // max hole wave energy

const double dMinEEnergy=0.005; // max electron wave energy
const double dMinHEnergy=0.005; // max hole wave energy


const double q=1.6e-19; // charge of a single electron

const double hbar=6.5875E-16;			//in unit of eV
const double hbarj=1.054E-34;		//in unit if Jole 
const double h_bar_sq_inv=1/(hbarj*hbarj)*0.5;	//		1/planck constant square
const double h_bar_sq=hbar*hbar;

const double pi=3.141592654;
const double PI=3.141592654;
const double epsilon0=8.854E-12;
const double me0=9.11E-31;	//free electron mass 
//const double Eso=0.35;	//spin orbit split-off energy
const double Eso=0.11;	//spin orbit split-off energy for InP

const double RB1=0.7;					// rb1 abd rb2 can always be set arbitrarily
const double RB2=1.0;

const double Meh=1.5;//1.5;			//polarization factor for heavy hole
const double Mel=0.5;//0.5;			//polarization factor for light hole	Meh+Mel=2


const double lsB1=0.1;
const double lsB2=0.4;
const double lsB3=0.8;

const int nW1=2000;
const int nW2=2000;
const int nW3=2000;

const double dStepPE1=lsB1/nW1;
const double dStepPE2=(lsB2-lsB1)/nW2;
const double dStepPE3=(lsB3-lsB2)/nW3;

const int MAXLEVEL=10;
const double ENGYTOLERANCE=1e-12; //ev
const double BOUNDEDLINEWIDTH=20e-3;//ev
const double maxError=0.01;

const int n=500;

struct TWODARRAY{
	double a[501][501];
};
struct ALFAMATRIX{
	double a[MAXLEVEL][MAXLEVEL];
};

class SolveEigen	//wave function solver, eigenvalue method
{
public:
	int n;
	double* e;
	double* d;
	double z[2501][2501];
	
	SolveEigen(int m){
		e=new double[m];
		d=new double[m];
		n=m-1;
	};

	void tqli();
	void eigsrt();

};


class CState;
class CTotalBoundedStates;
class CEngLevels;

class CQuantumWell{
public:
	CEngLevels* m_pElectronEng;
	CEngLevels* m_pLHEng;
//	CEngLevels* m_pHHEng;

	ALFAMATRIX m_dGuessAlfaLH;
	ALFAMATRIX m_dGuessAlfaHH;

private:
	int m_nNumberOfHHState;

	int m_nLineWidthFunc;
	double m_dELHBroadening;
	double m_dEHHBroadening;

	char m_strQWFile[200];
	char m_strDirectory[200];

	double m_dEg;

	double m_dRedHp;
	double m_dRedLp;

	double m_dRedH;
	double m_dRedL;

	double* m_pMEFFe;
	double* m_pMEFFlh;
	double* m_pMEFFhh;

	double* m_pWidth;
	double* m_Ratio;

	int m_nLayers;
	
	double m_dQWLength;
	double m_dBLength;
	int m_nPoints;
	double m_dStepSize;

	double* m_pCV;
	double* m_pVV;

	double m_dFKEg;

	double m_ab1;
	double m_ab2;
	int m_nHW;

	double* m_pAbsorp0;
	double* m_pAbsorp0_tm;


	double m_dHWStep;
	double rstep;

	double m_nr;
	double m_dielect;

	double m_Eref;
	double m_Emax;
	double m_Emin;
	double m_nFields;

	double m_HWHM_G;
	double m_HWHM_L;

	double m_BufferLength;
	double m_dFieldStep;

	double m_dLEMax;
	double m_dLEMin;

	double m_dLHMin;
	double m_dLHMax;

	double m_dHHMin;
	double m_dHHMax;

// compatible with Guanghai's code
	double* cptl;
	double* vptl;
	double* meff;
	double* mhhef;
	double* mlhef;
	double* mhhefp;
	double* mlhefp;

	int m_nThreads;

	bool LoadQWData(char* qwfile, char*);

	bool m_bFail;

	double m_dMaxHEnergy;
	double m_dMinHEnergy;

	double m_dMaxHHEnergy;
	double m_dMinHHEnergy;

	double m_dMaxEEnergy;
	double m_dMinEEnergy;

	double SeekMaxT_e(double* potential, double a, double b, double c, double fa, double fb, double fc, CTotalBoundedStates* pTotalStates , double dB1, double dB2);
	double SeekMaxT_lh(double* potential, double a, double b, double c, double fa, double fb, double fc, CTotalBoundedStates* pTotalStates, double dB1, double dB2);
	double SeekMaxT_hh(double* potential, double a, double b, double c, double fa, double fb, double fc, CTotalBoundedStates* pTotalStates, double dB1, double dB2);

	double TMM_e(double* potential, double incident_E, CCOMPLEX *wave_function, CCOMPLEX *dwave_function, double dB1, double dB2);
	double TMM_lh(double* potential, double incident_E, CCOMPLEX *wave_function, CCOMPLEX *dwave_function, double dB1, double dB2);
	double TMM_hh(double* potential, double incident_E, CCOMPLEX *wave_function, CCOMPLEX *dwave_function, double dB1, double dB2);

protected:
	void FindElectronWavesInRange(double dB1, double dB2, double* potentail, CCOMPLEX* wave, CCOMPLEX* dWave, double dEngMax, double dEngMin, double dEStep, CTotalBoundedStates* pTotalStates );
	void FindLightHoleWavesInRange(double dB1, double dB2, double* potentail, CCOMPLEX* wave, CCOMPLEX* dWave, double dEngMax, double dEngMin, double dEStep, CTotalBoundedStates* pTotalStates );
	void FindHeavyHoleWavesInRange(double dB1, double dB2, double* potentail, CCOMPLEX* wave, CCOMPLEX* dWave, double dEngMax, double dEngMin, double dEStep, CTotalBoundedStates* pTotalStates );
	void FindHeavyHoleWavesEigen(double E, CTotalBoundedStates* pTotalStates);

public:

	int GetnHW(){return m_nHW;};
	int GetnThreads(){return m_nThreads;};
	double GetEMin(){return m_Emin;};
	double GetEMax(){return m_Emax;};
	double GetEFieldStep(){return m_dFieldStep;};

//	CQuantumWell(){m_pAbsorp0=NULL;};
	CQuantumWell(char* qwfile, char* strDir){ 
		LoadQWData(qwfile, strDir);
		m_pAbsorp0=NULL;
		m_pElectronEng=NULL;
		m_pLHEng=NULL;
		m_nNumberOfHHState=0;
		//m_pHHEng=NULL;
		for(int k=0;k<MAXLEVEL;k++)
			for(int j=0;j<MAXLEVEL;j++)
			{
				m_dGuessAlfaHH.a[k][j]=0.0;
				m_dGuessAlfaLH.a[k][j]=0.0;
			}
	};

	virtual ~CQuantumWell();


	void FindElectronWaves(double E, double dEStep, CTotalBoundedStates* pTotalStates, CEngLevels* pEngLevels=NULL);
	void FindLightHoleWaves(double E, double dEStep, CTotalBoundedStates* pTotalStates, CEngLevels* pEngLevels=NULL);
	void FindHeavyHoleWaves(double E, double dEStep, CTotalBoundedStates* pTotalStates, CEngLevels* pEngLevels=NULL );

	void SaveEWaves(char* filename);
	void SaveLHWaves(char* filename);
	void SaveHHWaves(char* filename);

	void exciton_e_lh(int i, int j);
	void exciton_e_hh(int i, int j);
	bool fail() {return m_bFail;};

	bool SeekExcitonBondEng_e_lh(double initGuess, CState* pEState, CState* pHState, double& maxAlfa, double& maxf); 
	bool SeekExcitonBondEng_e_hh(double initGuess, CState* pEState, CState* pHState, double& maxAlfa, double& maxf); 

	void f012(double alfa, double mredlp, double* cwave_mag_sq, double* hwave_mag_sq, double& f0, double& f1, double& f2);

	void AbsorptionCoefCal(double E, double* absorp, double* absorp_tm, CEngLevels* pElectronLevels=NULL, CEngLevels* pLHLevels=NULL, CEngLevels* pHHLevels=NULL, ALFAMATRIX* pGuessAlfaLH=NULL, ALFAMATRIX* pGuessAlfaHH=NULL, bool bSaveStates=1, bool bSaveAbsorption=1);
	void RefAbsorptionCoefCal(bool bSaveStates=1, bool bSaveAbsorption=1);

	void SaveAbsorptionCoefs(double E, double* pBB, double* pExct, double* pAbs);
	void SaveAbsorptionCoefs_tm(double E, double* pBB, double* pExct, double* pAbs);


	void DeltaRIndexCal(double E, double* absorp, double* pIndex, bool bSaveIndex=1);
	void DeltaRIndexCal_tm(double E, double* absorp, double* pIndex_tm, bool bSaveIndex=1);
	void SaveDeltaRIndex(double E, double* pIndex, int nRB1, int nRB2);
	void SaveDeltaRIndex_tm(double E, double* pIndex, int nRB1, int nRB2);

	void SaveHWPoints(); // stores the photon evergy points for absorption and index calculations
	void CQuantumWell::eigsrt(double* diag, double* wave_matrix);
	void CQuantumWell::tqli(double* diag, double* sup_diag, double* wave_matrix);

	static unsigned __stdcall ThreadStaticSolveIndex(void * pParam);

};

class CState{
	double m_dLineWidth;
	double m_dEig;
	CCOMPLEX* m_pWave;
	int m_nPoints;
public:
	CState(double eig, CCOMPLEX* pWave,int nPoints,double dLineWidth){m_dEig=eig; m_pWave=pWave;m_nPoints=nPoints;m_dLineWidth=dLineWidth;};
	CState(){m_dEig=0.0;m_pWave=NULL;m_nPoints=0;m_dLineWidth=0;};
	~CState(){if(m_pWave!=NULL) delete m_pWave;};
	double eig(){return m_dEig;};
	CCOMPLEX* wave(){return m_pWave;};
	double GetLineWidth(){return m_dLineWidth;};
};

class CTotalBoundedStates{
	int m_nTotalStates;
	int m_nPoints;
public:
	CState* m_pStates[MAXLEVEL];
	CTotalBoundedStates(int nPoints){
		m_nPoints=nPoints;
		m_nTotalStates=0;
		for(int i=0;i<MAXLEVEL;i++) 
			m_pStates[i]=NULL;
	};
	void AddState(CState* pStates){m_pStates[m_nTotalStates]=pStates;m_nTotalStates++;};
	int NumberOfStates(){return m_nTotalStates;};
	~CTotalBoundedStates(){
		for(int i=0;i<MAXLEVEL;i++)
			if(m_pStates[i]!=NULL)
				delete m_pStates[i];
	};
	void SaveStates(double E, char* dir, int nStateType);
};

class CLineShape{

	double* F1;
	double* F2;
	double* F3;

	double* intgF1;
	double* intgF2;
	double* intgF3;

	int nW1;
	int nW2;
	int nW3;

	double lsB1;
	double lsB2;
	double lsB3;

public:
	CLineShape(double LineWidth, double B1, double B2, double B3, int n1, int n2, int n3, int nType=0);
	CLineShape(double LorenLineWidth, double GaussLineWidth, double B1, double B2, double B3, int n1, int n2, int n3);
	~CLineShape();

	void load_voigt_lineshape(double tL, double tG);
	void load_gaussian_lineshape(double tG);
	void load_lorentzian_lineshape(double tL);
	void load_modi_lorentzian_lineshape(double tL);

	double F_ev(double ev);
	double intg_F_ev_to_inf(double dPE);

};
	

double erfc(double x);
double erf(double x);

struct THREADPARAM{
	CQuantumWell* pQW;

	CEngLevels* m_pPrevElectronEng;
	CEngLevels* m_pPrevLHEng;
	CEngLevels* m_pPrevHHEng;

	ALFAMATRIX  m_dGuessAlfaLH;
	ALFAMATRIX  m_dGuessAlfaHH;

	double E;
};

class CEngLevels{
	double* m_pEngLevels;
	int m_nNumberOfStates;
public:
//	CEngLevels(double* pEngLevels, int nMaxLevel){m_pEngLevels=pEngLevels;m_nMaxLevel=nMaxLevel;};
	CEngLevels(CTotalBoundedStates* pTotalStates)
	{
		m_nNumberOfStates=pTotalStates->NumberOfStates();
		m_pEngLevels=new double[m_nNumberOfStates];
		for(int i=0;i<m_nNumberOfStates;i++)
			m_pEngLevels[i]=pTotalStates->m_pStates[i]->eig();
	};

	CEngLevels(CEngLevels* pEng)
	{
		m_nNumberOfStates=pEng->m_nNumberOfStates;
		m_pEngLevels=new double[m_nNumberOfStates];
		for(int i=0;i<m_nNumberOfStates;i++)
			m_pEngLevels[i]=pEng->m_pEngLevels[i];
	};


	void refresh(CTotalBoundedStates* pTotalStates)
	{
		if(m_pEngLevels!=NULL) delete[] m_pEngLevels;
		m_nNumberOfStates=pTotalStates->NumberOfStates();
		m_pEngLevels=new double[m_nNumberOfStates];
		for(int i=0;i<m_nNumberOfStates;i++)
			m_pEngLevels[i]=pTotalStates->m_pStates[i]->eig();
	};


	~CEngLevels(){
		if(m_pEngLevels!=NULL) delete[] m_pEngLevels;
	};
	int GetNumberOfStates(){return m_nNumberOfStates;};
	double* GetEngLevels(){return m_pEngLevels;};
};

double F_K_absorption(double HW, double F, double Eg, double m_reduce, double n);

#endif 

