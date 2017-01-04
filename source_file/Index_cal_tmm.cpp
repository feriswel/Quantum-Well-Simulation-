#include "qwdesign_tmm.h"

using namespace std;

/********************************dielectric constant calcualtion****************************/
void CQuantumWell::DeltaRIndexCal(double E, double* absorp, double* refract, bool bSaveIndex)
{
	double step=m_dHWStep;
	double* absorp0=m_pAbsorp0;

	double rb1=m_ab1+(m_ab2-m_ab1)*0.1;
	double rb2=m_ab2-(m_ab2-m_ab1)*0.1;

	int nRB1=int((rb1-m_ab1)/step);
	int nRB2=int((rb2-m_ab1)/step);

	rb1=m_ab1+double(nRB1)*step;
	rb2=m_ab1+double(nRB2)*step;

	int n2=m_nHW;
	
	for (int nHW=nRB1;nHW<nRB2;nHW++)
	{
		double HW=m_ab1+double(nHW)*step;
		double deltaN=0.0;
		int k=1;
		for(int k=1;(k<nHW) ||((k+nHW)<(n2-1));k++)
		{
			double p1=0.0;
			double p2=0.0;
			if(k<nHW)
			{
				double HW_p=m_ab1+double(nHW-k)*step;
				p1=(3E8*100*hbar/pi)*step*(absorp[nHW-k]-absorp0[nHW-k])/(HW_p*HW_p-HW*HW);
			}
			if((k+nHW)<(n2-1))
			{
				double HW_p=m_ab1+double(nHW+k)*step;
				p2=(3E8*100*hbar/pi)*step*(absorp[nHW+k]-absorp0[nHW+k])/(HW_p*HW_p-HW*HW);
			}
			deltaN+=(p1+p2);
		}
		refract[nHW]=deltaN;
	}

	if(bSaveIndex)
		SaveDeltaRIndex(E,refract, nRB1, nRB2);


}

void CQuantumWell::DeltaRIndexCal_tm(double E, double* absorp, double* refract, bool bSaveIndex)
{
	double step=m_dHWStep;
	double* absorp0=m_pAbsorp0_tm;

	double rb1=m_ab1+(m_ab2-m_ab1)*0.1;
	double rb2=m_ab2-(m_ab2-m_ab1)*0.1;

	int nRB1=int((rb1-m_ab1)/step);
	int nRB2=int((rb2-m_ab1)/step);

	rb1=m_ab1+double(nRB1)*step;
	rb2=m_ab1+double(nRB2)*step;

	int n2=m_nHW;
	
	for (int nHW=nRB1;nHW<nRB2;nHW++)
	{
		double HW=m_ab1+double(nHW)*step;
		double deltaN=0.0;
		int k=1;
		for(int k=1;(k<nHW) ||((k+nHW)<(n2-1));k++)
		{
			double p1=0.0;
			double p2=0.0;
			if(k<nHW)
			{
				double HW_p=m_ab1+double(nHW-k)*step;
				p1=(3E8*100*hbar/pi)*step*(absorp[nHW-k]-absorp0[nHW-k])/(HW_p*HW_p-HW*HW);
			}
			if((k+nHW)<(n2-1))
			{
				double HW_p=m_ab1+double(nHW+k)*step;
				p2=(3E8*100*hbar/pi)*step*(absorp[nHW+k]-absorp0[nHW+k])/(HW_p*HW_p-HW*HW);
			}
			deltaN+=(p1+p2);
		}
		refract[nHW]=deltaN;
	}

	if(bSaveIndex)
		SaveDeltaRIndex_tm(E,refract, nRB1, nRB2);


}

unsigned __stdcall CQuantumWell::ThreadStaticSolveIndex(void * pParam){

	struct THREADPARAM* pObj=(struct THREADPARAM*)pParam;
	CQuantumWell* qw=pObj->pQW;

	CEngLevels* pPrevElectronEng=(pObj->m_pPrevElectronEng);
	CEngLevels* pPrevLHEng=pObj->m_pPrevLHEng;
	CEngLevels* pPrevHHEng=pObj->m_pPrevHHEng;

	ALFAMATRIX* dGuessAlfaHH=&(pObj->m_dGuessAlfaHH);
	ALFAMATRIX* dGuessAlfaLH=&(pObj->m_dGuessAlfaLH);

	double* pAbs=new double[qw->GetnHW()];
	double* pIndex=new double[qw->GetnHW()];

	double* pAbs_tm=new double[qw->GetnHW()];
	double* pIndex_tm=new double[qw->GetnHW()];

		
	double dFieldStep=qw->m_dFieldStep*qw->m_nThreads;

	for(double E=pObj->E;E<=qw->m_Emax;E+=dFieldStep)
	{
		qw->AbsorptionCoefCal(E, pAbs,pAbs_tm, pPrevElectronEng, pPrevLHEng, pPrevHHEng, dGuessAlfaLH, dGuessAlfaHH);
		qw->DeltaRIndexCal(E,pAbs,pIndex);
		qw->DeltaRIndexCal_tm(E,pAbs_tm,pIndex_tm);

	}

	delete[] pAbs;
	delete[] pIndex;

	delete[] pAbs_tm;
	delete[] pIndex_tm;

	_endthreadex( 0 );

	return 0;
}
