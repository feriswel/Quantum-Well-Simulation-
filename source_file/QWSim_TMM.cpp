#include "qwdesign_tmm.h"

using namespace std;
/************************function*******************************************/
int main(int nParam, char* Params[])
{

	char strQWfile[200];
	char strDirectory[200];

	if(nParam==1)
	{
		cout<<"Insufficient parameters\n";
		cout<<"Usage:\n qw_cal qw_layer_file [output_director]\n";
		cout<<"Please try again\n";
		return 0;
	}
	if(nParam==2)
		strcpy_s(strQWfile,Params[1]);

	if(nParam==3)
		strcpy_s(strDirectory,Params[2]);
	else
		strcpy_s(strDirectory,".//output//");

	CQuantumWell* qw=new CQuantumWell(strQWfile, strDirectory);

	if(qw->fail()){
		cout<<"falied to load QW parameters\n";
		return 0;
	}

	double E;		//electric field in unit of V/cm
	E=0.0; // Yifei's modifications

	qw->RefAbsorptionCoefCal();
	qw->SaveHWPoints();


	//double* pAbs=new double[qw->GetnHW()];
	//qw->AbsorptionCoefCal(E, pAbs);
	//double* pIndex=new double[qw->GetnHW()];
	//qw->DeltaRIndexCal(E,pAbs,pIndex);
	//delete[] pAbs;
	//delete[] pIndex;

	HANDLE* hProcess=new HANDLE [qw->GetnThreads()];
	struct THREADPARAM* params=new struct THREADPARAM[qw->GetnThreads()];

	for(int i=0;i<qw->GetnThreads();i++){
		
		params[i].E=qw->GetEMin()+i*qw->GetEFieldStep();
		params[i].pQW=qw;
		
		params[i].m_pPrevElectronEng=new CEngLevels(qw->m_pElectronEng);
		params[i].m_pPrevLHEng=new CEngLevels(qw->m_pLHEng);
	//	params[i].m_pPrevHHEng=new CEngLevels(qw->m_pHHEng);
		
		for(int k=0;k<MAXLEVEL;k++)
			for(int j=0;j<MAXLEVEL;j++)
			{
				params[i].m_dGuessAlfaHH.a[k][j]=qw->m_dGuessAlfaHH.a[k][j];
				params[i].m_dGuessAlfaLH.a[k][j]=qw->m_dGuessAlfaLH.a[k][j];
			}


		unsigned ID;

		hProcess[i] = (HANDLE)_beginthreadex( NULL, // security
                      0,             // stack size
                      CQuantumWell::ThreadStaticSolveIndex,// entry-point-function
                      &(params[i]),           // arg list holding the "this" pointer
                      0, // so we can later call ResumeThread()
                      &ID );
	}

	for(int i=0;i<qw->GetnThreads();i++){
		WaitForSingleObject( hProcess[i], INFINITE );
		CloseHandle( hProcess[i]);
	}

	delete[] hProcess;

	// release memory hold by pPrevEngs...

	for(int i=0;i<qw->GetnThreads();i++){
		
		
		delete params[i].m_pPrevElectronEng;
//		delete params[i].m_pPrevHHEng;
		delete params[i].m_pPrevLHEng;
	}
		

	delete[] params;

	delete qw;
	return 1;
	
}

