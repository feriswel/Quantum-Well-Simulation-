#include "qwdesign_tmm.h"

using namespace std;



void CTotalBoundedStates::SaveStates(double E, char* dir,int nStateType)
{

	char filename[200];

	char index[20];
	int field=int(E);
	_itoa_s(field,index,20,10);
	
	ofstream ofstr;

	strcpy_s(filename,dir);

	if(nStateType==0)//electron state
		strcat_s(filename,"EWaves_");
	else
		if(nStateType==1)
			strcat_s(filename,"LHWaves_");
		else
			strcat_s(filename,"HHWaves_");


	strcat_s(filename,index);
	strcat_s(filename,".txt");

	ofstr.open(filename);

	if(ofstr.fail())
	{
		cout<<"failed to open file  "<<filename<<"\n";
		return;
	}

	int i;

	if(nStateType==0)//electron state
		ofstr<<"Electron states under E(kV/CM) -"<<E/1000<<endl;
	
	if(nStateType==1)//light hole state
		ofstr<<"Hole states under E(kV/CM) -"<<E/1000<<endl;



	ofstr<<"Eigen values (eV):  ";

	for(i=0;i<m_nTotalStates;i++)
	{
		ofstr<<m_pStates[i]->eig()<<" ";
	}
	ofstr<<endl;


	ofstr<<"Wavefunctions:  "<<endl;


	for(i=0;i<m_nPoints;i++)
	{
		for(int j=0;j<m_nTotalStates;j++)
		{
			CCOMPLEX* pWave=m_pStates[j]->wave();
			ofstr<<pWave[i].real()<<" "<<pWave[i].imag()<<" "<<abs(pWave[i])<<" ";
		}
		ofstr<<"\n";
	}
	ofstr.close();


}

void CQuantumWell::SaveAbsorptionCoefs(double E, double* im_diel_band, double* im_diel_exci, double* absorp)
{

	int n2=m_nHW;
	int hw;

	char filename[200];

	char index[20];
	int field=int(E);
	_itoa_s(field,index,20,10);
	
	ofstream ofstr;

	if(1){
		// save band to band transition

		strcpy_s(filename,m_strDirectory);
		strcat_s(filename,"Absp_BB_at_");
		strcat_s(filename,index);
		strcat_s(filename,".txt");
		
		ofstr.open(filename);
		ofstr<<"Absorption(Band to Band) at "<<E/1000<<"kV/CM:\n";

		for (hw=0; hw<n2; hw++)
		{
			ofstr<<im_diel_band[hw]<<"\t"<<endl;
		}

		ofstr.close();

		// save excition transition

		strcpy_s(filename,m_strDirectory);
		strcat_s(filename,"Absp_Exctn_at_");
		strcat_s(filename,index);
		strcat_s(filename,".txt");
		
		ofstr.open(filename);
		ofstr<<"Absorption(exciton) at "<<E/1000<<"kV/CM:\n";

		for (hw=0; hw<n2; hw++)
		{
			ofstr<<im_diel_exci[hw]<<"\t"<<endl;
		}

		ofstr.close();

		// save total absorption

		strcpy_s(filename,m_strDirectory);
		strcat_s(filename,"Absp_Totl_at_");
		strcat_s(filename,index);
		strcat_s(filename,".txt");
		
		ofstr.open(filename);
		ofstr<<"Absorption(total) at "<<E/1000<<"kV/CM:\n";

		for (hw=0; hw<n2; hw++)
		{
			ofstr<<absorp[hw]<<"\t"<<endl;
		}

		ofstr.close();
	}

}


void CQuantumWell::SaveAbsorptionCoefs_tm(double E, double* im_diel_band, double* im_diel_exci, double* absorp)
{

	int n2=m_nHW;
	int hw;

	char filename[200];

	char index[20];
	int field=int(E);
	_itoa_s(field,index,20,10);
	
	ofstream ofstr;

	if(1){
		// save band to band transition

		strcpy_s(filename,m_strDirectory);
		strcat_s(filename,"TM_Absp_BB_at_");
		strcat_s(filename,index);
		strcat_s(filename,".txt");
		
		ofstr.open(filename);
		ofstr<<"Absorption(Band to Band) at "<<E/1000<<"kV/CM:\n";

		for (hw=0; hw<n2; hw++)
		{
			ofstr<<im_diel_band[hw]<<"\t"<<endl;
		}

		ofstr.close();

		// save excition transition

		strcpy_s(filename,m_strDirectory);
		strcat_s(filename,"TM_Absp_Exctn_at_");
		strcat_s(filename,index);
		strcat_s(filename,".txt");
		
		ofstr.open(filename);
		ofstr<<"Absorption(exciton) at "<<E/1000<<"kV/CM:\n";

		for (hw=0; hw<n2; hw++)
		{
			ofstr<<im_diel_exci[hw]<<"\t"<<endl;
		}

		ofstr.close();

		// save total absorption

		strcpy_s(filename,m_strDirectory);
		strcat_s(filename,"TM_Absp_Totl_at_");
		strcat_s(filename,index);
		strcat_s(filename,".txt");
		
		ofstr.open(filename);
		ofstr<<"Absorption(total) at "<<E/1000<<"kV/CM:\n";

		for (hw=0; hw<n2; hw++)
		{
			ofstr<<absorp[hw]<<"\t"<<endl;
		}

		ofstr.close();
	}

}
void CQuantumWell::SaveDeltaRIndex(double E, double* refract, int nRB1, int nRB2)
{

	int n2=m_nHW;

	char filename[200];

	char index[20];
	int field=int(E);
	_itoa_s(field,index,20,10);
	
	ofstream ofstr;

	if(1){
		// save band to band transition

		strcpy_s(filename,m_strDirectory);
		strcat_s(filename,"IndexChange_");
		strcat_s(filename,index);
		strcat_s(filename,".txt");
		
		ofstr.open(filename);
		ofstr<<"Index Change at "<<E/1000<<"kV/CM:\n";

		for (int nHW=nRB1;nHW<nRB2;nHW++)
		{
			ofstr<<refract[nHW]<<endl;
		}

		ofstr.close();

	}

}

void CQuantumWell::SaveDeltaRIndex_tm(double E, double* refract, int nRB1, int nRB2)
{

	int n2=m_nHW;

	char filename[200];

	char index[20];
	int field=int(E);
	_itoa_s(field,index,20,10);
	
	ofstream ofstr;

	if(1){
		// save band to band transition

		strcpy_s(filename,m_strDirectory);
		strcat_s(filename,"TM_IndexChange_");
		strcat_s(filename,index);
		strcat_s(filename,".txt");
		
		ofstr.open(filename);
		ofstr<<"Index Change at "<<E/1000<<"kV/CM:\n";

		for (int nHW=nRB1;nHW<nRB2;nHW++)
		{
			ofstr<<refract[nHW]<<endl;
		}

		ofstr.close();

	}

}

void CQuantumWell::SaveHWPoints(){

	int n2=m_nHW;
	double step=m_dHWStep;

	char filename[200];

	strcpy_s(filename,m_strDirectory);
	strcat_s(filename,"PEs_Absp.txt");

	ofstream ofstr;
	ofstr.open(filename);

	ofstr<<"Photon energies for absorption calculation\n";
	for(int hw=0; hw<n2;hw++)
	{
		ofstr<<double(m_ab1+hw*step)<<"\n";
	}

	ofstr.close();


	double rb1=m_ab1+(m_ab2-m_ab1)*0.1;
	double rb2=m_ab2-(m_ab2-m_ab1)*0.1;

	int nRB1=int((rb1-m_ab1)/step);
	int nRB2=int((rb2-m_ab1)/step);


	strcpy_s(filename,m_strDirectory);
	strcat_s(filename,"PEs_Index.txt");

	ofstr.open(filename);
	//output wavelength(ev) for delta ref index calculation
//	ofstream fstrPEDeltaIndex("PEDeltaIndex.txt");

	ofstr<<"Photon energies for refractive calculation\n";
	for (int nHW=nRB1;nHW<nRB2;nHW++)
	{
		ofstr<<double(m_ab1+double(nHW)*step)<<"\n";
	}
	ofstr.close();
}

