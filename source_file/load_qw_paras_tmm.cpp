#include "qwdesign_tmm.h"

using namespace std;

bool CQuantumWell::LoadQWData(char* strQWfile, char* strDirectory)
{
	m_pMEFFe=NULL;
	m_pMEFFlh=NULL;
	m_pMEFFhh=NULL;
	m_pWidth=NULL;
	m_pCV=NULL;
	m_pVV=NULL;
	cptl=NULL;
	vptl=NULL;
	meff=NULL;
	mhhef=NULL;
	mlhef=NULL;
	mhhefp=NULL;
	mlhefp=NULL;
	m_pAbsorp0=NULL;

	m_bFail=0;
	
	strcpy_s(m_strQWFile,strQWfile);
	strcpy_s(m_strDirectory, strDirectory);

	ifstream fenergy(strQWfile);

	if(fenergy.fail())
	{
		cout<<"QW data file does not exist!!!!\n Please check the filename and re-run\n";
		m_bFail=1;
		return 0;
	}

	int i;

	istringstream strDataLine;
	char a[200];

	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_nPoints;
	cout<<"Total number of simulation points is: \n"<<m_nPoints<<endl;

	m_pCV=new double[m_nPoints];
	m_pVV=new double[m_nPoints];
	m_pMEFFe=new double[m_nPoints];
	m_pMEFFlh=new double[m_nPoints];
	m_pMEFFhh=new double[m_nPoints];


	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dQWLength;
	cout<<"total length of the MQW is (A)\n"<<m_dQWLength<<endl;
	//fenergy.seakg(1);

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dBLength;
	cout<<"barrier width and the two boundaies is (A)\n"<<m_dBLength<<endl;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_nLayers;
	cout<<"total number of the MQW is \n"<<m_nLayers<<endl;
	
	cout<<"band structure of the MQW is"<<endl;
	
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);

	
	m_pWidth=new double[m_nLayers];

	for (i=0; i<m_nLayers; i++)		//input the each layer width
	{
		strDataLine>>m_pWidth[i];
		cout<<m_pWidth[i]<<"\t";
	}
	cout<<endl;
	
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	
	cptl=new double[m_nLayers];
	vptl=new double[m_nLayers];

	for (i=0; i<m_nLayers; i++)		//input conduction band potential
	{
		strDataLine>>cptl[i];
		cout<<cptl[i]<<"\t";
	}

	//cptlmax=cptl[0];			//the barrier potential that defines the MQW for electrons is always the first barrier potential, whether there is E or not
	cout<<endl;					//actually, the cptmax, which is the max energy level, really depends on the electric field, can be either two sides

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);

	for (i=0; i<m_nLayers; i++)		//input valence band potential
	{
		strDataLine>>vptl[i];
		cout<<vptl[i]<<"\t";
	}
	
	cout<<endl;					//for commmon valence band structure, vptlmax=vptl[layer-1], the last one,
								//but for the engineered one, it is reversed,  vptlmax=vptl[0] the first one, the same as the conduction band
	
	
	//input the conduction band effective mass
	cout<<"electron effective mass:";
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);

	meff=new double[m_nLayers];
	for (i=0; i<m_nLayers; i++)
	{
		strDataLine>>meff[i];
		meff[i]=meff[i]*9.11E-31;
		cout<<meff[i]<<"     ";
	}
	
	//input the valence band heavy hole effective mass perpendicular to the layer
	cout<<endl<<"heavy hole effective mass perpendicular to the layer: \n";

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);

	mhhef=new double[m_nLayers];
	for (i=0; i<m_nLayers; i++)
	{
		strDataLine>>mhhef[i];
		mhhef[i]=mhhef[i]*9.11E-31;
		cout<<mhhef[i]<<"     ";
	}

	//input the valence band light hole effective mass perpendicular to the layer
	cout<<endl<<"light hole effective mass perpendicular to the layer:\n";
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);

	mlhef=new double[m_nLayers];
	for (i=0; i<m_nLayers; i++)
	{
		strDataLine>>mlhef[i];
		mlhef[i]=mlhef[i]*9.11E-31;
		cout<<mlhef[i]<<"     ";
	}

	cout<<endl<<"heavy hole effective mass parallel to the layer:\n";
	
	mhhefp=new double[m_nLayers]; 
	for (i=0; i<m_nLayers; i++)
	{
		mhhefp[i]=4*mhhef[i]*mlhef[i]/(3*mhhef[i]+mlhef[i]);
		cout<<mhhefp[i]<<"     ";
	}

	cout<<endl<<"light hole effective mass parallel to the layer:\n";

	mlhefp=new double[m_nLayers]; 
	for (i=0; i<m_nLayers; i++)
	{
		mlhefp[i]=4*mhhef[i]*mlhef[i]/(3*mlhef[i]+mhhef[i]);
		cout<<mlhefp[i]<<"     ";
	}

	m_dRedHp=meff[1]*mhhefp[1]/(meff[1]+mhhefp[1]); // yifei raise alarm!!!!!!!
	m_dRedLp=meff[1]*mlhefp[1]/(meff[1]+mlhefp[1]); // questionable code

	// reduced mass perpendicular to the qw (added by Yifei)
	m_dRedH=meff[0]*mhhef[0]/(meff[0]+mhhef[0]); 
	m_dRedL=meff[0]*mlhef[0]/(meff[0]+mlhef[0]); 


	// input Eg 
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dEg;
	cout<<"\ngap between minimum conduction band and maximum valence band:"<<m_dEg<<endl;

	//input hw
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_ab1>>m_ab2>>m_nHW;

	cout<<"Photon energy interested: \nfrom "<<m_ab1<<" to "<<m_ab2<<"with "<<m_nHW<<" points"<<endl;
	m_dHWStep=(m_ab2-m_ab1)/double(m_nHW);
	rstep=1.0/m_dHWStep;

	m_pAbsorp0=NULL;		//absorption coefficient without electric field

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dFKEg;

	//input nr
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_nr;
	cout<<"refractive index of the material in common cases is:"<<m_nr<<endl;

	//input dielect
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dielect;
	cout<<"dielectric constant of the material in common cases is:"<<m_dielect<<endl;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_Eref;
	cout<<"The reference E field is:"<<m_Eref<<endl;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_Emin>>m_Emax>>m_nFields;
	cout<<"The E field range is:"<<m_Emin<<" "<<m_Emax<<" "<<m_nFields<<endl;
	m_dFieldStep=(m_Emax-m_Emin)/double(m_nFields);
/*different components*/
	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	
	m_Ratio=new double[3];

	for (i=0; i<3; i++)		//input the each layer width
	{
		strDataLine>>m_Ratio[i];
	}
	cout<<"The ratio of components is:"<<m_Ratio[0]<<" "<<m_Ratio[1]<<" "<<m_Ratio[2]<<endl;
/*-----------------------------------------*/

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_HWHM_G>>m_HWHM_L;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_BufferLength;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_nThreads;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dMinEEnergy>>m_dMaxEEnergy;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dLEMin>>m_dLEMax;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dMinHEnergy>>m_dMaxHEnergy;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dLHMin>>m_dLHMax;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dMinHHEnergy>>m_dMaxHHEnergy;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dHHMin>>m_dHHMax;


	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_nLineWidthFunc;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dELHBroadening;

	fenergy.getline(a, sizeof(a)); 
	fenergy.getline(a, sizeof(a));
	fenergy.getline(a, sizeof(a));
	strDataLine.clear();
	strDataLine.str(a);
	strDataLine>>m_dEHHBroadening;

	fenergy.close();

	// load mess, potential, and ....
	int j=0;
	i=0;
	int k=0;

	m_dStepSize=m_dQWLength/m_nPoints;

	for (j=0; (j<m_nPoints) && (i<m_nLayers); j++){
	
		m_pCV[j]=cptl[i];
		m_pVV[j]=vptl[i];
		m_pMEFFe[j]=meff[i];
		m_pMEFFlh[j]=mlhef[i];
		m_pMEFFhh[j]=mhhef[i];

		if((j-k)>m_pWidth[i]/m_dStepSize)
		{
			k=j;
			i++;
		}
	}
	return 1;
		
}
CQuantumWell::~CQuantumWell(){

	if(m_pCV!=NULL)
		::delete[] m_pCV;
	if(m_pVV!=NULL)
		::delete[] m_pVV;
	if(m_pMEFFe!=NULL)
		::delete[] m_pMEFFe;
	if(m_pMEFFlh!=NULL)
		::delete[] m_pMEFFlh;
	if(m_pMEFFhh!=NULL)
		::delete[] m_pMEFFhh;
	if(m_pWidth!=NULL)
		::delete[] m_pWidth;
	if(cptl!=NULL)
		::delete[] cptl;
	if(vptl!=NULL)
		::delete[] vptl;
	if(meff!=NULL)
		::delete[] meff;
	if(mhhef!=NULL)
		::delete[] mhhef;
	if(mlhef!=NULL)
		::delete[] mlhef;
	if(mhhefp!=NULL)
		::delete[] mhhefp;
	if(mlhefp!=NULL)
		::delete[] mlhefp;
	if(m_pAbsorp0!=NULL)
		::delete[] m_pAbsorp0;
	if(m_pAbsorp0!=NULL)
		::delete[] m_pAbsorp0_tm;

}

