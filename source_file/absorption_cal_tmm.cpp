#include "qwdesign_tmm.h"

using namespace std;

/********************************dielectric constant calcualtion****************************/
void CQuantumWell::AbsorptionCoefCal(double E, double* absorp, double* absorp_tm, CEngLevels* pElectronLevels, CEngLevels* pLHLevels, CEngLevels* pHHLevels, ALFAMATRIX* pGuessAlfaLH, ALFAMATRIX* pGuessAlfaHH, bool bSaveStates, bool bSaveAbsorption)
{
	
	// calculate electron, and hole wavefunctions
	CTotalBoundedStates eStates(m_nPoints);
	CTotalBoundedStates lhStates(m_nPoints);
	CTotalBoundedStates hhStates(m_nPoints);


	FindElectronWaves(E,0.0005,&eStates, pElectronLevels);
	FindLightHoleWaves(E,0.0005,&lhStates, pLHLevels);
	//FindHeavyHoleWaves(E,0.00000001,&hhStates, pHHLevels);
	FindHeavyHoleWavesEigen(E,&hhStates);

	int nEState=eStates.NumberOfStates();
	int nLHState=lhStates.NumberOfStates();
	int nHHState=hhStates.NumberOfStates();


	int i, j;
	int hw;

	double HW;

	int n2=m_nHW;

	double* im_diel_band=new double[m_nHW];
	double* im_diel_exci=new double[m_nHW];

	double* absorp_bulk=new double[m_nHW];

	double* im_diel_band_tm=new double[m_nHW];
	double* im_diel_exci_tm=new double[m_nHW];

	

/****************************************band to band transition********************************************/

	//electron-light hole band to band transition

	double Mb_squared=me0*me0*m_dEg*q*(m_dEg+Eso)/(12*meff[1]*(m_dEg+2*Eso/3));
	Mb_squared*=2.0;
	
	for (hw=0; hw<n2; hw++){
		im_diel_band[hw]=0.0;
		im_diel_exci[hw]=0.0;

		im_diel_band_tm[hw]=0.0;
		im_diel_exci_tm[hw]=0.0;

	}

//		Electron - light hole interaction
	double coef=q*q*(m_dRedLp)*Mb_squared/(epsilon0*me0*me0*q*q*(m_dQWLength-m_dBLength)*1E-10);
	double coef_ex=2*Mb_squared*pi*q*q/(epsilon0*me0*me0*(1.0)/hbar*(1.0)/hbar*(m_dQWLength-m_dBLength)*1E-10);
	
	for (i=0; i<nEState; i++)
	{
		for (j=0; j<nLHState; j++)
		{
			double dLineWidth=m_HWHM_L; //linewidth due to phono scattering
			double dw1=eStates.m_pStates[i]->GetLineWidth();
			double dw2=lhStates.m_pStates[j]->GetLineWidth();

			dLineWidth=sqrt(dLineWidth*dLineWidth+dw1*dw1+dw2*dw2);
			double g_linewidth=m_HWHM_G+m_dELHBroadening*E/100e3;

			CLineShape* ls=NULL;
			if(m_nLineWidthFunc==2)
				ls=new CLineShape(g_linewidth,0.1,0.4,0.8,2000,2000,2000,2);
			else if(m_nLineWidthFunc==3)
				ls=new CLineShape(dLineWidth,g_linewidth, 0.1,0.4,0.8,2000,2000,2000);
			else
				ls=new CLineShape(dLineWidth,0.1,0.4,0.8,2000,2000,2000);

			double dEEnergy=eStates.m_pStates[i]->eig();
			double dLHEnergy=lhStates.m_pStates[j]->eig();

			CCOMPLEX* pCplxEWaves=eStates.m_pStates[i]->wave();
			CCOMPLEX* pCplxLHWaves=lhStates.m_pStates[j]->wave();

			CCOMPLEX cplexIntegWave=CCOMPLEX(0.0,0.0);

			for(int z=0;z<m_nPoints;z++)
				cplexIntegWave+=pCplxEWaves[z]*(pCplxLHWaves[z]);  // Yifei comment: integration does not consider the step size!!!!

			double dIntegWaveSq=norm(cplexIntegWave);

			for (hw=0, HW=m_ab1; hw<n2; hw++, HW+=m_dHWStep)
			{
					double dPE=HW-(m_dEg+dEEnergy+dLHEnergy); // energy difference between b-b and photon in eV
					double band_ij_contrib=coef*dIntegWaveSq*ls->intg_F_ev_to_inf(dPE)*Mel/(HW*HW);	
					im_diel_band[hw]+=band_ij_contrib;
			}

			// next we calculate the exciton bounding energy 
			double maxAlfa;
			double maxf;

//			SeekExcitonBondEng_e_lh(1e8, eStates.m_pStates[i], lhStates.m_pStates[j], maxAlfa, maxf );

			if(m_dGuessAlfaLH.a[i][j]<0.001)
			{
				SeekExcitonBondEng_e_lh(1e8, eStates.m_pStates[i], lhStates.m_pStates[j], maxAlfa, maxf );
				m_dGuessAlfaLH.a[i][j]=maxAlfa;
			}
			else
			{
				SeekExcitonBondEng_e_lh(pGuessAlfaLH->a[i][j], eStates.m_pStates[i], lhStates.m_pStates[j], maxAlfa, maxf );
				pGuessAlfaLH->a[i][j]=maxAlfa;
			}
			
			double exbd_min=maxf;
			double phiex0_sq=2.0*maxAlfa*maxAlfa/PI;

			//double coef_ex=2*Mb_squared*pi*q*q/(epsilon0*me0*me0*(ab1+hw*step)/hbar*(ab1+hw*step)/hbar*(length-blength)*1E-10);
			

			for (hw=0, HW=m_ab1; hw<n2; hw++, HW+=m_dHWStep)
			{
				double dPE;
				dPE=HW-(m_dEg+dEEnergy+dLHEnergy+exbd_min); // energy difference between b-b and photon in eV
				im_diel_exci[hw]+=coef_ex*phiex0_sq*dIntegWaveSq*Mel*ls->F_ev(dPE)/q/(HW*HW);
			}

			delete ls;

		}
	}
	
	for (hw=0; hw<n2; hw++){

		im_diel_band_tm[hw]=im_diel_band[hw]*4;
		im_diel_exci_tm[hw]=im_diel_exci[hw]*4;

	}


	// electro-heavy hole interaction

	coef=q*q*(m_dRedHp)*Mb_squared/(epsilon0*me0*me0*q*q*(m_dQWLength-m_dBLength)*1E-10);
	coef_ex=2*Mb_squared*pi*q*q/(epsilon0*me0*me0*(1.0)/hbar*(1.0)/hbar*(m_dQWLength-m_dBLength)*1E-10);

	for (i=0; i<nEState; i++)
	{
		for (j=0; j<nHHState; j++)
		{
			double dLineWidth=m_HWHM_L; //linewidth due to phono scattering
			double dw1=eStates.m_pStates[i]->GetLineWidth();
			double dw2=hhStates.m_pStates[j]->GetLineWidth();

			dLineWidth=sqrt(dLineWidth*dLineWidth+dw1*dw1+dw2*dw2);
			double g_linewidth=m_HWHM_G+0.02*m_dEHHBroadening/100e3;

			CLineShape* ls=NULL;
			if(m_nLineWidthFunc==2) // modified lorenzian function
				ls=new CLineShape(g_linewidth,0.1,0.4,0.8,2000,2000,2000,2);
			else if(m_nLineWidthFunc==3)
				ls=new CLineShape(dLineWidth,g_linewidth, 0.1,0.4,0.8,2000,2000,2000);
			else
				ls=new CLineShape(dLineWidth,0.1,0.4,0.8,2000,2000,2000);


			double dEEnergy=eStates.m_pStates[i]->eig();
			double dHHEnergy=hhStates.m_pStates[j]->eig();

			CCOMPLEX* pCplxEWaves=eStates.m_pStates[i]->wave();
			CCOMPLEX* pCplxHHWaves=hhStates.m_pStates[j]->wave();

			CCOMPLEX cplexIntegWave=CCOMPLEX(0.0,0.0);

			for(int z=0;z<m_nPoints;z++)
				cplexIntegWave+=pCplxEWaves[z]*(pCplxHHWaves[z]);  // Yifei comment: integration does not consider the step size!!!!

			double dIntegWaveSq=norm(cplexIntegWave);

			for (hw=0, HW=m_ab1; hw<n2; hw++, HW+=m_dHWStep)
			{
					double dPE=HW-(m_dEg+dEEnergy+dHHEnergy); // energy difference between b-b and photon in eV
					double band_ij_contrib=coef*dIntegWaveSq*ls->intg_F_ev_to_inf(dPE)*Meh/(HW*HW);	
					im_diel_band[hw]+=band_ij_contrib;
			}

			// next we calculate the exciton bounding energy 
			double maxAlfa;
			double maxf;

//			SeekExcitonBondEng_e_hh(1e8, eStates.m_pStates[i], hhStates.m_pStates[j], maxAlfa, maxf );
			if(m_dGuessAlfaHH.a[i][j]<0.001)
			{
				SeekExcitonBondEng_e_hh(1e8, eStates.m_pStates[i], hhStates.m_pStates[j], maxAlfa, maxf );
				m_dGuessAlfaHH.a[i][j]=maxAlfa;
			}
			else
			{
				SeekExcitonBondEng_e_hh(pGuessAlfaHH->a[i][j], eStates.m_pStates[i], hhStates.m_pStates[j], maxAlfa, maxf );
				pGuessAlfaHH->a[i][j]=maxAlfa;
			}
			
			double exbd_min=maxf;
			double phiex0_sq=2.0*maxAlfa*maxAlfa/PI;

			//double coef_ex=2*Mb_squared*pi*q*q/(epsilon0*me0*me0*(ab1+hw*step)/hbar*(ab1+hw*step)/hbar*(length-blength)*1E-10);
			

			for (hw=0, HW=m_ab1; hw<n2; hw++, HW+=m_dHWStep)
			{
				double dPE;
				dPE=HW-(m_dEg+dEEnergy+dHHEnergy+exbd_min); // energy difference between b-b and photon in eV
				im_diel_exci[hw]+=coef_ex*phiex0_sq*dIntegWaveSq*Meh*ls->F_ev(dPE)/q/(HW*HW);
			}

			delete ls;

		}
	}
	

	for (hw=0; hw<n2; hw++)				//the absorption coefficient
		im_diel_band[hw]=im_diel_band[hw]*(m_ab1+hw*m_dHWStep)/(3E8*m_nr*hbar)/100;			//100 is an coefficient to convert the unit from 1/m to 1/cm

	for (hw=0; hw<n2; hw++)				//the absorption coefficient		the unit is 1/cm
		im_diel_exci[hw]=im_diel_exci[hw]*(m_ab1+hw*m_dHWStep)/(3E8*m_nr*hbar)/100;

	
	// calculate light hole electron interaction for TM polarizations
	for (hw=0; hw<n2; hw++)				//the absorption coefficient
		im_diel_band_tm[hw]=im_diel_band_tm[hw]*(m_ab1+hw*m_dHWStep)/(3E8*m_nr*hbar)/100;			//100 is an coefficient to convert the unit from 1/m to 1/cm

	for (hw=0; hw<n2; hw++)				//the absorption coefficient		the unit is 1/cm
		im_diel_exci_tm[hw]=im_diel_exci_tm[hw]*(m_ab1+hw*m_dHWStep)/(3E8*m_nr*hbar)/100;


	// ************** calculate the bulk absorption coefficient

	for (hw=0, HW=m_ab1; hw<n2; hw++, HW+=m_dHWStep)
	{
		absorp_bulk[hw]=F_K_absorption(HW,E,m_dFKEg,m_dRedH,m_nr)+F_K_absorption(HW,E,m_dFKEg,m_dRedL,m_nr);
	}


//****************************total absorption coefficient****************************


	for (hw=0; hw<n2; hw++)
			absorp[hw]=m_Ratio[0]*im_diel_exci[hw]+m_Ratio[1]*im_diel_band[hw]+m_Ratio[2]*absorp_bulk[hw];

	for (hw=0; hw<n2; hw++)
			absorp_tm[hw]=m_Ratio[0]*im_diel_exci_tm[hw]+m_Ratio[1]*im_diel_band_tm[hw]+m_Ratio[2]*absorp_bulk[hw];

	// Save waves, absoprtion, .... if required

	if(bSaveStates){
		eStates.SaveStates(E, m_strDirectory,0);
		lhStates.SaveStates(E, m_strDirectory,1);
		hhStates.SaveStates(E, m_strDirectory,2);

	}

	if(bSaveAbsorption){
		SaveAbsorptionCoefs(E, im_diel_band, im_diel_exci, absorp);
		SaveAbsorptionCoefs_tm(E, im_diel_band_tm, im_diel_exci_tm, absorp_tm);

	}

	delete im_diel_band;
	delete im_diel_exci;
	delete im_diel_band_tm;
	delete im_diel_exci_tm;
	
	delete absorp_bulk;


}


void CQuantumWell::RefAbsorptionCoefCal(bool bSaveStates, bool bSaveAbsorption){
	if(m_pAbsorp0!=NULL)
		delete m_pAbsorp0;
	m_pAbsorp0=new double[m_nHW];
	m_pAbsorp0_tm=new double[m_nHW];

	AbsorptionCoefCal(m_Eref,m_pAbsorp0, m_pAbsorp0_tm);
}

