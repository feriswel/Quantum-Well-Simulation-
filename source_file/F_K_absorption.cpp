#include "airyf.h"
#include "qwdesign_tmm.h"

//const double me0=me0=9.11E-31;

double F_K_absorption(double HW, double F, double Eg, double m_reduce, double n)
{
	//return 0.0;
	if(fabs(F)<1)
		F=10; //A_j*sqrt(B_j*(HW-Eg))/pi;


	double temp=pow(2.0*m_reduce/me0,1.0/3.0);

	double B_j=1.1e5*temp;;
	double A_j=3.55e4*temp*temp*temp*temp/n/HW;

	double F_3=pow(F,1.0/3.0);


	double beta_j=B_j*(Eg-HW)/(F_3*F_3); 

	double ai;
	double aip;
	double bi;
	double bip;

	airy(beta_j,ai,aip,bi,bip);

	return A_j*F_3*((aip*aip)-beta_j*ai*ai);

}