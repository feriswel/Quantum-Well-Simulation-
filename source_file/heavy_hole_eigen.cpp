#include "qwdesign_tmm.h"

using namespace std;


void CQuantumWell::FindHeavyHoleWavesEigen(double E, CTotalBoundedStates* pTotalStates)
{
	double* potential=new double[m_nPoints];
//	CCOMPLEX* wave=new CCOMPLEX[m_nPoints];

	double dVMax=E*m_dQWLength*1E-8;
	double dStepSize_sq=m_dStepSize*m_dStepSize*1e-20;
	//double dEngMax=m_dMaxHEnergy-E*m_dLHMin*1e-8;
	
	double dEngMax=m_dMaxHHEnergy-E*m_dHHMax*1e-8;
	double dEngMin=m_dMinHHEnergy-E*m_dHHMin*1e-8;
	int i;
	int j;
	
	SolveEigen* hh_solve_eigen=new SolveEigen(m_nPoints);

	//hh_solve_eigen->n=500;


	for(i=0;i<m_nPoints;i++)
	{
		if(i*m_dStepSize<m_dBLength || i*m_dStepSize>(m_dQWLength-m_dBLength)) // code added for V9 for update potential shape
			potential[i]=0.5;
		else
			potential[i]=m_pVV[i]+dVMax/2-(dVMax/m_nPoints)*(i);		//initialize potential profile
	}


	hh_solve_eigen->d[0]=0;
	for (i=1;i<m_nPoints;i++)
	{
		hh_solve_eigen->d[i]=potential[i]+q*h_bar_sq/(m_pMEFFhh[i]*dStepSize_sq);	//initialize diagnoal elements
	}

	hh_solve_eigen->e[0]=hh_solve_eigen->e[1]=0;
	for (i=2; i<m_nPoints; i++)
	{
		hh_solve_eigen->e[i]=-q*h_bar_sq/(2*m_pMEFFhh[i]*dStepSize_sq);			//initialize super diagnoal elements
	}
	
	for (i=0; i<m_nPoints; i++)											//initialize wave function matrix in the form of identity matrix
	{
		for (j=0; j<m_nPoints; j++)
		{
			if (i==j)
				hh_solve_eigen->z[i][j]=1;
			else
				hh_solve_eigen->z[i][j]=0;
		}
	}
	
	hh_solve_eigen->tqli();         //solve eigenvalues and wave functions
	hh_solve_eigen->eigsrt();		//sort eigenvalues and wave functions
	
	CCOMPLEX* wave;

	int n=m_nPoints-1;

	if(m_nNumberOfHHState==0)
	{
	for(j=0;(hh_solve_eigen->d[n-j])<dEngMax;j++)     //store states with eigenvalue less than dEngMax(real bounded states)
	{
		wave=new CCOMPLEX[m_nPoints];

		if(hh_solve_eigen->z[n/2][n-j]>0)    //postprocess wave function
		{
			for(i=0;i<n+1;i++)
			{
				wave[i]=CCOMPLEX(hh_solve_eigen->z[i][n-j],0);
			}
		}
		else if(hh_solve_eigen->z[n/2][n-j]<0)
		{
			for(i=0;i<n+1;i++)
			{
				wave[i]=CCOMPLEX(-hh_solve_eigen->z[i][n-j],0);
			}
		}			//store states (eigenvalue, wave function, number of points in space, linewidth=0.00001)
		pTotalStates->AddState(new CState(hh_solve_eigen->d[n-j],wave,m_nPoints,0.00001));    ///0.00001 linewidth???arbitrary given?????
	}

	m_nNumberOfHHState=pTotalStates->NumberOfStates();
	}

	else 
	{
	for(j=0;j<m_nNumberOfHHState;j++)     //store states with eigenvalue less than dEngMax(real bounded states)
	{
		 wave=new CCOMPLEX[m_nPoints];

		if(hh_solve_eigen->z[n/2][n-j]>0)    //postprocess wave function
		{
			for(i=0;i<n+1;i++)
			{
				wave[i]=CCOMPLEX(hh_solve_eigen->z[i][n-j],0);
			}
		}
		else if(hh_solve_eigen->z[n/2][n-j]<0)
		{
			for(i=0;i<n+1;i++)
			{
				wave[i]=CCOMPLEX(-hh_solve_eigen->z[i][n-j],0);
			}
		}			//store states (eigenvalue, wave function, number of points in space, linewidth=0.00001)
		pTotalStates->AddState(new CState(hh_solve_eigen->d[n-j],wave,m_nPoints,0.00001));    ///0.00001 linewidth???arbitrary given?????
	}

	}

	//if(wave!=NULL)
	//	delete [] wave;
	//if(diag!=NULL)
	//	delete [] diag;
	//if(sup_diag!=NULL)
	//	delete [] sup_diag;
	if(potential!=NULL)
		delete [] potential;
	//if(wave_matrix!=NULL)
	//	delete wave_matrix;
	if(hh_solve_eigen!=NULL)
		delete hh_solve_eigen;
}




/**************************  Tridiagonal QL algorithm -- Implicit  **********************/
void SolveEigen::tqli()

{
int m, l, iter, i, k;
double s, r, p, g, f, dd, c, b;
//void erhand();

for (i = 2; i <= n; i++)
    e[i-1] = e[i];
e[n] = 0.0;
for (l = 1; l <= n; l++)
    {
    iter = 0;
    do
      {
      for (m = l; m <= n-1; m++)
          {
          dd = fabs(d[m]) + fabs(d[m+1]);
          if (fabs(e[m]) + dd == dd) break;
          }
          if (m != l)
             {
             if (iter++ == 30) printf("No convergence in TLQI.");
             g = (d[l+1] - d[l]) / (2.0 * e[l]);
             r = sqrt((g * g) + 1.0);
             g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
             s = c = 1.0;
             p = 0.0;
             for (i = m-1; i >= l; i--)
                 {
                 f = s * e[i];
                 b = c * e[i];


                 if (fabs(f) >= fabs(g))
                    {
                    c = g / f;
                    r = sqrt((c * c) + 1.0);
                    e[i+1] = f * r;
                    c *= (s = 1.0/r);
                    }
                 else
                    {
                    s = f / g;
                    r = sqrt((s * s) + 1.0);
                    e[i+1] = g * r;
                    s *= (c = 1.0/r);
                    }
                 g = d[i+1] - p;
                 r = (d[i] - g) * s + 2.0 * c * b;
                 p = s * r;
                 d[i+1] = g + p;
                 g = c * r - b;


                 for (k = 1; k <= n; k++)
                     {
                     f = z[k][i+1];
                     z[k][i+1] = s * z[k][i] + c * f;
                     z[k][i] = c * z[k][i] - s * f;
                     }
                 }
                 d[l] = d[l] - p;
                 e[l] = g;
                 e[m] = 0.0;
             }
          }  while (m != l);
      }
 }

///////////////////////////////sort eigen values and vectors//////////////////////
void SolveEigen::eigsrt()
{
        int k,j,i;
        double p;

        for (i=1;i<n;i++) {
                p=d[k=i];
                for (j=i+1;j<=n;j++)
                        if (d[j] >= p) p=d[k=j];
                if (k != i) {
                        d[k]=d[i];
                        d[i]=p;
                        for (j=1;j<=n;j++) {
                                p=z[j][i];
                                z[j][i]=z[j][k];
                                z[j][k]=p;
                        }
                }
        }
}