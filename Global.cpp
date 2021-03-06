#include "Global.h"
const double q=1.6e-19;
const double Vd=3.*q;
const double Vg=-1.*q;
const double kb=1.38e-23;
const double T=300.0;
const double Ndp=1e24;
const double Nd=1e23;
const double Ndm=1e21;
const double Vb=0.68*q;
const double ep0=8.8542e-12;
const double ep=8.9*ep0;
const double ep2=5.35*ep0;
const double Eg=3.2*q;
const double ni=1.5e16;
const double Vt=kb*T/q;
const double dE=0.001;
const double hbar=6.626e-34/(2*M_PI);

const double n=Nd;
const double qd=sqrt(q*n/(ep*(kb*T/q)));
const double lambda=1.0/qd;
const double wlo=0.0912*q/hbar;
const double n0=1.0/(exp(hbar*wlo/(kb*T))-1.0);
const double rho=6150.0;
const double epz=0.375;
const double vs=6560.0;
const double v0=8.3*q;
const double Dif=1e11*q;
const double wif=wlo;
const double nif=n0;
const double Ndis=1e13;
const double EG[]={3.39*q,5.29*q,5.49*q};
const double m=9.1e-31;
const double M[]={0.2*m,m,m};
const double Alpha[]={0.189/q,0,0};
const double C=(L*1e-6)*(L*1e-6)*ni*q/(Vt*ep);
const double W=0.5e-6;

const double dt=0.5e-15;
const double Vdf=3*q;
const double L=0.002;
const double Lx=0.8;
const double Ls=0.2;
const double Ld=0.6;
const double Lg1=0.3;
const double Lg2=0.5;
const double Ly=0.65;
const double er=1e-10;
const int Nx=401;
const int Ny=326;
int numdt=1000;
const int Nmax=40000;
