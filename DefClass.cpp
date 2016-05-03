#include <iostream>
#include <math.h>
#include <vector>
#include <algorithm>
//#include <iterator>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include "DefClass.h"
#include "Global.h"

Particle::Particle(){
}

void Particle::Init(){
    x_=.0; y_=.0; z_=.0; Ex_=.0; Ey_=.0; Ez_=.0;
    Iv_=1;
    m_=M[Iv_-1];
    alpha_=Alpha[Iv_-1];
    rest_dt_=0.0;

    En_=-1.5*kb*T*log(rd());
    k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
    phi_=2.*M_PI*rd();
    ctheta_=1.-2.*rd();
    stheta_=sqrt(1-ctheta_*ctheta_);
    kx_=k_*stheta_*cos(phi_);
    ky_=k_*stheta_*sin(phi_);
    kz_=k_*ctheta_;
}

void Particle::Into(){
    if (ky_>0.){
        ky_=-ky_;
        cart2sph();
    }
}

void Particle::cart2sph(){
    k_=sqrt(kx_*kx_+ky_*ky_+kz_*kz_);
    theta_=acos(kz_/k_);
    phi_=atan(ky_/kx_);
    ctheta_=cos(theta_);
    stheta_=sin(theta_);
}

void Particle::getEn(){
    if (Iv_==1){
        En_=(sqrt(1.+4.*alpha_*hbar*hbar*k_*k_/(2.*m_))-1.)/(2.*alpha_);
    }
    else{
        En_=hbar*hbar*k_*k_/(2.*m_);
    }
}

void Particle::Drift(double dt){
    double kx=kx_;
    double ky=ky_;

    kx_=kx-q*Ex_*dt/hbar;
    ky_=ky-q*Ey_*dt/hbar;

    x_=x_+hbar*dt/(2.*m_)*(kx+kx_)/(1e-6);
    y_=y_+hbar*dt/(2.*m_)*(ky+ky_)/(1e-6);
    Particle::cart2sph();
    Particle::getEn();
}

bool Particle::CheckBound(){
    if ((x_<Ls && y_>Ly)||(x_>Ld && y_>Ly))
        return false;
    else if (x_<0. && y_<0.){
        x_=-x_;
        y_=-y_;
        kx_=-kx_;
        ky_=-ky_;
        return true;
    }
    else if (x_>Lx && y_<0.){
        x_=2.*Lx-x_;
        y_=-y_;
        kx_=-kx_;
        ky_=-ky_;
        return true;
    }
    else if (x_<0.){
        x_=-x_;
        kx_=-kx_;
        return true;
    }
    else if (x_>Lx){
        x_=2.*Lx-x_;
        kx_=-kx_;
        return true;
    }
    else if (y_<0.){
        y_=-y_;
        ky_=-ky_;
        return true;
    }
    else if ((x_>Ls && y_>Ly) && (x_<Ld && y_>Ly)){
        y_=2.*Ly-y_;
        ky_=-ky_;
        return true;
    }
    else return true;
}

void Particle::Scatter(int mech,std::vector<double>& xtheta_){
    int i;
    double tphi,trd,tctheta,tstheta,xi,xtheta,temp;
    std::vector <double> y1_;
    switch (mech) {
        case 1:
            y1_.resize(10000);
            temp=rd();
            for (i=0;i<xtheta_.size();i++){
                y1_[i]=atan(sqrt(1.+4.*k_*k_*lambda*lambda)*tan(xtheta_[i]/2.));
                if (y1_[i]<0.) y1_[i]=y1_[i]+M_PI;
                y1_[i]=y1_[i]/M_PI+sin(xtheta_[i])*sqrt(1.+4.*k_*k_*lambda*lambda)/(2.*M_PI*(-1.-2.*k_*k_*lambda*lambda+2.*k_*k_*lambda*lambda*cos(xtheta_[i])));
                y1_[i]=fabs(y1_[i]-temp);
            }
            tphi=2.*M_PI*rd();
            xtheta=xtheta_[std::distance(y1_.begin(),std::min_element(y1_.begin(),y1_.end()))];
            Particle::Rotate(sin(xtheta),cos(xtheta),tphi);
            y1_.clear();
            break;
        case 2:
            tphi=2.*M_PI*rd();
            trd=rd();
            tctheta=1.-2.*trd/(1.+4.*k_*k_*lambda*lambda*(1.-trd));
            tstheta=sqrt(1.-tctheta*tctheta);
            Particle::Rotate(tstheta,tctheta,tphi);
            break;
        case 3:
            tphi=2.*M_PI*rd();
            xi=2.*sqrt(En_*(En_+hbar*wlo))/pow(sqrt(En_)-sqrt(En_+hbar*wlo),2.0);
            trd=rd();
            tctheta=((1.+xi)-pow(1.+2.*xi,trd))/xi;
            tstheta=sqrt(1.-tctheta*tctheta);
            En_=En_+hbar*wlo;
            k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
            Particle::Rotate(tstheta,tctheta,tphi);
            break;
        case 4:
            tphi=2.*M_PI*rd();
            if ((En_-hbar*wlo)<0.) {
//                std::cout << "Add energy" << std::endl;
                En_=En_+dE*q;
//                std::cout << En_ << std::endl;
            }
            xi=2.*sqrt(En_*(En_-hbar*wlo))/pow(sqrt(En_)-sqrt(En_-hbar*wlo),2.0);
            trd=rd();
            tctheta=((1.+xi)-pow(1.+2.*xi,trd))/xi;
            tstheta=sqrt(1.-tctheta*tctheta);
            En_=En_-hbar*wlo;
            k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
            Particle::Rotate(tstheta,tctheta,tphi);
            break;
        case 5:
            tphi=2.*M_PI*rd();
            trd=rd();
            tctheta=(2.*k_*k_+qd*qd-pow((4.*k_*k_+qd*qd)/(qd*qd),trd)*qd*qd)/(2.*k_*k_);
            tstheta=sqrt(1.-tctheta*tctheta);
            Particle::Rotate(tstheta,tctheta,tphi);
            break;
        case 6:
            phi_=2.*M_PI*rd();
            trd=rd();
            ctheta_=1.-2.*trd;
            stheta_=sqrt(1.-ctheta_*ctheta_);
            kx_=k_*stheta_*cos(phi_);
            ky_=k_*stheta_*sin(phi_);
            kz_=k_*ctheta_;
            break;
        case 7:
            phi_=2.*M_PI*rd();
            trd=rd();
            if (Iv_==1) {
                Particle::SetIv(2);
                En_=En_+hbar*wif-(EG[1]-EG[0]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            else if (Iv_==2){
                Particle::SetIv(1);
                En_=En_+hbar*wif-(EG[0]-EG[1]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            else {
                Particle::SetIv(1);
                En_=En_+hbar*wif-(EG[0]-EG[2]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            break;
        case 8:
            phi_=2.*M_PI*rd();
            trd=rd();
            if (Iv_==1) {
                Particle::SetIv(2);
                if ((En_-hbar*wif-(EG[1]-EG[0]))<0.) En_=En_+dE*q;
                En_=En_-hbar*wif-(EG[1]-EG[0]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            else if (Iv_==2){
                Particle::SetIv(1);
                if ((En_-hbar*wif-(EG[0]-EG[1]))<0.) En_=En_+dE*q;
                En_=En_-hbar*wif-(EG[0]-EG[1]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            else {
                Particle::SetIv(1);
                if ((En_-hbar*wif-(EG[0]-EG[2]))<0.) En_=En_+dE*q;
                En_=En_-hbar*wif-(EG[0]-EG[2]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            break;
        case 9:
            phi_=2.*M_PI*rd();
            trd=rd();
            if (Iv_==1) {
                Particle::SetIv(3);
                En_=En_+hbar*wif-(EG[2]-EG[0]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            else if (Iv_==2){
                Particle::SetIv(3);
                En_=En_+hbar*wif-(EG[2]-EG[1]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            else {
                Particle::SetIv(2);
                En_=En_+hbar*wif-(EG[1]-EG[2]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            break;
        case 10:
            phi_=2.*M_PI*rd();
            trd=rd();
            if (Iv_==1) {
                Particle::SetIv(3);
                if ((En_-hbar*wif-(EG[2]-EG[0]))<0.) En_=En_+dE*q;
                En_=En_-hbar*wif-(EG[2]-EG[0]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            else if (Iv_==2){
                Particle::SetIv(3);
                if ((En_-hbar*wif-(EG[2]-EG[1]))<0.) En_=En_+dE*q;
                En_=En_-hbar*wif-(EG[2]-EG[1]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            else {
                Particle::SetIv(2);
                if ((En_-hbar*wif-(EG[1]-EG[2]))<0.) En_=En_+dE*q;
                En_=En_-hbar*wif-(EG[1]-EG[2]);
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            break;
        case 11:
            phi_=2.*M_PI*rd();
            trd=rd();
            if (Iv_==1 || Iv_==2) std::cerr<< "scatter LM error" << std::endl;
            else {
                En_=En_+hbar*wif;
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            break;
        case 12:
            phi_=2.*M_PI*rd();
            trd=rd();
            if (Iv_==1 || Iv_==2) std::cerr<< "scatter LM error" << std::endl;
            else {
                if ((En_-hbar*wif)<0.) En_=En_+dE*q;
                En_=En_-hbar*wif;
                k_=sqrt(2.*m_*En_*(1.+alpha_*En_))/hbar;
                ctheta_=1.-2.*trd;
                stheta_=sqrt(1.-ctheta_*ctheta_);
                kx_=k_*stheta_*cos(phi_);
                ky_=k_*stheta_*sin(phi_);
                kz_=k_*ctheta_;
            }
            break;
    };
}

void Particle::SetIv(int Iv){
    Iv_=Iv;
    m_=M[Iv-1];
    alpha_=Alpha[Iv-1];
}
void Particle::Rotate(double tstheta,double tctheta,double tphi){
    double kxp=k_*tstheta*cos(tphi);
    double kyp=k_*tstheta*sin(tphi);
    double kzp=k_*tctheta;

    kx_=kxp*cos(phi_)*ctheta_-kyp*sin(phi_)+kzp*cos(phi_)*stheta_;
    ky_=kxp*sin(phi_)*ctheta_+kyp*cos(phi_)+kzp*sin(phi_)*stheta_;
    kz_=-kxp*stheta_+kzp*ctheta_;
    Particle::cart2sph();
//    kx_=k_*stheta_*cos(phi_);
//    ky_=k_*stheta_*sin(phi_);
//    kz_=k_*ctheta_;
}

bool Particle::Needdt(){
    if (rest_dt_/dt<1e-4) return true;
    else return false;
}

Device::Device(){
    X_.resize(Nx*Ny);
    Y_.resize(Nx*Ny);
    for (int i=0;i<Nx;i++){
        for (int j=0;j<Ny;j++){
            X_[i*Ny+j]=L*i;
            Y_[i*Ny+j]=L*j;
        }
    }

    nEqi_=false;
    p_.reserve(Nmax);
    psi_.resize(Nx*Ny);
    old_psi_.resize(Nx*Ny);
    tpsi_.resize(Nx*Ny);
    error_.resize(Nx*Ny);
    ND_.resize(Nx*Ny);
    f_.resize(Nx*Ny);
    c_.resize(Nx*Ny);
    tn_.resize(Nx*Ny);
    tp_.resize(Nx*Ny);
    meshn_.resize(Nx*Ny);
    Ex_.resize(Nx*Ny);
    Ey_.resize(Nx*Ny);

    std::vector <std::vector <double> > temp;
    temp.resize(10000,std::vector<double> (13,0.));
    S_Table_.resize(3,temp);
    temp.clear();

    double dx=2.*M_PI/10000.;
    for (int i=0;i<10000;i++){
        xtheta_.push_back(i*dx);
    }

    Contact_.resize(round((Ls+Lx-Ld)/L));

    Source_=0;
    Drain_=0;
    
    char filename[20];
    std::ofstream output;

    sprintf(filename,"Result.txt");
    out_.open(filename,std::ostream::out);
}

void Device::Add(double x,double dx,double y,double dy){
    Particle newPar; newPar.Init();
    newPar.x_=x+dx*rd();
    newPar.y_=y+dy*rd();
    p_.push_back(newPar);
}

void Device::Poisson_Init(){
    double tNdp=Ndp/ni;
    double tNd=Nd/ni;
    double tNdm=Ndm/ni;

    for (int i=0;i<psi_.size();i++){
        if (Y_[i]<0.2) {ND_[i]=tNdm;psi_[i]=log(tNdm);}
        else if ((X_[i]<=Ls || X_[i]>=Ld) && Y_[i]>=0.6) {ND_[i]=tNdp;psi_[i]=log(tNdp);}
        else {ND_[i]=tNd;psi_[i]=log(tNd);}
    }
    psi_s_=psi_[Ny-1];
    psi_d_=psi_[Nx*Ny-1];
    psi_g_=psi_[200*Ny-1]+(Eg/2.-Vb+Vg)/q/Vt;

    Equi_n=exp(psi_d_)*ni*(L*1e-6)*(L*1e-6)*W*0.5;
    std::cout << Equi_n << std::endl;
}

void Device::Poisson_SOR(){
    int neqi=0;
    if (nEqi_) neqi=1;
    bool condition=true;
    int Ind_error=1;

    double w=1.8;
    double terr=1e-5,err;

    int i;
    if (nEqi_){
        terr=1e-2;
        for (i=0;i<tn_.size();i++){
            tn_[i]=meshn_[i];
        }
    }

    while (condition){
        for (i=0;i<Nx*Ny;i++){
            tp_[i]=exp(-psi_[i]);
        }
        if (nEqi_==false){
            for (i=0;i<Nx*Ny;i++){
                tn_[i]=exp(psi_[i]);
            }
        }
        std::copy(psi_.begin(),psi_.end(),old_psi_.begin());
        std::copy(psi_.begin(),psi_.end(),tpsi_.begin());
        for (i=0;i<Nx*Ny;i++){
            f_[i]=-C*(tp_[i]-tn_[i]+ND_[i])-C*(tp_[i]+tn_[i])*psi_[i];
            c_[i]=4.+C*(tn_[i]+tp_[i]);
        }

        for (i=0;i<Nx*Ny;i++){
            if (Y_[i]>0. && Y_[i]<Ly && X_[i]>0. && X_[i]<Lx){
                tpsi_[i]=(tpsi_[i-1]+tpsi_[i+1]+tpsi_[i-Ny]+tpsi_[i+Ny]-f_[i])/c_[i];
            }
            else if (fabs(Y_[i])<er && X_[i]>0. && X_[i]<Lx){
                tpsi_[i]=(2.*tpsi_[i+1]+tpsi_[i-Ny]+tpsi_[i+Ny]-f_[i])/c_[i];
            }
            else if (fabs(Y_[i])<er && fabs(X_[i])<er){
                tpsi_[i]=(2.*tpsi_[i+1]+2.*tpsi_[i+Ny]-f_[i])/c_[i];
            }
            else if (fabs(Y_[i])<er && fabs(X_[i]-Lx)<er){
                tpsi_[i]=(2.*tpsi_[i+1]+2.*tpsi_[i-Ny]-f_[i])/c_[i];
            }
            else if (fabs(X_[i])<er && Y_[i]>0. && Y_[i]<Ly){
                tpsi_[i]=(tpsi_[i-1]+tpsi_[i+1]+2.*tpsi_[i+Ny]-f_[i])/c_[i];
            }
            else if (fabs(X_[i]-Lx)<er && Y_[i]>0. && Y_[i]<Ly){
                tpsi_[i]=(tpsi_[i-1]+tpsi_[i+1]+2.*tpsi_[i-Ny]-f_[i])/c_[i];
            }
            else if (fabs(Y_[i]-Ly)<er && ((X_[i]>Ls && X_[i]<Lg1) || (X_[i]>Lg2 && X_[i]<Ld))){
                tpsi_[i]=(2.*tpsi_[i-1]+tpsi_[i-Ny]+tpsi_[i+Ny]-f_[i])/c_[i];
            }
            else if (fabs(Y_[i]-Ly)<er && X_[i]<=Ls){
                tpsi_[i]=psi_s_;
            }
            else if (fabs(Y_[i]-Ly)<er && X_[i]>=Ld){
                tpsi_[i]=psi_d_;
            }
            else{
                tpsi_[i]=psi_g_;
            }
            tpsi_[i]=psi_[i]+w*(tpsi_[i]-psi_[i]);
        }
        std::copy(tpsi_.begin(),tpsi_.end(),psi_.begin());
        for (i=0;i<Nx*Ny;i++){
            error_[i]=fabs(psi_[i]-old_psi_[i])*Vt;
        }
        err=*std::max_element(error_.begin(),error_.end());
        condition=(err>terr);
        Ind_error++;
    }
}

void Device::MC_Init(){
    std::vector<double> n;
    n.resize(Nx*Ny);
    int i;
    for (i=0;i<Nx*Ny;i++){
        n[i]=tn_[i]*ni*(L*1e-6)*(L*1e-6)*W;
    }

    int tempn,j;
    for (i=0;i<Nx*Ny;i++){
        if (rd()<fmod(n[i],1.)) tempn=(int) ceil(n[i]);
        else tempn=(int) floor(n[i]);
        if (tempn>0){
            for (j=0;j<tempn;j++){
                Device::Add(X_[i]-L/2.,L,Y_[i]-L/2.,L);
            }
        }
    }

    for (i=0;p_.begin()+i<p_.end();i++){
        if (p_[i].x_<0. || p_[i].x_>Lx || p_[i].y_<0. || p_[i].y_>Ly){
            p_.erase(p_.begin()+i);
            i--;
        }
    }
    n.clear();

    Device::Read_S();
}

void Device::Read_S(){
    std::ifstream ifile("Dislocation.txt",std::ios::in);
    Device::Read_one(ifile,1);
    ifile.close();
    ifile.open("Coulomb.txt",std::ios::in);
    Device::Read_one(ifile,2);
    ifile.close();
    ifile.open("POP_Ab.txt",std::ios::in);
    Device::Read_one(ifile,3);
    ifile.close();
    ifile.open("POP_Em.txt",std::ios::in);
    Device::Read_one(ifile,4);
    ifile.close();
    ifile.open("Piezoelectric.txt",std::ios::in);
    Device::Read_one(ifile,5);
    ifile.close();
    ifile.open("AP.txt",std::ios::in);
    Device::Read_one(ifile,6);
    ifile.close();
    ifile.open("InterAb1.txt",std::ios::in);
    Device::Read_one(ifile,7);
    ifile.close();
    ifile.open("InterEm1.txt",std::ios::in);
    Device::Read_one(ifile,8);
    ifile.close();
    ifile.open("InterAb2.txt",std::ios::in);
    Device::Read_one(ifile,9);
    ifile.close();
    ifile.open("InterEm2.txt",std::ios::in);
    Device::Read_one(ifile,10);
    ifile.close();
    ifile.open("EInterAb.txt",std::ios::in);
    Device::Read_one(ifile,11);
    ifile.close();
    ifile.open("EInterEm.txt",std::ios::in);
    Device::Read_one(ifile,12);
    ifile.close();
//    ifile.open("Self.txt",std::ios::in);
//    Device::Read_one(ifile,13);
//    ifile.close();

    ifile.open("GammaM.txt",std::ios::in);
    if (!ifile.is_open()) std::cerr << "cannot open scatter file\n";
    GammaM_=new double [3];
    ifile>>GammaM_[0];
    ifile>>GammaM_[1];
    ifile>>GammaM_[2];
    GammaM_[0]=exp(GammaM_[0]);
    GammaM_[1]=exp(GammaM_[1]);
    GammaM_[2]=exp(GammaM_[2]);
    ifile.close();
}

void Device::Read_one(std::ifstream& ifile, int mech){
    if (!ifile.is_open()) std::cerr << "cannot open scatter file\n";
    int i;
    double temp[3];
    for (i=0;i<10000;i++){
        ifile>>S_Table_[0][i][mech-1]; ifile>>S_Table_[1][i][mech-1]; ifile>>S_Table_[2][i][mech-1];
    }
//    for (i=0;i<10;i++){
//        std::cout << S_Table_[0][9999-i][mech] << std::endl;
//    }
}

int Device::Scatter_Mech(int j){
    int EIdx=floor(p_[j].En_/q/dE);
    if (EIdx>10000) EIdx=9999;
    double trd=rd();
    int i;
    int Iv=p_[j].Iv_;
    if (Iv!=3){
        for (i=0;i<11;i++){
            if (i==0){
                if (trd<S_Table_[Iv-1][EIdx][0]) {
                    return 1;}
            }
            else if (i<10){
                if (trd>S_Table_[Iv-1][EIdx][i-1] && trd<S_Table_[Iv-1][EIdx][i]) {
                    return i+1;
                }
            }
            else return 13;
        }
    }
    else {
        for (i=0;i<13;i++){
            if (i==0){
                if (trd<S_Table_[Iv-1][EIdx][0]) return 1;
            }
            else if (i<12){
                if (trd>S_Table_[Iv-1][EIdx][i-1] && trd<S_Table_[Iv-1][EIdx][i]) return i+1;
            }
            else return 13;
        }
    }
}

void Device::MC_Field(){
    int i;
    for (i=0;i<Nx*Ny;i++){
        if (X_[i]<L*0.1){
            Ex_[i]=-(psi_[i+Ny]-psi_[i])*Vt/(L*1e-6);
        }
        else if (X_[i]<Lx-L*0.1){
            Ex_[i]=-0.5*(psi_[i+Ny]-psi_[i-Ny])*Vt/(L*1e-6);
        }
        else {
            Ex_[i]=-(psi_[i]-psi_[i-Ny])*Vt/(L*1e-6);
        }

        if (Y_[i]<L*0.1){
            Ey_[i]=-(psi_[i+1]-psi_[i])*Vt/(L*1e-6);
        }
        else if (Y_[i]<Ly-L*0.1){
            Ey_[i]=-0.5*(psi_[i+1]-psi_[i-1])*Vt/(L*1e-6);
        }
        else {
            Ey_[i]=-(psi_[i]-psi_[i-1])*Vt/(L*1e-6);
        }
    }
}

void Device::MC_Force(){
    int i,idx,idy;
    Device::MC_Field();
    for (i=0;i<p_.size();i++){
        idx=floor(p_[i].x_/L);
        idy=floor(p_[i].y_/L);
        p_[i].Ex_=0.25*(Ex_[idx*Ny+idy]+Ex_[(idx+1)*Ny+idy]+Ex_[idx*Ny+idy+1]+Ex_[(idx+1)*Ny+idy+1]);
        p_[i].Ey_=0.25*(Ey_[idx*Ny+idy]+Ey_[(idx+1)*Ny+idy]+Ey_[idx*Ny+idy+1]+Ey_[(idx+1)*Ny+idy+1]);
    }
}

void Device::ApplyVolt(){
    psi_d_=psi_d_+Vd/q/Vt;
    int i;
    for (i=0;i<psi_.size();i++){
        if (X_[i]>=Ld-L*0.1 && fabs(Y_[i]-Ly)<er) psi_[i]=psi_d_;
    }
}

void Device::Nequi_MC_Poisson(){
    int i,j;
    std::vector<double> dtau,dt2;
    dt2.reserve(Nmax);
    double dt3,dtp,dte2,dte;
    bool loop;
    int mech;
    Device::MC_Flight(dtau);
    Device::MC_Force();
    for (i=0;i<numdt;i++){
        dt2.resize(dtau.size());
        std::copy(dtau.begin(),dtau.end(),dt2.begin());

        for (j=0;j<dt2.size();j++){
            if (dt2[j]>dt) dt2[j]=dt;
            p_[j].Drift(dt2[j]);
            p_[i].rest_dt_=p_[i].rest_dt_-dt2[j];
        }

        for (j=0;j<p_.size();j++){
            if (dtau[j]>=dt) dtau[j]=dtau[j]-dt;
            else {
                loop=true;
                while (loop){
                    mech=Device::Scatter_Mech(j);
                    p_[j].Scatter(mech,xtheta_);
                    dt3=Device::MC_Flight_one(j);
                    dtp=dt-dtau[j];
                    if (dt3<=dtp) {
                        dt2[j]=dt3;
                    }
                    else {
                        dt2[j]=dtp;
                    }
                    p_[j].Drift(dt2[j]);
                    p_[j].rest_dt_=p_[j].rest_dt_-dt2[j];
                    dte2=dtau[j]+dt3;
                    dte=dte2;
                    if (dte<dt) {
                        dtau[j]=dte;
                        loop=true;
                    }
                    else {
                        dtau[j]=dte-dt;
                        loop=false;
                    }
                }
            }
        }

        for (j=0;j<p_.size();j++){
            if (!p_[j].CheckBound()){
                if (p_[j].x_<=Ls) Source_-=1;
                else Drain_-=1;
                p_.erase(p_.begin()+j);
                dtau.erase(dtau.begin()+j);
                dt2.erase(dt2.begin()+j);
                j--;
            }
        }
        Device::MC_CheckContact(dtau);

        Device::PtoM();
        nEqi_=true;
        Poisson_SOR();
        Device::MC_Force();
        if (fmod(i,1000)==0) Device::disp(i);
    }
}

void Device::PtoM(){
    int i,idx,idy;
    for (i=0;i<Nx*Ny;i++){
        meshn_[i]=0.0;
    }
    for (i=0;i<p_.size();i++){
        idx=floor(p_[i].x_/L);
        idy=floor(p_[i].y_/L);
        if (idx>0 && idx<Nx-2 && idy>0 && idy<Ny-2){
            meshn_[idx*Ny+idy]++;
            meshn_[idx*Ny+idy+1]++;
            meshn_[(idx+1)*Ny+idy]++;
            meshn_[(idx+1)*Ny+idy+1]++;
        }
        else if (idx==0 && idy>0 && idy<Ny-2){
            meshn_[idx*Ny+idy]+=2;
            meshn_[(idx+1)*Ny+idy]+=1;
            meshn_[idx*Ny+idy+1]+=2;
            meshn_[(idx+1)*Ny+idy+1]++;
//            std::cout << "ADD 2: " << X_[idx*Ny+idy] << " " << X_[idx*Ny+idy+1] << std::endl;
        }
        else if (idx==Nx-2 && idy>0 && idy<Ny-2){
            meshn_[idx*Ny+idy]++;
            meshn_[(idx+1)*Ny+idy]+=2;
            meshn_[idx*Ny+idy+1]+=1;
            meshn_[(idx+1)*Ny+idy+1]+=2;
//            std::cout << "ADD 2: " << X_[(idx+1)*Ny+idy] << " " << X_[(idx+1)*Ny+idy+1] << std::endl;
        }
        else if (idy==0 && idx>0 && idx<Nx-2){
            meshn_[idx*Ny+idy]+=2;
            meshn_[(idx+1)*Ny+idy]+=2;
            meshn_[idx*Ny+idy+1]+=1;
            meshn_[(idx+1)*Ny+idy+1]+=1;
//            std::cout << "ADD 2: Y= " << Y_[idx*Ny+idy] << " " << Y_[(idx+1)*Ny+idy] << std::endl;
        }
        else if (idy==Ny-2 && idx>0 && idx<Nx-2){
            meshn_[idx*Ny+idy]+=1;
            meshn_[(idx+1)*Ny+idy]+=1;
            meshn_[idx*Ny+idy+1]+=2;
            meshn_[(idx+1)*Ny+idy+1]+=2;
//            std::cout << "ADD 2: Y= " << Y_[idx*Ny+idy+1] << " " << Y_[(idx+1)*Ny+idy+1] << std::endl;
        }
        else if (idx==0 && idy==0){
            meshn_[idx*Ny+idy]+=4;
            meshn_[(idx+1)*Ny+idy]+=2;
            meshn_[idx*Ny+idy+1]+=2;
            meshn_[(idx+1)*Ny+idy+1]+=1;
//            std::cout << "ADD 4: x,y= " << X_[idx*Ny+idy] << " " << Y_[idx*Ny+idy] << std::endl;
        }
        else if (idx==0 && idy==Ny-2){
            meshn_[idx*Ny+idy]+=2;
            meshn_[(idx+1)*Ny+idy]+=1;
            meshn_[idx*Ny+idy+1]+=4;
            meshn_[(idx+1)*Ny+idy+1]+=2;
//            std::cout << "ADD 4: x,y= " << X_[idx*Ny+idy+1] << " " << Y_[idx*Ny+idy+1] << std::endl;
        }
        else if (idx==Nx-2 && idy==0){
            meshn_[idx*Ny+idy]+=2;
            meshn_[(idx+1)*Ny+idy]+=4;
            meshn_[idx*Ny+idy+1]+=1;
            meshn_[(idx+1)*Ny+idy+1]+=2;
//            std::cout << "ADD 4: x,y= " << X_[(idx+1)*Ny+idy] << " " << Y_[(idx+1)*Ny+idy] << std::endl;
        }
        else if (idx==Nx-2 && idy==Ny-2){
            meshn_[idx*Ny+idy]+=1;
            meshn_[(idx+1)*Ny+idy]+=2;
            meshn_[idx*Ny+idy+1]+=2;
            meshn_[(idx+1)*Ny+idy+1]+=4;
//            std::cout << "ADD 4: x,y= " << X_[(idx+1)*Ny+idy+1] << " " << Y_[(idx+1)*Ny+idy+1] << std::endl;
        }
    }

    for (i=0;i<Nx*Ny;i++){
        meshn_[i]=meshn_[i]/4.0/(ni*(L*1e-6)*(L*1e-6)*W);
    }
}

void Device::MC_CheckContact(std::vector<double>& dtau){
     int i,j;
     for (i=0;i<round((Ls+Lx-Ld)/L);i++){
         Contact_[i].resize(0);
     }

     int idx,idy,del,add;
     double rest;
     std::vector<int> ToDel;
     ToDel.resize(0);
     for (i=0;i<p_.size();i++){
         if (p_[i].x_<=Ls && (p_[i].y_-(Ly-0.5*L))>0){
             idx=floor(p_[i].x_/L);
             Contact_[idx].push_back(i);
         }
         else if ((p_[i].x_>=Ld) && (p_[i].y_-(Ly-0.5*L))>0){
             idx=floor((p_[i].x_-Ld)/L);
             Contact_[idx+100].push_back(i);
         }
     }

     for (i=0;i<round((Ls+Lx-Ld)/L);i++){
         if (Contact_[i].size()>Equi_n){
             del=floor(Contact_[i].size()-Equi_n);
             rest=Contact_[i].size()-del-Equi_n;
             if (rd()<rest) {
                 del++;
             }
             if (del>0){
                 std::random_shuffle(Contact_[i].begin(),Contact_[i].end());
                 for (j=0;j<del;j++){
                     ToDel.push_back(Contact_[i][j]);
                 }
             }
             if (i<round((Ls+Lx-Ld)/L)/2) Source_-=del;
             else Drain_-=del;
         }
         else if (Contact_[i].size()<Equi_n){
             add=floor(Equi_n-Contact_[i].size());
             rest=Equi_n-add-Contact_[i].size();
             if (rd()<rest){
                 add++;
             }
             if (add>0){
                 for (j=0;j<add;j++){
                     Device::MC_TempAdd(i*L+(Ld-Ls)*floor((i+1)/100.5),L,Ly-0.5*L,0.5*L);
                 }
             }
             if (i<round((Ls+Lx-Ld)/L)/2) Source_+=add;
             else Drain_+=add;
         }
     }

     std::sort(ToDel.begin(),ToDel.end());
     for (i=0;i<ToDel.size();i++){
         p_.erase(p_.begin()+ToDel[ToDel.size()-1-i]);
         dtau.erase(dtau.begin()+ToDel[ToDel.size()-1-i]);
     }
     Device::MC_Add(dtau);
     Device::OutSD();
}

void Device::MC_TempAdd(double x,double dx,double y,double dy){
    Particle newPar;
    newPar.Init();
    newPar.x_=x+dx*rd();
    newPar.y_=y+dy*rd();
    newPar.Into();
    temp_.push_back(newPar);
}

void Device::MC_Add(std::vector<double>& dtau){
    if (temp_.size()>0){
        for (int i=0;i<temp_.size();i++){
            p_.push_back(temp_[i]);
            dtau.push_back(Device::MC_Flight_one(p_.size()-1));
        }
    }
    temp_.resize(0);
}

void Device::MC_Flight(std::vector<double>& dtau){
    dtau.resize(p_.size());
    int i;
    for (i=0;i<p_.size();i++){
        dtau[i]=Device::MC_Flight_one(i);
    }
}

double Device::MC_Flight_one(int j){
    if (p_[j].Needdt()) {
        p_[j].rest_dt_=-log(rd())/GammaM_[p_[j].Iv_-1];
        return p_[j].rest_dt_;
    }
    else return p_[j].rest_dt_;
}

void Device::OutSD(){
    out_<< Source_ << "\t" << Drain_ << std::endl;
    Source_=0;
    Drain_=0;
}

void Device::disp(int Ind){
    char filename[20];
    std::ofstream output;

    if (Ind==0){
        sprintf(filename,"opsi.txt");
        output.open(filename,std::ostream::out);
        for (int i=0;i<Nx*Ny;i++){
            output<< psi_[i] << std::endl;
        }
        output.close();

        sprintf(filename,"omesh.txt");
        output.open(filename,std::ostream::out);
        for (int i=0;i<Nx*Ny;i++){
            output<< meshn_[i] << std::endl;
        }
        output.close();

        sprintf(filename,"out.txt");
        output.open(filename,std::ostream::out);
        for (int i=0;i<p_.size();i++){
            output<< p_[i].x_ << "\t" << p_[i].y_ << "\t" << p_[i].Ex_ << "\t" << p_[i].Ey_ << std::endl;
        }
        output.close();
    }
    else {
        sprintf(filename,"opsi%d.txt",Ind);
        output.open(filename,std::ostream::out);
        for (int i=0;i<Nx*Ny;i++){
            output<< psi_[i] << std::endl;
        }
        output.close();

        sprintf(filename,"omesh%d.txt",Ind);
        output.open(filename,std::ostream::out);
        for (int i=0;i<Nx*Ny;i++){
            output<< meshn_[i] << std::endl;
        }
        output.close();

        sprintf(filename,"out%d.txt",Ind);
        output.open(filename,std::ostream::out);
        for (int i=0;i<p_.size();i++){
            output<< p_[i].x_ << "\t" << p_[i].y_ << "\t" << p_[i].Ex_ << "\t" << p_[i].Ey_ << std::endl;
        }
        output.close();
    }
}

void Device::go(){
    Device::Poisson_Init();
    nEqi_=false;
    Device::Poisson_SOR();
    Device::MC_Init();
    Device::ApplyVolt();
    Device::Nequi_MC_Poisson();
    Device::disp(0);
}
