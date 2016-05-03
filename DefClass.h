#ifndef DEFCLASS_H
#define DEFCLASS_H
#include <vector>
#include <fstream>

class Particle{
    public:
        Particle();
        void Init();
        double x_,y_,z_,m_,alpha_,En_,Ex_,Ey_,Ez_,k_,phi_,ctheta_,stheta_,kx_,ky_,kz_,theta_,rest_dt_;
        int Iv_;
        void Into();
        void Drift(double);
        bool CheckBound();
        void getEn();
        void cart2sph();
        bool Needdt();
        void Scatter(int,std::vector<double>&);
    protected:
        void Rotate(double,double,double);
        void SetIv(int);
};

class Device{
    public:
        Device();
        void go();
    protected:
        bool nEqi_;
        double psi_s_,psi_d_,psi_g_,*GammaM_,Equi_n;
        std::vector<Particle> p_,temp_;
        std::vector<double> X_,Y_,psi_,tn_,tp_,meshn_,ND_,old_psi_,tpsi_,f_,c_,error_,xtheta_,Ex_,Ey_;
        std::vector< std::vector <std::vector <double> > > S_Table_;
        std::vector< std::vector <double> > Contact_;
        std::ofstream out_;
        int Source_,Drain_;
        void Add(double,double,double,double);
        void Poisson_Init();
        void Poisson_SOR();
        void MC_Init();
        void Read_S();
        void Read_one(std::ifstream&,int);
        int Scatter_Mech(int);
        void ApplyVolt();
        void MC_Field();
        void MC_Force();
        void Nequi_MC_Poisson();
        void MC_Flight(std::vector<double>&);
        double MC_Flight_one(int);
        void MC_CheckContact(std::vector<double>&);
        void MC_TempAdd(double,double,double,double);
        void MC_Add(std::vector<double>&);
        void PtoM();
        void OutSD();
        void disp(int);
};
#endif
