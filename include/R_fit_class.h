#include <iostream>
#include <fstream>

using namespace std;

#include "stdlib.h"
#include "TROOT.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "math.h"
#include "TFile.h"
#include "TObject.h"
#include "TKey.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TChain.h"
#include "TString.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEventList.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"

class R_fit_class
{
public:
    R_fit_class();

    void txt_data_read(string A);
    void txt_data_write(string A);
    void GE_model_gen(int id);
    void GE_noise_gen(int id);
    void z_trans();
    double z_calc(double Q2_in, double Tc_in, double T0_in);
    void GE_z_fit(int id);
    void GE_fit(int id);

    void set_npower(int np, int dnp, int fnorm);
    void set_para(double *parain, int npara);
    void set_Tc(double Tcin);
    void show_data(int nshow);
    void get_fit_result(double *fres);

private:
    int ndata;
    int npower;
    int dnpower;
    int f_norm;
    double Tc;
    double Q2_txt[1000],Q2[1000];
    double GE_txt[1000],GE[1000];
    double dGE_txt[1000],dGE[1000];
    double zero[1000];
    double z_arr[1000];
    double chi2fit,Rfit,Rfiterr;
    double param[100];

};
