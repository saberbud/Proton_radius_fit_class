#include "include/R_fit_class.h"
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

#include "TChain.h"
#include "TString.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEventList.h"
#include "TMath.h"

int main(Int_t argc, char *argv[])
{
  double fit_result[3];
  cout << "main" << endl;
  R_fit_class Rfc;

  string filename="xy_gen_Q2_GE.txt";
  Rfc.txt_data_read(filename);
  // Rfc.show_data(0);

  /*Rfc.GE_model_gen(0);
  Rfc.show_data(0);

  Rfc.GE_fit(1);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.GE_noise_gen(1);

  Rfc.GE_fit(1);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;
  */

  cout << endl;
  double parain[2];
  int npara=2;
  parain[0]=0.85;parain[1]=0.995;
  Rfc.set_para(parain, npara);
  Rfc.GE_model_gen(3);

  Rfc.GE_fit(1);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.GE_fit(3);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.GE_fit(5);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.set_Tc(2.0);
  Rfc.z_trans();
  Rfc.set_npower(2, 0, 1);
  Rfc.GE_z_fit(0);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  // return 0;

  cout << endl;
  cout << "Series tests" << endl;


  int id_gen,id_nois;
  id_gen=2;  //1: dipole, 2: monopole 3: Gaussian with float norm
  id_nois=2;  //Gaussian noise

  int nfits=8;
  int id_fit[100];
  /*id_fit[0]=5;  //Gaussian with float norm
  id_fit[1]=1;  //Dipole with float norm
  id_fit[2]=3;  //Monopole with float norm
  id_fit[3]=6;  //poly fit
  id_fit[4]=6;  //poly fit
  id_fit[5]=6;  //poly fit
  id_fit[6]=6;  //poly fit
  id_fit[7]=7;  //DH ratio fit
  */

  nfits=4;
  id_fit[0]=0;
  id_fit[1]=0;
  id_fit[2]=0;
  id_fit[3]=1;


  int npset;

  double fit_R[100],fit_Rerr[100],fit_chi2[100];

  TFile *outf = new TFile("out.root","RECREATE");
  TTree *tout = new TTree("T","T");
  tout->Branch("nfits", &nfits, "nfits/I");
  tout->Branch("fit_R", &fit_R, "fit_R[nfits]/D");
  tout->Branch("fit_Rerr", &fit_Rerr, "fit_Rerr[nfits]/D");
  tout->Branch("fit_chi2", &fit_chi2, "fit_chi2[nfits]/D");


  int nrp=3000;
  for(int k=0;k<nrp;k++){
    cout << endl;
    Rfc.GE_model_gen(id_gen);

    /*for(int j=0;j<nfits;j++){
      Rfc.GE_fit(id_fit[j]);
      Rfc.get_fit_result(fit_result);
      fit_R[j]=fit_result[0];
      fit_Rerr[j]=fit_result[1];
      fit_chi2[j]=fit_result[2];
      // cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;
      cout << "R= " << fit_R[j] << " , Rerr= " << fit_Rerr[j] << " , chi2= " << fit_chi2[j] << endl;
    }
    */

    Rfc.GE_noise_gen(id_nois);

    npset=0;
    for(int j=0;j<nfits;j++){
      // if(id_fit[j]==6){
      if(id_fit[j]==0){
        npset=npset+1;
        Rfc.set_npower(npset, 0, 1);
      }
      // Rfc.GE_fit(id_fit[j]);
      Rfc.GE_z_fit(id_fit[j]);
      Rfc.get_fit_result(fit_result);
      fit_R[j]=fit_result[0];
      fit_Rerr[j]=fit_result[1];
      fit_chi2[j]=fit_result[2];
      cout << "R= " << fit_R[j] << " , Rerr= " << fit_Rerr[j] << " , chi2= " << fit_chi2[j] << endl;
    }

    tout->Fill();
  }
  tout->Write();
  outf->Save();
  outf->Close();


  return 0;

  cout << endl;

  Rfc.GE_noise_gen(2);
  // Rfc.show_data(0);
  Rfc.GE_fit(5);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.GE_fit(1);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.GE_fit(3);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.set_npower(1, 0, 1);
  Rfc.GE_fit(6);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.set_npower(2, 0, 1);
  Rfc.GE_fit(6);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.set_npower(3, 0, 1);
  Rfc.GE_fit(6);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.set_npower(4, 0, 1);
  Rfc.GE_fit(6);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  return 0;

  for(int i=0;i<6;i++){
    Rfc.GE_fit(i);
    Rfc.get_fit_result(fit_result);
    cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;
  }

  Rfc.set_npower(2, 0, 1);
  Rfc.GE_fit(6);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.set_npower(2, 0, 0);
  Rfc.GE_fit(6);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  Rfc.GE_fit(7);
  Rfc.get_fit_result(fit_result);
  cout << "R= " << fit_result[0] << " , Rerr= " << fit_result[1] << " , chi2= " << fit_result[2] << endl;

  return 0;
}
















