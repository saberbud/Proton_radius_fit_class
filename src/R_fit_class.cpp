#include "R_fit_class.h"

R_fit_class::R_fit_class()
{
  //place holder
  for(int i=0;i<1000;i++){
    zero[i]=0.;
    Q2[i]=0.;
    GE[i]=0.;
    dGE[i]=0.;
    if(i<100)param[i]=0.;
  }

}

void R_fit_class::txt_data_read(string A)
{
  TString At=A;
  cout << "Read txt file: " << At << endl;

  ifstream infile(At);
  if(!infile.is_open()){
    cout << "File not exist." << endl;
    return;
  }
  infile >> ndata;

  for(int i=0;i<ndata;i++){
    infile >> Q2_txt[i] >> GE_txt[i] >> dGE_txt[i];
  }
  infile.close();
}

void R_fit_class::GE_model_gen(int id)
{
  cout << "GE_model_gen id= " << id << endl;
  if(id==0){  //Use input txt file values
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      GE[i]=GE_txt[i];
      dGE[i]=dGE_txt[i];
    }
  }else if(id==1){
    cout << "Generate, dipole: radius input= " << param[0] << " , norm= " << param[1] << endl;
    double p0=12./(param[0]*param[0]);
    double p1=param[1];
    TF1 *poly = new TF1("poly", "[1]/(1.+x/[0])/(1.+x/[0])", 0, 10.5);
    poly->SetParameter(0, p0);
    poly->SetParameter(1, p1);
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      GE[i]=poly->Eval(Q2[i]);
      GE[i]=GE[i]*p1;
      dGE[i]=dGE_txt[i];
    }

    delete poly;
  }else if(id==2){
    cout << "Generate, monopole: radius input= " << param[0] << " , norm= " << param[1] << endl;
    double p0=6./(param[0]*param[0]);
    double p1=param[1];
    TF1 *poly = new TF1("poly", "[1]/(1.+x/[0])", 0, 10.5);
    poly->SetParameter(0, p0);
    poly->SetParameter(1, p1);
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      GE[i]=poly->Eval(Q2[i]);
      GE[i]=GE[i]*p1;
      dGE[i]=dGE_txt[i];
    }

    delete poly;
  }else if(id==3){
    cout << "Generate, Gaussian: radius input= " << param[0] << " , norm= " << param[1] << endl;
    double p0=6./(param[0]*param[0]);
    double p1=param[1];
    TF1 *poly = new TF1("poly", "[1]*exp(-x/[0])", 0, 10.5);
    poly->SetParameter(0, p0);
    poly->SetParameter(1, p1);
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      GE[i]=poly->Eval(Q2[i]);
      GE[i]=GE[i]*p1;
      dGE[i]=dGE_txt[i];
    }

    delete poly;
  }


}

void R_fit_class::GE_noise_gen(int id)
{
  cout << "GE_noise_gen id= " << id << endl;

  TRandom3 * rdm = new TRandom3();
  rdm->SetSeed(0);

  double a0,a1;
  double temp;
  if(id==0){  //add uniform noise according to dGE
    cout << "Add uniform noise, type 0" << endl;
    for(int i=0;i<ndata;i++){
      a0=-1.*dGE_txt[i];
      a1=dGE_txt[i];
      temp=rdm->Uniform(a0, a1);

      GE[i]=GE[i]+temp;
      // cout << i << " | " << temp << endl;
    }
  }else if(id==1){
    cout << "Add BreitWigner noise, type 1" << endl;
    for(int i=0;i<ndata;i++){
      a0=0.;
      a1=dGE_txt[i];
      temp=rdm->BreitWigner(a0, a1);
      if(temp>(3.*a1))temp=0.;

      GE[i]=GE[i]+temp;
    }
  }else if(id==2){  
    cout << "Add Gaussian noise, type 2" << endl;
    for(int i=0;i<ndata;i++){
      a0=0.;
      a1=dGE_txt[i];
      temp=rdm->Gaus(a0, a1);

      GE[i]=GE[i]+temp;
    }
  }else if(id==3){
    cout << "Add Landau noise, type 3" << endl;
    for(int i=0;i<ndata;i++){
      a0=0.;
      a1=dGE_txt[i];
      temp=rdm->Landau(a0, a1);
      if(temp>(3.*a1))temp=0.;

      GE[i]=GE[i]+temp;
    }
  }else if(id==4){
    cout << "Add sine noise, type 4" << endl;
    for(int i=0;i<ndata;i++){
      a0=-1.5;
      a1=1.5;
      temp=rdm->Uniform(a0, a1);
      temp=sin(temp)*dGE_txt[i];

      GE[i]=GE[i]+temp;
    }
  }

  delete rdm;
}

void R_fit_class::GE_fit(int id)
{
  cout << "GE_fit id= " << id << endl;

  chi2fit=0.;Rfit=0.;Rfiterr=0.;

  TGraphErrors *gr=new TGraphErrors(ndata,Q2,GE,zero,dGE);

  double R_c,R_1,R_2,R_E;
  if(id==0){
    cout << "Fit: full dipole" << endl;
    TF1 *poly = new TF1("poly", "1./(1.+x/[0])/(1.+x/[0])", 0, 0.5);
    poly->SetParameter(0, 50.);
    gr->Fit("poly","0");
    R_c=sqrt(12./poly->GetParameter(0));
    R_1=sqrt(12./(poly->GetParameter(0)+poly->GetParError(0)));
    R_2=sqrt(12./(poly->GetParameter(0)-poly->GetParError(0)));
    R_E=(R_2-R_1)/2.;
    chi2fit=poly->GetChisquare();
    Rfit=R_c;
    Rfiterr=R_E;
    cout << "R= " << Rfit << " | Error= " << Rfiterr << " | chi2= " << chi2fit << endl;

    delete poly;
  }else if(id==1){
    cout << "Fit: full dipole with float norm" << endl;
    TF1 *poly = new TF1("poly", "[1]/(1.+x/[0])/(1.+x/[0])", 0, 0.5);
    poly->SetParameter(0, 50.);
    gr->Fit("poly","0");  
    R_c=sqrt(12./poly->GetParameter(0));
    R_1=sqrt(12./(poly->GetParameter(0)+poly->GetParError(0)));
    R_2=sqrt(12./(poly->GetParameter(0)-poly->GetParError(0)));
    R_E=(R_2-R_1)/2.;
    chi2fit=poly->GetChisquare();
    Rfit=R_c;
    Rfiterr=R_E;
    cout << "R= " << Rfit << " | Error= " << Rfiterr << " | chi2= " << chi2fit << endl;

    delete poly;
  }else if(id==2){
    cout << "Fit: full monopole" << endl;
    TF1 *poly = new TF1("poly", "1./(1.+x/[0])", 0, 0.5);
    poly->SetParameter(0, 50.);
    gr->Fit("poly","0");
    R_c=sqrt(6./poly->GetParameter(0));
    R_1=sqrt(6./(poly->GetParameter(0)+poly->GetParError(0)));
    R_2=sqrt(6./(poly->GetParameter(0)-poly->GetParError(0)));
    R_E=(R_2-R_1)/2.;
    chi2fit=poly->GetChisquare();
    Rfit=R_c;
    Rfiterr=R_E;
    cout << "R= " << Rfit << " | Error= " << Rfiterr << " | chi2= " << chi2fit << endl;

    delete poly;
  }else if(id==3){
    cout << "Fit: full monopole with float norm" << endl;
    TF1 *poly = new TF1("poly", "[1]/(1.+x/[0])", 0, 0.5);
    poly->SetParameter(0, 50.);
    gr->Fit("poly","0");
    R_c=sqrt(6./poly->GetParameter(0));
    R_1=sqrt(6./(poly->GetParameter(0)+poly->GetParError(0)));
    R_2=sqrt(6./(poly->GetParameter(0)-poly->GetParError(0)));
    R_E=(R_2-R_1)/2.;
    chi2fit=poly->GetChisquare();
    Rfit=R_c;
    Rfiterr=R_E;
    cout << "R= " << Rfit << " | Error= " << Rfiterr << " | chi2= " << chi2fit << endl;

    delete poly;
  }else if(id==4){
    cout << "Fit: full Gaussian" << endl;
    TF1 *poly = new TF1("poly", "exp(-x/[0])", 0, 0.5);
    poly->SetParameter(0, 50.);
    gr->Fit("poly","0");
    R_c=sqrt(6./poly->GetParameter(0));
    R_1=sqrt(6./(poly->GetParameter(0)+poly->GetParError(0)));
    R_2=sqrt(6./(poly->GetParameter(0)-poly->GetParError(0)));
    R_E=(R_2-R_1)/2.;
    chi2fit=poly->GetChisquare();
    Rfit=R_c;
    Rfiterr=R_E;
    cout << "R= " << Rfit << " | Error= " << Rfiterr << " | chi2= " << chi2fit << endl;

    delete poly;
  }else if(id==5){
    cout << "Fit: full Gaussian with float norm" << endl;
    TF1 *poly = new TF1("poly", "[1]*exp(-x/[0])", 0, 0.5);
    poly->SetParameter(0, 50.);
    gr->Fit("poly","0");
    R_c=sqrt(6./poly->GetParameter(0));
    R_1=sqrt(6./(poly->GetParameter(0)+poly->GetParError(0)));
    R_2=sqrt(6./(poly->GetParameter(0)-poly->GetParError(0)));
    R_E=(R_2-R_1)/2.;
    chi2fit=poly->GetChisquare();
    Rfit=R_c;
    Rfiterr=R_E;
    cout << "R= " << Rfit << " | Error= " << Rfiterr << " | chi2= " << chi2fit << endl;

    delete poly;
  }else if(id==6){
    cout << "Fit: multiple para poly, Q2 power: " << npower << " , float norm: " << f_norm << endl;

    TString f_express="1. - [0]*[0]*x/6.";
    if(f_norm>0)f_express="("+f_express;
    for(int p=1;p<npower;p++){
      f_express=f_express+Form("+[%d]",p);
      for(int k=0;k<=p;k++){
        f_express=f_express+"*x";
      }
    }
    if(f_norm>0)f_express=f_express+Form(")*[%d]",npower);
    cout << f_express << endl;

    TF1 *poly = new TF1("poly", f_express, 0, 0.5);
    poly->SetParameter(0, 0.9);
    gr->Fit("poly","0");
    R_c=poly->GetParameter(0);
    R_1=poly->GetParameter(0)+poly->GetParError(0);
    R_2=poly->GetParameter(0)-poly->GetParError(0);
    R_E=(R_1-R_2)/2.;
    chi2fit=poly->GetChisquare();
    Rfit=R_c;
    Rfiterr=R_E;
    cout << "R= " << Rfit << " | Error= " << Rfiterr << " | chi2= " << chi2fit << endl;

    delete poly;
  }else if(id==7){
    cout << "Fit: poly ratio 3 para R" << endl;
    TF1 *poly = new TF1("poly", "[2]*(1.-[0]*[0]*x/6.+[1]*x)/(1.+[1]*x)", 0, 0.5);
    poly->SetParameter(0, 0.9);
    gr->Fit("poly","0");
    R_c=poly->GetParameter(0);
    R_1=poly->GetParameter(0)-poly->GetParError(0);
    R_2=poly->GetParameter(0)+poly->GetParError(0);
    R_E=(R_2-R_1)/2.;
    chi2fit=poly->GetChisquare();
    Rfit=R_c;
    Rfiterr=R_E;
    cout << "R= " << Rfit << " | Error= " << Rfiterr << " | chi2= " << chi2fit << endl;

    delete poly;
  }



  delete gr;
}

void R_fit_class::set_npower(int np, int dnp, int fnorm)
{
  npower=np;
  dnpower=dnp;
  f_norm=fnorm;

  cout << "Set npower= " << npower << " , dnpower= " << dnpower << " , f_norm= " << f_norm << endl;
}

void R_fit_class::set_para(double *parain, int npara)
{
  cout << "Set parameter: " << endl;
  for(int i=0;i<npara;i++){
    param[i]=parain[i];
    cout << i << " " << param[i] << endl;
  }

}

void R_fit_class::get_fit_result(double *fres)
{
  fres[0]=Rfit;
  fres[1]=Rfiterr;
  fres[2]=chi2fit;
}

void R_fit_class::show_data(int nshow)
{
  if(nshow<1)nshow=ndata;
  cout << "show data: " << nshow << endl;
  cout << endl;

  for(int i=0;i<ndata;i++){
    cout << Q2_txt[i] << " " << zero[i] << " " << GE_txt[i] << " " << dGE_txt[i] << " | " << Q2[i] << " " << GE[i] << " " << dGE[i] << endl; 
  }

}
