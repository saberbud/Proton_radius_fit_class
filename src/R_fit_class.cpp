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

void R_fit_class::txt_data_write(string A)
{
  ofstream outfile(A);
  outfile << ndata << endl;
  for(int i=0;i<ndata;i++){
    outfile << Q2[i] << " " << GE[i] << " " << dGE[i] << endl;
  }
  outfile.close();

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
  }else if(id==4){
    cout << "Generate, Kelly: radius input= " << 0.863 << endl;
    double A=1./(25.68*4.*0.938*0.938);
    TF1 *poly = new TF1("poly", "((1.+[0]*x)/(1.+[1]*x+[2]*x*x+[3]*x*x*x))", 0, 10.5);
    double p0=-0.24*A;
    poly->SetParameter(0, p0);
    p0=10.98*A;
    poly->SetParameter(1, p0);
    p0=12.82*A*A;
    poly->SetParameter(2, p0);
    p0=21.97*A*A*A;
    poly->SetParameter(3, p0);
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      GE[i]=poly->Eval(Q2[i]);
      dGE[i]=dGE_txt[i];
    }

    delete poly;
  }else if(id==5){
    cout << "Generate, Arrington: radius input= " << 0.8779 << endl;
    double A=1./(25.68*4.*0.938*0.938);
    TF1 *poly = new TF1("poly", "((1.+[0]*x+[1]*x*x+[2]*x*x*x)/(1.+[3]*x+[4]*x*x+[5]*x*x*x+[6]*x*x*x*x+[7]*x*x*x*x*x))", 0, 10.5);
    double p0=2.90966*A;
    poly->SetParameter(0, p0);
    p0=-1.11542229*A*A;
    poly->SetParameter(1, p0);
    p0=3.866171e-2*A*A*A;
    poly->SetParameter(2, p0);
    p0=14.5187212*A;
    poly->SetParameter(3, p0);
    p0=40.88333*A*A;
    poly->SetParameter(4, p0);
    p0=99.999998*A*A*A;
    poly->SetParameter(5, p0);
    p0=4.579e-5*A*A*A*A;
    poly->SetParameter(6, p0);
    p0=10.3580447*A*A*A*A*A;
    poly->SetParameter(7, p0);
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      GE[i]=poly->Eval(Q2[i]);
      dGE[i]=dGE_txt[i];
    }

    delete poly;
  }else if(id==6){
    cout << "Generate, Arrington inverse poly: radius input= " << 0.8682 << endl;
    double A=1./25.68;
    TF1 *poly = new TF1("poly", "(1./(1.+[0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x+[4]*x*x*x*x*x+[5]*x*x*x*x*x*x))", 0, 10.5);
    double p0=3.226*A;
    poly->SetParameter(0, p0);
    p0=1.508*A*A;
    poly->SetParameter(1, p0);
    p0=-0.3773*A*A*A;
    poly->SetParameter(2, p0);
    p0=0.611*A*A*A*A;
    poly->SetParameter(3, p0);
    p0=-0.1853*A*A*A*A*A;
    poly->SetParameter(4, p0);
    p0=1.596e-2*A*A*A*A*A*A;
    poly->SetParameter(5, p0);
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      GE[i]=poly->Eval(Q2[i]);
      dGE[i]=dGE_txt[i];
    }

    delete poly;
  }else if(id==7){
    cout << "Generate, Arrington-Sick-2007: radius input= " << 0.8965 << endl;
    TString f_express="1. + [0]*x";
    int tnp=5;
    for(int p=1;p<tnp;p++){
      f_express=f_express+Form("/(1.+[%d]*x",p);
    }

    for(int p=1;p<tnp;p++){
      f_express=f_express+")";
    }

    f_express="1./("+f_express+")";
    cout << f_express << endl;

    double A=1./25.68;
    double p0;
    TF1 *poly = new TF1("poly", f_express, 0, 10.5);
    p0=3.440*A;
    poly->SetParameter(0, p0);
    p0=-0.178*A;
    poly->SetParameter(1, p0);
    p0=-1.212*A;
    poly->SetParameter(2, p0);
    p0=1.176*A;
    poly->SetParameter(3, p0);
    p0=-0.284*A;
    poly->SetParameter(4, p0);

    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      GE[i]=poly->Eval(Q2[i]);
      dGE[i]=dGE_txt[i];
    }

    delete poly;
  }else if(id==8){
    // cout << "Generate, YZJA: radius input= " << 0.8792 << endl;
    cout << "Generate, YZJA: radius input= " << 0.85 << endl;

    double tz;
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      tz=z_calc(Q2[i], 2.00252188772, -17.99);
      // GE[i]=0.239163298067 - 1.109858574410*tz + 1.444380813060*pow(tz,2) + 0.479569465603*pow(tz,3) - 2.286894741870*pow(tz,4) + 1.126632984980*pow(tz,5) + 1.250619843540*pow(tz,6) - 3.631020471590*pow(tz,7) + 4.082217023790*pow(tz,8) + 0.504097346499*pow(tz,9) - 5.085120460510*pow(tz,10) + 3.967742543950*pow(tz,11) - 0.981529071103*pow(tz,12);  //0.8792 fm
      GE[i]=0.239448638275 - 1.11264843899*tz + 1.44819766508*pow(tz,2) + 0.514648365159*pow(tz,3) - 2.36672103495*pow(tz,4) + 0.926165750483*pow(tz,5) + 2.05049172945*pow(tz,6) - 4.11073019989*pow(tz,7) + 2.64932410946*pow(tz,8) + 3.51485719222*pow(tz,9) - 7.5760640139*pow(tz,10) + 4.96350589461*pow(tz,11) - 1.14047565701*pow(tz,12);  //0.85 fm
      dGE[i]=dGE_txt[i];
    }

  }else if(id==9){
    cout << "Generate, Bernauer: radius input= " << 0.8868 << endl;

    double tz;
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      tz=Q2_txt[i]/25.68;
      GE[i]=1. - 3.36591660*tz + 1.45487683e+01*pow(tz,2) - 8.87959239e+01*pow(tz,3) + 4.61097705e+02*pow(tz,4) - 1.67562381e+03*pow(tz,5) + 4.07646487e+03*pow(tz,6) - 6.45411460e+03*pow(tz,7) + 6.34035079e+03*pow(tz,8) - 3.49373923e+03*pow(tz,9) + 8.22601568e+02*pow(tz,10);  //0.8868 fm
      dGE[i]=dGE_txt[i];
    }

  }else if(id==10){
    cout << "Generate, YZJA: radius input= " << 0.8792 << endl;

    double tz;
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      tz=z_calc(Q2[i], 2.00252188772, -17.99);
      GE[i]=0.239163298067 - 1.109858574410*tz + 1.444380813060*pow(tz,2) + 0.479569465603*pow(tz,3) - 2.286894741870*pow(tz,4) + 1.126632984980*pow(tz,5) + 1.250619843540*pow(tz,6) - 3.631020471590*pow(tz,7) + 4.082217023790*pow(tz,8) + 0.504097346499*pow(tz,9) - 5.085120460510*pow(tz,10) + 3.967742543950*pow(tz,11) - 0.981529071103*pow(tz,12);  //0.8792 fm
      dGE[i]=dGE_txt[i];
    }

  }else if(id==11){
    cout << "Generate, Bernauer orig: radius input= " << 0.8872 << endl;

    double tz;
    for(int i=0;i<ndata;i++){
      Q2[i]=Q2_txt[i];
      tz=Q2_txt[i]/25.68;
      GE[i]=1. - 3.3686*tz + 14.5606*pow(tz,2) - 88.1912*pow(tz,3) + 453.6244*pow(tz,4) - 1638.7911*pow(tz,5) + 3980.7174*pow(tz,6) - 6312.6333*pow(tz,7) + 6222.3646*pow(tz,8) - 3443.2251*pow(tz,9) + 814.4112*pow(tz,10);  //0.8868 fm
      dGE[i]=dGE_txt[i];
    }

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
  }else if(id==5){
    cout << "Add scaling noise, type 5" << endl;
    a0=1.;
    a1=param[0];
    temp=rdm->Gaus(a0, a1);
    cout << "Scaling: " << temp << endl;
    for(int i=0;i<ndata;i++){
      GE[i]=GE[i]*temp;
    }
  }else if(id==6){
    cout << "Add mannual shift, type 6" << endl;
    a0=param[0];
    a1=param[1];
    for(int i=0;i<ndata;i++){
      if(i<a0){
        GE[i]=GE[i]+a1*dGE_txt[i];
      }else{
        GE[i]=GE[i]-a1*dGE_txt[i];
      }
    }
  }

  delete rdm;
}

void R_fit_class::z_trans()
{
  double num,den;
  for(int i=0;i<ndata;i++){
    num=sqrt(Tc+Q2[i])-sqrt(Tc);
    den=sqrt(Tc+Q2[i])+sqrt(Tc);
    z_arr[i]=num/den;
  }
}

double R_fit_class::z_calc(double Q2_in, double Tc_in, double T0_in)
{
  double res=0.;
  double num=sqrt(Tc_in+Q2_in)-sqrt(Tc_in-T0_in);
  double den=sqrt(Tc_in+Q2_in)+sqrt(Tc_in-T0_in);
  res=num/den;

  return res;
}

void R_fit_class::GE_z_fit(int id)
{
  cout << "GE_z_fit id= " << id << endl;

  chi2fit=0.;Rfit=0.;Rfiterr=0.;

  TGraphErrors *gr=new TGraphErrors(ndata,z_arr,GE,zero,dGE);

  double R_c,R_1,R_2,R_E;
  if(id==0){
    cout << "Fit: multiple para poly, z power: " << npower << " , float norm: " << f_norm << endl;

    TString f_express=Form("1. - [0]*[0]*x*%.4f*2./3.",Tc);
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
  }else if(id==1){
    cout << "Fit: poly z ratio 3 para R" << endl;
    TString f_express=Form("[2]*(1.-[0]*[0]*x*%.4f*2./3.+[1]*x)/(1.+[1]*x)",Tc);
    cout << f_express << endl;
    TF1 *poly = new TF1("poly", f_express, 0, 0.5);
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
  }else if(id==8){
    cout << "Fit: poly Q4 with fixed Q6, Q8" << endl;
    TString f_express=Form("[2]*(1.-[0]*[0]*x/6.+[1]*x*x-%.9f*x*x*x+%.9f*x*x*x*x)",param[0],param[1]);
    cout << f_express << endl;
    TF1 *poly = new TF1("poly", f_express, 0, 0.5);
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
  }else if(id==9){
    cout << "Fit: continued fraction, order: " << npower << " , float norm: " << f_norm << endl;
    TString f_express="1. + [0]*[0]*x/6.";
    for(int p=1;p<npower;p++){
      f_express=f_express+Form("/(1.+[%d]*x",p);
    }

    for(int p=1;p<npower;p++){
      f_express=f_express+")";
    }

    if(f_norm>0){
      f_express=Form("[%d]/(",npower)+f_express+")";
    }else{
      f_express="1./("+f_express+")";
    }
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
  }else if(id==10){
    cout << "Fit: poly ratio np=" << npower << " ,dnpower= " << dnpower << " , float norm: " << f_norm << endl;

    int itemp;
    TString fn_express="1.-[0]*[0]*x/6.+[1]*x";
    for(int p=1;p<npower;p++){
      itemp=p+1;
      fn_express=fn_express+Form("+[%d]",itemp);
      for(int k=0;k<=p;k++){
        fn_express=fn_express+"*x";
      }
    }
    // cout << fn_express << endl;

    TString fd_express="1.+[1]*x";
    for(int p=1;p<dnpower;p++){
      itemp=p+npower;
      fd_express=fd_express+Form("+[%d]",itemp);
      for(int k=0;k<=p;k++){
        fd_express=fd_express+"*x";
      }
    }
    // cout << fd_express << endl;

    TString f_express="("+fn_express+")/("+fd_express+")";
    if(f_norm>0)f_express=Form("[%d]*",npower+dnpower)+f_express;
    cout << f_express << endl;

    TF1 *poly = new TF1("poly", f_express, 0, 0.5);
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

void R_fit_class::set_Tc(double Tcin)
{
  Tc=Tcin;

  cout << "Set Tc= " << Tc << endl;
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
