// Instructions: Call SmearSpectrum() after loading macro into interpreter (i.e. .L Li9Spec.C)
// Note: Fermi function corrections are negligible. Disabled by default.

#include <TF1.h>
#include <TAxis.h>

#define Z 4			// Be
#define m_e 0.511
#define N_LVL 4

double Spec(double *pT, double *par)
{
  double T = *pT;
  double norm = par[0];
  double Q = par[1];

  double p = sqrt(T*T + 2*m_e*T);
  double E = T + m_e;

  return T>Q ? 0 : norm * p * E * (Q-T)*(Q-T);
}

// 1985 J. Phys. G: Nucl. Phys. 11 359
double FermiSpec(double *pT, double *par)
{
  double W = (*pT+m_e)/m_e;
  double ferm = sqrt(1.232 + 0.0302/(W-1));
  return ferm * Spec(pT, par);
}

double FullSpec(double *pT, double *par)
{
  double sum = 0;

  for (int i = 0; i < N_LVL; ++i)
    sum += Spec(pT, par+2*i);

  return sum;
}

double FullFermiSpec(double *pT, double *par)
{
  double sum = 0;

  for (int i = 0; i < N_LVL; ++i)
    sum += FermiSpec(pT, par+2*i);

  return sum;
}

void SetQ(TF1 *f, double Q)
{
  f->SetParameters(1, Q);
  f->SetParameter(0, 1/f->Integral(0, Q));
  f->SetRange(0, Q);
}

TF1* GetSpectrum(bool fermi=false)
{
  double (*fun)(double*, double*) = fermi ? FermiSpec : Spec;
  TF1 spec("spec", fun, 0, 0, 2);

  double Q[N_LVL] = {11.176, 10.826, 1.976, 1.446};
  double br[N_LVL] = {30, 16, 1.1, 2.7};

  double sum_br = 0;
  for (int i = 0; i < N_LVL; ++i) sum_br += br[i];
  for (int i = 0; i < N_LVL; ++i) br[i] /= sum_br;

  double (*fullfun)(double*, double*) = fermi ? FullFermiSpec : FullSpec;
  TF1 *full = new TF1(fermi ? "ffullspec" : "fullspec",
                      fullfun, 0, Q[0], 2*N_LVL);

  for (int i = 0; i < N_LVL; ++i) {
    spec.SetParameters(1, Q[i]);
    full->SetParameter(2*i, br[i] / spec.Integral(0, Q[i]));
    full->SetParameter(2*i+1, Q[i]);
  }

  full->SetTitle("Theoretical Li9 beta spectrum, four branches, no alpha/n");
  full->GetXaxis()->SetTitle("E/MeV");
  full->Draw();

  return full;
}


void SmearSpectrum(){

  Double_t e_min = 0;
  Double_t e_max = 15;
  const Int_t nsteps = 3000;

  //rebinned histogram for fit 
  const Int_t n_evis_bins = 37;
  Double_t evis_bins[38]; // Single bins between 0.7 and 1.0 MeV. 0.2 MeV bins from 1.0 to 8.0 MeV. Single bin between 8.0 and 12 MeV. total 37 bins
  evis_bins[0] = 0.7;
  for (Int_t i = 0; i < 36; i++){
    evis_bins[i+1] = 0.2 *i + 1.0;
  }
  evis_bins[37] = 12.0;

  //rebinned histogram for toy 
  const Int_t n_evis_bins_toy = 240;
  TH1D *h_li9_smeared_toy = new TH1D("h_li9_smeared_toy","h_li9_smeared_toy",n_evis_bins_toy,0,12);
  Double_t evis_bins_toy[241]; // The toy spectra currently spans 0-12MeV in 240 bins
  for (Int_t i = 0; i <= n_evis_bins_toy; i++){
    evis_bins_toy[i] = 0.05 *i;
  }
  
  Double_t e_step = (e_max - e_min)/(Double_t)nsteps;

  Double_t li9_smeared[nsteps];
  
  for (Int_t i = 0; i < nsteps; i++){
    li9_smeared[i] = 0;
  }

  TH1D * h_li9_smeared = new TH1D("h_li9_smeared","h_li9_smeared",nsteps,e_min,e_max);
  
  TF1 * true_spec = GetSpectrum();

  for (Int_t i = 0; i < nsteps; i++){
    //  for (Int_t i = 2000; i < 2001; i++){
    Double_t e = (i+0.5) * e_step;
    Double_t amp = true_spec->Eval(e);
    // if (i == 2000){
    //   amp = 1;
    // }else{
    //   amp = 0;
    // }
    
    Double_t sigma = e * 0.01 * (7.5/sqrt(e) + 0.9);
    //    cout << e << " " << amp << " " << sigma << endl;
    for (Int_t j = 0; j < nsteps; j++){
      Double_t x = (j+0.5) * e_step;
      Double_t w = amp * e_step * exp(-pow((x-e),2)/pow(sigma,2)/2.)/sqrt(2 * TMath::Pi())/sigma;
      li9_smeared[j] += w;
      //      cout << x << " " << w << endl;
    }
    //    cout << "-----------------------------------" << endl;
    cout << ".";
    cout.flush();
  }

  //When filling h_li9_smeared also fill spectrum for toy; if do same rebin as for fit spectra get ugly results
  for (Int_t i = 0; i < nsteps; i++){
    h_li9_smeared->SetBinContent(i+1,li9_smeared[i]);
  }
  h_li9_smeared->Draw();
  true_spec->SetLineStyle(2);
  true_spec->Draw("same");

  //li9 spectrum for fit 
  TH1D * h_li9_smeared_rebin = (TH1D*)h_li9_smeared->Rebin(n_evis_bins,"h_li9_smeared_rebin",&evis_bins[0]);
  
  //li9 spectrum for toy 
  TH1D * h_li9_smeared_toy = (TH1D*)h_li9_smeared->Rebin(n_evis_bins_toy,"h_li9_smeared_toy",&evis_bins_toy[0]);
  h_li9_smeared_toy->GetXaxis()->SetTitle("E_{prompt} (MeV)");
  h_li9_smeared_toy->GetYaxis()->SetTitle("arbitrary units");
  h_li9_smeared_toy->SetTitle("");
  
  /*for(int ibin=1;ibin<=h_li9_smeared->GetXaxis()->GetNbins();++ibin){
    h_li9_smeared_toy->Fill(h_li9_smeared->GetBinCenter(ibin),h_li9_smeared->GetBinContent(ibin));
    }*/
  //-->set bin content to zero for bins with energy below 0.7 MeV (or else when doing shape distorsions the total number of events is not the same before and after, since the histos for fit start at 0.7)
  for(int ibin=0;ibin<h_li9_smeared_toy->GetXaxis()->GetNbins();++ibin){
    double x = h_li9_smeared_toy->GetBinCenter(ibin+1);
    if(x<0.7) h_li9_smeared_toy->SetBinContent(ibin+1,0);
  }

  //saving to file
  TFile * fout = new TFile("li9_spectrum.root","recreate");
  true_spec->Write();
  h_li9_smeared->Write();
  h_li9_smeared_rebin->Write();
  h_li9_smeared_toy->Write();
  fout->Close();
  
}
