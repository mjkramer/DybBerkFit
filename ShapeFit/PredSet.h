#include <iostream>
#include "TObject.h"
#include "FluxCalculator.h"
#include "TH1.h"

using namespace std;

class PredSet : public TObject {

 public:
  PredSet();
  ~PredSet();
  void SetPred(int istage, int idet1, int idet2, TH1F* val);
  TH1F* GetPred(int istage, int idet1, int idet2);
  void PrintToScreen(int istage);

  
 private:
  TH1F* PredEvts[Nstage][Ndetectors][Ndetectors];
  TH1F *h_comb;

 public:
  ClassDef(PredSet,1);
 
};

void PredSet::SetPred(int istage, int idet1, int idet2, TH1F* val){
  if(idet1<0 || idet1>=Ndetectors){
    std::cout << "WARNING!! idet1 indices in PredSet are out of range" << std::endl;
  } 
  if(idet2<0 || idet2>=Ndetectors){
    std::cout << "WARNING!! idet2 indices in PredSet are out of range" << std::endl;
  }
  val->Copy(*PredEvts[istage][idet1][idet2]);
  
}

TH1F* PredSet::GetPred(int istage, int idet1, int idet2){
  if(idet1<0 || idet1>=Ndetectors){
    std::cout << "WARNING!! idet1 indices in PredSet are out of range" << std::endl;
  } 
  if(idet2<0 || idet2>=Ndetectors){
    std::cout << "WARNING!! idet2 indices in PredSet are out of range" << std::endl;
  }
  return PredEvts[istage][idet1][idet2]; 
  
}

void PredSet::PrintToScreen(int istage){
  for(int idet1=4;idet1<8;++idet1){
    for(int idet2=0;idet2<4;++idet2){
      cout << "AD" << idet1+1 << " from AD" << idet2+1 << ": " << PredEvts[istage][idet1][idet2]->Integral() << endl;
    }
  }
}
