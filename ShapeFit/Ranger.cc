#include "Ranger.h"
Ranger::Ranger(){
  nsteps=0;
  min=-1;
  max=-1;
  logScale = false;
  useMap = false;
};

void Ranger::setLogScale(){
  min = log10(min);
  max = log10(max);
  logScale = true;
}

double Ranger::returnVal(int index){
  double val=min;
  if(nsteps>1){
    val=(max-min)*index*1./(nsteps-1)+min;
  }
  if (logScale) val = pow(10,val);
  return val;   
};
int Ranger::findIndex(double val){
  
  double best=min;
  if (logScale) val = log10(val);
  int ibest=0;
  if(max!=min){
    best=(val-min)*(nsteps-1)*1./(max-min);
      ibest=(int)floor(best+0.5);
  }
  return ibest;
  
};

