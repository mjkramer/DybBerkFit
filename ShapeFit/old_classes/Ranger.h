#ifndef RANGER_H
#define RANGER_H
#include <map>
#include "TObject.h"

class Ranger : public TObject {

 public:
  Ranger();
  ~Ranger(){};
  //Note: index runs from 0 to (nsteps-1); 
  double returnVal(int index);
  int findIndex(double val);
  void setLogScale();
  bool getLogScale(){return logScale;}
  
 public:
  int nsteps;
  double min;
  double max;

 private:
  bool logScale;
  bool useMap;
  
 public:
  ClassDef(Ranger,1);

};
#endif
