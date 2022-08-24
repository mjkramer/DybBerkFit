#pragma once

// Generic Data Set (for comparison with a physics model)
//
// Created by: dandwyer@caltech.edu 2012/02/01

#include "TMap.h"

class TH1F;
class TGraph;
class TObject;

class DataSet {
public:
  DataSet();
  virtual ~DataSet();
  // Add methods to get and set common data types
  double getDouble(const char* name);
  const char* getString(const char* name);
  TH1F* getHistogram(const char* name);
  TGraph* getGraph(const char* name);
  TObject* getObject(const char* name);
  // DataTable* getTable(const char* name);
  bool isDouble(const char* name);

  void setDouble(const char* name, double value);
  void setString(const char* name, const char* value);
  void setHistogram(const char* name, TH1F* value);
  void setGraph(const char* name, TGraph* value);
  void setObject(const char* name, TObject* value);
  // void setTable(const char* name, DataTable* value);

  // Load data set from a text file
  int load(const char* filename);
  // Print contents of data set to screen
  void dump();

private:
  void maybeSetNominal(const char* varname, double value);
  TMap m_data; // Container for data
  bool m_warnOnReplace = true;
};
