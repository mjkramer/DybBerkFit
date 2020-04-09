#include "DataSet.h"

#include "TGraph.h"
#include "TH1F.h"
#include "TObjString.h"
#include "TObject.h"
#include "TParameter.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

DataSet::DataSet() : m_data()
{
  ;
}

DataSet::~DataSet()
{
  ;
}

double DataSet::getDouble(const char* name)
{
  TParameter<double>* par =
      dynamic_cast<TParameter<double>*>(this->getObject(name));
  if (!par) {
    cerr << "Failed to find \"" << name << "\" (double)" << endl;
    return 0;
  }
  return par->GetVal();
}


const char* DataSet::getString(const char* name)
{
  TObjString* par = dynamic_cast<TObjString*>(this->getObject(name));
  if (!par) {
    cerr << "Failed to find \"" << name << "\" (string)" << endl;
    return 0;
  }
  return par->String().Data();
}

TH1F* DataSet::getHistogram(const char* name)
{
  return dynamic_cast<TH1F*>(this->getObject(name));
}

TGraph* DataSet::getGraph(const char* name)
{
  return dynamic_cast<TGraph*>(this->getObject(name));
}

TObject* DataSet::getObject(const char* name)
{
  return m_data.GetValue(name);
}

/*
  DataTable* DataSet::getTable(const char* name)
  {
  }
*/

void DataSet::setDouble(const char* name, double value)
{
  // Create new data object
  TParameter<double>* par = new TParameter<double>(name, value);
  this->setObject(name, par);
}


void DataSet::setString(const char* name, const char* value)
{
  // Create new data object
  TObjString* par = new TObjString(value);
  this->setObject(name, par);
}

void DataSet::setHistogram(const char* name, TH1F* value)
{
  this->setObject(name, value);
}

void DataSet::setGraph(const char* name, TGraph* value)
{
  this->setObject(name, value);
}

void DataSet::setObject(const char* name, TObject* value)
{
  TObject* pair = m_data.FindObject(name);
  if (pair) {
    // This data object exists, so replace
    cout << "Warning: Replacing object: " << name << endl;
    m_data.Remove(((TPair*)pair)->Key());
    // FIXME: Consider deleting object if there are no other owners...
  }
  // Add new data
  m_data.Add(new TObjString(name), value);
}

/*
  void DataSet::setTable(const char* name, DataTable* value)
  {
  }
*/

int DataSet::load(const char* filename)
{
  // Load a data set from a text file
  ifstream fileData(filename);
  if (!fileData.is_open() || !fileData.good()) {
    std::cout << "DataSet::load: "
              << "Error: Failed to open data file " << filename << std::endl;
    return -1;
  }
  string line;
  while (true) {
    getline(fileData, line);
    if (!fileData.good())
      break;
    if (line.empty())
      continue;
    istringstream lineStr(line);
    if (lineStr.peek() == '#') {
      // Skip lines starting with '#'
      continue;
    }
    string name;
    if (!(lineStr >> name)) {
      std::cout << "DataSet::load: "
                << "Error: Failed to read name from line " << line << std::endl;
      return -1;
    }
    string svalue;
    if (!(lineStr >> svalue)) {
      std::cout << "DataSet::load: "
                << "Error: Failed to read value from line " << line
                << std::endl;
      return -1;
    }
    bool isNumber = false;
    double dvalue = 0;
    istringstream testStream(svalue);
    if (testStream >> dvalue)
      isNumber = true;

    if (isNumber) {
      this->setDouble(name.c_str(), dvalue);
    } else {
      this->setString(name.c_str(), svalue.c_str());
    }
  }
  return 0;
}

void DataSet::dump()
{
  m_data.Print();
  return;
}
