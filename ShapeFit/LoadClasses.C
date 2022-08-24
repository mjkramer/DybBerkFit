void LoadFile(const char* fname)
{
  const char* suffix = gSystem->Getenv("LBNL_FIT_DEBUG") ? "g" : "";
  const char* line = Form(".L %s+%s", fname, suffix);
  gROOT->ProcessLine(line);
}

void LoadClasses()
{
  gROOT->Macro("EnableOpenMP.C");

  LoadFile("Utils.cc");
  LoadFile("Paths.cc");
  LoadFile("Binning.cc");
  LoadFile("Config.cc");
  LoadFile("OscCalc.cc");
  LoadFile("FluxCalculator.cc");
  LoadFile("PredSet.cc");
  LoadFile("TimePeriodData.cc");
  LoadFile("Predictor.cc");
  LoadFile("Ranger.cc");
  // LoadFile("ExtrapTable.cc");
  LoadFile("OscProbTable.cc");
  LoadFile("DataSet.cc");
}
