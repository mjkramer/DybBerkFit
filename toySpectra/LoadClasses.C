// Compile and load fitter classes
{
  gROOT->Macro("../ShapeFit/EnableOpenMP.C");

  TString include = ".include ";
  TString load = ".L ";

  TString prefix;
  TString suffix = gSystem->Getenv("LBNL_FIT_DEBUG") ? "g" : "";

  prefix = "../ShapeFit";
  gROOT->ProcessLine((include + prefix).Data());
  gROOT->ProcessLine((load + prefix + "/Paths.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/Utils.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/OscCalc.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/FluxCalculator.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/PredSet.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/TimePeriodData.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/Predictor.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/Binning.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/Config.cc+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/DataSet.cc+" + suffix).Data());

  // prefix = "./quickFit";
  // gROOT->ProcessLine( (include + prefix).Data() );
  // gROOT->ProcessLine( (load + prefix + "/FitSummary.C+" + suffix).Data() );
  // gROOT->ProcessLine( (load + prefix + "/PhysicsModel.C+" + suffix).Data() );
  // gROOT->ProcessLine( (load + prefix + "/PoissonModel.C+" + suffix).Data() );
  // gROOT->ProcessLine( (load + prefix + "/QuickFitter.C+" + suffix).Data() );

  prefix = "./reactor";
  gROOT->ProcessLine((include + prefix).Data());
  gROOT->ProcessLine((load + prefix + "/IsotopeTable.C+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/CrossSectionTable.C+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/CoreSpectrum.C+" + suffix).Data());
  gROOT->ProcessLine((load + prefix + "/Spectrum.C+" + suffix).Data());
  // gROOT->ProcessLine( (load + prefix + "/ChiSquare.C+" + suffix).Data() );
  // gROOT->ProcessLine( (load + prefix + "/ChiSquareEnergyUnc.C+" +
  // suffix).Data() );
}
