// Compile and load fitter classes
{
  TString include = ".include ";
  TString load = ".L ";

  TString suffix = gSystem->Getenv("LBNL_FIT_DEBUG") ? "g" : "";

  TString prefix = "../ShapeFit";
  gROOT->ProcessLine( (include + prefix).Data() ); 
  gROOT->ProcessLine( (load + prefix + "/OscCalc.cc+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/FluxCalculator.cc+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/PredSet.cc+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/TimePeriodData.cc+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/Predictor.cc+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/Binning.cc+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/DataSet.cc+" + suffix).Data() );

  // TString prefix = "./quickFit";  
  // gROOT->ProcessLine( (include + prefix).Data() ); 
  //gROOT->ProcessLine( (load + prefix + "/FitSummary.C+" + suffix).Data() );
  //gROOT->ProcessLine( (load + prefix + "/PhysicsModel.C+" + suffix).Data() );
  //gROOT->ProcessLine( (load + prefix + "/PoissonModel.C+" + suffix).Data() );
  //gROOT->ProcessLine( (load + prefix + "/QuickFitter.C+" + suffix).Data() );

  TString prefix = "./reactor";
  gROOT->ProcessLine( (include + prefix).Data() ); 
  gROOT->ProcessLine( (load + prefix + "/IsotopeTable.C+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/CrossSectionTable.C+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/CoreSpectrum.C+" + suffix).Data() );
  gROOT->ProcessLine( (load + prefix + "/Spectrum.C+" + suffix).Data() );
  //gROOT->ProcessLine( (load + prefix + "/ChiSquare.C+" + suffix).Data() );
  //gROOT->ProcessLine( (load + prefix + "/ChiSquareEnergyUnc.C+" + suffix).Data() );

  
 
}
