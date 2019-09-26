// Compile and load fitter classes
{
  TString include = ".include ";
  TString load = ".L ";

  TString prefix = "../ShapeFit";
  gROOT->ProcessLine( (include + prefix).Data() ); 
  gROOT->ProcessLine( (load + prefix + "/OscCalc.cc+").Data() );
  gROOT->ProcessLine( (load + prefix + "/FluxCalculator.cc+").Data() );
  gROOT->ProcessLine( (load + prefix + "/PredSet.cc+").Data() );
  gROOT->ProcessLine( (load + prefix + "/TimePeriodData.cc+").Data() );
  gROOT->ProcessLine( (load + prefix + "/Predictor.cc+").Data() );
  gROOT->ProcessLine( (load + prefix + "/Binning.cc+").Data() );
  gROOT->ProcessLine( (load + prefix + "/DataSet.cc+").Data() );

  // TString prefix = "./quickFit";  
  // gROOT->ProcessLine( (include + prefix).Data() ); 
  //gROOT->ProcessLine( (load + prefix + "/FitSummary.C+").Data() );
  //gROOT->ProcessLine( (load + prefix + "/PhysicsModel.C+").Data() );
  //gROOT->ProcessLine( (load + prefix + "/PoissonModel.C+").Data() );
  //gROOT->ProcessLine( (load + prefix + "/QuickFitter.C+").Data() );

  TString prefix = "./reactor";
  gROOT->ProcessLine( (include + prefix).Data() ); 
  gROOT->ProcessLine( (load + prefix + "/IsotopeTable.C+").Data() );
  gROOT->ProcessLine( (load + prefix + "/CrossSectionTable.C+").Data() );
  gROOT->ProcessLine( (load + prefix + "/CoreSpectrum.C+").Data() );
  gROOT->ProcessLine( (load + prefix + "/Spectrum.C+").Data() );
  //gROOT->ProcessLine( (load + prefix + "/ChiSquare.C+").Data() );
  //gROOT->ProcessLine( (load + prefix + "/ChiSquareEnergyUnc.C+").Data() );

  
 
}
