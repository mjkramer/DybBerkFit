void run_fit_shape_2d_P17B_IHEPdata_dm221(Int_t input = -1){

  gROOT->ProcessLine(".x LoadClasses.C");
  gROOT->ProcessLine(".L fit_shape_2d_P17B_IHEPdata_dm221.C+");

    double s22t13_min=0.06;
    double s22t13_max=0.10;
    const int s22t13_points=31;
    
    double dm221_min=0.0;
    double dm221_max=0.0003; //30x10^-5 eV^2
    const int dm221_points=31;
    
    double s2t_array[s22t13_points*dm221_points];
    double dm2_array[s22t13_points*dm221_points];
    int counter=0;
    
    for(int is2t=0; is2t<s22t13_points; ++is2t){
        for(int idm2=0; idm2<dm221_points; ++idm2){
            s2t_array[counter]=s22t13_min+is2t*(s22t13_max-s22t13_min)/(s22t13_points-1);
            dm2_array[counter]=dm221_min+idm2*(dm221_max-dm221_min)/(dm221_points-1);
            counter++;
        }
    }
    
  Char_t outputfilename[256];
 sprintf(outputfilename,"./fit_result_files/dm221_fit/fit_dmee_constrained_t2k_dm2ee_s22t12_%.4f_dm221_%.6f.root",s2t_array[input],dm2_array[input]);
 
    cout<<"Running input="<<input<<endl;
    cout<<"s22t12 is "<<s2t_array[input]<<endl;
    cout<<"dm221 is "<<dm2_array[input]<<endl;
    cout<<"outputfile is "<<outputfilename<<endl;
    if(dm2_array[input]==0) dm2_array[input]=0.00000001;
 fit_shape_2d_P17B_IHEPdata_dm221(s2t_array[input],dm2_array[input],outputfilename);
 
}
