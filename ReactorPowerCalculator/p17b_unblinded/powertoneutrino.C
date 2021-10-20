{
  const int N = 500;
  int maxperiod = 0, period, core, start[N], end[N];
  double anitnu_rate_per_period[N][6], uncert, u235[N][6], u238[N][6], pu239[N][6], pu241[N][6];
  string line;

  ifstream infile("weekly_nuRate_p17b.txt");
  
  while(1) {
    if (infile.peek()=='#') {
      getline(infile, line);
      continue;
    } else {
      infile >> period >> core >> start[period] >> end[period]
	     >> anitnu_rate_per_period[period][core-1] >> uncert >> u235[period][core-1]
	     >> u238[period][core-1] >> pu239[period][core-1]
	     >> pu241[period][core-1];
      maxperiod = period;
      if (!infile.good()) break;
    }
  }
  infile.close();
  infile.clear();
  
  maxperiod += 1;

  double energyPerFission[4] = {202.36, 205.99, 211.12, 214.26};
  
  double eNu[2000], dNdE[4][2000];

  infile.open("fissionIsotopeSpectra_Huber_Fit.txt");
  int i = 0, maxI;

  while(1) {
    if (infile.peek()=='#') {
      getline(infile, line);
      continue;
    } else {
      infile >> eNu[i] >> dNdE[0][i] >> dNdE[1][i] >> dNdE[2][i] >> dNdE[3][i];
      i++;
      maxI = i;
    }
    if (!infile.good()) break;
  }
  
  infile.close();
  infile.clear();

  maxI -= 1;
  
  i=0;
  maxI=0;
  double dxdE[2000];

  infile.open("inverseBetaCrossSection_v2.txt");
  while(1) {
    if (infile.peek()=='#') {
      getline(infile, line);
      continue;
    } else {
      infile >> eNu[i] >> dxdE[i];
      i++;
      maxI = i;
    }
    if (!infile.good()) break;
  }
  infile.close();
  infile.clear();

  maxI -= 1;

  double antinuPerFission[4] = {0,0,0,0};
  double ibdPerFission[4] = {0, 0, 0, 0};

  double dEnu = 0.01;

  for (int j=0; j<4; j++) {
    for (int i=0; i<maxI; i++) {
      antinuPerFission[j] += dNdE[j][i] * dEnu;
      ibdPerFission[j] += dNdE[j][i] * dxdE[i] * dEnu;
    }
    //cout << j << " " << antinuPerFission[j] << endl;
  }


  cout << "#  This is the unblinded weekly reactor neutrino rate for P17B on Apr 26, 2018" << endl;
  cout << "#week reactor timestart timeend nusRate error u235 u238 pu239 pu241" << endl;
  
  double power[N][6];
  double ibdRatePerPeriod[N][6];

  for (int i=0; i<maxperiod; i++) {
    for (int j=0; j<6; j++) {

      double avg_energyPerFission = 
	energyPerFission[0] * u235[i][j]
	+ energyPerFission[1] * u238[i][j]
	+ energyPerFission[2] * pu239[i][j]
	+ energyPerFission[3] * pu241[i][j];
      
      double avg_antinuPerFission = 
	antinuPerFission[0] * u235[i][j] 
	+ antinuPerFission[1] * u238[i][j]
	+ antinuPerFission[2] * pu239[i][j]
	+ antinuPerFission[3] * pu241[i][j]; 

      double avg_ibdPerFission = 
	ibdPerFission[0] * u235[i][j]
	+ ibdPerFission[1] * u238[i][j]
	+ ibdPerFission[2] * pu239[i][j]
	+ ibdPerFission[3] * pu241[i][j]; 
	
      /*antinuRatePerPeriod[i][j] = power[i][j]
	* 2.895 * 6.2415e21 / avg_energyPerFission * avg_antinuPerFission;*/
        
        
        power[i][j]=anitnu_rate_per_period[i][j]*avg_energyPerFission/(2.895 * 6.2415e21*avg_antinuPerFission);
      
     /* ibdRatePerPeriod[i][j] = power[i][j]
	* 2.895 * 6.2415e21 / avg_energyPerFission * avg_ibdPerFission;*/
      
      cout << i << " " << j+1 << " " << start[i] << " " << end[i]
	   << " " << power[i][j]  << " "
	   << " 0 " << u235[i][j] << " " 
	   << u238[i][j] << " " << pu239[i][j] << " " << pu241[i][j] << endl;
    }    
  }
    
    FILE* outputfile=fopen("./WeeklyAvg_P17B_by_Beda.txt","w");
    for (int i=0; i<maxperiod; i++) {
        for (int j=0; j<6; j++) {
            fprintf(outputfile,"%d %d %d %d %f 0 %f %f %f %f\n",i, j+1,start[i],end[i],power[i][j],u235[i][j],u238[i][j],pu239[i][j] , pu241[i][j]);
        }
    }
}
