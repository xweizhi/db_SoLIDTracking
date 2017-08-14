//generate the db_gemc.dat file for GEM digitization
//just use Root cint 

#include <iostream>
#include <fstream>
#include <iomanip>
#include "TVector2.h"
#include "TMath.h"

using namespace std;

int main()
{
  const int NPLANE = 6;
  const int NSECTOR = 30;
  const double DPHI = 12;
  const double DEPTH = 0.05;
  const double PITCH = 0.0004;
  const double UANGLE = 84.0;
  const double VANGLE = 96.0;
  const double FRAMEWIDTH = 0.0;
  const double Z0[NPLANE]     = {-1.75, -1.50, -1.19, -0.68,  0.05,  0.92};
  const double RIN[NPLANE]    = { 0.36,  0.21,  0.25,  0.32,  0.42,  0.55};
  const double ROUT[NPLANE]   = { 0.87,  0.98,  1.12,  1.35,  1.00,  1.23};
  const double PHIOFF[NPLANE] = {0.};
  const double PHISTART = -6;
  ofstream db_file;
  
  int count = 1;
  db_file.open ("db_gemc.dat");
  db_file<<"950101\n";
  for (int i=0; i<NPLANE; i++){
    double phi0 = PHISTART;
    for (int j=0; j<NSECTOR; j++){
      phi0 = PHISTART + j*DPHI + PHIOFF[i]; //same definition in GEMC
      db_file<<"gemc.gem"<<count<<".r0 = "<<setprecision(6)<<fixed<<RIN[i]<<endl;
      db_file<<"gemc.gem"<<count<<".r1 = "<<setprecision(6)<<fixed<<ROUT[i]<<endl;
      db_file<<"gemc.gem"<<count<<".phi0 = "<<setprecision(6)<<fixed<<phi0<<endl;
      db_file<<"gemc.gem"<<count<<".dphi = "<<setprecision(0)<<fixed<<DPHI<<endl;
      db_file<<"gemc.gem"<<count<<".z0 = "<<setprecision(6)<<fixed<<Z0[i]<<endl;
      db_file<<"gemc.gem"<<count<<".depth = "<<setprecision(2)<<fixed<<DEPTH<<endl;
      db_file<<"gemc.gem"<<count<<".frame_width = "<<setprecision(6)<<fixed<<FRAMEWIDTH<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"x.stripangle = "<<setprecision(1)<<fixed<<UANGLE<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"x.pitch = "<<setprecision(4)<<fixed<<PITCH<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"y.stripangle = "<<setprecision(1)<<fixed<<VANGLE<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"y.pitch = "<<setprecision(4)<<fixed<<PITCH<<endl;
      db_file<<"gemc.gem"<<count<<".n_HV_sector_off = "<<"0"<<endl;
      db_file<<endl;
      count++;
    }
  }
  db_file.close();
  
  return 1;
  
}
