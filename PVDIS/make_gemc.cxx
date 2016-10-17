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
  const double toDEG = 180./TMath::Pi();
  const double toRAD = TMath::Pi()/180.;
  const int NPLANE = 5;
  const int NSECTOR = 30;
  const double DPHI = 12;
  const double DEPTH = 0.05;
  const double PITCH = 0.0004;
  const double UANGLE = 84.0;
  const double VANGLE = 96.0;
  const double Z0[NPLANE]     = {1.575, 1.855, 1.900, 3.060,  3.150 };
  const double RIN[NPLANE]    = {0.48,  0.59,  0.65,  1.05,   1.09  };
  const double ROUT[NPLANE]   = {1.22,  1.43,  1.43,  2.30,   2.37  };
  const double PHIOFF[NPLANE] = {0.5,   0.0,   0.0,   -0.5,   -0.5  };
  const double PHISTART = -96;
  ofstream db_file;
  
  int count = 1;
  db_file.open ("db_gemc.dat");
  db_file<<"950101\n";
  for (int i=0; i<NPLANE; i++){
    double phi0 = PHISTART;
    for (int j=0; j<NSECTOR; j++){
      phi0 = TVector2::Phi_0_2pi((-1.*(PHISTART + j*DPHI) + PHIOFF[i])*toRAD)*toDEG; //same definition in GEMC
      db_file<<"gemc.gem"<<count<<".r0 = "<<setprecision(6)<<fixed<<RIN[i]<<endl;
      db_file<<"gemc.gem"<<count<<".r1 = "<<setprecision(6)<<fixed<<ROUT[i]<<endl;
      db_file<<"gemc.gem"<<count<<".phi0 = "<<setprecision(6)<<fixed<<phi0<<endl;
      db_file<<"gemc.gem"<<count<<".dphi = "<<setprecision(0)<<fixed<<DPHI<<endl;
      db_file<<"gemc.gem"<<count<<".z0 = "<<setprecision(6)<<fixed<<Z0[i]<<endl;
      db_file<<"gemc.gem"<<count<<".depth = "<<setprecision(2)<<fixed<<DEPTH<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"x.stripangle = "<<setprecision(1)<<fixed<<UANGLE<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"x.pitch = "<<setprecision(4)<<fixed<<PITCH<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"y.stripangle = "<<setprecision(1)<<fixed<<VANGLE<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"y.pitch = "<<setprecision(4)<<fixed<<PITCH<<endl;
      db_file<<endl;
      count++;
    }
  }
  db_file.close();
  
  return 1;
  
}
