//generate the db_gemc.dat file for GEM digitization

#include <iostream>
#include <fstream>
#include <iomanip>
#include "TVector2.h"
#include "TMath.h"
#include <vector>

using namespace std;

struct HVSector{
    Double_t xmin;
    Double_t xmax;
    Double_t ymin;
    Double_t ymax;
    HVSector(const Double_t& x1, const Double_t& x2, const Double_t& y1, const Double_t y2)
    : xmin(x1), xmax(x2), ymin(y1), ymax(y2) {} 
};

void InitHVSector(vector< vector<HVSector> >& v);

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
  const double FRAMEWIDTH = 0.011;
  const double Z0[NPLANE]     = {1.575, 1.855, 1.900, 3.060,  3.150 };
  const double RIN[NPLANE]    = {0.48,  0.59,  0.65,  1.05,   1.09  };
  const double ROUT[NPLANE]   = {1.22,  1.43,  1.43,  2.30,   2.37  };
  const double PHIOFF[NPLANE] = {3.2,   2.2,   2.2,   -0.8,   -0.8  };
  const double PHISTART = -6;
  const double ROTDIR = 1.;
  ofstream db_file;
  vector< vector<HVSector> > HVSECTOR;
  InitHVSector(HVSECTOR);
  
  int count = 1;
  db_file.open ("db_gemc.dat");
  db_file<<"950101\n";
  for (int i=0; i<NPLANE; i++){
    double phi0 = PHISTART;
    for (int j=0; j<NSECTOR; j++){
      phi0 = TVector2::Phi_0_2pi((ROTDIR*(PHISTART + j*DPHI) - PHIOFF[i])*toRAD)*toDEG; //same definition in GEMC
      db_file<<"gemc.gem"<<count<<".r0 = "<<setprecision(6)<<fixed<<RIN[i]<<endl;
      db_file<<"gemc.gem"<<count<<".r1 = "<<setprecision(6)<<fixed<<ROUT[i]<<endl;
      db_file<<"gemc.gem"<<count<<".phi0 = "<<setprecision(6)<<fixed<<phi0<<endl;
      db_file<<"gemc.gem"<<count<<".dphi = "<<setprecision(0)<<fixed<<DPHI<<endl;
      db_file<<"gemc.gem"<<count<<".z0 = "<<setprecision(6)<<fixed<<Z0[i]<<endl;
      db_file<<"gemc.gem"<<count<<".depth = "<<setprecision(2)<<fixed<<DEPTH<<endl;
      db_file<<"gemc.gem"<<count<<".frame_width = "<<setprecision(6)<<fixed<<FRAMEWIDTH<<endl;
      db_file<<"gemc.gem"<<count<<".n_HV_sector_off = "<<HVSECTOR[i].size()<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"x.stripangle = "<<setprecision(1)<<fixed<<UANGLE<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"x.pitch = "<<setprecision(4)<<fixed<<PITCH<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"y.stripangle = "<<setprecision(1)<<fixed<<VANGLE<<endl;
      db_file<<"gemc.gem"<<count<<".gem"<<count<<"y.pitch = "<<setprecision(4)<<fixed<<PITCH<<endl;
      
      for (unsigned int k=0; k<HVSECTOR[i].size(); k++){
        db_file<<"gemc.gem"<<count<<".HV"<<k+1<<".bound = "<<HVSECTOR[i][k].xmin<<" "<<HVSECTOR[i][k].xmax<<" "<<HVSECTOR[i][k].ymin<<" "<<HVSECTOR[i][k].ymax<<endl;
      }
      
      db_file<<endl;
      count++;
    }
  }
  db_file.close();
  
  return 1;
  
}
//_____________________________________________________________________________________________________
void InitHVSector(vector< vector<HVSector> >& v)
{
    vector<HVSector> thisV; thisV.clear();
    //tracker 1
    thisV.push_back(HVSector(0.46, 0.52, -0.5, 0.5));
    thisV.push_back(HVSector(0.46, 0.65, 0.02, 0.2));
    thisV.push_back(HVSector(0.46, 1.23, -0.3, -0.05));
    v.push_back(thisV);
    thisV.clear();
    
    //tracker 2
    thisV.push_back(HVSector(0.58, 0.62, -0.5, 0.5));
    thisV.push_back(HVSector(0.58, 0.75, 0.02, 0.2));
    thisV.push_back(HVSector(0.58, 1.44, -0.3, -0.06)); 
    v.push_back(thisV);
    thisV.clear();
    
    //tracker 3
    thisV.push_back(HVSector(0.64, 0.8, 0.02, 0.1));
    thisV.push_back(HVSector(0.64, 1.44, -0.3, -0.06));
    v.push_back(thisV);
    thisV.clear();
    
    //tracker 4
    thisV.push_back(HVSector(1.04, 1.50, 0.06, 0.2));
    thisV.push_back(HVSector(1.04, 2.31, -0.3, -0.13));
    v.push_back(thisV);
    thisV.clear();
    
    //tracker 5
    thisV.push_back(HVSector(1.08, 1.55, 0.06, 0.2));
    thisV.push_back(HVSector(1.08, 2.38, -0.3, -0.13));
    v.push_back(thisV);
    thisV.clear();
}
