//*-- Author :    Ole Hansen  08-March-2013
//*-- Modified by: Weizhi Xiong 16-Oct-2016
// dbconvert.cxx
//
// Utility to convert a libsolgem-style database file to the format expected
// by the SoLIDTracking

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "TMath.h"
#include "TVector2.h"
#include "TRotation.h"
#include "VarDef.h"
#include <cstdlib>
#include <cctype>
#include <algorithm>
#include <functional>
#include <map>
#include <utility>
#include <iomanip>
#include <cassert>
#include <unistd.h>

using namespace std;

#define ALL(c) (c).begin(), (c).end()

#define INFILE_DEFAULT  "db_gemc.dat"
#define OUTFILE_DEFAULT "db_solid.trackersystem.dat"

static bool do_debug = false;
static bool do_ECal  = true;
static string infile = INFILE_DEFAULT;
static string outfile = OUTFILE_DEFAULT;
static const char* prgname = "";
static const string spc = "                   ";

//_____________________________________________________________________________
static string find_key( ifstream& inp, const string& key )
{
  static const string empty("");
  string line;
  string::size_type keylen = key.size();
  inp.seekg(0); // could probably be more efficient, but it's fast enough
  while( getline(inp,line) ) {
    if( line.size() <= keylen )
      continue;
    if( line.compare(0,keylen,key) == 0 ) {
      if( keylen < line.size() ) {
	string::size_type pos = line.find_first_not_of(" \t=", keylen);
	if( pos != string::npos )
	  return line.substr(pos);
      }
      break;
    }
  }
  return empty;
}
//______________________________________________________________________________
static int load_db( ifstream& inp, DBRequest* request, const string& prefix )
{
  DBRequest* item = request;
  while( item->name ) {
    ostringstream sn(prefix, ios_base::ate);
    sn << item->name;
    const string& key = sn.str();
    string val = find_key(inp,key);
    if( !val.empty() ) {
      istringstream sv(val);
      sv >> *((double*)item->var);
      if( !sv ) {
	cerr << "Error converting key/value = " << key << "/" << val << endl;
	return 1;
      }
      if( do_debug )
	cout << "Read: " << key << " = " << *((double*)item->var) << endl;
    } else {
      cerr << "key \"" << key << "\" not found" << endl;
      return 2;
    }
    ++item;
  }
  return 0;
}
//_________________________________________________________________________________
static
void print_req( const DBRequest* req )
{
  const DBRequest* item = req;
  while( item->name ) {
    cout << " " << item->name << " = " << *((double*)item->var) << endl;
    ++item;
  }
}
//_________________________________________________________________________________
// Parameters for one GEM chamber (two readout coordinates)
struct ValueSet_t {
  // Values read from source file
  double rmin, rmax, phi0, dphi, z, dz;
  double xangle, xpitch, yangle, ypitch;
  // Computed/converted quantities
  double phi, phioff;
  double angle[2], start[2], pitch[2];  // angle is that of the axis!
  int nstrips[2];                       // number of u/v strips
  int iplane, isector;                  // plane index (1-5), sector number

};
//__________________________________________________________________________________
// Helper functors for STL algos
struct BySectorThenPlane
  : public std::binary_function< const ValueSet_t&, const ValueSet_t&, bool >
{
  bool operator() ( const ValueSet_t& L, const ValueSet_t& R ) const
  {
    return (L.isector != R.isector) ?
      (L.isector < R.isector) : (L.iplane < R.iplane);
  }
};
//___________________________________________________________________________________
void usage()
{
  // Print usage message and exit
  cerr << "Usage: " << prgname << "[-hdm] [-o outfile] [infile]" << endl;
  cerr << " Convert libsolgem database <infile> to TreeSearch-SoLID database"
       << " <outfile>" << endl;
  cerr << " -h: Print this help message" << endl;
  cerr << " -d: Output extensive debug information" << endl;
  cerr << " -o <outfile>: Write output to <outfile>. Default: "
       << OUTFILE_DEFAULT << endl;
  cerr << " <infile>: Read input from <infile>. Default: "
       << INFILE_DEFAULT << endl;
  exit(255);
}
//______________________________________________________________________________________
void getargs( int argc, const char** argv )
{
  // Get command line parameters

  while (argc-- > 1) {
    const char *opt = *++argv;
    if (*opt == '-') {
      while (*++opt != '\0') {
	switch (*opt) {
	case 'h':
	  usage();
	  break;
	case 'd':
	  do_debug = true;
	  break;
	case 'o':
	  if (!*++opt) {
	    if (argc-- < 1)
	      usage();
	    opt = *++argv;
	  }
	  outfile = opt;
	  opt = "?";
	  break;
	default:
	  usage();
	}
      }
    } else {
      infile = *argv;
    }
  }
}
//________________________________________________________________________________________
static inline
void write_module( ostream& outp, int ic, int slot_hi, int model, int nchan )
{
  outp << spc << ic << "    ";
  if( ic < 10 )  outp << " ";
  outp << 0 << "       " << slot_hi << "      ";
  if( slot_hi < 10 ) outp << " ";
  outp << model << "   " << nchan;
}
//_________________________________________________________________________________________
static inline
void write_cslh( ostream& outp, int cr, int sl, int lo, int hi )
{
  outp << cr << " ";
  if( cr < 10 )  outp << " ";
  outp << sl << " ";
  if( sl < 10 ) outp << " ";
  outp << lo << " " << hi;
}
//__________________________________________________________________________________________
int main( int argc, const char** argv )
{
  prgname = basename(argv[0]);
  getargs(argc,argv);

  const int    nsystem   = 30; //PVDIS has 30 independent tracker system
  const int    nplanes   = 5;  //each system has 5 GEM plane
  const int    nchamber  = 1;  //each plane has 1 GEM chamber
  const int    nproj     = 2;  //each chamber has two readout plane (or projection)
  const double Uangle    = 84; //12 stereo angle symmetric readout strips
  const double Vangle    = 96; //these angles are not the direction the strips point to
                               //but vertical to that direction
  const double z_shift   = -0.31825/100.; //GEMC gives z position of a chamber as the center
                                          //of the GEM, we want the center of the first gas layer
  const double phi_start = 0;//starting phi angle of the 0-th chamber   
  
  double strip_angle[2] = { Uangle, Vangle};
  
  //phi angle offset of the GEM chamber on each plane
  //double phi_offset[nplanes] = {3.0,   2.0,   2.0,   0.,   0. }; //CLEO
  double phi_offset[nplanes] = {3.2, 2.2, 2.2, -0.8, -0.8};
  
  //nominal z position of each GEM plane
  double plane_z[nplanes]    = {1.575, 1.855, 1.900, 3.060,  3.150};
  for (int i=0; i<nplanes; i++) plane_z[i] += z_shift;
  
  const string prefix = "gemc.";
  const string out_prefix = "solid.trackersystem.";
  //  const string out_prefix = "${DET}.";
  const string allsystems_prefix = out_prefix + "${allsystems}.";
  const string alltrackers_prefix = allsystems_prefix + "${alltrackers}.";
  const string allchambers_prefix = alltrackers_prefix + "${allchambers}.";
  const string allreadouts_prefix = allchambers_prefix + "${allreadouts}.";
  const char* proj_name[2] = { "u", "v" };
  const string dashes =
    "#-----------------------------------------------------------------";
    
  vector<ValueSet_t> values;
  
  ifstream inp(infile.c_str());
  if( !inp ) {
    cerr << "Error opening " << infile << endl;
    exit(1);
  }

  int max_nstrips = 0;
  map<int,int> nstrip_map;

  int nplanes_eff = nplanes;
  int nchamber_eff = nsystem;
  
  for( int ip = 0; ip < nplanes_eff; ++ip ) {
    for( int is = 0; is < nchamber_eff; ++is ) {
      ValueSet_t vals;
      
      DBRequest request[] = {
	{ "r0",    &vals.rmin },
	{ "r1",    &vals.rmax },
	{ "phi0",  &vals.phi0 },
	{ "dphi",  &vals.dphi },
	{ "z0",    &vals.z },
	{ "depth", &vals.dz },
	{ 0 }
      };
      DBRequest plane_request[] = {
	{ "x.stripangle", &vals.xangle },
	{ "x.pitch",      &vals.xpitch },
	{ "y.stripangle", &vals.yangle },
	{ "y.pitch",      &vals.ypitch },
	{ 0 }
      };
      if( ip < nplanes ) {
	// Regular GEM planes
	int idx = is + nchamber_eff*ip; // linear index of this plane/sector combo
	ostringstream sector_prefix(prefix, ios_base::ate);
	sector_prefix << "gem" << idx+1 << ".";
	int err = load_db( inp, request, sector_prefix.str() );
	if( err )
	  exit(2);

	if( vals.rmin <= 0 or vals.rmax <= 0 or vals.rmax <= vals.rmin ) {
	  cerr << "Invalid radii r0 = " << vals.rmin
	       << ", r1 = " << vals.rmax << endl;
	  exit(3);
	}
	if( vals.dphi < 0 or vals.dphi >= 90.0 ) {
	  cerr << "Invalid opening angle dphi = " << vals.dphi << endl;
	  exit(3);
	}
	if( vals.dz < 0 ) {
	  cerr << "Invalid  z = " << vals.z << endl;
	  exit(3);
	}
	// Override the database z value - it is not accurate enough
	if( fabs(vals.z-plane_z[ip]) > 1e-2 ) { // max 1cm tolerance
	  cerr << "Input z = " << vals.z <<  "differs from \"correct\" z = "
	       << plane_z[ip] << "by more than 1cm. Update dbconvert."
	       << endl;
	  exit(3);
	}
	vals.z = plane_z[ip];

	ostringstream plane_prefix(sector_prefix.str(), ios_base::ate);
	plane_prefix << "gem" << idx+1; // sic, the same thing again
	err = load_db( inp, plane_request, plane_prefix.str() );
	if( err )
	  exit(2);

	if( vals.xpitch <= 0 or vals.ypitch <= 0 ) {
	  cerr << "Illegal strip pitch xpitch = " << vals.xpitch
	       << ", ypitch = " << vals.ypitch << endl;
	  exit(3);
	}
	vals.phioff = -1.*phi_offset[ip];
      }

      // Convert parameters from libsolgem conventions to ours
      double phi2 = 0.5*vals.dphi;  // half opening angle
      vals.phi    = vals.phi0 + phi2 - vals.phioff;

      //keep phi from 0 to 360 in data base
      vals.phi = TVector2::Phi_0_2pi(vals.phi*TMath::DegToRad())*TMath::RadToDeg();

      // Calculate strip start positions in the same way as in
      // TSolGEMPlane::ReadGeometry
      double torad = TMath::DegToRad(), phi2rad = phi2 * torad;
      double xs = 0.5 * ( vals.rmax - vals.rmin * TMath::Cos(phi2rad) );
      double ys = vals.rmax * TMath::Sin(phi2rad);
      if( do_debug )
	cout << " xs/ys = " << xs << "/" << ys << endl;
      TRotation plane_to_xstrip, plane_to_ystrip;
      plane_to_xstrip.RotateZ(-vals.xangle*torad);
      plane_to_ystrip.RotateZ(-vals.yangle*torad);
      TVector3 TR(xs,ys,0), BR(xs,-ys,0), TL(-xs,ys,0), BL(-xs,-ys,0);
      TVector3 C[4] = { TR, BR, TL, BL };
      int sminx = 1e9, sminy = 1e9, smaxx = -1e9, smaxy = -1e9;
      for( int i = 0; i < 4; ++i ) {
	if( do_debug ) {
	  cout << " i = " <<  i << " ";
	  C[i].Print();
	}
	TVector3 vx = plane_to_xstrip * C[i];
	TVector3 vy = plane_to_ystrip * C[i];
	int sx = (int) (vx.X() / vals.xpitch);
	int sy = (int) (vy.X() / vals.ypitch);
	sminx = min(sx,sminx);
	smaxx = max(sx,smaxx);
	sminy = min(sy,sminy);
	smaxy = max(sy,smaxy);
      }
      vals.nstrips[0] = smaxx - sminx + 1;
      vals.nstrips[1] = smaxy - sminy + 1;
      vals.pitch[0] = vals.xpitch;
      vals.pitch[1] = vals.ypitch;
      for( int ij = 0; ij < nproj; ij++ ) {
	vals.start[ij] = -vals.nstrips[ij] * vals.pitch[ij] * 0.5;
	// For the tracking code, the strip positions should be the _center_
	// of the strips, not the lower edge
	vals.start[ij] += 0.5 * vals.pitch[ij];
      }
      // Correct angles for offset. The results are the projection axis angles.
      // For a given projection, they must be the same for all planes.
      vals.angle[0] = vals.xangle + vals.phioff;
      vals.angle[1] = vals.yangle + vals.phioff;

      // Save the number of strips of each readout for later use when writing
      // the detector maps
      if( ip < nplanes ) {
	for( int ij = 0; ij < nproj; ij++ ) {
	  max_nstrips = max(max_nstrips,vals.nstrips[ij]);
	  int jx = ij + nproj*( ip + nplanes*is );
	  pair<map<int,int>::iterator,bool> itx =
	    nstrip_map.insert( make_pair(jx,vals.nstrips[ij]) );
	  if( !itx.second ) {
	    cerr << "Duplicate index " << jx << " for sector/plane/proj = "
		 << is+1 << "/" << ip+1 << "/" << ij << endl;
	    cerr << "Bug - should never happen." << endl;
	    exit(6);
	  }
	}
      }

      // Remember who we are. Counting starts at 1, as the database is
      // meant for humans ...
      vals.iplane  = ip+1;
      vals.isector = is+1;

      // Save results
      values.push_back( vals );

      if( do_debug ) {
	// Display results for debugging
	//      cout << sector_prefix.str() << endl;
	print_req( request );
	print_req( plane_request );
	cout << " phi/offset = " << vals.phi << "/" << vals.phioff << endl;

	cout << " " << proj_name[0] << "/" << proj_name[1] << " n/ang/start/pitch = ";
	for( int i=0; i < 2; ++i ) {
	  cout << vals.nstrips[i] << "/" << vals.angle[i] << "/" << vals.start[i]
	       << "/" << vals.pitch[i];
	  if( i == 0 )
	    cout << "  ";
	}
	cout << endl;
      } // do_debug

    } // all sectors
  }   // all planes
  
  inp.close();
  
  vector<double> proj_angle(nproj,-1e10);
  vector<double> sect_phi(nchamber_eff,-1e10);
  for( vector<ValueSet_t>::size_type i = 0; i < values.size(); ++i ) {
    ValueSet_t& v = values[i];
    
    int is = v.isector-1;
    if( sect_phi[is] < -1e9 ) {
      sect_phi[is] = v.phi;
    } else if( fabs(sect_phi[is] - v.phi) > 1e-10 ) {
      cerr << "Error: inconsistent sector angle = " << v.phi
	   << " in sector " << is+1 << " at plane " << v.iplane
	   << endl
	   << "Expected " << sect_phi[is] << endl;
      exit(4);
    }
  }
  
   //==== Write output ====
  ofstream outp( outfile.c_str(), ios_base::out|ios_base::trunc );
  if( !outp ) {
    cerr << " Error opening output file " << outfile << endl;
    exit(5);
  }

  // Header
  //TDatime now;
  outp << "# -*- mode: Text -*-" << endl
       << "#" << endl
       << "# Database for SoLIDTracking" << endl
       << "#" << endl
       << "# Converted from " << infile << endl
       // << "# on " << now.AsString() << " by " << basename(argv[0]) << endl
       << endl;
       
  
  //************************parameter common for the tracker system***********************************// 
  // including parameters used for pattern recognition and track fit, since they are done here//
  outp << dashes<<endl;
  outp << "# parameters common for all tracker systems"<<endl;
  outp << dashes<<endl;
  outp << allsystems_prefix << "detconf = 1" << endl;
  outp << allsystems_prefix << "MCdata = 1" << endl;
  outp << allsystems_prefix << "ntracker = "<<nplanes << endl;
  outp << allsystems_prefix << "do_rawdecode = "<<1<<endl;
  outp << allsystems_prefix << "do_coarsetrack = "<<1<<endl;
  outp << allsystems_prefix << "do_finetrack = "<<0<<endl;
  outp << allsystems_prefix << "do_chi2 = "<<0<<endl;
  outp << allsystems_prefix << "chi2_cut = "<<1e4<<endl;
  outp << allsystems_prefix << "max_miss_hit = "<<1<<endl;
  int model = 6425;    // Dummy model for virtual APV25
  int nchan = 1600;    // Number of channels per module (arbitrary)
  int MAXSLOT = 30;    // Max slots per crate (from THaCrateMap.h)
  int modules_per_readout = max_nstrips/nchan+1;
  int modules_per_chamber = nproj*modules_per_readout; // Modules needed per chamber
  int chambers_per_crate = (MAXSLOT/modules_per_chamber/nplanes)*nplanes;
  cout<<max_nstrips<<" "<<chambers_per_crate<<" "<<modules_per_chamber<<" "<<MAXSLOT<<" "<<nplanes<<endl;
  int slot_hi = chambers_per_crate*modules_per_chamber-1;
  if( do_debug ) {
    cout << "Crate map: modules_per_readout = " << modules_per_readout << endl;
    cout << "           modules_per_chamber = " << modules_per_chamber << endl;
    cout << "           chambers_per_crate  = " << chambers_per_crate << endl;
    cout << "           sectors_per_crate   = " << chambers_per_crate/nplanes
	 << endl;
  }
  outp << endl;
  outp << "# \"Crate map\". Specifies the overall DAQ module configuration." << endl;
  outp << "# The map can be common to all sectors. It's just a lookup table" << endl;
  outp << "# for (crate,slot) -> (model,nchan)" << endl;
  outp << "#" << endl;
  outp << "# Each row is:     crate slot_lo slot_hi model# nchan" << endl;
  outp << "#" << endl;
  outp << allsystems_prefix << "cratemap = \\" << endl;

  // Remember that each sector must have its own virtual hardware.
  // However, all sectors can share the same crate map - it's just a catalog
  // of which modules are in which slot, regardless of whether or not we use
  // the slot.
  int maxcrates = TMath::CeilNint( (double)nplanes*nchamber_eff /
				   (double)chambers_per_crate );
  for( int ic = 0; ic < maxcrates; ++ ic ) {
    write_module( outp, ic, slot_hi, model, nchan );
   
    if( ic+1 != maxcrates or do_ECal) outp << " \\";
    outp << endl;
  }
  //for FAEC and LAEC
  write_module( outp, maxcrates, 3, model, nchamber_eff );
  outp<<endl;
//*****************************************************************************************//
////***********************parameters for ecal********************************//
  outp << endl;
  outp << dashes<<endl;
  outp << "# parameters for ecal"<<endl;
  outp << dashes<<endl;
  outp << allsystems_prefix << "ecal.laec_z = " <<0<< endl;
  outp << allsystems_prefix << "ecal.faec_z = " <<3.23<< endl;
  outp << allsystems_prefix << "ecal.ec_pos_reso = " <<0.01<< endl;
  outp << allsystems_prefix << "ecal.ec_energy_reso = " <<0.1<< endl;
  outp << allsystems_prefix << "ecal.laec_detmap_pos = ";
  write_cslh( outp, maxcrates, 0, 0, 29 ); outp<<endl;
  outp << allsystems_prefix << "ecal.laec_detmap_edp = ";
  write_cslh( outp, maxcrates, 1, 0, 29 ); outp<<endl;
  outp << allsystems_prefix << "ecal.faec_detmap_pos = ";
  write_cslh( outp, maxcrates, 2, 0, 29 ); outp<<endl;
  outp << allsystems_prefix << "ecal.faec_detmap_edp = ";
  write_cslh( outp, maxcrates, 3, 0, 29 ); outp<<endl;
//*****************************************************************************************//  

//***********************parameters common for all trackers********************************//
  outp << endl;
  outp << dashes<<endl;
  outp << "# parameters common for all trackers"<<endl;
  outp << dashes<<endl;
  outp << alltrackers_prefix << "nchamber = " <<nchamber<< endl;
  outp << alltrackers_prefix << "combine_hits = " <<1<< endl;
  outp << alltrackers_prefix << "kill_cross_talk = "<<1<<endl;
  outp << alltrackers_prefix << "cross_talk_thres = "<<0.1<<endl;
  outp << alltrackers_prefix << "cross_strip_apart = "<<32<<endl;
//******************************************************************************************//
  
//**********************parameters common for all chambers**********************************//
//include also parameters used for hit coordinate matching since it is done here//
  outp << endl;
  outp << dashes<<endl;
  outp << "# parameters common for all chambers"<<endl;
  outp << dashes<<endl;
  outp << allchambers_prefix << "do_3d_amcorr = " <<1<< endl;
  outp << allchambers_prefix << "3d_amcorr_cut = " <<1<< endl;
  outp << allchambers_prefix << "nreadout = " <<2<< endl;
  outp << allchambers_prefix << "acc_file_name = " <<"pos_acceptance.root"<< endl;
//******************************************************************************************//
//***********************parameters common for all readouts*********************************//
//include also parameters used for readout decoding and hit clustering since they are done here//
  outp << endl;
  outp << dashes<<endl;
  outp << "# parameters common for all readouts"<<endl;
  outp << dashes<<endl;
  
  outp << allreadouts_prefix << "xp_res = " << 6e-05<< endl;
  outp << allreadouts_prefix << "maxclustsiz = "<< 4 <<endl;
  outp << allreadouts_prefix << "adc_min = " <<95 <<endl;
  outp << allreadouts_prefix << "split_frac = "<<0.2 <<endl;
  outp << allreadouts_prefix << "maxhits = " <<1000<<endl;
  outp << allreadouts_prefix << "maxsamp = "<<1<<endl;
  outp << allreadouts_prefix << "adc_sigma = "<<0.2<<endl;
  outp << allreadouts_prefix << "do_noise = "<<0<<endl;
  outp << allreadouts_prefix << "check_pulse_shape = "<<1<<endl;
  outp << allreadouts_prefix << "do_histos = "<<0<<endl;
  outp << allreadouts_prefix << "strip_pitch = "<<0.0004<<endl;
  outp << allreadouts_prefix << "dz = "<<0.05<<endl;
  outp << allreadouts_prefix << "dphi = "<<12.<<endl;
  outp << allreadouts_prefix << "deconmode = "<<0<<endl; 
//********************************************************************************************//  

//******************************tracker system specific parameters****************************//
  outp << endl;
  outp << dashes<<endl;
  outp << "# tracker system specific parameters"<<endl;
  outp << "# for SIDIS there is only 1 tracker system"<<endl;
  outp << "# for PVDIS there are 30, each has a rotation with respect to the lab frame"<<endl;
  outp << dashes<<endl;
  
  double tmp_phi_angle = phi_start;
  for (int i=0; i< nsystem; i++){
    tmp_phi_angle = TVector2::Phi_0_2pi(tmp_phi_angle*TMath::DegToRad())*TMath::RadToDeg();
    outp << out_prefix <<i<<"."<<"phi = "<<tmp_phi_angle<<endl;
    tmp_phi_angle += 12.;
  }
//********************************************************************************************//

//******************************tracker specific parameters***********************************//
  outp << endl;
  outp << dashes<<endl;
  outp << "# tracker specific parameters"<<endl;
  outp << "# for PVDIS, the concept of tracker and chamber is the same"<<endl;
  outp << "# so I will leave the phi_offset to chamber specific parameter"<<endl;
  outp << dashes<<endl;
  
  for (int j=0; j<nplanes; j++){
    outp<<allsystems_prefix<<j<<"."<<"tracker_z = "<<plane_z[j]<<endl;
  }
  
//*********************************************************************************************//

//******************************chamber specific parameters************************************//
//include phi offset, z offset, size of the chamber and so on//
  outp << endl;
  outp << dashes<<endl;
  outp << "# chamber specific parameters"<<endl;
  outp << dashes<<endl;
  
  for (int i=0; i<nsystem; i++){
    for (int j=0; j<nplanes; j++){
      for (int k=0; k<nchamber; k++){
        outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<"rmin = "<<values.at(j*nsystem+i).rmin<<endl;
        outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<"rmax = "<<values.at(j*nsystem+i).rmax<<endl;
        outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<"dz = "<<0<<endl;
        outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<"phi = "<<0<<endl;
        outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<"phi_cover = "<<12.<<endl;
        outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<"phi_offset = "<<phi_offset[j]<<endl;
        outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<"reference = "<<0.<<" "<<0.<<endl;
        outp << endl;
      }
    }
  }  

//***********************************************************************************************//

//********************readout plane specific parameters******************************************//
//include number of strips, starting position and so on//
  outp << endl;
  outp << dashes<<endl;
  outp << "# readout plane specific parameters"<<endl;
  outp << dashes<<endl;
  
    for (int i=0; i<nsystem; i++){
      for (int j=0; j<nplanes; j++){
        for (int k=0; k<nchamber; k++){
          for (int n=0; n<nproj; n++){
            outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<n<<"."<<"start = "<<values.at(j*nsystem+i).start[n]<<endl;
            outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<n<<"."<<"strip_angle = "<<strip_angle[n]<<endl;
            outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<n<<"."<<"nstrips = "<<values.at(j*nsystem+i).nstrips[n]<<endl;
            outp << out_prefix<<i<<"."<<j<<"."<<k<<"."<<n<<"."<<"detmap = ";
            if( modules_per_readout > 1 ) outp << "\\" << endl;
            // Find the actual number of channels used by this particular
            // readout. We don't strictly need this, but setting a limit here
            // provides extra protection against decoder bugs.
            int jx = n + nproj*( j + nplanes*i );
            map<int,int>::iterator it = nstrip_map.find(jx);
            if( it == nstrip_map.end() ) {
            cerr << "Error retrieving nstrips for sector/plane/proj = "
               << k << "/" << j << "/" << n << endl;
            cerr << "Bug - should never happen." << endl;
            exit(6);
          }
          int the_nstrips = (*it).second;
          // Write detector map for this readout
          for( int im = 0; im < modules_per_readout; ++im ) {
          int ix = im + modules_per_readout*( n + nproj*( j + nplanes*i ));
          int cr = ix / (chambers_per_crate*modules_per_chamber);
          int sl = ix - cr*chambers_per_crate*modules_per_chamber;
          int lo = im*nchan;
          int hi = min((im+1)*nchan,the_nstrips)-1;
          if( modules_per_readout > 1 ) outp << spc;
          write_cslh( outp, cr, sl, lo, hi );
          if( 2*(im+1) < modules_per_chamber ) outp << " \\";
            outp<<endl;
            outp<<endl;
            }
          }
        }
      }
    }
  
//***********************************************************************************************//
 outp.close();
  
}






