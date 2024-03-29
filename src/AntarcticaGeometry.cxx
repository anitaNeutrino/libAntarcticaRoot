#include "AntarcticaGeometry.h" 
#include "Geoid.h" 
#include "RampdemReader.h" 
#include "Adu5Pat.h" 
#include "FFTtools.h"
#include "TRandom.h" 
#include "TGraph2D.h" 
#include "TString.h" 
#include "TFile.h" 

#ifdef USE_GEOGRAPHIC_LIB
#include "GeographicLib/GeodesicLine.hpp" 
#include "GeographicLib/Geodesic.hpp" 
#endif

ClassImp(AntarcticCoord); 
ClassImp(AntarcticSegmentationScheme); 
ClassImp(StereographicGrid); 

template <int coarseness>
class SurfaceWrapper
{

  public:
    double compute(double E, double N) 
    {
//      printf("%g,%g\n", E,N); 
      return map->Interpolate(E,N); 
    }

    TProfile2D *map; 
    SurfaceWrapper(RampdemReader::dataSet set) 
    {
      map = RampdemReader::getMap(set, coarseness); 
    }

    ~SurfaceWrapper()
    {
//      delete map; 

    }
}; 

template <int coarseness> 
static SurfaceWrapper<coarseness> & getSurface(RampdemReader::dataSet set) 
{
  if (set == RampdemReader::surface)
  {
    static SurfaceWrapper<coarseness> w(set); 
    return w; 
  }

  //the only other thing that makes sense is rampdem
  static SurfaceWrapper<coarseness> w(RampdemReader::rampdem); 
  return w; 
}





void AntarcticCoord::asString(TString * s) const
{

  if (type == WGS84) 
  {
    s->Form("%g N %g E %g m",x,y,z); 
  }
  
  if (type == STEREOGRAPHIC)
  {
    s->Form("%g m(E), %g m(N), %g m",x,y,z); 
  }

  if (type == CARTESIAN)
  {
    s->Form("(%g,%g,%g) m",x,y,z); 
  }

}


static void cart2stereo(double*,double*,double*) __attribute__((optimize("fast-math"),optimize("O3"))); 
static void stereo2cart(double*,double*,double*) __attribute__((optimize("fast-math"),optimize("O3"))); 

// for ECE -> lat/lon 
static const double a = 6378137; 
static const double b = 6356752.31424518; 
static const double binv = 1./6356752.31424518; 
static const double e = sqrt((a*a-b*b)/(a*a)); 
static const double ep = e * a/b; 


static void cart2stereo(double *x, double * y, double *z)
{
   //Turns out this is an important conversion for performance. Will hackily put it in right here.
   //I'm using a slightly different conversion from cartesian to WGS84, and also avoiding all trig functions, which just cause havoc. 
   // This is based on  http://www.microem.ru/pages/u_blox/tech/dataconvert/GPS.G1-X-00006.pdf

   double X = *y; //silly 
   double Y = *x; //silly 
   double Z = *z;  

   double H = sqrt(X*X+Y*Y); 
   double Hinv = 1./H; 

   double sin_lon = Y*Hinv; 
   double cos_lon = Y!=0 ? (X/Y * sin_lon) : ( (X > 0) - (X < 0));  

   double tan_theta = Z*Hinv*(a*binv); 
   double cos_theta = pow(1+ tan_theta* tan_theta,-0.5); 
   double sin_theta = tan_theta * cos_theta; 
      
   //this might be able to be simplified 
   double num = Z + (ep *ep*b) * (sin_theta * sin_theta * sin_theta); 
   double denom = H - (e *e*a) * (cos_theta * cos_theta * cos_theta); 

   double H2inv = pow(num*num + denom*denom,-0.5); 
   double sin_lat =  num * H2inv * ( (denom > 0) - ( denom < 0) ); //do I need this sign change here? I suspect denom is always well greater than zero... 
   double cos_lat =  denom / num * sin_lat; 


   double N = a * pow(1. - (e*e) * sin_lat * sin_lat,-0.5); 

   *z = H / cos_lat - N; 

   ///ok now, we have to go to stereographic. We already reversed the sign of latitude

   using namespace Geoid;
   double R = scale_factor * c_0 * pow( (1 + eccentricity * sin_lat ) / (1 - eccentricity * sin_lat), eccentricity/2) ; 
      
   //use some trig identities here... 
   double tan_lat_over_2 = sin_lat / (1 + cos_lat); 
   R *= (1 - tan_lat_over_2) / (1 + tan_lat_over_2); 

   *x = R * sin_lon; 
   *y = R * cos_lon; 

}



static void stereo2cart(double *x, double *y, double *z) 
{
  double E = *x; 
  double N = *y; 

  if (!E && !N) //special case the south pole
  {
    *z +=b; 
    return; 
  }

  double alt = *z; 

  double H2 = E*E + N *N; 
  double H = sqrt(H2); 
  double Hinv =  1./H; 
  double sin_lon = E*Hinv;
  double cos_lon = N*Hinv; 


  // Here I'm just taking the equations from EastingNorthingToLonLat and avoiding trig using
  // identities:  
  // sin(pi/2-x) = cos(x) ,
  // cos(pi/2-x) = sin(x)  , 
  // cos(2x) = (1 - tan^2 x) / (1 + tan^2 x) ,
  // sin(2x) = 2 tan x / ( 1 + tan^2 x) 
  using namespace Geoid; 
  double tan2_factor = H2 /(scale_factor * scale_factor * c_0 * c_0);  
  double sin_iso_lat = (1 - tan2_factor) / (1 + tan2_factor);  
  double cos_iso_lat = 2 * H/(scale_factor * c_0) / (1 + tan2_factor); 

  // and now the fun begins
  double s2 = sin_iso_lat * sin_iso_lat;
  double sc = sin_iso_lat * cos_iso_lat; 
  double s3c = s2 * sc; 
  double s5c = s2 * s3c; 
  double s7c = s2 * s5c; 

  double sin2_iso_lat = 2 * sc; 
  double sin4_iso_lat = 4 * sc - 8 * s3c; 
  double sin6_iso_lat = 6 * sc - 32 * s3c + 32 *s5c; 
  double sin8_iso_lat = 8 * sc - 80 * s3c + 192* s5c - 128* s7c; 

  //Well, I guess I will use trig functions here... since otherwise I'll have a very very nasty sum 
  //There is probably a better way to do this calculation, but it's quite tricky 

  double correction =  a_bar * sin2_iso_lat + b_bar * sin4_iso_lat + c_bar * sin6_iso_lat + d_bar * sin8_iso_lat; 

  // sadly, one trig call 
  double sin_corr = sin(correction); 
  double cos_corr = cos(correction); 
  double sin_lat = fabs(sin_iso_lat * cos_corr + cos_iso_lat * sin_corr);  //force positivity 
  double cos_lat = cos_iso_lat * cos_corr - sin_iso_lat * sin_corr; 

  //printf("    lat: %g lon: %g\n", atan2(sin_lat, cos_lat) * TMath::RadToDeg(), atan2(sin_lon,cos_lon) * TMath::RadToDeg());  

  //copied and pasted from getCartesianCoords, with a few intermediate values stored
  
  const double C1 = (1-Geoid::FLATTENING_FACTOR) * (1-Geoid::FLATTENING_FACTOR); 
  
  double C2 = pow(cos_lat * cos_lat + C1 * sin_lat*sin_lat,-0.5);
  double Q2 = C1* C2; 
  
  double C3 = Geoid::R_EARTH * C2 + alt; 

  //same silly convention as before here 
  *x = C3 * cos_lat * sin_lon; 
  *y = C3 * cos_lat * cos_lon; 
  *z = (Geoid::R_EARTH * Q2 + alt) * sin_lat; 
}



void AntarcticCoord::convert(CoordType t) 
{

  if (type == WGS84) 
  {
    if (t == CARTESIAN) 
    {
      double cartesian[3]; 
      Geoid::getCartesianCoords(x,y,z, cartesian); 
      x = cartesian[0]; 
      y = cartesian[1]; 
      z = cartesian[2]; 
    }

    else if ( t==STEREOGRAPHIC) 
    {
      RampdemReader::LonLatToEastingNorthing(y,x,x,y); 
    }

  }

  else if (type == STEREOGRAPHIC)
  {


    if (t == WGS84) 
    {
      RampdemReader::EastingNorthingToLonLat(x,y,y,x); 
    }

    else if (t == CARTESIAN) 
    {
//      convert(WGS84); 
//      convert(CARTESIAN); 
      stereo2cart(&x,&y,&z); 
      
    }
  }

  else if (type == CARTESIAN) 
  {

    if (t == WGS84) 
    {
      double cartesian[3]; 
      cartesian[0] = x; 
      cartesian[1] = y; 
      cartesian[2] = z; 
      Geoid::getLatLonAltFromCartesian(cartesian, x,y,z); 
    }
    else if ( t== STEREOGRAPHIC)
    {

//      convert(WGS84); 
//      convert(STEREOGRAPHIC); 

     cart2stereo(&x,&y,&z); 

     }
  }

  else 
  {
    fprintf(stderr,"Invalid CoordType: %d!", type); 
    return; 
  }


  type = t; 

}



AntarcticSegmentationScheme * AntarcticSegmentationScheme::factory(const char * key) 
{

  if (strstr(key,"stereographic_grid"))
  {
    int nx,ny,maxE,maxN; 
    ny = -1; 
    nx = 0; 
    maxE = 3330000; 
    maxN = 3330000; 
    sscanf(key,"stereographic_grid_%d_%d_%d_%d", &nx, &ny, &maxN, &maxE); 
    if (!nx) return 0; 
    if (ny < 0) ny = nx; 
    return new StereographicGrid(nx,ny,maxE,maxN); 
  }

  return 0; 
}

const int nsamples = 64; 

void AntarcticSegmentationScheme::DrawI(const char * opt, const int * idata, const double * range) const 
{
  int N = NSegments(); 
  double * data = 0;
  if (idata) 
  {
    data = new double[N]; 
    for (int i = 0; i < N; i++) { data[i] = idata[i]; } 
    Draw(opt,data,range); 
    delete []  data; 
  }
}

void AntarcticSegmentationScheme::Draw(const char * opt, const double * data, const double * range) const 
{
  int N = NSegments(); 

  TGraph2D * g  = new TGraph2D(N*nsamples); 
  g->SetBit(kCanDelete); 

  AntarcticCoord samples[nsamples]; 

  for (int i =0; i < N; i++) 
  {
    /* do not need altitude for this! */ 
    sampleSegment(i, nsamples, samples, true,false); 

    for (int j =0;  j<  nsamples; j++) 
    {
      samples[j].to(AntarcticCoord::STEREOGRAPHIC); 
      g->SetPoint(nsamples*i+j,samples[j].x,samples[j].y,data ? data[i] : i); 
    }
  }

  if (range) 
  {
    g->GetXaxis()->SetRangeUser(range[0],range[1]); 
    g->GetYaxis()->SetRangeUser(range[2],range[3]); 
  }

  g->Draw(opt); 
}



StereographicGrid::StereographicGrid(int NX, int NY, double MAX_E, double MAX_N)
  : AntarcticSegmentationScheme(), nx(NX), ny(NY), max_E(MAX_E), max_N(MAX_N) 
{

  dx = (2*max_E)/nx; 
  dy = (2*max_N)/ny; 
}

int StereographicGrid::getSegmentIndex(const AntarcticCoord & coord) const
{
  AntarcticCoord stereo = coord.as(AntarcticCoord::STEREOGRAPHIC); 

  if (stereo.x >  max_E || stereo.x < -max_E || stereo.y > max_N || stereo.y < -max_N) return -1; 


  int xbin = nx/2 + stereo.x / dx; 
  int ybin = ny/2 - stereo.y / dy; 

  return xbin + nx * ybin; 

}

void StereographicGrid::getSegmentCenter(int idx, AntarcticCoord * fillme, bool fillalt) const 
{

  int xbin = idx % nx; 
  int ybin = idx / nx; 
  double x = (xbin +0.5) *dx - max_E; 
  double y = max_N - (ybin +0.5) *dy ; 
  double z = fillalt ? getSurface<4> (dataset).compute(x,y): 0; 
  fillme->set(AntarcticCoord::STEREOGRAPHIC,x,y,z); 
}


AntarcticCoord * StereographicGrid::sampleSegment(int idx, int N, AntarcticCoord * fill, bool random, bool fillalt) const 
{

  if (!fill) fill = new AntarcticCoord[N]; 

  int xbin = idx % nx; 
  int ybin = idx / nx; 
  double lx = (xbin) *dx - max_E; 
  double ly = max_N - (ybin) *dy ; 


  if (random) 
  {
    for (int i = 0; i < N; i++)
    {
      double x = gRandom->Uniform(lx, lx+dx); 
      double y = gRandom->Uniform(ly, ly-dy); 
      double z = fillalt ? getSurface<4> (dataset).compute(x,y): 0; 
      fill[i].set(AntarcticCoord::STEREOGRAPHIC,x,y,z); 
    }
  }
  else
  {
    int grid = sqrt(N) + 0.5 ; 

    for (int i = 0; i < N; i++)
    {
      //TODO check if this is what i want! 
      double x = lx + (dx / (grid-1)) * ((i % grid)); 
      double y = ly - (dy / (grid-1)) * ((i / grid)); 
      double z = fillalt ? getSurface<4> (dataset).compute(x,y): 0; 
      fill[i].set(AntarcticCoord::STEREOGRAPHIC,x,y,z); 
    }
  }

  return fill; 
}



void StereographicGrid::Draw(const char * opt, const double * data, const double * range) const
{
  if(gDirectory->Get("tmp")){
    delete gDirectory->Get("tmp");
  }
  TH2 * h = 0; 

  TString stropt(opt); 
  stropt.ToUpper(); 

  bool do_map = false;
  if (stropt.Contains("MAP"))
  {
    stropt.ReplaceAll("MAP",""); 
    h = new TH2DAntarctica("tmp","Stereographic Grid", nx, ny); 
    do_map = true;
    AntarcticCoord c;
    for (int i = 0; i < NSegments(); i++) 
    {
       getSegmentCenter(i,&c,false);
       c.to(AntarcticCoord::WGS84);
       h->Fill(c.y, c.x, data ? data[i]: i); 
       // std::cout << c.x<<" "<<c.y<<" "<<data[i]<<std::endl;
    }
  }
  else 
  {
    h = new TH2D("tmp","Stereographic Grid", nx, -max_E, max_E, ny, -max_N, max_N); 
    AntarcticCoord c;
    for (int i = 0; i < NSegments(); i++) 
    {
       getSegmentCenter(i,&c,false);
       h->Fill(c.x, c.y, data ? data[i]: i); 
       // std::cout << c.x<<" "<<c.y<<" "<<data[i]<<std::endl;
    }
  }

  

  h->SetStats(0); 

  if (range) 
  {

    h->GetXaxis()->SetRangeUser(range[0], range[1]); 
    h->GetYaxis()->SetRangeUser(range[2], range[3]); 

  }
  // h->DrawCopy(opt);
  // std::cout<<opt<< " "<< stropt<<std::endl; 
  h->Draw(stropt); 
  // delete h; 
}

void StereographicGrid::asString(TString * str) const
{
  str->Form("stereographic_grid_%d_%d_%d_%d", nx,ny,(int) max_E,(int) max_N); 
}



int StereographicGrid::getNeighbors(int segment, std::vector<int> * neighbors) const
{

  int N = 0; 

  bool left_edge = (segment % nx ) == 0; 
  bool right_edge = (segment % nx ) == nx-1 ; 
  bool top_edge = (segment / ny ) == 0; 
  bool bottom_edge = (segment / ny ) == ny-1 ; 

  int start_x = left_edge ? 0 : -1; 
  int end_x = right_edge ? 0 : 1; 
  int start_y = top_edge ? 0 : -1; 
  int end_y = bottom_edge ? 0 : 1; 

  for (int i = start_x; i <= end_x; i++) 
  {
    for (int j = start_y; j <= end_y; j++)
    {
      if (i == 0 && j == 0) continue;
      if (i * j != 0) continue; // added by peng, only looking for 4 neighbours instead of 8 neighbours. 

      N++; 
      if (neighbors) neighbors->push_back(segment + i + j *nx); 
    }
  }

  return N; 
}










CartesianSurfaceMap::CartesianSurfaceMap(double resolution , RampdemReader::dataSet d) 
{
  // figure out the corners of the stereographic projection in cartesian. 
 
  AntarcticCoord sa(AntarcticCoord::STEREOGRAPHIC); 
  AntarcticCoord sb(AntarcticCoord::STEREOGRAPHIC); 
  AntarcticCoord sc(AntarcticCoord::STEREOGRAPHIC); 
  AntarcticCoord sd(AntarcticCoord::STEREOGRAPHIC); 
  
  double xu,xl, yl,yu; 
  RampdemReader::getMapCoordinates(xl,yl,xu,yu, d); 

  sa.x = xl; 
  sa.y = yl; 
  sb.x = xu; 
  sb.y = yl; 
  sc.x = xu; 
  sc.y = yu; 
  sd.x = xl; 
  sd.y = yu; 


  sa.to(AntarcticCoord::CARTESIAN); 
  sb.to(AntarcticCoord::CARTESIAN); 
  sc.to(AntarcticCoord::CARTESIAN); 
  sd.to(AntarcticCoord::CARTESIAN); 

  double xmin = TMath::Min (  TMath::Min( sa.x,sb.x), TMath::Min(sc.x,sd.x)); 
  double xmax = TMath::Max (  TMath::Max( sa.x,sb.x), TMath::Max(sc.x,sd.x)); 
  double ymin = TMath::Min (  TMath::Min( sa.y,sb.y), TMath::Min(sc.y,sd.y)); 
  double ymax = TMath::Max (  TMath::Max( sa.y,sb.y), TMath::Max(sc.y,sd.y)); 


  int nbinsx = (xmax - xmin)/ resolution; 
  int nbinsy = (ymax - ymin)/ resolution; 
  map = new TH2D("cartmap","Cartesian Map", nbinsx, xmin, xmax, nbinsy, ymin,ymax); 
  map->SetDirectory(0); 

  for (int i = 1; i <= nbinsx; i++) 
  {
    for (int j = 1; j <= nbinsy; j++)
    {
      // get cartesian point on geoid... 
      double x = map->GetXaxis()->GetBinCenter(i); 
      double y =  map->GetYaxis()->GetBinCenter(j); 

      static const double a = 6378137; 
      static const double b = 6356752.31424518; 
      double z = b * sqrt( 1  -  (1./(a*a))*(x*x+y*y)); 
      AntarcticCoord c(AntarcticCoord::CARTESIAN, x,y,z); 

      double R = c.v().Mag(); 

      c.to(AntarcticCoord::STEREOGRAPHIC); 
      double alt = 0; 
      if (c.x > xl && c.x < xu && c.y > yl && c.y < yu) 
      {
        alt = getSurface<4>(d).map->Interpolate(c.x, c.y); 
      }

      map->SetBinContent(i,j, R + alt); 
    }
  }
}


double CartesianSurfaceMap::CartesianSurfaceMap::surface(double x, double y) const 
{
  return map->Interpolate(x,y); 
}

double CartesianSurfaceMap::CartesianSurfaceMap::z(double x, double y) const 
{
  double r = surface(x,y); 
  return sqrt( r*r - x*x - y*y); 
}



double CartesianSurfaceMap::metersAboveIce(double x, double y, double z) const 
{
  return sqrt(x*x + y*y+z*z) - surface(x,y); 
}



#ifndef USE_GEOGRAPHIC_LIB
static int nagged_about_surface = 0; 
static double hav(double phi) { return ((1-cos(phi)/2)); }
#endif

double AntarcticCoord::surfaceDistance(const AntarcticCoord & a, const AntarcticCoord & b) 
{
  AntarcticCoord A = a.as(WGS84); 
  AntarcticCoord B = b.as(WGS84); 
  return surfaceDistance(A.x, B.x, A.y, B.y); 


}
double AntarcticCoord::surfaceDistance(double lat0 , double lat1, double lon0, double lon1) 
{
  if (lat0 < -90 || lat0 < -90 || lon0 < -180 || lon1 < -180) return -999; 
#ifndef USE_GEOGRAPHIC_LIB

  if (!nagged_about_surface)
  {
    nagged_about_surface++; 
    fprintf(stderr,"WARNING: compiled without geographic lib support. Using haversine distance between points\n"); 
  }
  return 2 * Geoid::R_EARTH * asin( hav(lat1-lat0) + cos(lat0) * cos(lat1) * hav(lon1-lon0)); 
#else

  double distance; 
  GeographicLib::Geodesic::WGS84().Inverse( lat0,lon0,lat1,lon1, distance); 
  return distance; 
#endif

}


