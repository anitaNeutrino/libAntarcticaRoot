#include "GeoidModel.h"
#include "TMath.h"
#include "RampdemReader.h"
#include <iostream>


inline void convert_cartesian_coordinates_between_wgs84_and_anita_conventions(double& x, double& y, double& z){
  // The cartesian coordinate system in WGS84 has:
  // The origin at the center of mass of the Earth.
  // The +ve z-axis running through the north pole.
  // The +ve x-axis running through the prime merian (0 longitude) at the equator
  // The +ve y-axis is picked such that the coordinate system is right handed.

  // ANITA, for historical reasons, does things a little differently.
  // The +ve z-axis runs through the south pole.
  // Then the x and y axes are swapped, which maintains a right handed coordinate system
  // So, the +ve y-axis comes out of the Earth at the equator at 0 longitude.
  // This does have the nice property that it's still right handed,
  // the +ve x aligns with easting, +ve y aligns with northing,
  // and quantities like elevation have more +ve z for higher altitude in Antarctica.

  // Still, it's a bit confusing, but we're well over a decade in now...

  z = -z;
  double oldX = x;
  x = y;
  y = oldX;
}


void GeoidModel::getCartesianCoords(Double_t lat, Double_t lon, Double_t alt, Double_t p[3]){

  // see page 71 onwards of https://web.archive.org/web/20120118224152/http://mercator.myzen.co.uk/mercator.pdf  
  lat *= TMath::DegToRad();
  lon *= TMath::DegToRad();
  //calculate x,y,z coordinates

  Double_t C2 = pow(TMath::Cos(lat)*TMath::Cos(lat)+(1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR)*TMath::Sin(lat)*TMath::Sin(lat),-0.5);
  Double_t Q2 = (1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR)*C2;
  p[0]=(R_EARTH*C2+alt)*TMath::Cos(lat)*TMath::Cos(lon);
  p[1]=(R_EARTH*C2+alt)*TMath::Cos(lat)*TMath::Sin(lon);
  p[2]=(R_EARTH*Q2+alt)*TMath::Sin(lat);

  convert_cartesian_coordinates_between_wgs84_and_anita_conventions(p[0], p[1],  p[2]);
}

void GeoidModel::getLatLonAltFromCartesian(const Double_t p[3], Double_t &lat, Double_t &lon, Double_t &alt){


  // see page 71 onwards of https://web.archive.org/web/20120118224152/http://mercator.myzen.co.uk/mercator.pdf  
  Double_t x=p[0];
  Double_t y=p[1];
  Double_t z=p[2];

  convert_cartesian_coordinates_between_wgs84_and_anita_conventions(x, y, z);  

  static Double_t cosaeSq=(1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR);
  
  const Double_t lonVal   = TMath::ATan2(y,x);
  const Double_t xySq     = TMath::Sqrt(x*x+y*y);
  const Double_t tanPsit  = z/xySq;
  Double_t latGuess = TMath::ATan(tanPsit/cosaeSq);
  Double_t nextLat  = latGuess;
  Double_t geomBot  = R_EARTH*R_EARTH*xySq;

  const double deltaLatCloseEnough = 1e-6; // this corresponds to < 1m at the equator
  do {
    latGuess=nextLat;
    Double_t N      = R_EARTH/TMath::Sqrt(cos(latGuess)*cos(latGuess)+cosaeSq*sin(latGuess)*sin(latGuess));
    Double_t top    = (R_EARTH*R_EARTH*z + (1-cosaeSq)*cosaeSq*TMath::Power(N*TMath::Sin(latGuess),3));
    Double_t bottom = geomBot-(1-cosaeSq)*TMath::Power(N*TMath::Cos(latGuess),3);
    nextLat = TMath::ATan(top/bottom);
    // std::cout << latGuess << "\t" << nextLat << "\n";
  } while(TMath::Abs(nextLat-latGuess) > deltaLatCloseEnough);
  latGuess=nextLat;

  Double_t N = R_EARTH/TMath::Sqrt(cos(latGuess)*cos(latGuess)+cosaeSq*sin(latGuess)*sin(latGuess));
  Double_t height=(xySq/TMath::Cos(nextLat))-N;
  
  lat = latGuess*TMath::RadToDeg();
  lon = lonVal*TMath::RadToDeg();
  alt = height;
}






Double_t GeoidModel::getDistanceToCentreOfEarth(Double_t lat)
{
  Vector v;
  v.SetLonLatAlt(0, lat, 0);
  return v.Mag();
}


Double_t GeoidModel::getGeoidRadiusAtLatitude(Double_t latitude) {
  Vector v;
  v.SetLonLatAlt(0, latitude, 0);
  return getGeoidRadiusAtCosTheta(v.CosTheta());
}






void GeoidModel::Vector::updateCartesianFromGeoid() {
  // always called after any lon/lat/alt has been updated
  GeoidModel::getCartesianCoords(fLatitude, fLongitude, fAltitude, fCartAtLastGeoidCalc);
  SetXYZ(fCartAtLastGeoidCalc[0], fCartAtLastGeoidCalc[1], fCartAtLastGeoidCalc[2]);
}



void GeoidModel::Vector::updateGeoidFromCartesian() const {
  // called when Longitude(), Latitude(), Altitude() is requested
  if(X() != fCartAtLastGeoidCalc[0] ||
     Y() != fCartAtLastGeoidCalc[1] ||
     Z() != fCartAtLastGeoidCalc[2]){
    // Then we've moved, so must recalculate lon, lat alt;
    GetXYZ(fCartAtLastGeoidCalc);
    GeoidModel::getLatLonAltFromCartesian(fCartAtLastGeoidCalc, fLatitude, fLongitude, fAltitude);
  }  
}


void GeoidModel::Vector::updateAnglesFromCartesian() const {

  bool xDirty = X() != fCartAtLastAngleCal[0];
  bool yDirty = Y() != fCartAtLastAngleCal[1];
  
  if(xDirty || yDirty){
    fPhi = TVector3::Phi();
    fCartAtLastAngleCal[0] = X();
    fCartAtLastAngleCal[1] = Y();
  }

  if(xDirty || yDirty || Z() != fCartAtLastAngleCal[2]){
    fTheta = TVector3::Theta();
    // if x or y was dirty, already stored them in fCartAtLastAngleCal
    fCartAtLastAngleCal[2] = Z();
  }
}



void GeoidModel::Vector::updateEastingNorthingFromLonLat() const {

  if(fLongitude != fLonLatAtLastEastNorthCalc[0] ||
     fLatitude != fLonLatAtLastEastNorthCalc[1]){

    fLonLatAtLastEastNorthCalc[0] = Longitude();
    fLonLatAtLastEastNorthCalc[1] = Latitude();
    
    RampdemReader::LonLatToEastingNorthing(fLonLatAtLastEastNorthCalc[0],
					   fLonLatAtLastEastNorthCalc[1],
					   fEasting, fNorthing);
    
  }
}


void GeoidModel::Vector::updateLonLatFromEastingNorthing(bool mustRecalcuateAltitudeFirst) {
  // Since we are going to eventually update cartesian from lon/lat/alt
  // we must make sure altitude is up to date...
  //
  double lon, lat;
  double alt = mustRecalcuateAltitudeFirst ? Altitude() : fAltitude;
  
  RampdemReader::EastingNorthingToLonLat(fEasting, fNorthing, lon, lat);
  SetLonLatAlt(lon, lat, alt);
}
