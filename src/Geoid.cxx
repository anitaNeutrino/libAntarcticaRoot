#include "Geoid.h"
#include "TMath.h"
#include <iostream>



void Geoid::getCartesianCoords(Double_t lat, Double_t lon, Double_t alt, Double_t p[3]){

  // see page 71 onwards of https://web.archive.org/web/20120118224152/http://mercator.myzen.co.uk/mercator.pdf  
  lat *= TMath::DegToRad();
  lon *= TMath::DegToRad();
  //calculate x,y,z coordinates

  // Double_t C2 = pow(TMath::Cos(lat)*TMath::Cos(lat)+(1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR)*TMath::Sin(lat)*TMath::Sin(lat),-0.5);
  Double_t C2 = 1./TMath::Sqrt(TMath::Cos(lat)*TMath::Cos(lat)+(1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR)*TMath::Sin(lat)*TMath::Sin(lat));
  Double_t Q2 = (1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR)*C2;

  // Swapping x/y and inverting z
  p[1]=(R_EARTH*C2+alt)*TMath::Cos(lat)*TMath::Cos(lon);
  p[0]=(R_EARTH*C2+alt)*TMath::Cos(lat)*TMath::Sin(lon);
  p[2]=-(R_EARTH*Q2+alt)*TMath::Sin(lat);

}

void Geoid::getLatLonAltFromCartesian(const Double_t p[3], Double_t &lat, Double_t &lon, Double_t &alt){


// swapping x,y and inverting z
  Double_t x=p[1]; 
  Double_t y=p[0];
  Double_t z=-p[2];

  // see page 71 onwards of https://web.archive.org/web/20120118224152/http://mercator.myzen.co.uk/mercator.pdf  

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






Double_t Geoid::getDistanceToCentreOfEarth(Double_t lat)
{
  Position v;
  v.SetLonLatAlt(0, lat, 0);
  return v.Mag();
}


Double_t Geoid::getGeoidRadiusAtLatitude(Double_t latitude) {
  Position v;
  v.SetLonLatAlt(0, latitude, 0);
  return getGeoidRadiusAtCosTheta(v.CosTheta());
}






void Geoid::Position::updateCartesianFromGeoid() {
  // always called after any lon/lat/alt has been updated
  Geoid::getCartesianCoords(fLatitude, fLongitude, fAltitude, fCartAtLastGeoidCalc);
  SetXYZ(fCartAtLastGeoidCalc[0], fCartAtLastGeoidCalc[1], fCartAtLastGeoidCalc[2]);
}



void Geoid::Position::updateGeoidFromCartesian() const {
  // called when Longitude(), Latitude(), Altitude() is requested
  if(X() != fCartAtLastGeoidCalc[0] ||
     Y() != fCartAtLastGeoidCalc[1] ||
     Z() != fCartAtLastGeoidCalc[2]){

    // Then we've moved, so must recalculate lon, lat alt;
    GetXYZ(fCartAtLastGeoidCalc);
    Geoid::getLatLonAltFromCartesian(fCartAtLastGeoidCalc, fLatitude, fLongitude, fAltitude);
  }  
}


void Geoid::Position::updateAnglesFromCartesian() const {

  bool xDirty = X() != fCartAtLastAngleCalc[0];
  bool yDirty = Y() != fCartAtLastAngleCalc[1];
  
  if(xDirty || yDirty){
    fPhi = TVector3::Phi();
    fCartAtLastAngleCalc[0] = X();
    fCartAtLastAngleCalc[1] = Y();
  }

  if(xDirty || yDirty || Z() != fCartAtLastAngleCalc[2]){
    fTheta = TVector3::Theta();
    // if x or y was dirty, already stored them in fCartAtLastAngleCalc
    fCartAtLastAngleCalc[2] = Z();
  }
}



void Geoid::Position::updateEastingNorthingFromLonLat() const {

  if(fLongitude != fLonLatAtLastEastNorthCalc[0] ||
     fLatitude != fLonLatAtLastEastNorthCalc[1]){

    fLonLatAtLastEastNorthCalc[0] = Longitude();
    fLonLatAtLastEastNorthCalc[1] = Latitude();
    
    LonLatToEastingNorthing(fLonLatAtLastEastNorthCalc[0],
			    fLonLatAtLastEastNorthCalc[1],
			    fEasting, fNorthing);
    
  }
}


void Geoid::Position::updateLonLatFromEastingNorthing(bool mustRecalcuateAltitudeFirst) {
  // Since we are going to eventually update cartesian from lon/lat/alt
  // we must make sure altitude is up to date...
  //
  double lon, lat;
  double alt = mustRecalcuateAltitudeFirst ? Altitude() : fAltitude;
  
  EastingNorthingToLonLat(fEasting, fNorthing, lon, lat);
  SetLonLatAlt(lon, lat, alt);
}
