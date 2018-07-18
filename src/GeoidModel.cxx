#include "GeoidModel.h"
#include "TMath.h"
#include "RampdemReader.h"
#include <iostream>

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
}

void GeoidModel::getLatLonAltFromCartesian(const Double_t p[3], Double_t &lat, Double_t &lon, Double_t &alt){

  // see page 71 onwards of https://web.archive.org/web/20120118224152/http://mercator.myzen.co.uk/mercator.pdf
  Double_t x=p[0];
  Double_t y=p[1];
  Double_t z=p[2];

  static Double_t cosaeSq=(1-FLATTENING_FACTOR)*(1-FLATTENING_FACTOR);
  
  const Double_t lonVal   = TMath::ATan2(y,x);
  const Double_t xySq     = TMath::Sqrt(x*x+y*y);
  const Double_t tanPsit  = z/xySq;
  Double_t latGuess = TMath::ATan(tanPsit/cosaeSq);
  Double_t nextLat  = latGuess;
  Double_t geomBot  = R_EARTH*R_EARTH*xySq;

  const double epsilon = 1e-6; // this corresponds to < 1m at the equator
  do {
    latGuess=nextLat;
    Double_t N      = R_EARTH/TMath::Sqrt(cos(latGuess)*cos(latGuess)+cosaeSq*sin(latGuess)*sin(latGuess));
    Double_t top    = (R_EARTH*R_EARTH*z + (1-cosaeSq)*cosaeSq*TMath::Power(N*TMath::Sin(latGuess),3));
    Double_t bottom = geomBot-(1-cosaeSq)*TMath::Power(N*TMath::Cos(latGuess),3);
    nextLat = TMath::ATan(top/bottom);
    // std::cout << latGuess << "\t" << nextLat << "\n";
  } while(TMath::Abs(nextLat-latGuess)>epsilon);
  latGuess=nextLat;

  Double_t N = R_EARTH/TMath::Sqrt(cos(latGuess)*cos(latGuess)+cosaeSq*sin(latGuess)*sin(latGuess));
  Double_t height=(xySq/TMath::Cos(nextLat))-N;
  
  lat = latGuess*TMath::RadToDeg();  lon = lonVal*TMath::RadToDeg();
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
  

void GeoidModel::Vector::setCartesianFromGeoid() {
  GeoidModel::getCartesianCoords(fLatitude, fLongitude, fAltitude, fCartAtLastGeoidCalc);
  SetXYZ(fCartAtLastGeoidCalc[0], fCartAtLastGeoidCalc[1], fCartAtLastGeoidCalc[2]);
}



void GeoidModel::Vector::setGeoidFromCartesian() const {

  if(X() != fCartAtLastGeoidCalc[0] ||
     Y() != fCartAtLastGeoidCalc[1] ||
     Z() != fCartAtLastGeoidCalc[2]){
    // Then we've moved, so must recalculate lon, lat alt;
    GetXYZ(fCartAtLastGeoidCalc);
    GeoidModel::getLatLonAltFromCartesian(fCartAtLastGeoidCalc, fLatitude, fLongitude, fAltitude);
  }  
}


void GeoidModel::Vector::setAnglesFromCartesian() const {

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



void GeoidModel::Vector::setEastingNorthingFromLonLat() const {

  if(fLongitude != fLonLatAtLastEastNorthCalc[0] ||
     fLatitude != fLonLatAtLastEastNorthCalc[1]){

    fLonLatAtLastEastNorthCalc[0] = Longitude();
    fLonLatAtLastEastNorthCalc[1] = Latitude();
    
    RampdemReader::LonLatToEastingNorthing(fLonLatAtLastEastNorthCalc[0],
					   fLonLatAtLastEastNorthCalc[1],
					   fEasting, fNorthing);
    
  }
}
