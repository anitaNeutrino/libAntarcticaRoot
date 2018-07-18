#ifndef GEOID_MODEL_H
#define GEOID_MODEL_H

#include "TVector3.h"
#include "TMath.h"

/**
 * @namespace GeoidModel
 * @brief Get positions, radii, latitudes, longitudes, and other goodies when modelling the Earth
 */
namespace GeoidModel {

  static constexpr double R_EARTH = 6.378137E6;
  static constexpr double GEOID_MAX = 6.378137E6; // parameters of geoid model
  static constexpr double GEOID_MIN = 6.356752E6; // parameters of geoid model
  static constexpr double FLATTENING_FACTOR = (1./298.257223563);

  inline Double_t getGeoidRadiusAtCosTheta(Double_t cosTheta);
  Double_t getGeoidRadiusAtLatitude(Double_t lat);
  inline Double_t getGeoidRadiusAtTheta(Double_t theta);

  void getCartesianCoords(Double_t lat, Double_t lon, Double_t alt, Double_t p[3]);
  void getLatLonAltFromCartesian(const Double_t p[3], Double_t &lat, Double_t &lon, Double_t &alt);
  Double_t getDistanceToCentreOfEarth(Double_t lat);



  /**
   * @class Vector
   * 
   * @brief The ultimate method to represent a position a Geoid modelled Earth
   * 
   * Very much in the spirit of https://xkcd.com/927/
   * 
   * This class is supposed to combine the best features of a TVector, an icemc::Vector, 
   * and the proper GeoidModel transformations.
   * 
   * Features:
   * Lazy, behind the scenes conversion of x,y,z to lon, lat, alt, easting, northing
   * 
   */
  class Vector : public TVector3 {

  public:

    Vector() : TVector3() {};
    Vector(TVector3& v) : TVector3(v) {};
    Vector(Double_t x, Double_t y, Double_t z) : TVector3(x, y, z) {};
    template <class T> Vector(const T& t);
    template <class T> Vector(const T* t);

    /**
     * Longitude/Latitude/Altitude Getter functions
     */
    inline Double_t Latitude() const;
    inline Double_t Longitude() const;
    inline Double_t Altitude() const;
    inline Double_t Easting() const;
    inline Double_t Northing() const;

    Double_t GetGeoidRadius() const;
    inline void GetLonLatAlt(Double_t& lon, Double_t& lat, Double_t& alt) const;
    template <class T> void GetLonLatAlt(T& t) const;
    template <class T> void GetLonLatAlt(T* t) const;

    inline void SetLongitude(double longitude);
    inline void SetLatitude(double latitude);
    inline void SetAltitude(double altitude);
    inline void SetLonLatAlt(double lon, double lat, double alt);

    inline Double_t Theta() const;
    inline Double_t Phi() const;

  private:

    /**
     * The fX, fY, fZ in the base TVector class are the definitive 
     * representation of the position. In other words, any non-const 
     * functions ALWAYS update the fX, fY, fZ with setCartesianFromGeoid().
     * The fLongitude, fLatitude, fAltitude, and fTheta, fPhi, 
     * and fEasting, fNorthing get calculated lazily as needed.
     * The results are cached in the mutable (non-ROOT persistent) variables
     * (marked by comments after variable beginning with //!).
     */

    // can't make const as the cartesian x,y,z are represented in TVector3
    // and accessors aren't virtual and so can't be overloaded. 
    void setCartesianFromGeoid();
    void setGeoidFromCartesian() const;
    void setAnglesFromCartesian() const;

    void setEastingNorthingFromLonLat() const;
    // void setGeoidFromEastingNorthing() const;

    mutable Double_t fLongitude = 0; //! cached longitude is not stored!
    mutable Double_t fLatitude  = 0; //! cached latitude is not stored!
    mutable Double_t fAltitude  = 0; //! cached altitude is not stored!
    mutable Double_t fTheta     = 0; //! cached angle is not stored!
    mutable Double_t fPhi       = 0; //! cached angle is not stored!
    mutable Double_t fEasting   = 0; //! cached polar coordinate is not stored!
    mutable Double_t fNorthing  = 0; //! cached polar coordinate is not stored!

    mutable Double_t fCartAtLastGeoidCalc[3]       = {-1, -1, -1};    //! fX, fY, fZ when the geoid was last updated. Is not stored!
    mutable Double_t fCartAtLastAngleCal[3]        = {0, 0, 0};       //! fX, fY, fZ when the angles were last updated. Is not stored!
    mutable Double_t fLonLatAtLastEastNorthCalc[2] = {-9999, -9999};  //! Longitude(), Latitude() when the Easting/Northing were last updated. Is not stored!    
  };






















  template <class T> Vector::Vector(const T& t){
    SetLonLatAlt(t.longitude,  t.latitude, t.altitude);
    setCartesianFromGeoid();
  };

  template <class T> Vector::Vector(const T* t) : Vector(*t){;}


  inline void Vector::SetLongitude(double lon) {
    if(lon > 180){
      lon -= 360;
    }
    fLongitude = lon;
    setCartesianFromGeoid();
  }
  inline void Vector::SetLatitude(double lat){
    fLatitude = lat;
    setCartesianFromGeoid();
  }
  inline void Vector::SetAltitude(double alt) {
    fAltitude = alt;
    setCartesianFromGeoid();
  }
  inline void Vector::SetLonLatAlt(double lon, double lat, double alt) {
    fLatitude = lat;
    fAltitude = alt;
    SetLongitude(lon);
  }


  // inline void Vector::SetEasting(double easting) {
  //   fEasting = easting;
  //   setGeoidFromEastingNorthing();
  // }
  // inline void Vector::SetNorthing(double easting) {
  //   fEasting = easting;
  //   setGeoidFromEastingNorthing();
  // }
  // inline void Vector::SetEastingNorthingAlt(double easting, double northing, double alt) {
  //   fEasting = easting;
  //   fNorthing = northing;
  //   fAltitude = alt;
  // }



  inline Double_t Vector::Latitude() const {
    setGeoidFromCartesian();
    return fLatitude;
  }
  inline Double_t Vector::Longitude() const {
    setGeoidFromCartesian();
    return fLongitude;
  }
  inline Double_t Vector::Altitude() const {
    setGeoidFromCartesian();
    return fAltitude;
  }

  inline Double_t Vector::Theta() const {
    setAnglesFromCartesian();
    return fTheta;
  }
  inline Double_t Vector::Phi() const {
    setAnglesFromCartesian();
    return fPhi;
  }


  inline Double_t Vector::Easting() const {
    setEastingNorthingFromLonLat();
    return fEasting;
  }

  inline Double_t Vector::Northing() const {
    setEastingNorthingFromLonLat();
    return fNorthing;
  }

  inline Double_t Vector::GetGeoidRadius() const {
    return getGeoidRadiusAtCosTheta(CosTheta());
  }



  void Vector::GetLonLatAlt(Double_t& lon, Double_t& lat, Double_t& alt) const {
    // call functions rather than direct access so cache is updated!
    lon = Longitude();
    lat = Latitude();
    alt = Altitude();
  }

  template <class T>
  void Vector::GetLonLatAlt(T& t) const {
    GetLonLatAlt(t.longitude,  t.latitude,  t.altitude);
  }

  template <class T>
  void Vector::GetLonLatAlt(T* t) const {
    GetLonLatAlt(*t);
  }












  inline Double_t getGeoidRadiusAtCosTheta(Double_t cosTheta) {
    return GEOID_MIN*GEOID_MAX/TMath::Sqrt(GEOID_MIN*GEOID_MIN-(GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*cosTheta*cosTheta);
  }
  inline Double_t getGeoidRadiusAtTheta(Double_t theta) {
    return getGeoidRadiusAtCosTheta(TMath::Cos(theta));
  }



  
}


#endif
