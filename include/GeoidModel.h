#ifndef GEOID_MODEL_H
#define GEOID_MODEL_H

#include "TVector3.h"
#include "TMath.h"
#include <iostream>

/**
 * @namespace GeoidModel
 * @brief Get positions, radii, latitudes, longitudes, and other goodies when modelling the Earth
 */
namespace GeoidModel {

  /**
   * Constants
   */
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



  enum class Pole {North,South}; // for chosing solutions for Geoid z as a function of x,y


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

    Vector(Double_t x=0, Double_t y=0) : TVector3(x, y, 0) {
      moveToGeoidZ();
    }
    Vector(Pole pole) : TVector3(0, 0, 0) {
      moveToGeoidZ(pole);
    }

    Vector(Double_t x, Double_t y, Double_t z) : TVector3(x, y, z) {};
    Vector(TVector3& v) : TVector3(v) {};

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

    inline void SetEasting(double easting);
    inline void SetNorthing(double northing);
    inline void SetEastingNorthing(double easting, double northing);
    inline void SetEastingNorthingAlt(double easting, double northing, double alt);

    inline Double_t Theta() const;
    inline Double_t Phi() const;


    /** 
     * Find the value of z on Geoid surface given the values X(), Y()
     * Note: There are two solutions here, pole chooses which.
     * 
     * @param pole solution of z to find (Pole::North or Pole::South), default is Pole::South (near Antarctica).
     * 
     * @return The value of Z
     */
    inline double surfaceZ(Pole pole = Pole::South);

    /** 
     * Change z so that we are on the surface at this X(),  Y()
     * @see surfaceZ(Pole pole)
     * 
     * @param signZ sign of Z() to move to the surface, default is towards Pole::South (near Antarctica).
     */
    inline void moveToGeoidZ(Pole pole = Pole::South);





  private:

    /**
     * A note on Cartesian coordinates: We don't use WGS84 convention!
     * 
     * The cartesian coordinate system in WGS84 has:
     * The origin at the center of mass of the Earth.
     * The +ve z-axis running through the north pole.
     * The +ve x-axis running through the prime merian (0 longitude) at the equator
     * The +ve y-axis is picked such that the coordinate system is right handed.
     *
     * ANITA and icemc, for historical reasons, does things a little differently.
     * The origin is in the same place but:
     * The +ve z-axis runs through the south pole.
     * Then the x-axis and y-axis are swapped relative to WGS84.
     * That x/y swap maintains a right handed coordinate system.
     * i.e. the +ve y-axis comes out of the Earth at the equator at 0 longitude.
     * The ANITA/icemc coordinate system is still right handed, and has the property
     * that +ve x-axis aligns with easting, +ve y-axis aligns with northing,
     * and quantities like elevation are more +ve in z for higher altitude in Antarctica.
     * Which does make plots of things in Cartesian coordinates a bit easier to look at.
     */
    inline void swap_wgs84_anita(double& x, double& y, double& z){
      z = -z;
      double oldX = x;
      x = y;
      y = oldX;
    }

    inline int signOfZ(Pole pole){
      switch(pole){
      case Pole::North: return -1;
      default:
      case Pole::South: return 1;
      }
    }


    /**
     * How the class actually works:
     *
     * The fX, fY, fZ in the base TVector class are the definitive 
     * representation of the position. In other words, any non-const 
     * functions ALWAYS update the fX, fY, fZ immediately
     * (fX, fY, fZ are the variable names in TVector3.)
     * However, the fLongitude, fLatitude, fAltitude, and fTheta, fPhi,
     * and fEasting, fNorthing are not updated, and are recalculated on
     * demand. When calculated they are stored and if fX, fY, fZ didn't change
     * The last answer will be returned again.
     *
     * The model for the object constness is that only changing the vector position
     * is non-const. Therefore the derived quantities (which are calculated and
     * stored on demand) are mutable.
     *
     * The derived coordiante caches are non-ROOT persistent,
     * (marked by comments after variable beginning with //!).
     * i.e. they won't be stored in a TTree.
     */

    void updateGeoidFromCartesian() const;
    void updateAnglesFromCartesian() const;
    void updateEastingNorthingFromLonLat() const;

    // can't make const as the cartesian x,y,z are represented in TVector3
    // and accessors aren't virtual and so can't be overloaded. 
    void updateLonLatFromEastingNorthing(bool mustRecalcuateAltitudeFirst);
    void updateCartesianFromGeoid();

    mutable Double_t fLongitude = 0; //! cached longitude is not stored!
    mutable Double_t fLatitude  = 0; //! cached latitude is not stored!
    mutable Double_t fAltitude  = 0; //! cached altitude is not stored!
    mutable Double_t fTheta     = 0; //! cached angle is not stored!
    mutable Double_t fPhi       = 0; //! cached angle is not stored!
    mutable Double_t fEasting   = 0; //! cached polar coordinate is not stored!
    mutable Double_t fNorthing  = 0; //! cached polar coordinate is not stored!

    mutable Double_t fCartAtLastGeoidCalc[3]       = {-1, -1, -1};    //! fX, fY, fZ when the geoid was last updated. Is not stored!
    mutable Double_t fCartAtLastAngleCal[3]        = {0, 0, 0};       //! fX, fY, fZ when the angles were last updated. Is not stored!
    mutable Double_t fLonLatAtLastEastNorthCalc[2] = {-9999, -9999};  //! Longitude(), Latitude() when the Easting/Northing were last calculated. Is not stored!
  };






















  template <class T> Vector::Vector(const T& t){
    SetLonLatAlt(t.longitude,  t.latitude, t.altitude);
    updateCartesianFromGeoid();
  };

  template <class T> Vector::Vector(const T* t) : Vector(*t){;}


  inline void Vector::SetLongitude(double lon) {
    if(lon > 180){
      lon -= 360;
    }
    fLongitude = lon;
    updateCartesianFromGeoid();
  }
  inline void Vector::SetLatitude(double lat){
    fLatitude = lat;
    updateCartesianFromGeoid();
  }
  inline void Vector::SetAltitude(double alt) {
    fAltitude = alt;
    updateCartesianFromGeoid();
  }
  inline void Vector::SetLonLatAlt(double lon, double lat, double alt) {
    fLatitude = lat;
    fAltitude = alt;
    SetLongitude(lon);
  }


  inline void Vector::SetEasting(double easting) {
    fNorthing = Northing();
    fEasting = easting;
    updateLonLatFromEastingNorthing(false);
  }
  inline void Vector::SetNorthing(double northing) {    
    fEasting = Easting();
    fNorthing = northing;
    updateLonLatFromEastingNorthing(false);
  }
  inline void Vector::SetEastingNorthing(double easting, double northing) {
    fEasting = easting;
    fNorthing = northing;
    updateLonLatFromEastingNorthing(true);
  }
  inline void Vector::SetEastingNorthingAlt(double easting, double northing, double alt) {
    fEasting = easting;
    fNorthing = northing;
    fAltitude = alt;
    updateLonLatFromEastingNorthing(false);
  }



  inline Double_t Vector::Latitude() const {
    updateGeoidFromCartesian();
    return fLatitude;
  }
  inline Double_t Vector::Longitude() const {
    updateGeoidFromCartesian();
    return fLongitude;
  }
  inline Double_t Vector::Altitude() const {
    updateGeoidFromCartesian();
    return fAltitude;
  }

  inline Double_t Vector::Theta() const {
    updateAnglesFromCartesian();
    return fTheta;
  }
  inline Double_t Vector::Phi() const {
    updateAnglesFromCartesian();
    return fPhi;
  }


  inline Double_t Vector::Easting() const {
    updateEastingNorthingFromLonLat();
    return fEasting;
  }

  inline Double_t Vector::Northing() const {
    updateEastingNorthingFromLonLat();
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




  inline double Vector::surfaceZ(Pole pole){
    // ellipse defined by: p^{2}/(geoid_max^{2}) + z^{2}/(geoid_min^{2}) = 1
    const double pSq = X()*X() + Y()*Y(); //lateral width of the geoid
    const double zSq = GEOID_MIN*GEOID_MIN*(1 - pSq/(GEOID_MAX*GEOID_MAX));
    if(zSq < 0){
      std::cerr <<  "Error in" << __PRETTY_FUNCTION__ << " can't find z if outside Geoid in x/y plane! "
		<< " x = " << X() << ", y = " << Y() << ", GEOID_MAX = " << GEOID_MAX << std::endl;
      return TMath::QuietNaN();
    }
    else {
      return signOfZ(pole)*TMath::Sqrt(zSq);
    }
  }

  inline void Vector::moveToGeoidZ(Pole pole){
    SetZ(surfaceZ(pole));
  }








  inline Double_t getGeoidRadiusAtCosTheta(Double_t cosTheta) {
    return GEOID_MIN*GEOID_MAX/TMath::Sqrt(GEOID_MIN*GEOID_MIN-(GEOID_MIN*GEOID_MIN-GEOID_MAX*GEOID_MAX)*cosTheta*cosTheta);
  }
  inline Double_t getGeoidRadiusAtTheta(Double_t theta) {
    return getGeoidRadiusAtCosTheta(TMath::Cos(theta));
  }



  
}


#endif
