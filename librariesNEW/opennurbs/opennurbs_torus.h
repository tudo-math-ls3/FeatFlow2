/* $Header: /src4/opennurbs/opennurbs_torus.h 2     12/13/04 1:35p Dalelear $ */
/* $NoKeywords: $ */
/*
//
// Copyright (c) 1993-2001 Robert McNeel & Associates. All rights reserved.
// Rhinoceros is a registered trademark of Robert McNeel & Assoicates.
//
// THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY.
// ALL IMPLIED WARRANTIES OF FITNESS FOR ANY PARTICULAR PURPOSE AND OF
// MERCHANTABILITY ARE HEREBY DISCLAIMED.
//				
// For complete openNURBS copyright information see <http://www.opennurbs.org>.
//
////////////////////////////////////////////////////////////////
*/

#if !defined(ON_TORUS_INC_)
#define ON_TORUS_INC_

class ON_RevSurface;
class ON_TextLog;

/*
Description:
  The torus is defined by a major circle and minor radius.  The
  torus is parameterized by (major_angle,minor_angle).  The angles
  are specified in radians.  The domain of both parameters is (0,2pi).
*/
class ON_CLASS ON_Torus
{

public:
  // for expert users

  ON_Plane plane; // major circle plane
  double major_radius;  // > minor_radius
  double minor_radius;  // > 0

public:

  ON_Torus();
  ON_Torus( const ON_Plane& major__plane, double major__radius, double minor__radius );
  ON_Torus( const ON_Circle& major__circle, double minor__radius );
  ~ON_Torus();

  BOOL IsValid( ON_TextLog* text_log = NULL ) const;

  BOOL Create( const ON_Plane& major__plane, double major__radius, double minor__radius );
  BOOL Create( const ON_Circle& major__circle, double minor__radius);

  /*
  Description:
    Get the circle that is the isocurve on the torus
    at a specified minor angle.
  Parameteters:
    minor_angle_radians - [in]
  Returns:
    A circle with normal major_circle.plane.zaxis that starts
    at PointAt( 0.0, minor_angle_radians ).
  See Also:
    ON_Torus::MajorCircleRadians
    ON_Torus::MajorCircleDegrees
    ON_Torus::MinorCircleRadians
    ON_Torus::MinorCircleDegrees
  */
  ON_Circle MajorCircleRadians(double minor_angle_radians ) const;

  /*
  Description:
    Get the circle that is the isocurve on the torus
    at a specified minor angle.
  Parameteters:
    minor_angle_degrees - [in]
  Returns:
    A circle with normal major_circle.plane.zaxis that starts
    at PointAt( 0.0, minor_angle_degrees*ON_PI/180.0 ).
  See Also:
    ON_Torus::MajorCircleRadians
    ON_Torus::MajorCircleDegrees
    ON_Torus::MinorCircleRadians
    ON_Torus::MinorCircleDegrees
  */
  ON_Circle MajorCircleDegrees(double minor_angle_degrees) const;

  /*
  Description:
    Get the minor circle that is the isocurve on the torus
    at a specified major angle.
  Parameteters:
    major_angle_radians - [in]
  Returns:
    A circle with radius = minor_radis, 
    center = major_circle.PointAt(major_angle_radians), and
    starting point PointAt( major_angle_radians, 0.0 ).
  See Also:
    ON_Torus::MajorCircleRadians
    ON_Torus::MajorCircleDegrees
    ON_Torus::MinorCircleRadians
    ON_Torus::MinorCircleDegrees
  */
  ON_Circle MinorCircleRadians(double major_angle_radians) const;

  /*
  Description:
    Get the minor circle that is the isocurve on the torus
    at a specified major angle.
  Parameteters:
    major_angle_degrees - [in]
  Returns:
    A circle with radius = minor_radis, 
    center = major_circle.PointAt(major_angle_degrees*ON_PI/180.0), and
    starting point PointAt( major_angle_degrees*ON_PI/180.0, 0.0 ).
  See Also:
    ON_Torus::MajorCircleRadians
    ON_Torus::MajorCircleDegrees
    ON_Torus::MinorCircleRadians
    ON_Torus::MinorCircleDegrees
  */
  ON_Circle MinorCircleDegrees(double major_angle_degrees) const;

  ON_3dPoint Center() const;
  ON_3dVector Axis() const;
  double MajorRadius() const;
  double MinorRadius() const;

  ON_3dPoint PointAt(
    double major_angle_radians, 
    double minor_angle_radians
    ) const;

  ON_3dVector NormalAt(
    double major_angle_radians, 
    double minor_angle_radians
    ) const;

  // returns parameters of point on torus that is closest to test_point.
  BOOL ClosestPointTo( 
         ON_3dPoint test_point, 
         double* major_angle_radians, 
         double* minor_angle_radians
         ) const;

  // returns point on torus that is closest to test_point
  ON_3dPoint ClosestPointTo( 
         ON_3dPoint test_point
         ) const;

  // rotate torus about its origin
  BOOL Rotate(
        double sin_angle,               // sin(angle)
        double cos_angle,               // cos(angle)
        const ON_3dVector& axis_of_rotation // axis of rotation
        );

  BOOL Rotate(
        double angle_radians,               // angle in radians
        const ON_3dVector& axis_of_rotation // axis of rotation
        );

  // rotate torus about a point and axis
  BOOL Rotate(
        double sin_angle,               // sin(angle)
        double cos_angle,               // cos(angle)
        const ON_3dVector& axis_of_rotation, // axis of rotation
        const ON_3dPoint& center_of_rotation  // center of rotation
        );

  BOOL Rotate(
        double angle_radians,               // angle in radians
        const ON_3dVector& axis_of_rotation, // axis of rotation
        const ON_3dPoint& center_of_rotation  // center of rotation
        );

  BOOL Translate(
        const ON_3dVector&
        );

  BOOL Transform( const ON_Xform& );

  // parameterization of NURBS surface does not match torus's transcendental paramaterization
  int GetNurbForm( ON_NurbsSurface& ) const; // returns 0=failure, 2=success

  /*
  Description:
    Creates a surface of revolution definition of the torus.
  Parameters:
    srf - [in] if not NULL, then this srf is used.
  Result:
    A surface of revolution or NULL if the torus is not valid.
  */
  ON_RevSurface* RevSurfaceForm( ON_RevSurface* srf = NULL ) const;
};

#endif
