/* $Header: /src4/opennurbs/opennurbs_light.cpp 6     9/01/06 4:35p Dalelear $ */
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

#include "opennurbs.h"

ON_OBJECT_IMPLEMENT(ON_Light,ON_Geometry,"85A08513-F383-11d3-BFE7-0010830122F0");

void ON_Light::Default()
{
  m_light_name.Destroy();
  m_bOn = 1;
  m_intensity = 1.0;
  m_watts = 0.0;
  m_style =  ON::camera_directional_light;
  m_ambient = ON_Color(0,0,0);
  m_diffuse = ON_Color(255,255,255);
  m_specular = ON_Color(255,255,255);
  m_direction = ON_3dVector(0.0,0.0,-1.0);
  m_location = ON_3dPoint(0.0,0.0,0.0);
  m_length = ON_3dVector(0.0,0.0,0.0);
  m_width = ON_3dVector(0.0,0.0,0.0);
  m_spot_angle = 180.0;
  m_spot_exponent = 0.0;
  m_hotspot = 1.0;
  m_attenuation = ON_3dVector(1.0,0.0,0.0);
  m_shadow_intensity = 1.0;
  m_light_index = 0;
  memset(&m_light_id,0,sizeof(m_light_id));
}

ON_Light::ON_Light() 
{
  Default();
}

ON_Light::~ON_Light()
{
}

BOOL ON_Light::IsValid( ON_TextLog* text_log ) const
{
  int s = Style();
  if ( s <= ON::unknown_light_style || s >= ON::light_style_count ) {
    ON_ERROR("ON_Light::IsValid(): illegal light style.");
    return FALSE;
  }
  return TRUE;
}

void ON_Light::Dump( ON_TextLog& dump ) const
{
  BOOL bDumpDir = FALSE;
  BOOL bDumpLength = FALSE;
  BOOL bDumpWidth = FALSE;

  const char* sStyle = "unknown";
  switch(Style())
  {
  //case ON::view_directional_light:
  //  sStyle = "view_directional_light";
  //  bDumpDir = TRUE;
  //  break;
  //case ON::view_point_light:
  //  sStyle = "view_point_light";
  //  break;
  //case ON::view_spot_light:
  //  sStyle = "view_spot_light";
  //  bDumpDir = TRUE;
  //  break;
  case ON::camera_directional_light:
    sStyle = "camera_directional_light";
    bDumpDir = TRUE;
    break;
  case ON::camera_point_light:
    sStyle = "camera_point_light";
    break;
  case ON::camera_spot_light:
    sStyle = "camera_spot_light";
    bDumpDir = TRUE;
    break;
  case ON::world_directional_light:
    sStyle = "world_directional_light";
    bDumpDir = TRUE;
    break;
  case ON::world_point_light:
    sStyle = "world_point_light";
    break;
  case ON::world_spot_light:
    sStyle = "world_spot_light";
    bDumpDir = TRUE;
    break;
  case ON::world_linear_light:
    sStyle = "linear_light";
    bDumpDir = TRUE;
    bDumpLength = TRUE;
    break;
  case ON::world_rectangular_light:
    sStyle = "rectangular_light";
    bDumpDir = TRUE;
    bDumpLength = TRUE;
    bDumpWidth = TRUE;
    break;
  case ON::ambient_light:
    sStyle = "ambient_light";
    break;
  case ON::unknown_light_style:
    sStyle = "unknown";
    break;
  default:
    sStyle = "unknown";
    break;
  }
  dump.Print("index = %d  style = %s\n",LightIndex(),sStyle);

  dump.Print("location = "); dump.Print(Location()); dump.Print("\n");
  if ( bDumpDir )
    dump.Print("direction = "); dump.Print(Direction()); dump.Print("\n");
  if ( bDumpLength )
    dump.Print("length = "); dump.Print(Length()); dump.Print("\n");
  if ( bDumpWidth )
    dump.Print("width = "); dump.Print(Width()); dump.Print("\n");

  dump.Print("intensity = %g%%\n",Intensity()*100.0);

  dump.Print("ambient rgb = "); dump.PrintRGB(Ambient()); dump.Print("\n");
  dump.Print("diffuse rgb = "); dump.PrintRGB(Diffuse()); dump.Print("\n");
  dump.Print("specular rgb = "); dump.PrintRGB(Specular()); dump.Print("\n");

  dump.Print("spot angle = %g degrees\n",SpotAngleDegrees());
}

BOOL ON_Light::Write(
       ON_BinaryArchive& file
     ) const
{
  int i;
  BOOL rc = file.Write3dmChunkVersion(1,2);
  // version 1.0 fields
  if ( rc ) rc = file.WriteInt( m_bOn );
  i = m_style;
  if ( rc ) rc = file.WriteInt( i );
  if ( rc ) rc = file.WriteDouble( m_intensity );
  if ( rc ) rc = file.WriteDouble( m_watts );
  if ( rc ) rc = file.WriteColor( m_ambient );
  if ( rc ) rc = file.WriteColor( m_diffuse );
  if ( rc ) rc = file.WriteColor( m_specular );
  if ( rc ) rc = file.WriteVector( m_direction );
  if ( rc ) rc = file.WritePoint( m_location );
  if ( rc ) rc = file.WriteDouble( m_spot_angle );
  if ( rc ) rc = file.WriteDouble( m_spot_exponent );
  if ( rc ) rc = file.WriteVector( m_attenuation );
  if ( rc ) rc = file.WriteDouble( m_shadow_intensity );
  if ( rc ) rc = file.WriteInt( m_light_index );
  if ( rc ) rc = file.WriteUuid( m_light_id );
  if ( rc ) rc = file.WriteString( m_light_name );
  // version 1.1 added support for linear and rectangular
  if ( rc ) rc = file.WriteVector( m_length );
  if ( rc ) rc = file.WriteVector( m_width );
  // version 1.2 added m_hotspot support
  if ( rc ) rc = file.WriteDouble( m_hotspot );
  return rc;
}

BOOL ON_Light::Read(
       ON_BinaryArchive&  file
     )
{
  Default();
  int major_version = 0;
  int minor_version = 0;
  BOOL rc = file.Read3dmChunkVersion(&major_version,&minor_version);
  if ( rc && major_version == 1 ) {
    int i;
    // version 1.0 fields
    if ( rc ) rc = file.ReadInt( &i );
    if ( rc ) Enable(i);
    if ( rc ) rc = file.ReadInt( &i );
    if ( rc ) SetStyle(ON::LightStyle(i));
    if ( rc ) rc = file.ReadDouble( &m_intensity );
    if ( rc ) rc = file.ReadDouble( &m_watts );
    if ( rc ) rc = file.ReadColor( m_ambient );
    if ( rc ) rc = file.ReadColor( m_diffuse );
    if ( rc ) rc = file.ReadColor( m_specular );
    if ( rc ) rc = file.ReadVector( m_direction );
    if ( rc ) rc = file.ReadPoint( m_location );
    if ( rc ) rc = file.ReadDouble( &m_spot_angle );
    if ( rc ) rc = file.ReadDouble( &m_spot_exponent );
    if ( rc ) rc = file.ReadVector( m_attenuation );
    if ( rc ) rc = file.ReadDouble( &m_shadow_intensity );
    if ( rc ) rc = file.ReadInt( &m_light_index );
    if ( rc ) rc = file.ReadUuid( m_light_id );
    if ( rc ) rc = file.ReadString( m_light_name );

    if ( minor_version < 2 ) {
      // set hotspot from 1.0 or 1.1 m_spot_exponent
      double h = 1.0 - m_spot_exponent/128.0;
      if ( h < 0.0 ) 
        h = 0.0;
      else if ( h > 1.0 )
        h = 1.0;
      m_hotspot = h;
      m_spot_exponent = 0.0;
    }

    if ( minor_version >= 1 ) {
      // version 1.1 fields
      if ( rc ) rc = file.ReadVector( m_length );
      if ( rc ) rc = file.ReadVector( m_width );
      if ( minor_version >= 2 ) {
        // version 1.2 fields
        if ( rc ) rc = file.ReadDouble( &m_hotspot );
      }
    }
  }
  return rc;
}

ON::object_type ON_Light::ObjectType() const
{
  return ON::light_object;
}

int ON_Light::Dimension() const
{
  return 3;
}

BOOL ON_Light::GetBBox( // returns TRUE if successful
       double* boxmin,    // boxmin[dim]
       double* boxmax,    // boxmax[dim]
       BOOL bGrowBox
       ) const
{
  bool rc = true;
  ON_3dPointArray points(16);

  switch(m_style)
  {
  case ON::camera_directional_light:
  case ON::world_directional_light:
    points.Append(m_location);
    points.Append(m_location+m_direction);
    break;

  case ON::camera_point_light:
  case ON::world_point_light:
    points.Append(m_location);
    break;

  case ON::camera_spot_light:
  case ON::world_spot_light:
    if ( m_spot_angle > 0.0 && m_spot_angle < 90.0 )
    {
      double r = m_direction.Length()*tan(ON_PI*m_spot_angle/180.0);
      ON_Circle c(ON_Plane(m_location+m_direction,m_direction),r);
      ON_BoundingBox cbox = c.BoundingBox();
      cbox.GetCorners( points );
    }
    else
    {
      points.Append(m_location+m_direction);
    }
    points.Append(m_location);
    break;

  case ON::ambient_light:
    points.Append(m_location);
    rc = false;
    break;
  
  case ON::world_linear_light:
    points.Append(m_location);
    points.Append(m_location+m_length);
    break;

  case ON::world_rectangular_light:
    points.Append(m_location);
    points.Append(m_location+m_length);
    points.Append(m_location+m_width);
    points.Append(m_location+m_width+m_length);
    break;

  default:
    rc = false;
    break;
  }

  if ( rc && points.Count() > 0 )
  {
     rc = ON_GetPointListBoundingBox( 3, 0, points.Count(), 3, 
                                      (double*)points.Array(), 
                                      boxmin, boxmax, 
                                      bGrowBox?true:false )
        ? true 
        : false;
  }

  return rc;
}

BOOL ON_Light::Transform( 
       const ON_Xform& xform
       )
{
  ON_3dVector v;
  double vlen;
  TransformUserData(xform);
  m_location = xform*m_location;
  
  v = xform*m_direction;
  vlen = v.Length();
  if ( vlen > 0.0 ) {
    m_direction = v;
  }
  
  v = xform*m_length;
  vlen = v.Length();
  if ( vlen > 0.0 ) {
    m_length = v;
  }
  
  v = xform*m_width;
  vlen = v.Length();
  if ( vlen > 0.0 ) {
    m_width = v;
  }
  return TRUE;
}

BOOL ON_Light::Enable(BOOL b)
{
  BOOL oldb = m_bOn;
  m_bOn = (b)?TRUE:FALSE;
  return oldb;
}

BOOL ON_Light::IsEnabled() const
{
  return m_bOn;
}

ON::light_style ON_Light::Style() const
{
  return m_style;
}

const BOOL ON_Light::IsPointLight() const
{
  BOOL rc;
  switch(m_style)
  {
  //case ON::view_point_light:
  case ON::camera_point_light:
  case ON::world_point_light:
    rc = TRUE;
    break;
  default:
    rc = FALSE;
    break;
  }
  return rc;
}

const BOOL ON_Light::IsDirectionalLight() const
{
  BOOL rc;
  switch(m_style)
  {
  //case ON::view_directional_light:
  case ON::camera_directional_light:
  case ON::world_directional_light:
    rc = TRUE;
    break;
  default:
    rc = FALSE;
    break;
  }
  return rc;
}

const BOOL ON_Light::IsSpotLight() const
{
  BOOL rc;
  switch(m_style)
  {
  //case ON::view_spot_light:
  case ON::camera_spot_light:
  case ON::world_spot_light:
    rc = TRUE;
    break;
  default:
    rc = FALSE;
    break;
  }
  return rc;
}

const BOOL ON_Light::IsLinearLight() const
{
  BOOL rc;
  switch(m_style)
  {
  //case ON::view_linear_light:
  //case ON::camera_linear_light:
  case ON::world_linear_light:
    rc = TRUE;
    break;
  default:
    rc = FALSE;
    break;
  }
  return rc;
}

const BOOL ON_Light::IsRectangularLight() const
{
  BOOL rc;
  switch(m_style)
  {
  //case ON::view_rectangular_light:
  //case ON::camera_rectangular_light:
  case ON::world_rectangular_light:
    rc = TRUE;
    break;
  default:
    rc = FALSE;
    break;
  }
  return rc;
}


void ON_Light::SetStyle(ON::light_style s )
{
  m_style = s;
}

void ON_Light::SetLightName( const char* s )
{
  m_light_name = s;
  m_light_name.TrimLeftAndRight();
}

void ON_Light::SetLightName( const wchar_t* s )
{
  m_light_name = s;
  m_light_name.TrimLeftAndRight();
}

const ON_wString& ON_Light::LightName() const
{
  return m_light_name;
}

//void ON_Light::SetLightUuid( const ON_UUID& uuid )
//{
//  // m_light_id is set and maintained by Rhino - if your
//  // plug-in is messing with this field, fix the plugin
//  m_light_id = uuid;
//}

//const ON_UUID& ON_Light::LightUuid() const
//{
//  // m_light_id is set and maintained by Rhino - if your
//  // plug-in is messing with this field, fix the plugin
//  return m_light_id;
//}


void ON_Light::SetAttenuation(double a,double b,double c)
{
  m_attenuation = ON_3dVector(a,b,c);;
}

void ON_Light::SetAttenuation(const ON_3dVector& att )
{
  m_attenuation = att;
}

ON_3dVector ON_Light::Attenuation() const
{
  return m_attenuation;
}

double ON_Light::Attenuation(double d) const
{
  // computes 1/(a[0] + d*a[1] + d^2*a[2]) where d = argument
  // returns 0 if a[0] + d*a[1] + d^2*a[2] <= 0
  double a = m_attenuation.x + d*(m_attenuation.y + d*m_attenuation.z);
  if ( a > 0.0 )
    a = 1.0/a;
  else
    a = 0.0;
  return a;
}

/////////////////////////////////////////////////////////
//
// spot light parameters (ignored for non-spot lights)
//
// angle = 0 to 90 degrees
// exponent = 0 to 128 (0=uniform, 128=high focus)
//
void ON_Light::SetSpotAngleRadians( double a )
{
  a *= 180.0/ON_PI;
  if ( a > 90.0 )
    m_spot_angle = 90.0;
  else if ( a > 0.0 )
    m_spot_angle = a;
}

double ON_Light::SpotAngleRadians() const
{
  return m_spot_angle*ON_PI/180.0;
}

void ON_Light::SetSpotAngleDegrees( double a )
{
  if ( a >= 90.0 )
    m_spot_angle = 90.0;
  else if ( a > 0.0 )
    m_spot_angle = a;
}

double ON_Light::SpotAngleDegrees() const
{
  return m_spot_angle;
}

// The spot exponent "e", hot spot "h", and spotlight cone angle "a"
// are mutually constrained by the formula
//   cos(h*angle)^e = hotspot_min
// where hotspot_min = value of spotlight exponential attenuation factor 
// at the hot spot radius.  hotspot_min must be >0, < 1, and should be >= 1/2;
//static double log_hotspot_min = log(0.5);
static double log_hotspot_min = log(0.70710678118654752440084436210485);

void ON_Light::SetSpotExponent( double e )
{
  // cos(h)^e = 0.5
  if ( e < 0.0 || !ON_IsValid(e) )
    m_spot_exponent = 0.0;
  else
    m_spot_exponent = e;
  m_hotspot = ON_UNSET_VALUE; // indicates hotspot should be computed from m_spot_exponent
}

void ON_Light::SetHotSpot( double h )
{
  if ( h == ON_UNSET_VALUE || !ON_IsValid(h) )
    m_hotspot = ON_UNSET_VALUE;
  else if ( h <= 0.0 ) 
    m_hotspot = 0.0;
  else if ( h >= 1.0 )
    m_hotspot = 1.0;
  else
    m_hotspot = h;
}

double ON_Light::SpotExponent() const
{
  double e = m_spot_exponent;
  if ( m_hotspot >= 0.0 && m_hotspot <= 1.0 ) {
    // spotlight is using hot spot interface
    double h = m_hotspot;
    if ( h < 0.015 )
      h = 0.015;
    if ( h >= 1.0 || m_spot_angle <= 0.0 || m_spot_angle > 90.0)
      e = 0.0;
    else {
      // compute SpotExponent() from  cos(h*angle)^e = hotspot_min
      double a, c;
      a = h*SpotAngleRadians();
      c = cos(a);
      if ( c <= 0.0 )
        e = 1.0;
      else {
        e = log_hotspot_min/log(c);
        if ( e < 0.0 )
          e = 0.0;
      }
    }
  }
  return e;
}

double ON_Light::HotSpot() const
{
  double h = m_hotspot;
  if ( h < 0.0 || h > 1.0 ) {
    // spotlight is using spot exponent interface
    if ( m_spot_exponent >= 65536.0 )
      h = 0.0;
    else if ( m_spot_exponent <= 0.0 || m_spot_angle <= 0.0 || m_spot_angle > 90.0 )
      h = 1.0;
    else {
      // compute HotSpot() from cos(h*angle)^e = hotspot_min
      double x, a, cos_ha;
      x = log_hotspot_min/m_spot_exponent; // note that x < 0.0
      if ( x < -690.0 ) {
        // prevent underflow.  cos_ha is close to zero so
        h = 1.0;
      }
      else 
      {
        cos_ha = exp(x);
        if (!ON_IsValid(cos_ha))  cos_ha =  0.0;
        else if ( cos_ha >  1.0 ) cos_ha =  1.0;
        else if ( cos_ha < -1.0 ) cos_ha = -1.0;
        a = SpotAngleRadians();
        h = acos(cos_ha)/a;
        if ( h < 0.0 )
          h = 0.0;
        else if ( h > 1.0 ) {
          // happens for smaller e
          h = 1.0;
        }
      }
    }
  }
  return h;
}

void ON_Light::SetLength( const ON_3dVector& v )
{
  m_length = v;
}

ON_3dVector ON_Light::Length() const
{
  return m_length;
}

void ON_Light::SetWidth( const ON_3dVector& v )
{
  m_width = v;
}

ON_3dVector ON_Light::Width() const
{
  return m_width;
}


void ON_Light::SetShadowIntensity(double s )
{
  if ( s < 0.0 )
    s = 0.0;
  else if (s > 1.0 )
    s = 1.0;
  m_shadow_intensity = s;
}

double ON_Light::ShadowIntensity() const
{
  return m_shadow_intensity;
}

void ON_Light::SetLightIndex( int i )
{
  m_light_index = i;
}

int ON_Light::LightIndex() const
{
  return m_light_index;
}

void ON_Light::SetAmbient(  ON_Color c )
{
  m_ambient = c;
}

void ON_Light::SetDiffuse(  ON_Color c )
{
  m_diffuse = c;
}

void ON_Light::SetSpecular( ON_Color c )
{
  m_specular = c;;
}

ON_Color ON_Light::Ambient() const
{
  return m_ambient;
}

ON_Color ON_Light::Diffuse() const
{
  return m_diffuse;
}

ON_Color ON_Light::Specular() const
{
  return m_specular;
}

ON::coordinate_system ON_Light::CoordinateSystem() const // determined by style
{
  ON::coordinate_system cs = ON::world_cs;
  switch( m_style ) {
  case ON::unknown_light_style:
    cs = ON::world_cs;
    break;
  //case ON::view_directional_light:
  //case ON::view_point_light:
  //case ON::view_spot_light:
  //  cs = ON::clip_cs;
  //  break;
  case ON::camera_directional_light:
  case ON::camera_point_light:
  case ON::camera_spot_light:
    cs = ON::camera_cs;
    break;
  case ON::world_directional_light:
  case ON::world_point_light:
  case ON::world_spot_light:
  case ON::world_linear_light:
  case ON::world_rectangular_light:
  case ON::ambient_light:
    cs = ON::world_cs;
    break;
  default:
    cs = ON::world_cs;
    break;
  }
  return cs;
}

BOOL ON_Light::GetLightXform( 
         const ON_Viewport& vp,
         ON::coordinate_system dest_cs, 
         ON_Xform& xform 
         ) const
{
  ON::coordinate_system src_cs = CoordinateSystem();
  return vp.GetXform( src_cs, dest_cs, xform );
}



void ON_Light::SetLocation( const ON_3dPoint& loc )
{
  m_location = loc;
}

void ON_Light::SetDirection( const ON_3dVector& dir )
{
  m_direction = dir;
}

ON_3dPoint ON_Light::Location() const
{
  return m_location;
}

ON_3dVector ON_Light::Direction() const
{
  return m_direction;
}


double ON_Light::Intensity() const
{
  // 0.0 = 0%  1.0 = 100%
  return m_intensity;
}

void ON_Light::SetIntensity(double v)
{
  if ( v <= 0.0 )
    m_intensity = 0.0;
  else if (v >= 1.0 )
    m_intensity = 1.0;
  else
    m_intensity = v;
}

double ON_Light::PowerWatts() const
{
  return m_watts;
}

double ON_Light::PowerLumens() const
{
  return m_watts/683.0;
}

double ON_Light::PowerCandela() const
{
  return m_watts/683.0;
}

void ON_Light::SetPowerWatts( double watts )
{
  m_watts = (watts > 0.0) ? watts : 0.0;
}

void ON_Light::SetPowerLumens( double lumens )
{
  m_watts = (lumens > 0.0) ? lumens*683.0 : 0.0;
}

void ON_Light::SetPowerCandela( double candela )
{
  m_watts = (candela > 0.0) ? candela*683.0 : 0.0;
}



