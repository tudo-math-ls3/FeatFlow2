#include "opennurbs.h"
//  class ON_HatchLine
/////////////////////////////////////////////////////////////////

ON_HatchLine::ON_HatchLine()
: m_angle( 0.0), m_base( 0.0,0.0), m_offset( 0.0, 1.0)
{
}

ON_HatchLine::ON_HatchLine(double angle, 
                           const ON_2dPoint& base, 
                           const ON_2dVector& offset,
                           const ON_SimpleArray<double> dashes)
: m_angle( angle), m_base( base), m_offset( offset), m_dashes( dashes)
{
}

bool ON_HatchLine::operator==(const ON_HatchLine& src) const
{
  return( m_angle == src.m_angle && 
          m_base == src.m_base &&
          m_offset == src.m_offset && 
          m_dashes == src.m_dashes);
}

bool ON_HatchLine::operator!=(const ON_HatchLine& src) const
{
  return !operator==( src);
}

BOOL ON_HatchLine::IsValid( ON_TextLog* text_log) const
{
  bool rc = m_angle >= 0.0;
  if( !rc)
  {
    if( text_log)
      text_log->Print( "Angle ( %lf) must be >= 0.0\n", m_angle);
    return false;
  }
  rc = m_angle < ON_PI * 2.0;
  if( !rc)
  {
    if( text_log)
      text_log->Print( "Angle ( %lf) must be < 2*Pi.\n", m_angle);
    return false;
  }
  rc = m_base != ON_2dPoint( ON_UNSET_VALUE, ON_UNSET_VALUE);
  if( !rc)
  {
    if( text_log)
      text_log->Print( "Base is not a valid point.\n");
    return false;
  }
  rc = m_offset.x != ON_UNSET_VALUE;
  if( !rc)
  {
    if( text_log)
      text_log->Print( "Offset is not a valid vector.\n");
    return false;
  }
  rc = m_offset.y > ON_SQRT_EPSILON;
  if( !rc)
  {
    if( text_log)
      text_log->Print( "Offset.y ( %lf) must be > 0.0", m_offset.y);
    return false;
  }
  return true;
}

void ON_HatchLine::Dump( ON_TextLog& dump) const
{
  dump.Print( "ON_HatchLine: angle = %lf radians ( %lf degrees) ", 
    Angle(), ON_RADIANS_TO_DEGREES * Angle());
  dump.Print( " base = ");
  dump.Print( m_base);
  dump.Print( " offset = ");
  dump.Print( m_offset);
  int count = m_dashes.Count();
  dump.Print( "\nDash count = %d: ", count);
  for( int i = 0; i < count; i++)
  {
    dump.Print( "%lf", Dash( i));
    if( i < count-1)
      dump.Print( ", ");
  }
  dump.Print( "\n");
}

BOOL ON_HatchLine::Write( ON_BinaryArchive& ar) const
{
  BOOL rc = ar.Write3dmChunkVersion(1,1);

  if (rc) rc = ar.WriteDouble( m_angle);
  if (rc) rc = ar.WritePoint( m_base);
  if (rc) rc = ar.WriteVector( m_offset);
  if (rc) rc = ar.WriteArray( m_dashes);

  return rc;
}

BOOL ON_HatchLine::Read( ON_BinaryArchive& ar)
{
  m_angle = 0.0;
  m_base.Set( 0.0, 0.0);
  m_offset.Set( 0.0, 1.0);
  m_dashes.Empty();
  int major_version = 0;
  int minor_version = 0;
  BOOL rc = ar.Read3dmChunkVersion( &major_version, &minor_version);
  if ( major_version == 1 ) 
  {
    if ( rc) rc = ar.ReadDouble( &m_angle);
    if ( rc) rc = ar.ReadPoint( m_base);
    if ( rc) rc = ar.ReadVector( m_offset);
    if ( rc) rc = ar.ReadArray( m_dashes);
  }
  return rc;
}

// ON_HatchLine Interface
double ON_HatchLine::Angle() const
{
  return m_angle;
}

void ON_HatchLine::SetAngle( double angle)
{
  m_angle = angle;
  double twopi = ON_PI * 2.0;

  // clamp between [0  2pi)
  while( m_angle < 0.0)
    m_angle += twopi;
  while( m_angle > twopi)
    m_angle -= twopi;
}

ON_2dPoint ON_HatchLine::Base() const
{
  return m_base;
}

void ON_HatchLine::SetBase( const ON_2dPoint& base)
{
  m_base = base;
}

ON_2dVector ON_HatchLine::Offset() const
{
  return m_offset;
}

void ON_HatchLine::SetOffset( const ON_2dVector& offset)
{
  m_offset = offset;
}

int ON_HatchLine::DashCount() const
{
  return m_dashes.Count();
}

double ON_HatchLine::Dash( int index) const
{
  if( index >= 0 && index < m_dashes.Count())
    return m_dashes[index];
  return 0.0;
}

void ON_HatchLine::AppendDash( double dash)
{
//  if( fabs( dash) > ON_SQRT_EPSILON)
    m_dashes.Append( dash);
}

void ON_HatchLine::SetPattern( const ON_SimpleArray<double>& dashes)
{
  m_dashes = dashes;
}

void ON_HatchLine::GetLineData( double& angle, 
                                ON_2dPoint& base, 
                                ON_2dVector& offset, 
                                ON_SimpleArray<double>& dashes) const
{
  angle = m_angle;
  base = m_base;
  offset = m_offset;  dashes = m_dashes;
}

double ON_HatchLine::GetPatternLength() const
{
  int i;
  double length = 0.0;
  for( i = 0; i < m_dashes.Count(); i++)
    length += fabs( m_dashes[i]);

  return length;
}


//  class ON_HatchPattern
/////////////////////////////////////////////////////////////////
ON_OBJECT_IMPLEMENT( ON_HatchPattern, ON_Object, "064E7C91-35F6-4734-A446-79FF7CD659E1" );

ON_HatchPattern::ON_HatchPattern()
{
  m_hatchpattern_index = -1;
  memset(&m_hatchpattern_id,0,sizeof(m_hatchpattern_id));
}

ON_HatchPattern::~ON_HatchPattern()
{
}

BOOL ON_HatchPattern::IsValid( ON_TextLog* text_log) const
{
  eFillType type = FillType();
  BOOL rc = true;
  if( type != ftSolid && type != ftLines && type != ftGradient)
  {
    if( text_log)
      text_log->Print( "Type field not set correctly.\n");
    rc = false;
  }
  if( type == ftLines)
  {
    int count = m_lines.Count();
    if( count < 1)
    {
      if( text_log)
        text_log->Print( "Line type patetern with no lines.\n");
      return false;
    }
    for( int i = 0; i < count; i++)
    {
      if( !m_lines[i].IsValid())
      {
        if( text_log)
          text_log->Print( "Line[%d] is not valid.\n", i);
        return false;
      }
    }
    return true;
  }
  return rc;
}

void ON_HatchPattern::Dump( ON_TextLog& dump) const
{
  dump.Print( "Hatch pattern - ");
  switch( m_type)
  {
  case ftSolid:
    dump.Print( "fill type: Solid");
    break;
  case ftLines:
    dump.Print( "fill type: Lines");
    break;
  case ftGradient:
    dump.Print( "fill type: Gradient");
    break;
  case ftLast:
    // no action, but this keeps gcc happy
    break;
  }
  dump.Print( "\n");

  const wchar_t* s = m_hatchpattern_name;
  if ( 0 == s )
    s = L"";
  dump.Print( L"Name: %s\n", s);

  s =  m_description;
  if ( 0 == s )
    s = L"";
  dump.Print( L"Description: %s\n", s);

  if( m_type == ftLines)
  {
    int count = m_lines.Count();
    dump.Print( "Line count = %d\n", count);
    for( int i = 0; i < count; i++)
    {
      m_lines[i].Dump( dump);
    }
    dump.Print( L"\n");
  }
}

BOOL ON_HatchPattern::Write( ON_BinaryArchive& ar) const
{
  BOOL rc = ar.Write3dmChunkVersion(1,2);

  if (rc) rc = ar.WriteInt( m_hatchpattern_index);
  if (rc) rc = ar.WriteInt( m_type);
  if (rc) rc = ar.WriteString( m_hatchpattern_name);
  if (rc) rc = ar.WriteString( m_description);
  if( rc)
  {
    if( m_type == ftLines)
    {
      int i, count = m_lines.Count();
      if ( count < 0 )
        count = 0;
      rc = ar.WriteInt( count );
      for( i = 0; i < count && rc; i++)
        rc = m_lines[i].Write( ar);
    }
  }
  // version 1.2 field
  if (rc) rc = ar.WriteUuid(m_hatchpattern_id);

  return rc;
}

BOOL ON_HatchPattern::Read( ON_BinaryArchive& ar)
{
  m_hatchpattern_index = -1;
  memset(&m_hatchpattern_id,0,sizeof(m_hatchpattern_id));
  m_type = ftSolid;
  m_hatchpattern_name.Empty();
  m_description.Empty();
  m_lines.Empty();
  int i;

  int major_version = 0;
  int minor_version = 0;
  BOOL rc = ar.Read3dmChunkVersion( &major_version, &minor_version);
  if ( major_version == 1 ) 
  {
    if( rc) rc = ar.ReadInt( &m_hatchpattern_index);
    if( rc) rc = ar.ReadInt( &i);
    if( rc) 
    {
      switch( i)
      {
      case 0:  m_type = ftSolid;    break;
      case 1:  m_type = ftLines;    break;
      case 2:  m_type = ftGradient; break;
      default: rc = false;          break;
      }
    }
    if( rc) rc = ar.ReadString( m_hatchpattern_name);
    if( rc) rc = ar.ReadString( m_description);
    if( rc)
    {
      if( m_type == ftLines)
      {
        m_lines.Empty();
        int count = 0;
        rc = ar.ReadInt( &count);
        if( rc && count > 0 ) 
        {
          m_lines.SetCapacity( count);
          int i;
          for( i = 0; rc && i < count; i++)
          {
            ON_HatchLine& line = m_lines.AppendNew();
            rc = line.Read( ar);
          }
        }
      }
    }
    if ( minor_version >= 2 )
    {
      rc = ar.ReadUuid(m_hatchpattern_id);
    }
  }
  return rc;
}

ON_HatchPattern::eFillType ON_HatchPattern::FillType() const
{
  if( m_type >= ftSolid && m_type < ftLast)
    return m_type;

  return ftLast;
}

void ON_HatchPattern::SetFillType( eFillType type)
{
  m_type = type;
}

void ON_HatchPattern::SetName( const wchar_t* pName)
{
  m_hatchpattern_name = pName;
  m_hatchpattern_name.TrimLeftAndRight();
}

void ON_HatchPattern::SetName( const char* pName)
{
  m_hatchpattern_name = pName;
  m_hatchpattern_name.TrimLeftAndRight();
}

void ON_HatchPattern::GetName( ON_wString& string) const
{
  string = m_hatchpattern_name;
}

const wchar_t* ON_HatchPattern::Name() const
{
  return m_hatchpattern_name;
}


void ON_HatchPattern::SetDescription( const wchar_t* pDescription)
{
  m_description = pDescription;
}

void ON_HatchPattern::SetDescription( const char* pDescription)
{
  m_description = pDescription;
}

void ON_HatchPattern::GetDescription( ON_wString& string) const
{
  string = m_description;
}

const wchar_t* ON_HatchPattern::Description() const
{
  return m_description;
}


void ON_HatchPattern::SetIndex( int i)
{
  m_hatchpattern_index = i;
}

int ON_HatchPattern::Index() const
{
  return m_hatchpattern_index;
}


//  Line HatchPattern functions

int ON_HatchPattern::HatchLineCount() const
{
  return m_lines.Count();
}

int ON_HatchPattern::AddHatchLine( const ON_HatchLine& line)
{
  m_lines.Append( line);
  return m_lines.Count()-1;
}

const ON_HatchLine* ON_HatchPattern::HatchLine( int index) const
{
  if( index >= 0 && index < m_lines.Count())
    return &m_lines[index];
  else
    return NULL;
}

bool ON_HatchPattern::RemoveHatchLine( int index)
{
  if( index >= 0 && index < m_lines.Count())
  {
    m_lines.Remove( index);
    return true;
  }
  return false;
}

void ON_HatchPattern::RemoveAllHatchLines()
{
  m_lines.Empty();
}

int ON_HatchPattern::SetHatchLines( const ON_ClassArray<ON_HatchLine> lines)
{
  m_lines = lines;
  return m_lines.Count();
}





//  class ON_HatchLoop
/////////////////////////////////////////////////////////////////

ON_HatchLoop::ON_HatchLoop()
: m_type( ON_HatchLoop::ltOuter), m_p2dCurve( NULL)
{
}

ON_HatchLoop::ON_HatchLoop( ON_Curve* pCurve2d, eLoopType type)
: m_type( type), m_p2dCurve( pCurve2d)
{
}

ON_HatchLoop::ON_HatchLoop( const ON_HatchLoop& src)
: m_type( src.m_type), m_p2dCurve( NULL)
{ 
  if( src.m_p2dCurve)
    m_p2dCurve = src.m_p2dCurve->DuplicateCurve();
}

ON_HatchLoop::~ON_HatchLoop()
{
  delete m_p2dCurve;
}

ON_HatchLoop& ON_HatchLoop::operator=( const ON_HatchLoop& src)
{
  if( this != &src)
  {
    if( m_p2dCurve)
      delete m_p2dCurve;
    m_p2dCurve = src.m_p2dCurve->DuplicateCurve();

    m_type = src.m_type;
  }
  return *this;
}

BOOL ON_HatchLoop::IsValid( ON_TextLog* text_log) const
{
  BOOL rc = m_p2dCurve != NULL;
  if( !rc)
  {
    if( text_log)
      text_log->Print( "2d loop curve is NULL\n");
  }
  if( rc)
  {
    rc = m_p2dCurve->IsValid( text_log);
    if( !rc)
    {
      if( text_log)
        text_log->Print( "Loop curve is not valid\n");
    }
  }

  if( rc)
  {
    ON_BoundingBox box;
    m_p2dCurve->GetBoundingBox( box);
    rc = ( box.Max().z == box.Min().z && box.Max().z == 0.0);
    if( !rc)
    {
      if( text_log)
        text_log->Print( "2d loop curve has non-zero z coordinates\n");
    }
  }

  if( rc && m_type != ltOuter && m_type != ltInner)
  {
    if( text_log)
      text_log->Print( "Loop type is invalid.\n");
    rc = false;
  }

  return rc;
}

void ON_HatchLoop::Dump( ON_TextLog& dump) const
{
  if( m_type == ltOuter)
    dump.Print( "Loop type: Outer\n");
  if( m_type == ltInner)
    dump.Print( "Loop type: Inner\n");
  dump.Print( "2d curve: %p\n", m_p2dCurve);
}

BOOL ON_HatchLoop::Write( ON_BinaryArchive& ar) const
{
  BOOL rc = ar.Write3dmChunkVersion(1,1);
  if( rc) rc = ar.WriteInt( m_type);
  if( rc) rc = ar.WriteObject( m_p2dCurve);
  return rc;
}

BOOL ON_HatchLoop::Read( ON_BinaryArchive& ar)
{
  m_type = ltOuter;
  delete m_p2dCurve;
  m_p2dCurve = NULL;
  int major_version = 0;
  int minor_version = 0;
  BOOL rc = ar.Read3dmChunkVersion( &major_version, &minor_version);
  if ( major_version == 1 ) 
  {
    int type;
    if( rc) rc = ar.ReadInt( &type);
    if( rc) 
    {
      switch( type)
      {
      case ltOuter:  m_type = ltOuter; break;
      case ltInner:  m_type = ltInner; break;
      default: rc = false; break;
      }
    }
    if( rc)
    {
      ON_Object* pObj = NULL;
      rc = ar.ReadObject( &pObj);
      if( pObj)
      {
        m_p2dCurve = ON_Curve::Cast( pObj);
        if( !m_p2dCurve) // read something, but it wasn't right
        {
          rc = false;
          delete pObj;
        }
      }
    }
  }
  return rc;
}

const ON_Curve* ON_HatchLoop::Curve() const
{
  return m_p2dCurve;
}

bool ON_HatchLoop::SetCurve( const ON_Curve& curve)
{
  ON_Curve* pC = curve.DuplicateCurve();
  if( pC)
  {
    if( m_p2dCurve)
      delete m_p2dCurve;
    m_p2dCurve = pC;
  }
  return true;
}
ON_HatchLoop::eLoopType ON_HatchLoop::Type() const
{
  return m_type;
}

void ON_HatchLoop::SetType( eLoopType type)
{
  m_type = type;
}

//  class ON_Hatch
/////////////////////////////////////////////////////////////////
ON_OBJECT_IMPLEMENT( ON_Hatch, ON_Geometry, "0559733B-5332-49d1-A936-0532AC76ADE5");


ON_Hatch::ON_Hatch()
: m_pattern_scale( 1.0),
  m_pattern_rotation( 0.0),
  m_pattern_index( -1)
{
}

ON_Hatch::ON_Hatch( const ON_Hatch& src)
:  ON_Geometry(src),
   m_plane( src.m_plane), 
   m_pattern_scale( src.m_pattern_scale),
   m_pattern_rotation( src.m_pattern_rotation),
   m_pattern_index( src.m_pattern_index)
{
  m_loops.Reserve( src.m_loops.Count());
  for( int i = 0; i < src.m_loops.Count(); i++)
  {
    ON_HatchLoop* pL = new ON_HatchLoop( *src.m_loops[i]);
    m_loops.Append( pL);
  }
}

ON_Hatch& ON_Hatch::operator=( const ON_Hatch& src)
{
  if( this != &src)
  {
    // Nov 3 2004 Dale Lear:
    //   Delete existing loops so we don't leak the memory;
    int i;
    for ( i = 0; i < m_loops.Count(); i++ )
    {
      ON_HatchLoop* pL = m_loops[i];
      if ( pL )
      {
        m_loops[i] = 0;
        delete pL;
      }
    }
    m_loops.SetCount(0);

    ON_Geometry::operator =(src);

    m_plane = src.m_plane;
    m_pattern_index = src.m_pattern_index;
    m_pattern_scale = src.m_pattern_scale;
    m_pattern_rotation = src.m_pattern_rotation;
    m_loops.Reserve( src.m_loops.Count());
    for( i = 0; i < src.m_loops.Count(); i++)
    {
      ON_HatchLoop* pL = new ON_HatchLoop( *src.m_loops[i]);
      m_loops.Append( pL);
    }
  }
  return *this;
}

ON_Hatch::~ON_Hatch()
{
  int i;
  for ( i = 0; i < m_loops.Count(); i++ )
  {
    ON_HatchLoop* pL = m_loops[i];
    if ( pL )
    {
      m_loops[i] = 0;
      delete pL;
    }
  }
}


ON_Hatch* ON_Hatch::DuplicateHatch() const
{
  return Duplicate();
}

BOOL ON_Hatch::IsValid( ON_TextLog* text_log) const
{
  BOOL rc = m_plane.IsValid();
  if( !rc)
  {
    if( text_log)
      text_log->Print( "Plane is not valid\n");
    return false;
  }
  int i;
  int count = m_loops.Count();
  for( i = 0; i < count; i++)
  {
    rc = m_loops[i]->IsValid( text_log);
    if( !rc)
    {
      if( text_log)
        text_log->Print( "Loop[%d] is not valid\n", i);
      return false;
    }
  }
  return true;
}

void ON_Hatch::Dump( ON_TextLog& dump) const
{
  dump.Print( "Hatch: Solid fill");
  int count = m_loops.Count();
  dump.Print( "Loop count = %d\n", count);
  for( int i = 0; i < count; i++)
    m_loops[i]->Dump( dump);
}

BOOL ON_Hatch::Write( ON_BinaryArchive& ar) const
{
  BOOL rc = ar.Write3dmChunkVersion(1,1);
  if (rc) rc = ar.WritePlane( m_plane);
  if (rc) rc = ar.WriteDouble( m_pattern_scale);
  if (rc) rc = ar.WriteDouble( m_pattern_rotation);
  if (rc) rc = ar.WriteInt( m_pattern_index);
  if (rc)
  {
    int i, count = m_loops.Count();
    if( count < 0 )
      count = 0;
    BOOL rc = ar.WriteInt( count);
    for( i = 0; i < count && rc; i++)
      rc = m_loops[i]->Write( ar);
  }
  return rc;
}

BOOL ON_Hatch::Read( ON_BinaryArchive& ar)
{
  m_plane.CreateFromNormal( ON_origin, ON_zaxis);
  m_pattern_scale = 1.0;
  m_pattern_rotation = 0.0;
  m_pattern_index = -1;
  m_loops.Empty();
  int major_version = 0;
  int minor_version = 0;
  BOOL rc = ar.Read3dmChunkVersion( &major_version, &minor_version);
  if ( major_version == 1 ) 
  {
    if( rc) rc = ar.ReadPlane( m_plane);
    if( rc) rc = ar.ReadDouble( &m_pattern_scale);
    if( rc) rc = ar.ReadDouble( &m_pattern_rotation);
    if( rc) rc = ar.ReadInt( &m_pattern_index);
    if( rc)
    {
      m_loops.Empty();
      int i, count = 0;
      rc = ar.ReadInt( &count);
      if( rc && count > 0)
      {
        m_loops.SetCapacity( count );
        for( i = 0; rc && i < count; i++)
        {
          ON_HatchLoop*& pLoop = m_loops.AppendNew();
          pLoop = new ON_HatchLoop;
          if( pLoop)
            rc = pLoop->Read( ar);
          else
            rc = false;
        }
      }
    }
  }
  return rc;
}

ON::object_type ON_Hatch::ObjectType() const
{
  return ON::hatch_object;
}

int ON_Hatch::Dimension() const
{
  return 3;
}

// Copy the 2d curve, make it 3d, and transform it 
// to the 3d plane position
ON_Curve* ON_Hatch::LoopCurve3d( int index) const
{
  int count = m_loops.Count();
  ON_Curve* pC = NULL;

  if( index >= 0 && index < count)
  {
    if( m_loops[index]->Curve())
    {
      pC = m_loops[index]->Curve()->DuplicateCurve();
      if( pC)
      {
        pC->ChangeDimension( 3);

        ON_Xform xf;
        xf.Rotation( ON_xy_plane, m_plane);

        pC->Transform( xf);
      }
    }
  }
  return pC;
}


int ON_Hatch::PatternIndex() const
{
  return m_pattern_index;
}

void ON_Hatch::SetPatternIndex( int index)
{
  m_pattern_index = index;
}


BOOL ON_Hatch::GetBBox( double* bmin, double* bmax, BOOL bGrowBox) const
{
  int i;
  int count = m_loops.Count();
  BOOL rc = true;
  ON_Curve* pC;
  for( i = 0; rc && i < count; i++)
  {
    pC = LoopCurve3d( i);
    if( pC)
    {
      rc = pC->GetBBox( bmin, bmax, i?true:bGrowBox);
      delete pC;
    }
  }
  return rc;
}

bool ON_Hatch::GetTightBoundingBox( ON_BoundingBox& tight_bbox, int bGrowBox, const ON_Xform* xform) const
{
  int i;
  int count = m_loops.Count();
  ON_CurveArray curves(count);
  for( i = 0; i < count; i++)
  {
    curves.Append( LoopCurve3d(i) );
  }
  return curves.GetTightBoundingBox(tight_bbox,bGrowBox,xform);
}

BOOL ON_Hatch::Transform( const ON_Xform& xform)
{
  if( fabs( fabs( xform.Determinant()) - 1.0) > 1.0e-4)
  {
    // xform has a scale component
    ON_Plane tmp( m_plane);
    tmp.Transform( xform);
    ON_Xform A, B, T;
    A.Rotation( ON_xy_plane, m_plane);
    B.Rotation( tmp, ON_xy_plane);
    T = B * xform * A;

    // kill translation and z-scaling
    T[0][2] = T[0][3] = 0.0;
    T[1][2] = T[1][3] = 0.0;
    T[2][0] = T[2][1] = 0.0; T[2][2] = 1.0; T[2][3] = 0.0; 
    T[3][0] = T[3][1] = T[3][2] = 0.0; T[3][3] = 1.0;

    for( int i = 0; i < LoopCount(); i++)
      m_loops[i]->m_p2dCurve->Transform( T);
  }
  int rc = m_plane.Transform( xform);

  return rc;
}

bool ON_Hatch::Create( const ON_Plane& plane,
                       const ON_SimpleArray<const ON_Curve*> loops, 
                       int pattern_index, 
                       double pattern_rotation, 
                       double pattern_scale)
{
  if( loops.Count() < 1)
    return false;
  if( pattern_index < 0)
    return false;
  SetPlane( plane);
  for( int i = 0; i < loops.Count(); i++)
  {
    ON_HatchLoop* pLoop = new ON_HatchLoop;
    pLoop->SetCurve( *loops[i]);
    pLoop->SetType( i?ON_HatchLoop::ltInner:ON_HatchLoop::ltOuter);
    AddLoop( pLoop);
  }
  SetPatternIndex( pattern_index);
  SetPatternRotation( pattern_rotation);
  SetPatternScale( pattern_scale);
  return true;
}

const ON_Plane& ON_Hatch::Plane() const
{
  return m_plane;
}

void ON_Hatch::SetPlane( const ON_Plane& plane)
{
  m_plane = plane;
}

double ON_Hatch::PatternRotation() const
{
  return m_pattern_rotation;
}

void ON_Hatch::SetPatternRotation( double rotation)
{
  m_pattern_rotation = rotation;
}

double ON_Hatch::PatternScale() const
{
  return m_pattern_scale;
}

void ON_Hatch::SetPatternScale( double scale)
{
  if( scale > ON_SQRT_EPSILON)
    m_pattern_scale = scale;
}

int ON_Hatch::LoopCount() const
{
  return m_loops.Count();
}

void ON_Hatch::AddLoop( ON_HatchLoop* pLoop)
{
  m_loops.Append( pLoop);
}

const ON_HatchLoop* ON_Hatch::Loop( int index) const
{
  if( index >= 0 && index < m_loops.Count())
    return m_loops[index];
  
  return NULL;
}

