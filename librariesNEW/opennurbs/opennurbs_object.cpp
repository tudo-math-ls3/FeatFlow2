#include "opennurbs.h"


#if defined(ON_DLL_EXPORTS)

////////////////////////////////////////////////////////////////////////////////
//
// When openNURBS is used as a MS Windows DLL, it is possible 
// for new/delete to allocate memory in one executable and delete
// it in another.  Because MS Windows has incompatible memory 
// managers in its plethora of C libraries and the choice of which
// C library actually gets used depends on the code generation 
// options you choose,  we get lots of support questions asking
// about hard to trace crashes.
//
// If you are using openNURBS as a Windows DLL, you are sure you know
// what you are doing, and you promise never to ask for support, then
// feel free to delete these overrides.
//
//
#pragma message( " --- OpenNURBS overriding ON_Object new and delete" )

void* ON_Object::operator new(size_t sz)
{
  // all classes derived from ON_Object use this operator new
  // ON_Object new
  return onmalloc(sz);
}

void ON_Object::operator delete(void* p)
{
  // ON_Object delete
  onfree(p);
}

void* ON_Object::operator new[] (size_t sz)
{
  // ON_Object array new
  return onmalloc(sz);
}

void ON_Object::operator delete[] (void* p)
{
  // ON_Object array delete
  onfree(p);
}

void* ON_Object::operator new(size_t, void* p)
{
  // ON_Object placement new
  return p;
}

void ON_Object::operator delete(void*, void*)
{
  // ON_Object placement delete
  return;
}

//
//
////////////////////////////////////////////////////////////////////////////////

#endif

const ON_UUID ON_nil_uuid = {0,0,0,{0,0,0,0,0,0,0,0}};
const ON_UUID ON_max_uuid = {0xFFFFFFFF,0xFFFF,0xFFFF,{0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF}};

void ON::End()
{
}

void ON::Begin()
{
  // By explicityly calling all the ON_Object::Cast overrides,
  // we can insure that the class definitions get linked in
  // by making a single call to ON::Begin().  Otherwise, when
  // openNURBS is used as a static library, a smart linker will
  // skip linking some class definitions that are needed for the
  // reading code to work right.

  const ON_Object* p=0;
  ON_Annotation::Cast(p);
  ON_Annotation2::Cast(p);
  ON_Bitmap::Cast(p);
  ON_Curve::Cast(p);
  ON_Geometry::Cast(p);
  ON_Object::Cast(p);
  ON_Surface::Cast(p);
  ON_UserData::Cast(p);
  ON_LinearDimension::Cast(p);
  ON_RadialDimension::Cast(p);
  ON_AngularDimension::Cast(p);
  ON_TextEntity::Cast(p);
  ON_Leader::Cast(p);
  ON_AnnotationTextDot::Cast(p);
  ON_AnnotationArrow::Cast(p);
  ON_LinearDimension2::Cast(p);
  ON_RadialDimension2::Cast(p);
  ON_AngularDimension2::Cast(p);
  ON_TextEntity2::Cast(p);
  ON_Leader2::Cast(p);
  ON_TextDot::Cast(p);
  ON_ArcCurve::Cast(p);
  // OBSOLETE // ON_OpenGLBitmap::Cast(p);
  ON_WindowsBitmap::Cast(p);
  ON_BrepVertex::Cast(p);
  ON_BrepEdge::Cast(p);
  ON_BrepTrim::Cast(p);
  ON_BrepLoop::Cast(p);
  ON_BrepFace::Cast(p);
  ON_Brep::Cast(p);
  ON_CurveOnSurface::Cast(p);
  ON_CurveProxy::Cast(p);
  ON_DimStyle::Cast(p);
  ON_Font::Cast(p);
  ON_Group::Cast(p);
  ON_Layer::Cast(p);
  ON_Light::Cast(p);
  ON_LineCurve::Cast(p);
  ON_Material::Cast(p);
  ON_Mesh::Cast(p);
  ON_NurbsCurve::Cast(p);
  ON_NurbsSurface::Cast(p);
  ON_PlaneSurface::Cast(p);
  ON_PointCloud::Cast(p);
  ON_Point::Cast(p);
  ON_PointGrid::Cast(p);
  ON_PolyCurve::Cast(p);
  ON_PolylineCurve::Cast(p);
  ON_RevSurface::Cast(p);
  ON_SumSurface::Cast(p);
  ON_SurfaceProxy::Cast(p);
  ON_UnknownUserData::Cast(p);
  ON_Viewport::Cast(p);
}


ON_ClassId* ON_ClassId::m_p0 = 0; // static pointer to first id in list
int ON_ClassId::m_mark0 = 0;

int ON_ClassId::IncrementMark()
{
  m_mark0++;
  return m_mark0;
}

int ON_ClassId::Purge( int mark_value )
{
  // Fundamental openNURBS class ids have a mark value of 0 and cannot be purged.
  int purge_count = 0;
  if ( mark_value > 0 ) {
    ON_ClassId* prev = 0;
    ON_ClassId* next = 0;
    ON_ClassId* p;
    for ( p = m_p0; p; p = next )
    {
      next = p->m_pNext;
      if ( p->m_mark == mark_value ) {
        purge_count++;
        if ( prev )
          prev->m_pNext = next;
        else
          m_p0 = next;
        p->m_pNext = 0;
      }
      else
        prev = p;
    }
  }
  return purge_count;
}

//////////////////////////////////////////////////////////////////////////////

static BOOL g_bDisableDemotion = FALSE;

static void IntToString( int i, char s[7] )
{
  // avoid format printing during early start up
  int j;
  int digit;
  char sdig[10];
  sdig[0] = '0';
  sdig[1] = '1';
  sdig[2] = '2';
  sdig[3] = '3';
  sdig[4] = '4';
  sdig[5] = '5';
  sdig[6] = '6';
  sdig[7] = '7';
  sdig[8] = '8';
  sdig[9] = '9';
  for ( digit = 5; digit > 0; digit-- )
  {
    j = i%10;
    if ( j > 9 || j < 0 )
    {
      s[digit] = '-';
    }
    else
    {
      s[digit] = sdig[j];
    }
    i /= 10;
  }
  s[0] = '-';
  s[6] = 0;
}

ON_ClassId::ON_ClassId( const char* sClassName, 
                        const char* sBaseClassName, 
                        ON_Object* (*create)(),
                        const char* sUUID // UUID in registry format from guidgen
                        ) 
                        : m_pNext(0),
                          m_pBaseClassId(0),
                          m_create(create),
                          m_mark(m_mark0)
{
  memset( m_sClassName, 0, sizeof(m_sClassName) );
  memset( m_sBaseClassName, 0, sizeof(m_sBaseClassName) );
  m_uuid = ON_UuidFromString(sUUID);
  if ( sClassName ) {
    strncpy( m_sClassName, sClassName, sizeof(m_sClassName)-1 );
  }
  if ( sBaseClassName ) {
    strncpy( m_sBaseClassName, sBaseClassName, sizeof(m_sBaseClassName)-1 );
  }
  m_pBaseClassId = ClassId( m_sBaseClassName );

  if ( !m_sClassName[0] ) {
    ON_ERROR("ON_ClassId::ON_ClassId() - missing class name");
    return;
  }

  const ON_ClassId* duplicate_class = ClassId( m_sClassName );
  // The m_mark0 > 2 test prevents opennurbs and Rhino from
  // having two ON_Object derived classes that have the same
  // name.  Plug-ins are free to use any name.
  if ( 0 != duplicate_class && m_mark0 > 2 )
  {
    char s[7];
    int ver;
    ON_WARNING("ON_ClassId::ON_ClassId() - class name already in use.  Will append number to make it unique.");
    for ( ver = 1; ver < 10000 && 0 != duplicate_class; ver++ )
    {
      IntToString(ver,s);
      s[6] = 0;
      strncpy( m_sClassName, sClassName, sizeof(m_sClassName)-1 );
      strncat( m_sClassName, s, sizeof(m_sClassName)-1 );
      duplicate_class = ClassId( m_sClassName );
    }
  }

  if ( 0 != duplicate_class )
  {
    // Do NOT permit core classes to have duplicate names.
    ON_ERROR("ON_ClassId::ON_ClassId() - class name already in use.");
    return;
  }

  if (    m_sClassName[0] != 'O'
       || m_sClassName[1] != 'N'
       || m_sClassName[2] != '_'
       || m_sClassName[3] != 'O'
       || m_sClassName[4] != 'b'
       || m_sClassName[5] != 'j'
       || m_sClassName[6] != 'e'
       || m_sClassName[7] != 'c'
       || m_sClassName[8] != 't'
       || m_sClassName[9] != 0 ) {
    if ( !m_sBaseClassName[0] ) 
    {
      ON_ERROR("ON_ClassId::ON_ClassId() - missing baseclass name.");
      return;
    }
  }

  g_bDisableDemotion = TRUE;
  if ( ClassId( m_uuid ) ) 
  {
    g_bDisableDemotion = FALSE;
    ON_ERROR("ON_ClassId::ON_ClassId() - class uuid already in use.");
    return;
  }
  g_bDisableDemotion = FALSE;

  if ( ON_UuidIsNil( m_uuid ) ) {
    ON_ERROR("ON_ClassId::ON_ClassId() - class uuid is nill.");
    return;
  }

  // see if any derived classes need to be updated because their static
  // members got initialized first
  if ( m_sClassName[0] ) 
  {
    this->m_pNext = m_p0;
    m_p0 = this;

    ON_ClassId* p = this->m_pNext;
    while ( p ) {
      if ( !p->m_pBaseClassId && p->m_sBaseClassName ) {
        if ( !strcmp( m_sClassName, p->m_sBaseClassName ) )
          p->m_pBaseClassId = this;
      }
      p = p->m_pNext;
    }
  }
}

ON_ClassId::~ON_ClassId()
{}

ON_Object* ON_ClassId::Create() const
{
  return m_create ? m_create() : 0;
}

const ON_ClassId* ON_ClassId::ClassId( const char* sClassName )
{
  // static member function
  // search list of class ids for one with a matching class name
  ON_ClassId* p;
  const char* s0;
  const char* s1;
  if ( !sClassName || !sClassName[0] || sClassName[0] == '0' )
    return NULL;
  for(p = m_p0; p; p = p->m_pNext) {
    // avoid strcmp() because it crashes on NULL strings
    s0 = sClassName;
    s1 = p->m_sClassName;
    if ( s0 && s1 && *s0 ) {
      while ( *s0 && *s0 == *s1 )
        {s0++; s1++;}
      if ( !(*s0) && !(*s1) )
        break;
    }
    else {
      break;
    }
  }
  return p;
}


const ON_ClassId* ON_ClassId::ClassId( ON_UUID uuid )
{
  // static member function
  // search list of class ids for one with a matching typecode
  const ON_ClassId* p;
  for(p = m_p0; p; p = p->m_pNext) 
  {
    if ( !ON_UuidCompare(&p->m_uuid,&uuid) )
      break;
  }

  if ( !p && !g_bDisableDemotion) 
  {
    // enable OpenNURBS toolkit to read files that contain old uuids even when
    // old class definitions are not loaded.

    // 5EAF1119-0B51-11d4-BFFE-0010830122F0 = TL_NurbsCurve
    ON_UUID nc0 = {0x5EAF1119,0x0B51,0x11d4,{0xBF,0xFE,  0x00,0x10,0x83,0x01,0x22,0xF0}};

    // 76A709D5-1550-11d4-8000-0010830122F0 = old nurbs curve
    ON_UUID nc1 = {0x76A709D5,0x1550,0x11d4,{0x80,0x00,  0x00,0x10,0x83,0x01,0x22,0xF0}};

    // 4760C817-0BE3-11d4-BFFE-0010830122F0 = TL_NurbsSurface
    ON_UUID ns0 = {0x4760C817,0x0BE3,0x11d4,{0xBF,0xFE,  0x00,0x10,0x83,0x01,0x22,0xF0}};

    // FA4FD4B5-1613-11d4-8000-0010830122F0 = old nurbs surface
    ON_UUID ns1 = {0xFA4FD4B5,0x1613,0x11d4,{0x80,0x00,  0x00,0x10,0x83,0x01,0x22,0xF0}};

    // EF638317-154B-11d4-8000-0010830122F0 = old poly curve
    ON_UUID pc0 = {0xEF638317,0x154B,0x11d4,{0x80,0x00,  0x00,0x10,0x83,0x01,0x22,0xF0}};

    // 0705FDEF-3E2A-11d4-800E-0010830122F0 = old trimmed surface
    ON_UUID br0 = {0x0705FDEF,0x3E2A,0x11d4,{0x80,0x0E,  0x00,0x10,0x83,0x01,0x22,0xF0}};

    // 2D4CFEDB-3E2A-11d4-800E-0010830122F0 = old b-rep
    ON_UUID br1 = {0x2D4CFEDB,0x3E2A,0x11d4,{0x80,0x0E,  0x00,0x10,0x83,0x01,0x22,0xF0}};

    // F06FC243-A32A-4608-9DD8-A7D2C4CE2A36 = TL_Brep
    ON_UUID br2 = {0xF06FC243,0xA32A,0x4608,{0x9D,0xD8, 0xA7,0xD2,0xC4,0xCE,0x2A,0x36}};

    // 0A8401B6-4D34-4b99-8615-1B4E723DC4E5 = TL_RevSurface
    ON_UUID revsrf = { 0xa8401b6, 0x4d34, 0x4b99, { 0x86, 0x15, 0x1b, 0x4e, 0x72, 0x3d, 0xc4, 0xe5 } };

    // 665F6331-2A66-4cce-81D0-B5EEBD9B5417 = TL_SumSurface
    ON_UUID sumsrf = { 0x665f6331, 0x2a66, 0x4cce, { 0x81, 0xd0, 0xb5, 0xee, 0xbd, 0x9b, 0x54, 0x17 } };

    if      ( !ON_UuidCompare( &uuid, &nc0 ) || !ON_UuidCompare( &uuid, &nc1 ) )
      p = &ON_NurbsCurve::m_ON_NurbsCurve_class_id;
    else if ( !ON_UuidCompare( &uuid, &ns0 ) || !ON_UuidCompare( &uuid, &ns1 ) )
      p = &ON_NurbsSurface::m_ON_NurbsSurface_class_id;
    else if ( !ON_UuidCompare( &uuid, &pc0 ) )
      p = &ON_PolyCurve::m_ON_PolyCurve_class_id;
    else if ( !ON_UuidCompare( &uuid, &br0 ) || !ON_UuidCompare( &uuid, &br1 ) || !ON_UuidCompare( &uuid, &br2 ) )
      p = &ON_Brep::m_ON_Brep_class_id;
    else if ( !ON_UuidCompare( &uuid, &revsrf ) )
      p = &ON_RevSurface::m_ON_RevSurface_class_id;
    else if ( !ON_UuidCompare( &uuid, &sumsrf ) )
      p = &ON_SumSurface::m_ON_SumSurface_class_id;
    else
      p = 0; // <- does nothing but it's a good place for debugger breakpoint
  }
  return p;
}

class ON__ClassIdDumpNode
{
public:
  ON__ClassIdDumpNode();
  ~ON__ClassIdDumpNode();
  const ON_ClassId* m_class_id;
  class ON__ClassIdDumpNode* m_parent_node;
  int m_depth;
  ON_SimpleArray<class ON__ClassIdDumpNode*> m_child_nodes;
  int CompareClassUuid( const class ON__ClassIdDumpNode& ) const;
  int CompareClassName( const class ON__ClassIdDumpNode& ) const;
  bool Dump( int depth, ON_TextLog& text_log );
};

ON__ClassIdDumpNode::ON__ClassIdDumpNode() 
{
  m_class_id=0;
  m_parent_node=0;
  m_depth=0;
};

ON__ClassIdDumpNode::~ON__ClassIdDumpNode() 
{
}

int ON__ClassIdDumpNode::CompareClassUuid( const class ON__ClassIdDumpNode& other ) const
{
  int rc = 0;
  const ON_ClassId* a = m_class_id;
  const ON_ClassId* b = other.m_class_id;
  if ( a != b )
  {
    if ( 0 == a )
    {
      rc = -1;
    }
    else if ( 0 == b )
      rc = 1;
    else
    {
      rc = ON_UuidCompare(a->Uuid(),b->Uuid());
      if ( 0 == rc )
      {
        rc = CompareClassName(other);
      }
    }
  }
  return rc;
}

int ON__ClassIdDumpNode::CompareClassName( const class ON__ClassIdDumpNode& other ) const
{
  int rc = 0;
  const ON_ClassId* a = m_class_id;
  const ON_ClassId* b = other.m_class_id;
  if ( a != b )
  {
    if ( 0 == a )
    {
      rc = -1;
    }
    else if ( 0 == b )
      rc = 1;
    else
    {
      const char* a_name = a->ClassName();
      const char* b_name = b->ClassName();
      if ( 0 == a_name )
      {
        if ( 0 == b_name )
        {
          rc = b->Mark() - a->Mark();
          if ( 0 == rc )
            rc = ON_UuidCompare(a->Uuid(),b->Uuid());
        }
        else
          rc = -1;
      }
      else if ( 0 == b_name )
      {
        rc = 1;
      }
      else
      {
        rc = on_stricmp(a_name,b_name);
        if ( 0 == rc )
        {
          rc = strcmp(a_name,b_name);
          if ( 0 == rc )
          {
            rc = b->Mark() - a->Mark();
            if ( 0 == rc )
              rc = ON_UuidCompare(a->Uuid(),b->Uuid());
          }
        }
      }
    }
  }
  return rc;
}


static int ON__ClassIdDumpNode_CompareUuid( const ON__ClassIdDumpNode* a, const ON__ClassIdDumpNode* b )
{
  int rc = 0;
  if ( 0 == a )
  {
    rc = (0 == b) ? 0 : -1;
  }
  else if ( 0 == b )
  {
    rc = 1;
  }
  else
  {
    rc = a->CompareClassUuid(*b);
  }
  return rc;
}

static int ON__ClassIdDumpNode_CompareName( ON__ClassIdDumpNode *const* a, ON__ClassIdDumpNode *const* b )
{
  int rc = 0;
  if ( 0 == a )
  {
    rc = (0 == b) ? 0 : -1;
  }
  else if ( 0 == b )
  {
    rc = 1;
  }
  else
  {
    rc = (*a)->CompareClassName(*(*b));
  }
  return rc;
}

bool ON__ClassIdDumpNode::Dump( int depth, ON_TextLog& text_log )
{
  bool rc = true;
  if ( 0 == m_class_id || m_depth != 0 || depth < 1)
    rc = false;
  else
  {
    m_depth = depth;
    const char* class_name = m_class_id->ClassName();
    if ( 0 == class_name )
    {
      class_name = "!!ERROR!!";
      rc = false;
    }
    text_log.Print("%s::ClassId: ",m_class_id->ClassName());
    text_log.Print( "mark=%d ",m_class_id->Mark() );
    text_log.Print( m_class_id->Uuid() );
    text_log.Print("  (%08x)\n",m_class_id);
    int i, count = m_child_nodes.Count();
    if ( count > 0 )
    {
      // dump children names alphabetically
      m_child_nodes.HeapSort( ON__ClassIdDumpNode_CompareName );

      text_log.PushIndent();
      for ( i = 0; i < count; i++ )
      {
        ON__ClassIdDumpNode* child_node = m_child_nodes[i];
        if ( 0 == child_node )
          rc = false;
        else
        {
          if ( !child_node->Dump(depth+1,text_log) )
            rc = false;
        }
      }
      text_log.PopIndent();
    }
  }
  return rc;
}

void ON_ClassId::Dump( ON_TextLog& dump )
{
  int i, j, count = 0;
  const ON_ClassId* p;
  for(p = m_p0; p && count < 1000000; p = p->m_pNext) 
  {
    count++;
  }
  if ( 0 != p )
  {
    dump.Print("ON_ClassId::m_p0 list is damaged.\n");
  }
  else
  {
    ON__ClassIdDumpNode tmp_node;
    ON_ClassArray<ON__ClassIdDumpNode> nodes(count);
    for(p = m_p0; p; p = p->m_pNext) 
    {
      ON__ClassIdDumpNode& node = nodes.AppendNew();
      node.m_class_id = p;
    }

    // sort nodes by class id's uuid
    nodes.HeapSort(ON__ClassIdDumpNode_CompareUuid);

    // fill in m_parent_node and m_child_nodes[]
    for ( i = 0; i < count; i++ )
    {
      ON__ClassIdDumpNode& node = nodes[i];
      p = node.m_class_id;
      if ( 0 != p )
      {
        tmp_node.m_class_id = p->BaseClass();
        j = nodes.BinarySearch(&tmp_node,ON__ClassIdDumpNode_CompareUuid);
        if ( j >= 0 && i != j)
        {
          ON__ClassIdDumpNode& base_node = nodes[j];
          node.m_parent_node = &base_node;
          base_node.m_child_nodes.Append(&node);
        }
      }      
    }

    // print class tree
	tmp_node.m_class_id = &ON_Object::m_ON_Object_class_id;
    i = nodes.BinarySearch(&tmp_node,ON__ClassIdDumpNode_CompareUuid);
    bool rc = false;
    if ( i >= 0 )
    {
      // recursively dump class tree
      rc = nodes[i].Dump(1,dump);
      for ( i = 0; i < count && rc; i++ )
      {
        if ( nodes[i].m_depth <= 0 )
          rc = false;
      }
    }

    if (!rc)
    {
      // should never happen
      for(p = m_p0; p; p = p->m_pNext) 
      {
        dump.Print("%s::ClassId: ",p->m_sClassName);
        dump.Print( "mark=%d ",p->m_mark );
        dump.Print( p->m_uuid );
        dump.Print("  (%08x)\n",p);
      }
    }
  }
}

const char* ON_ClassId::ClassName() const
{
  return m_sClassName;
}

const char* ON_ClassId::BaseClassName() const
{
  return m_sBaseClassName;
}

ON_UUID ON_ClassId::Uuid() const
{
  return m_uuid;
}

int ON_ClassId::Mark() const
{
  return m_mark;
}


const ON_ClassId* ON_ClassId::BaseClass() const
{
  return m_pBaseClassId;
}

BOOL ON_ClassId::IsDerivedFrom( const ON_ClassId* pBaseClassId ) const
{
  // determine if this is derived from pBaseClassId
  BOOL b = FALSE;
  if ( pBaseClassId ) {
    const ON_ClassId* p = this;
    for(;p;) {
      if ( p == pBaseClassId ) {
        b = TRUE;
        break;
      }
      p = p->m_pBaseClassId;
    }
  }
  return b;
}

//////////////////////////////////////////////////////////////////////////////

ON_VIRTUAL_OBJECT_IMPLEMENT(ON_Object,0,"60B5DBBD-E660-11d3-BFE4-0010830122F0");


ON_Object::ON_Object() : m_userdata_list(0)
{}

ON_Object::ON_Object(const ON_Object& src) : m_userdata_list(0)
{
  CopyUserData(src);
}

ON_Object& ON_Object::operator=(const ON_Object& src)
{
  if ( this !=&src ) {
    PurgeUserData();
    CopyUserData(src);
  }
  return *this;
}

ON_Object::~ON_Object()
{
  PurgeUserData();
}

void ON_Object::PurgeUserData()
{
  if ( m_userdata_list ) {
    ON_UserData* p = m_userdata_list;
    ON_UserData* next;
    while(p) {
      next = p->m_userdata_next;
      p->m_userdata_owner = 0;
      p->m_userdata_next = 0;
      delete p;
      p = next;
    }
    m_userdata_list = 0;
  }
}

BOOL ON_Object::AttachUserData( ON_UserData* p )
{
  BOOL rc = FALSE;
  if ( p 
       && NULL == p->m_userdata_owner
       && ON_UuidCompare( &ON_nil_uuid, &p->m_userdata_uuid) 
       && NULL == GetUserData( p->m_userdata_uuid )
       ) {
    if ( p->IsUnknownUserData() ) {
      // make sure we have valid user data - the first beta release of Rhino 2.0 
      // created empty user data.
      ON_UnknownUserData* uud = ON_UnknownUserData::Cast(p);
      if (uud)
        rc = uud->IsValid();
      if ( !rc ) {
        ON_ERROR("ON_Object::AttachUserData() - attempt to attach invalid UnknownUserData.");
      }
    }
    else
      rc = TRUE;
    if (rc) {
      p->m_userdata_owner = this;
      p->m_userdata_next = m_userdata_list;
      m_userdata_list = p;
    }
  }
  return rc;
}

BOOL ON_Object::DetachUserData( ON_UserData* p )
{
  BOOL rc = FALSE;
  if ( p && p->m_userdata_owner == this ) {
    ON_UserData* prev = 0;
    ON_UserData* ud = m_userdata_list;
    while ( ud ) {
      if ( ud == p ) {
        if ( prev )
          prev->m_userdata_next = ud->m_userdata_next;
        else
          m_userdata_list = ud->m_userdata_next;
        ud->m_userdata_owner = 0;
        ud->m_userdata_next = 0;
        rc = TRUE;
        break;
      }
      prev = ud;
      ud = ud->m_userdata_next;
    }
  }
  return rc;
}


ON_UserData* ON_Object::GetUserData( const ON_UUID& userdata_uuid ) const
{
  ON_UserData* prev = NULL;
  ON_UserData* p;
  for ( p = m_userdata_list; p; prev = p, p = p->m_userdata_next ) 
  {
    if ( !ON_UuidCompare( &p->m_userdata_uuid, &userdata_uuid ) ) 
    {
      if ( p->IsUnknownUserData() ) 
      {
        // See if we can convert this unknown user data into something useful.
        // Unknown user data is created when a 3dm archive is read and
        // the definition of the specific user data class is not loaded.
        // If something is getting around to asking for a specific kind
        // of user data, the class definition has probably be dynamically
        // loaded.
        ON_UnknownUserData* uud = ON_UnknownUserData::Cast(p);
        if ( uud ) {
          ON_UserData* realp = uud->Convert();
          if ( realp ) 
          {
            // replace unknown user data with the real thing
            if ( prev )
              prev->m_userdata_next = realp;
            else if ( p == m_userdata_list ) 
            {
              // little white lie to attach the "real" user
              // data to the object in place of the unknown
              // user data.
              ON_Object* pNotConst = const_cast<ON_Object*>(this);
              pNotConst->m_userdata_list = realp;
              realp->m_userdata_owner = pNotConst; // Dale Lear added 22 Jan 2004 to fix I/O bug 
            }
            realp->m_userdata_next = p->m_userdata_next;
            p->m_userdata_next = 0;
            p->m_userdata_owner = 0;
            delete p;
            p = realp;
          }
        }
      }
      break;
    }
  }
  return p; 
}

ON_UserData* ON_Object::FirstUserData() const
{
  return m_userdata_list;
}


void ON_Object::TransformUserData( const ON_Xform& x )
{
  ON_UserData *p, *next;
  for ( p = m_userdata_list; p; p = next ) {
    next = p->m_userdata_next;
    if ( !p->Transform(x) )
      delete p;
  }
}

void ON_Object::CopyUserData( const ON_Object& src )
{
  const ON_UserData* p;
  for ( p = src.m_userdata_list; p; p = p->m_userdata_next ) {
    if ( p->m_userdata_copycount ) {
      ON_Object* o = p->Duplicate();
      if ( o ) {
        if ( !AttachUserData(ON_UserData::Cast(o)) )
          delete o;
      }
    }
  }
}

void ON_Object::MoveUserData( ON_Object& src )
{
  ON_UserData *p, *next;

  if ( 0 == m_userdata_list )
  {
    // quick and simple when the "this" doesn't
    // have any user data.
    if ( 0 != src.m_userdata_list )
    {
      m_userdata_list = src.m_userdata_list;
      src.m_userdata_list = 0;
      for ( p = m_userdata_list; p; p = p->m_userdata_next )
      {
        p->m_userdata_owner = this;
      }
    }
  }
  else
  {    
    // Carefully move userdata an item at a time to
    // avoid conflicts with existing items on "this".

    // delete source user data that already exists on this
    for ( p = src.m_userdata_list; p; p = next ) {
      next = p->m_userdata_next;
      if ( GetUserData( p->m_userdata_uuid ) )
        delete p;
    }

    // append source user data to this user data
    next = src.m_userdata_list;
    src.m_userdata_list = 0;
    for ( p = next; p; p = p->m_userdata_next ) {
      p->m_userdata_owner = this;
    }

    if ( !m_userdata_list ) 
      m_userdata_list = next;
    else 
    {
      p = m_userdata_list;
      while ( p->m_userdata_next )
        p = p->m_userdata_next;
      p->m_userdata_next = next;
    }
  }
}


void ON_Object::MemoryRelocate()
{
  // When the memory location of an ON_Object changes,
  // the back-pointers on its user data need to be updated.
  ON_UserData* ud;
  for ( ud = m_userdata_list; ud; ud = ud->m_userdata_next )
  {
    ud->m_userdata_owner = this;
  }
}

BOOL ON_Object::IsKindOf( const ON_ClassId* pBaseClassId ) const
{
  BOOL b = FALSE;
  const ON_ClassId* p = ClassId();
  if ( p )
    b = p->IsDerivedFrom( pBaseClassId );
  return b;
}


ON::object_type ON_Object::ObjectType() const
{
  // virtual function that is generally overridden
  return ON::unknown_object_type;
}

ON_UUID ON_Object::ModelObjectId() const
{
  return ON_nil_uuid;
}

ON_UUID ON_Material::ModelObjectId() const
{
  return m_material_id;
}

ON_UUID ON_Layer::ModelObjectId() const
{
  return m_layer_id;
}

ON_UUID ON_Font::ModelObjectId() const
{
  return m_font_id;
}

ON_UUID ON_DimStyle::ModelObjectId() const
{
  return m_dimstyle_id;
}

ON_UUID ON_HatchPattern::ModelObjectId() const
{
  return m_hatchpattern_id;
}

ON_UUID ON_Linetype::ModelObjectId() const
{
  return m_linetype_id;
}

ON_UUID ON_Bitmap::ModelObjectId() const
{
  return m_bitmap_id;
}

ON_UUID ON_Light::ModelObjectId() const
{
  return m_light_id;
}

ON_UUID ON_TextureMapping::ModelObjectId() const
{
  return m_mapping_id;
}

ON_UUID ON_InstanceDefinition::ModelObjectId() const
{
  return m_uuid;
}

unsigned int ON_Object::SizeOf() const
{
  unsigned int sz = sizeof(*this);
  const ON_UserData* ud = m_userdata_list;
  while ( ud )
  {
    sz += ud->SizeOf();
    ud = ud->m_userdata_next;
  }
  return sz;
}

ON__UINT32 ON_Object::DataCRC(ON__UINT32 current_remainder) const
{
  // do not include user data.
  return current_remainder;
}

void ON_Object::Dump( ON_TextLog& dump ) const
{
  const ON_ClassId* p = ClassId();
  if ( p ) 
  {
    const char* class_name = p->ClassName();
    if ( 0 == class_name ) 
      class_name = "unknown";
    dump.Print("class name: %s\n",class_name);
    dump.Print("class uuid: ");
    dump.Print(p->Uuid());
    dump.Print("\n");
  }
  else 
  {
    dump.Print("ON_Object::ClassId() FAILED\n");
  }
}

BOOL ON_Object::Write(
       ON_BinaryArchive&
     ) const
{
  // default Write() does nothing.
  return FALSE;

  // object derived from ON_Object should have a Write() that looks
  // something like

  /*
  BOOL rc = file.Write3dmChunkVersion(1,0);
  if (rc) {
    // TODO
  }
  return rc;
  */

}

BOOL ON_Object::Read(
       ON_BinaryArchive&
     )
{
  // default Read() does nothing.
  return FALSE;

  // object derived from ON_Object should have a Read() that looks
  // something like

  /*
  int major_version = 0;
  int minor_version = 0;
  BOOL rc = file.Read3dmChunkVersion(&major_version,&minor_version);
  if (rc && major_version==1) {
    // common to all 1.x versions
    // TODO
  }
  return rc;
  */

}

void ON_Object::DestroyRuntimeCache( bool bDelete )
{
}

void ON_Curve::DestroyRuntimeCache( bool bDelete )
{
  if ( m_ctree ) 
  {
#if defined(OPENNURBS_PLUS_INC_)
    if ( bDelete )
      delete m_ctree;
#endif
    m_ctree = 0;
  }
}


void ON_CurveProxy::DestroyRuntimeCache( bool bDelete )
{
  if ( m_ctree ) 
  {
#if defined(OPENNURBS_PLUS_INC_)
    if ( bDelete )
      delete m_ctree;
#endif
    m_ctree = 0;
  }
  if ( 0 != m_real_curve && m_real_curve != this )
  {
    ON_Curve* curve = const_cast<ON_Curve*>(m_real_curve);
    if ( 0 != curve )
      curve->DestroyRuntimeCache( bDelete );
  }
}


void ON_Surface::DestroyRuntimeCache( bool bDelete )
{
  if ( m_stree ) 
  {
#if defined(OPENNURBS_PLUS_INC_)
    if ( bDelete )
      delete m_stree;
#endif
    m_stree = 0;
  }
}


void ON_SurfaceProxy::DestroyRuntimeCache( bool bDelete )
{
  if ( m_stree ) 
  {
#if defined(OPENNURBS_PLUS_INC_)
    if ( bDelete )
      delete m_stree;
#endif
    m_stree = 0;
  }
  if ( 0 != m_surface && m_surface != this )
  {
    ON_Surface* surface = const_cast<ON_Surface*>(m_surface);
    if ( 0 != surface )
      surface->DestroyRuntimeCache( bDelete );
  }
}

void ON_Brep::DestroyRuntimeCache( bool bDelete )
{
  int i, count;

  count = m_C2.Count();
  for ( i = 0; i < count; i++ )
  {
    if ( m_C2[i] )
      m_C2[i]->DestroyRuntimeCache(bDelete);
  }

  count = m_C3.Count();
  for ( i = 0; i < count; i++ )
  {
    if ( m_C3[i] )
      m_C3[i]->DestroyRuntimeCache(bDelete);
  }

  count = m_S.Count();
  for ( i = 0; i < count; i++ )
  {
    if ( m_S[i] )
      m_S[i]->DestroyRuntimeCache(bDelete);
  }

  count = m_T.Count();
  for ( i = 0; i < count; i++ )
  {
    m_T[i].DestroyRuntimeCache(bDelete);
  }

  count = m_E.Count();
  for ( i = 0; i < count; i++ )
  {
    m_E[i].DestroyRuntimeCache(bDelete);
  }

  count = m_F.Count();
  for ( i = 0; i < count; i++ )
  {
    m_F[i].DestroyRuntimeCache(bDelete);
  }

  // 15 August 2003 Dale Lear:
  //    I added the line to destroy the face's m_bbox.
  //    Since m_bbox is private, it will be recalculated
  //    when it is needed.  (We hope.)  The fact the face
  //    m_bbox is private and recalculated as needed makes
  //    it different than the m_pbox info on trims and loops.
  m_bbox.Destroy();
}

struct ON_Vtable* ON_ClassVtable(void* p)
{
  void* vtable_ptr = *((void**)p);
  struct ON_Vtable* the_vtable = (ON_Vtable*)vtable_ptr;
  return the_vtable;
}
