/* $Header: /src4/opennurbs/opennurbs_uuid.h 8     9/06/06 11:31a Dalelear $ */
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

#if !defined(OPENNURBS_UUID_INC_)
#define OPENNURBS_UUID_INC_

// ON_UUID is a 16 byte universally unique identifier
#if defined(UUID_DEFINED)
typedef UUID ON_UUID;
#elif defined(GUID_DEFINED)
typedef GUID ON_UUID;
#else

#define ON_UUID_DECLARED_AS_CLASS
// For uuids, it is critical that the DataN fields have
// exactly the sizes specified below.  For that reason,
// the ON__UINTnn typedefs are used.
class ON_CLASS ON_UUID
{
public:
  ON__UINT32     Data1;    // 32 bit unsigned integer
  ON__UINT16     Data2;    // 16 bit unsigned integer
  ON__UINT16     Data3;    // 16 bit unsigned integer
  unsigned char  Data4[8]; 

  bool operator==(const ON_UUID& other) const;
  bool operator!=(const ON_UUID& other) const;
};

#endif

ON_BEGIN_EXTERNC

// All bits are zero in ON_nil_uuid and
// ON_UuidCompare( ON_nil_uuid, U ) < 0 if U != ON_nil_uuid.
extern ON_EXTERN_DECL const ON_UUID ON_nil_uuid;

// All bits are one in ON_max_uuid and
// ON_UuidCompare( U, ON_max_uuid ) < 0 if U != ON_max_uuid.
extern ON_EXTERN_DECL const ON_UUID ON_max_uuid;

ON_END_EXTERNC

#if defined(ON_CPLUSPLUS)

/*
Description:
  Creates a new uuid.(&a,&b) compares two uuids.
Parameters:
  new_uuid - [out]
Returns:
  True if successful.
Remarks:
  Only works on Windows.
*/
ON_DECL 
bool ON_CreateUuid( ON_UUID& uuid );

/*
Description:
  This class is used by ON_UuidIndexList.  It is used when
  uuids are used to search for items that can be found by
  an integer index.
*/
class ON_CLASS ON_UuidIndex
{
public:
  ON_UuidIndex();

  static 
  int CompareIdAndIndex( const ON_UuidIndex* a, const ON_UuidIndex* b );

  static 
  int CompareId( const ON_UuidIndex* a, const ON_UuidIndex* b );

  static 
  int CompareIndex( const ON_UuidIndex* a, const ON_UuidIndex* b );

  ON_UUID m_uuid;
  int m_i;
};

/*
Description:
  ON_UuidCompare(&a,&b) compares two uuids.
Parameters:
  a - [in]
  b - [in]
Returns:
  @untitled table
  -1    a < b
   0    a == b
  +1    a > b
Remarks:
  A NULL pointer is considered < a non-NULL pointer.
*/
ON_DECL 
int ON_UuidCompare( 
        const ON_UUID* a, 
        const ON_UUID* b 
        );

/*
Description:
  ON_UuidCompare(a,b) compares two uuids.
Parameters:
  a - [in]
  b - [in]
Returns:
  @untitled table
  -1    a < b
   0    a == b
  +1    a > b
*/
ON_DECL 
int ON_UuidCompare( 
        const ON_UUID& a, 
        const ON_UUID& b
        );

/*
Description:
  Test uuid to see if it is nil (identically zero).
Parameters:
  uuid - [in]
Returns:
  true if uuid is nil.
*/
ON_DECL
bool ON_UuidIsNil( 
        const ON_UUID& uuid 
        );

/*
Description:
  Test uuid to see if it is not nil (not identically zero).
Parameters:
  uuid - [in]
Returns:
  true if uuid is not nil (non zero)
*/
ON_DECL
bool ON_UuidIsNotNil( 
        const ON_UUID& uuid 
        );

/*
Description:
  Converts a string like
    "{85A08515-f383-11d3-BFE7-0010830122F0}" 
  into a uuid.
  The brackets are optional and are ignored.
  Hyphens can appear anywhere or be missing.
  The hex digits can be upper or lower case.
Parameters:
  s - [in]
Returns:
  uuid.  
  If the string is not a uuid, then ON_nil_uuid is returnd.
*/
ON_DECL 
ON_UUID ON_UuidFromString( const char* s );

/*
Description:
  Converts a string like
    "{85A08515-f383-11d3-BFE7-0010830122F0}" 
  into a uuid.
  The brackets are optional and are ignored.
  Hyphens can appear anywhere or be missing.
  The hex digits can be upper or lower case.
Parameters:
  s - [in]
Returns:
  uuid.  
  If the string is not a uuid, then ON_nil_uuid is returnd.
*/
ON_DECL 
ON_UUID ON_UuidFromString( const wchar_t* s );

/*
Description:
  Converts a uuid to a null termintated ASCII string like 
     "85a08515-f383-11d3-bfe7-0010830122f0". 
Parameters:
  uuid - [in]
  s - [out]  The s[] char array must have length >= 37.  
             The returned char array will have a 36 
             character uuid in s[0..35] and a null in s[36].
Returns:
  The pointer to the array is returned.
*/
ON_DECL 
char* ON_UuidToString( const ON_UUID& uuid, char* s );


/*
Description:
  Converts a uuid to a null termintated UNICODE string like 
     "85a08515-f383-11d3-bfe7-0010830122f0". 
Parameters:
  uuid - [in]
  s - [out]  The s[] wchar_t array must have length >= 37.  
             The returned char array will have a 36 
             character uuid in s[0..35] and a null in s[36].
Returns:
  The pointer to the array is returned.
*/
ON_DECL 
wchar_t* ON_UuidToString( const ON_UUID& uuid, wchar_t* s );

class ON_String;

/*
Description:
  Converts a uuid to a null termintated string like 
     "85a08515-f383-11d3-bfe7-0010830122f0". 
Parameters:
  uuid - [in]
  s - [out]
Returns:
  The pointer to the array is returned.
*/
ON_DECL 
const char* ON_UuidToString( const ON_UUID& uuid, ON_String& s);

class ON_wString;

/*
Description:
  Converts a uuid to a null termintated string like 
     "85a08515-f383-11d3-bfe7-0010830122f0". 
Parameters:
  uuid - [in]
  s - [out]
Returns:
  The pointer to the array is returned.
*/
ON_DECL 
const wchar_t* ON_UuidToString( const ON_UUID& uuid, ON_wString& s);

#endif

#endif
