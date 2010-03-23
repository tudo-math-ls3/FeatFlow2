/* $Header: /src4/opennurbs/opennurbs_dll.cpp 5     8/30/05 2:11p Dalelear $ */
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

// opennurbs_dll.cpp : Defines the entry point for the Windows DLL application.
//
#include "opennurbs.h"

#if defined(ON_OS_WINDOWS) && defined(ON_DLL_EXPORTS)

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
  static BOOL bRunning = 0;
  if ( !bRunning ) {
    ON_OpenNURBS_Begin();
    bRunning = TRUE;
  }

  switch( ul_reason_for_call ) {

  case DLL_PROCESS_ATTACH:
    ::OutputDebugStringA("OpenNURBS DllMain() ul_reason_for_call = DLL_PROCESS_ATTACH\n");
    ON_ClassId::IncrementMark(); // make sure each DLL that each process that 
                                 // uses OpenNURBS has a unique mark.
    break;

  case DLL_THREAD_ATTACH:
    ::OutputDebugStringA("OpenNURBS DllMain() ul_reason_for_call = DLL_THREAD_ATTACH\n");
    break;

  case DLL_THREAD_DETACH:
    ::OutputDebugStringA("OpenNURBS DllMain() ul_reason_for_call = DLL_THREAD_DETACH\n");
    break;

  case DLL_PROCESS_DETACH:
    ::OutputDebugStringA("OpenNURBS DllMain() ul_reason_for_call = DLL_PROCESS_DETACH\n");
    break;

  default:
    ::OutputDebugStringA("OpenNURBS DllMain() ul_reason_for_call = ?\n");
    break;
  }

  return TRUE;
}


///////////////////////////////////////////////////////////////////////////////
//
// For testing crash handling in opennurbs.dll
//

typedef int (*CRASHTEST__FUNCTION__POINTER__)(int);

static
void CrashTestHelper_GetNullIntPrt(int ** pp )
{
  memset(pp,0,sizeof(*pp));
}

static
void CrashTestHelper_GetBogusIntPtr(int ** pp )
{
  // caca[] is 16 bytes to insure it at least as large as any
  // pointer that will appear in my lifetime.
  unsigned int caca[4] = {0xCACAF00D,0xCACAF00D,0xCACAF00D,0xCACAF00D};
  memcpy(pp,&caca,sizeof(*pp));
}

static
void CrashTestHelper_GetNullFuncPtr(CRASHTEST__FUNCTION__POINTER__* pp )
{
  memset(pp,0,sizeof(*pp));
}

static
void CrashTestHelper_GetBogusFuncPtr(CRASHTEST__FUNCTION__POINTER__* pp )
{
  // caca[] is 16 bytes to insure it at least as large as any
  // pointer that will appear in my lifetime.
  unsigned int caca[4] = {0xCACAF00D,0xCACAF00D,0xCACAF00D,0xCACAF00D};
  memcpy(pp,&caca[0],sizeof(*pp));
}

static
int CrashTestHelper_DerefNullIntPtr( int crash_type, int* stack_ptr )
{
  int* ptr;
  CrashTestHelper_GetNullIntPrt(&ptr); // sets ptr = NULL
  *stack_ptr = *ptr;    // dereferences NULL pointer and crashes
  return crash_type;
}

static
int CrashTestHelper_DerefBogusIntPtr( int crash_type, int* stack_ptr )
{
  int* ptr;
  CrashTestHelper_GetBogusIntPtr(&ptr); // sets ptr = 0xCACAF00D
  *stack_ptr = *ptr;     // dereferences bogus pointer and crashes
  return crash_type;
}

static
int CrashTestHelper_CallNullFuncPtr( int crash_type, int* stack_ptr )
{
  CRASHTEST__FUNCTION__POINTER__ fptr;
  CrashTestHelper_GetNullFuncPtr(&fptr); // sets ptr = NULL
  *stack_ptr = fptr(crash_type); // dereferences NULL function pointer and crashes
  return crash_type;
}

static
int CrashTestHelper_CallBoguslFuncPtr( int crash_type, int* stack_ptr )
{
  CRASHTEST__FUNCTION__POINTER__ fptr;
  CrashTestHelper_GetNullFuncPtr(&fptr); // sets ptr = NULL
  *stack_ptr = fptr(crash_type); // dereferences bogus function pointer and crashes
  return crash_type;
}

static
bool CrashTestHelper_DivideByZero( const char* zero )
{
  int iz = 0, iy = 0;
  double z;
  sscanf( zero, "%g", &z );
  double y = 13.0/z;
  if ( ON_IsValid(y) )
  {
    iz = (int)z;
    iy = 17/iz;
  }
  return (123456 != iz);
}

static
bool CrashTestHelper_LogNegativeNumber( const char* minus_one )
{
  double z;
  sscanf( minus_one, "%g", &z );
  double y = log(z);
  return ON_IsValid(y);
}

/*
Description:
  Create a condition that should result in a crash.
  The intended use is for testing application crash handling.
Parameters:
  crash_type - [in]
    0: dereference NULL data pointer
    1: dereference bogus data pointer (0xCACAF00D)
    2: dereference NULL function pointer
    3: dereference bogus function pointer (0xCACAF00D)
    4: call abort()
    5: call exit(99);
    6: divide by zero
    7: log(negative number) - should call math error handler
Returns:
  Value of crash_type parameter.
*/
int ON_CrashTest( int crash_type, ON_TextLog& text_log )
{
  // Note: The code and calls are intentionally obtuse
  //       so that it is difficult for an optimizer to
  //       not perform the calculation.
  int stack_int = 0;
  int rc = crash_type;

  switch( crash_type )
  {
  case 0: // dereference NULL pointer
    rc = CrashTestHelper_DerefNullIntPtr( crash_type, &stack_int );
    break;

  case 1: // dereference bogus pointer
    rc = CrashTestHelper_DerefBogusIntPtr( crash_type, &stack_int );
    break;

  case 2: // corrupt stack so return sets IP to bogus value
    rc = CrashTestHelper_CallNullFuncPtr( crash_type, &stack_int );
    break;

  case 3: // divide by zero
    rc = CrashTestHelper_CallBoguslFuncPtr( crash_type, &stack_int );
    break;

  case 4:
    abort();
    text_log.Print("abort() didn't crash.\n");
    break;

  case 5:
    exit(99);
    text_log.Print("exit(99) didn't crash.\n");
    break;

  case 6:
    if ( CrashTestHelper_DivideByZero( "0.0" ) )
    {
      text_log.Print("Divide by 0.0 didn't crash - exception must have been handled or ignored.\n");
    }
    break;

  case 7:
    if ( CrashTestHelper_LogNegativeNumber( "-1.0" ) )
    {
      text_log.Print("log(negative number) didn't crash - exception must have been handled or ignored.\n");
    }
    break;

  default:
    break;
  }

  return rc;
}

#endif
