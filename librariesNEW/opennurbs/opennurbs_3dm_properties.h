/* $Header: /src3/opennurbs/opennurbs_3dm_properties.h 1     7/25/01 2:14p Dalelear $ */
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

#if !defined(OPENNURBS_3DM_PROPERTIES_INC_)
#define OPENNURBS_3DM_PROPERTIES_INC_

//////////////////////////////////////////////////////////////////////////////////////////

class ON_CLASS ON_3dmRevisionHistory
{
public:
  ON_3dmRevisionHistory();
  ON_3dmRevisionHistory( const ON_3dmRevisionHistory& );
  ~ON_3dmRevisionHistory();
  ON_3dmRevisionHistory& operator=( const ON_3dmRevisionHistory& );

  void Default();
  BOOL IsValid() const;
  int NewRevision(); // returns updated revision count

  BOOL Read( ON_BinaryArchive& );
  BOOL Write( ON_BinaryArchive& ) const;

  void Dump( ON_TextLog& ) const;

  ON_wString m_sCreatedBy;
  ON_wString m_sLastEditedBy;
  struct tm m_create_time;
  struct tm m_last_edit_time;
  int       m_revision_count;
};

//////////////////////////////////////////////////////////////////////////////////////////

class ON_CLASS ON_3dmNotes
{
public:
  ON_3dmNotes();
  ON_3dmNotes( const ON_3dmNotes& );
  ~ON_3dmNotes();
  ON_3dmNotes& operator=(const ON_3dmNotes&);

  void Default();
  BOOL IsValid() const;

  BOOL Read( ON_BinaryArchive& );
  BOOL Write( ON_BinaryArchive& ) const;

  void Dump(ON_TextLog&) const;

  ////////////////////////////////////////////////////////////////
  //
  // Interface - this information is serialized.  Applications
  // may want to derive a runtime class that has additional
  // window and font information.
  ON_wString m_notes; // UNICODE
  BOOL m_bVisible;    // TRUE if notes window is showing
  BOOL m_bHTML;       // TRUE if notes are in HTML

  // last window position
  int m_window_left;
  int m_window_top;
  int m_window_right;
  int m_window_bottom;
};

//////////////////////////////////////////////////////////////////////////////////////////

class ON_CLASS ON_3dmApplication
{
  // application that created the 3dm file
public:
  ON_3dmApplication();
  ON_3dmApplication( const ON_3dmApplication& );
  ~ON_3dmApplication();
  ON_3dmApplication& operator=(const ON_3dmApplication&);

  void Default();
  BOOL IsValid() const;

  BOOL Read( ON_BinaryArchive& );
  BOOL Write( ON_BinaryArchive& ) const;

  void Dump( ON_TextLog& ) const;

  ON_wString m_application_name;    // short name like "Rhino 2.0"
  ON_wString m_application_URL;     // URL
  ON_wString m_application_details; // whatever you want
};

//////////////////////////////////////////////////////////////////////////////////////////

class ON_CLASS ON_3dmProperties
{
public:
  ON_3dmProperties();
  ~ON_3dmProperties();
  ON_3dmProperties(const ON_3dmProperties&);
  ON_3dmProperties& operator=(const ON_3dmProperties&);

  void Default();

  BOOL Read(ON_BinaryArchive&);
  BOOL Write(ON_BinaryArchive&) const;

  void Dump( ON_TextLog& ) const;

  ON_3dmRevisionHistory  m_RevisionHistory;
  ON_3dmNotes            m_Notes;
  ON_WindowsBitmap       m_PreviewImage;     // preview image of model
  ON_3dmApplication      m_Application;      // application that created 3DM file
};

//////////////////////////////////////////////////////////////////////////////////////////

#endif
