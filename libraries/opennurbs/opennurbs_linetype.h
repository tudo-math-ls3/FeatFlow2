#if !defined(OPENNURBS_LINETYPE_INC_)
#define OPENNURBS_LINETYPE_INC_






//////////////////////////////////////////////////////////////////////
// class ON_Linetype

class ON_CLASS ON_Linetype : public ON_Object
{
  ON_OBJECT_DECLARE(ON_Linetype);

public:

  /*
  Description:
    Sets index = -1.
  */
  ON_Linetype();

  ~ON_Linetype();


  /*
  Description:
    Sets index = -1 and emptys name and segment list.
  */
  void Default();

  /*
    Description:
      Tests that name is set and there is at least one non-zero length segment
  */
  BOOL IsValid( ON_TextLog* text_log = NULL ) const;

  void Dump( ON_TextLog& ) const; // for debugging

  /*
    Description:
      Write to file
  */
  BOOL Write(
         ON_BinaryArchive&  // serialize definition to binary archive
       ) const;

  /*
    Description:
      Read from file
  */
  BOOL Read(
         ON_BinaryArchive&  // restore definition from binary archive
       );

  // virtual
  ON_UUID ModelObjectId() const;


  //////////////////////////////////////////////////////////////////////
  //
  // Interface

  /*
    Unique name for each linetype
  */
  bool SetLinetypeName( const char*);
  bool SetLinetypeName( const wchar_t*);
	const wchar_t* LinetypeName() const;

  /*
    Index of each linetype
    This index is used by geometry objects to 
    reference a specific linetype
  */
  bool SetLinetypeIndex( int);
  int LinetypeIndex() const;

  /*
    Description:
      Returns the total length of one repeat of the pattern
  */
  double PatternLength() const;


  /*
    Description:
      Returns the number of segments in the pattern
  */
  int SegmentCount() const;

  /*
  Description:
    Adds a segment to the pattern
  Returns:
    Index of the added segment.
  */
  int AppendSegment( const ON_LinetypeSegment& segment);

  /*
  Description:
    Removes a segment in the linetype.
  Parameters:
    index - [in]
      Zero based index of the segment to remove.
  Returns:
    True if the segment index was removed.
  */
  bool RemoveSegment( int index );

  /*
    Description:
      Sets the segment at index to match segment
  */
  bool SetSegment( int index, const ON_LinetypeSegment& segment);

  /*
    Description:
      Sets the length and type of the segment at index
  */
  bool SetSegment( int index, double length, ON_LinetypeSegment::eSegType type);

  /*
    Description:
      Returns a copy of the segment at index
  */
  ON_LinetypeSegment Segment( int index) const;

  /*
    Description:
      Expert user function to get access to the segment array
      for rapid calculations.
  */
  ON_SimpleArray<ON_LinetypeSegment>& Segments();
  const ON_SimpleArray<ON_LinetypeSegment>& Segments() const;

public:
  int m_linetype_index;
  ON_UUID m_linetype_id;    // Set by Rhino - unique id of this linetype
  ON_wString m_linetype_name;

private:
  ON_SimpleArray<ON_LinetypeSegment> m_segments;
};

#endif

