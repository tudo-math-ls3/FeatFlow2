/* $Header: /src4/opennurbs/opennurbs_version.h 41    9/07/06 9:38a Dalelear $ */
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

#if !defined(OPENNURBS_VERSION_DEFINITION)
#error Do NOT include opennurbs_version.h in your code.  Use ON::Version() instead.
#endif

// OpenNURBS core developers:
//   If you change OpenNURBS definitions, set 
//   OPENNURBS_VERSION = YYYYMMDDn, where "n" starts at 0.

// OpenNURBS users:
//   Do not change OPENNURBS_VERSION or the OpenNURBS code
//   that reads 3DM files not work correctly.

// version 200210210: 
//   First V3 opennurbs source code release

// version 200211190: 
//   First Rhino 3.0 SDK opennurbs source code release

// version 200303050: 
//   Rhino 3.0 saves layers layer dialog sorted order

// version 200306060: 
//   Rhino 3.0 SR2

// version 200310100: 
//   Rhino 3.0 SR3

// version 200310290: 
//   Opennurbs source code release corresponding to Rhino 3.0 SR3.

// version 200311100
//   Beginning of Rhino 4.0

// version 200401070
//   Added user data support to ON_3dmObjectAttributes

// version 200404080
//   Added ON_Hatch and related classes

// version 200405030
//   Added Hatch archiving

// version 200405180
//   OpenNURBS.org source code update released

// version 200410010
//   ON_COMPONENT_INDEX added.

// version 200410190
//   Rhino V4 WIP began writing V4 files.

// version 200501310
//   Added ON_Viewport::m_viewport_id and 
//   ON_3dmObjectAttributes::m_dmref[] to I/O support
//   to 3dm binary archives.

// version 200503150
//   files incorrectly written

// version 200503170
//   ON_Linetype table

// version 200503250
//   ON_PlugInRef list added to ON_3dmSettings

// version 200505020
//   ON_Material expansion 1

// version 200505030
//   ON_Material expansion 2

// version 200505110
//   Added ON_TextureMapping support

// version 200508110
//   Annotation overhaul

// version 200508300
//   MBCS to UNICODE

// version 200509020
//   ON_3dmObjectAttrbutes IO version 1.6

// version 200509080
//   Remove prototype ON_LayerSet that was never used

// version 200509090
//   Remove obsolete ON_MeshParameters m_bWeld and m_combine_angle

// version 200510040
//   Changes so the code works with Win32 and Win64

// version 200510140
//   New ON_MeshParameters settings

// version 200510140
//   Added ON_Linetype::m_linetype_id
//         ON_Group::m_group_id
//         ON_Font::m_font_id
//         ON_Dimstyle::m_dimstyle_id
//         ON_HatchPattern::m_hatchpattern_id

// version 200511010
//    Fixed hatchpattern record to use ReadObject/WriteObject

// version 200511110
//    Added texture mapping table to 3dm archive

// version 200511150
//    Added IO support for V4 material definition to 3dm archive

// version 200512070
//    Added IO support for point cloud normals and colors

// version 200601180
//    Added history record table

// version 200602080
//    Added ON_3dmSettings.m_PageUnitsAndTolerances

// version 200603070
//    Added 
//      ON_3dmSettings.m_model_basepoint, 
//      ON_InstanceDefinition.m_source_unit_system
//      ON_InstanceDefinition.m_source_bRelativePath

// version 200603100
//    ON_CheckSum change
//

// version 200604190
//    Added ON_3dmSettings.m_bSaveMaterialBitmapsInFile
//    and started saving more bitmaps in files.
//

// version 200605260
//    Changes to texture mapping
//

// version 200606060
//    Added IO support for m_rendering_attributes on
//    layers and object attributes.
//

// version 200607130
//    Added ON_ObjRefEvaluationParameter 
//    and enhanced ON_ObjRef_IRefID
//

// version 200609070
//    First openNURBS V4 release
//

#define OPENNURBS_VERSION 200609070
