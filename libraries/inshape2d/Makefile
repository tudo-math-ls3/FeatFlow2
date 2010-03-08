###############################################################
# 
#
# In order to get full wide character (UNICODE) support, you
# need to define _GNU_SOURCE

AR = ar qvl 
RANLIB = ranlib

# C compiler and flags
CC = gcc
CFLAGS = -g -Wall -c -D_GNU_SOURCE -DMY_ZCALLOC -DZ_PREFIX -I.  
#CFLAGS = -O -Wall -D_GNU_SOURCE -DMY_ZCALLOC -DZ_PREFIX -I.  

# C++ compiler and flags
CCC = g++
CCFLAGS = -g -w -c -D_GNU_SOURCE -I. -I./src -I./include -I../opennurbs
#CCFLAGS = -O -Wall -D_GNU_SOURCE -I.
LDFLAGS = -L. 

###############################################################
LINK = $(CCC)
LINKFLAGS =

###############################################################

#%.o: %.cpp
#	$(CCC) $(CCFLAGS) -c $< -o $@

.cpp.o :
	$(CCC) $(CCFLAGS) -c $*.cpp -o $*.o

# If your make doesn't have a rule for .c to .o, uncomment the
# next two lines.
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@


###############################################################
# Name of the library
OPENNURBS_LIB_NAME = inshape2d
OPENNURBS_LIB_FILE = lib$(OPENNURBS_LIB_NAME).a


ON_INC= ../opennurbs/opennurbs.h \
		../opennurbs/opennurbs_3dm.h \
		../opennurbs/opennurbs_3dm_attributes.h \
		../opennurbs/opennurbs_3dm_properties.h \
		../opennurbs/opennurbs_3dm_settings.h \
		../opennurbs/opennurbs_annotation.h \
		../opennurbs/opennurbs_annotation2.h \
		../opennurbs/opennurbs_arc.h \
		../opennurbs/opennurbs_arccurve.h \
		../opennurbs/opennurbs_archive.h \
		../opennurbs/opennurbs_array.h \
		../opennurbs/opennurbs_array_defs.h \
		../opennurbs/opennurbs_bezier.h \
		../opennurbs/opennurbs_bitmap.h \
		../opennurbs/opennurbs_bounding_box.h \
		../opennurbs/opennurbs_brep.h \
		../opennurbs/opennurbs_circle.h \
		../opennurbs/opennurbs_color.h \
		../opennurbs/opennurbs_cone.h \
		../opennurbs/opennurbs_crc.h \
		../opennurbs/opennurbs_curve.h \
		../opennurbs/opennurbs_curveonsurface.h \
		../opennurbs/opennurbs_curveproxy.h \
		../opennurbs/opennurbs_cylinder.h \
		../opennurbs/opennurbs_defines.h \
		../opennurbs/opennurbs_detail.h \
		../opennurbs/opennurbs_dimstyle.h \
		../opennurbs/opennurbs_ellipse.h \
		../opennurbs/opennurbs_error.h \
		../opennurbs/opennurbs_evaluate_nurbs.h \
		../opennurbs/opennurbs_extensions.h \
		../opennurbs/opennurbs_font.h \
		../opennurbs/opennurbs_fpoint.h \
		../opennurbs/opennurbs_geometry.h \
		../opennurbs/opennurbs_gl.h \
		../opennurbs/opennurbs_group.h \
		../opennurbs/opennurbs_hatch.h \
		../opennurbs/opennurbs_instance.h \
		../opennurbs/opennurbs_intersect.h \
		../opennurbs/opennurbs_knot.h \
		../opennurbs/opennurbs_layer.h \
		../opennurbs/opennurbs_light.h \
		../opennurbs/opennurbs_line.h \
		../opennurbs/opennurbs_linecurve.h \
		../opennurbs/opennurbs_linestyle.h \
		../opennurbs/opennurbs_linetype.h \
		../opennurbs/opennurbs_mapchan.h \
		../opennurbs/opennurbs_massprop.h \
		../opennurbs/opennurbs_material.h \
		../opennurbs/opennurbs_math.h \
		../opennurbs/opennurbs_matrix.h \
		../opennurbs/opennurbs_memory.h \
		../opennurbs/opennurbs_mesh.h \
		../opennurbs/opennurbs_nurbscurve.h \
		../opennurbs/opennurbs_nurbssurface.h \
		../opennurbs/opennurbs_object.h \
		../opennurbs/opennurbs_object_history.h \
		../opennurbs/opennurbs_objref.h \
		../opennurbs/opennurbs_offsetsurface.h \
		../opennurbs/opennurbs_optimize.h \
		../opennurbs/opennurbs_plane.h \
		../opennurbs/opennurbs_planesurface.h \
		../opennurbs/opennurbs_pluginlist.h \
		../opennurbs/opennurbs_point.h \
		../opennurbs/opennurbs_pointcloud.h \
		../opennurbs/opennurbs_pointgeometry.h \
		../opennurbs/opennurbs_pointgrid.h \
		../opennurbs/opennurbs_polycurve.h \
		../opennurbs/opennurbs_polyline.h \
		../opennurbs/opennurbs_polylinecurve.h \
		../opennurbs/opennurbs_rendering.h \
		../opennurbs/opennurbs_revsurface.h \
		../opennurbs/opennurbs_sphere.h \
		../opennurbs/opennurbs_string.h \
		../opennurbs/opennurbs_sumsurface.h \
		../opennurbs/opennurbs_surface.h \
		../opennurbs/opennurbs_surfaceproxy.h \
		../opennurbs/opennurbs_system.h \
		../opennurbs/opennurbs_textlog.h \
		../opennurbs/opennurbs_texture.h \
		../opennurbs/opennurbs_texture_mapping.h \
		../opennurbs/opennurbs_torus.h \
		../opennurbs/opennurbs_userdata.h \
		../opennurbs/opennurbs_uuid.h \
		../opennurbs/opennurbs_version.h \
		../opennurbs/opennurbs_viewport.h \
		../opennurbs/opennurbs_workspace.h \
		../opennurbs/opennurbs_x.h \
		../opennurbs/opennurbs_xform.h \
		../opennurbs/opennurbs_zlib.h

ON_SRC= ../opennurbs/opennurbs_3dm_attributes.cpp \
		../opennurbs/opennurbs_3dm_properties.cpp \
		../opennurbs/opennurbs_3dm_settings.cpp \
		../opennurbs/opennurbs_annotation.cpp \
		../opennurbs/opennurbs_annotation2.cpp \
		../opennurbs/opennurbs_arc.cpp \
		../opennurbs/opennurbs_arccurve.cpp \
		../opennurbs/opennurbs_archive.cpp \
		../opennurbs/opennurbs_array.cpp \
		../opennurbs/opennurbs_basic.cpp \
		../opennurbs/opennurbs_bezier.cpp \
		../opennurbs/opennurbs_beziervolume.cpp \
		../opennurbs/opennurbs_bitmap.cpp \
		../opennurbs/opennurbs_bounding_box.cpp \
		../opennurbs/opennurbs_brep.cpp \
		../opennurbs/opennurbs_brep_changesrf.cpp \
		../opennurbs/opennurbs_brep_extrude.cpp \
		../opennurbs/opennurbs_brep_io.cpp \
		../opennurbs/opennurbs_brep_isvalid.cpp \
		../opennurbs/opennurbs_brep_kinky.cpp \
		../opennurbs/opennurbs_brep_tools.cpp \
		../opennurbs/opennurbs_brep_v2valid.cpp \
		../opennurbs/opennurbs_circle.cpp \
		../opennurbs/opennurbs_color.cpp \
		../opennurbs/opennurbs_cone.cpp \
		../opennurbs/opennurbs_crc.cpp \
		../opennurbs/opennurbs_curve.cpp \
		../opennurbs/opennurbs_curveonsurface.cpp \
		../opennurbs/opennurbs_curveproxy.cpp \
		../opennurbs/opennurbs_cylinder.cpp \
		../opennurbs/opennurbs_defines.cpp \
		../opennurbs/opennurbs_detail.cpp \
		../opennurbs/opennurbs_dimstyle.cpp \
		../opennurbs/opennurbs_ellipse.cpp \
		../opennurbs/opennurbs_error.cpp \
		../opennurbs/opennurbs_error_message.cpp \
		../opennurbs/opennurbs_evaluate_nurbs.cpp \
		../opennurbs/opennurbs_extensions.cpp \
		../opennurbs/opennurbs_font.cpp \
		../opennurbs/opennurbs_geometry.cpp \
		../opennurbs/opennurbs_gl.cpp \
		../opennurbs/opennurbs_group.cpp \
		../opennurbs/opennurbs_hatch.cpp \
		../opennurbs/opennurbs_instance.cpp \
		../opennurbs/opennurbs_intersect.cpp \
		../opennurbs/opennurbs_knot.cpp \
		../opennurbs/opennurbs_layer.cpp \
		../opennurbs/opennurbs_light.cpp \
		../opennurbs/opennurbs_line.cpp \
		../opennurbs/opennurbs_linecurve.cpp \
		../opennurbs/opennurbs_linetype.cpp \
		../opennurbs/opennurbs_massprop.cpp \
		../opennurbs/opennurbs_material.cpp \
		../opennurbs/opennurbs_math.cpp \
		../opennurbs/opennurbs_matrix.cpp \
		../opennurbs/opennurbs_memory.c \
		../opennurbs/opennurbs_memory_util.c \
		../opennurbs/opennurbs_mesh.cpp \
		../opennurbs/opennurbs_mesh_tools.cpp \
		../opennurbs/opennurbs_morph.cpp \
		../opennurbs/opennurbs_nurbscurve.cpp \
		../opennurbs/opennurbs_nurbssurface.cpp \
		../opennurbs/opennurbs_nurbsvolume.cpp \
		../opennurbs/opennurbs_object.cpp \
		../opennurbs/opennurbs_object_history.cpp \
		../opennurbs/opennurbs_objref.cpp \
		../opennurbs/opennurbs_offsetsurface.cpp \
		../opennurbs/opennurbs_optimize.cpp \
		../opennurbs/opennurbs_plane.cpp \
		../opennurbs/opennurbs_planesurface.cpp \
		../opennurbs/opennurbs_pluginlist.cpp \
		../opennurbs/opennurbs_point.cpp \
		../opennurbs/opennurbs_pointcloud.cpp \
		../opennurbs/opennurbs_pointgeometry.cpp \
		../opennurbs/opennurbs_pointgrid.cpp \
		../opennurbs/opennurbs_polycurve.cpp \
		../opennurbs/opennurbs_polyline.cpp \
		../opennurbs/opennurbs_polylinecurve.cpp \
		../opennurbs/opennurbs_revsurface.cpp \
		../opennurbs/opennurbs_sphere.cpp \
		../opennurbs/opennurbs_string.cpp \
		../opennurbs/opennurbs_sum.cpp \
		../opennurbs/opennurbs_sumsurface.cpp \
		../opennurbs/opennurbs_surface.cpp \
		../opennurbs/opennurbs_surfaceproxy.cpp \
		../opennurbs/opennurbs_textlog.cpp \
		../opennurbs/opennurbs_torus.cpp \
		../opennurbs/opennurbs_userdata.cpp \
		../opennurbs/opennurbs_uuid.cpp \
		../opennurbs/opennurbs_viewport.cpp \
		../opennurbs/opennurbs_workspace.cpp \
		../opennurbs/opennurbs_wstring.cpp \
		../opennurbs/opennurbs_x.cpp \
		../opennurbs/opennurbs_xform.cpp \
		../opennurbs/opennurbs_zlib.cpp \
		../opennurbs/opennurbs_zlib_memory.cpp

ON_OBJ=opennurbs_3dm_attributes.o \
		opennurbs_3dm_properties.o \
		opennurbs_3dm_settings.o \
		opennurbs_annotation.o \
		opennurbs_annotation2.o \
		opennurbs_arc.o \
		opennurbs_arccurve.o \
		opennurbs_archive.o \
		opennurbs_array.o \
		opennurbs_basic.o \
		opennurbs_bezier.o \
		opennurbs_beziervolume.o \
		opennurbs_bitmap.o \
		opennurbs_bounding_box.o \
		opennurbs_brep.o \
		opennurbs_brep_changesrf.o \
		opennurbs_brep_extrude.o \
		opennurbs_brep_io.o \
		opennurbs_brep_isvalid.o \
		opennurbs_brep_kinky.o \
		opennurbs_brep_tools.o \
		opennurbs_brep_v2valid.o \
		opennurbs_circle.o \
		opennurbs_color.o \
		opennurbs_cone.o \
		opennurbs_crc.o \
		opennurbs_curve.o \
		opennurbs_curveonsurface.o \
		opennurbs_curveproxy.o \
		opennurbs_cylinder.o \
		opennurbs_defines.o \
		opennurbs_detail.o \
		opennurbs_dimstyle.o \
		opennurbs_ellipse.o \
		opennurbs_error.o \
		opennurbs_error_message.o \
		opennurbs_evaluate_nurbs.o \
		opennurbs_extensions.o \
		opennurbs_font.o \
		opennurbs_geometry.o \
		opennurbs_group.o \
		opennurbs_hatch.o \
		opennurbs_instance.o \
		opennurbs_intersect.o \
		opennurbs_knot.o \
		opennurbs_layer.o \
		opennurbs_light.o \
		opennurbs_line.o \
		opennurbs_linecurve.o \
		opennurbs_linetype.o \
		opennurbs_massprop.o \
		opennurbs_material.o \
		opennurbs_math.o \
		opennurbs_matrix.o \
		opennurbs_memory.o \
		opennurbs_memory_util.o \
		opennurbs_mesh.o \
		opennurbs_mesh_tools.o \
		opennurbs_morph.o \
		opennurbs_nurbscurve.o \
		opennurbs_nurbssurface.o \
		opennurbs_nurbsvolume.o \
		opennurbs_object.o \
		opennurbs_object_history.o \
		opennurbs_objref.o \
		opennurbs_offsetsurface.o \
		opennurbs_optimize.o \
		opennurbs_plane.o \
		opennurbs_planesurface.o \
		opennurbs_pluginlist.o \
		opennurbs_point.o \
		opennurbs_pointcloud.o \
		opennurbs_pointgeometry.o \
		opennurbs_pointgrid.o \
		opennurbs_polycurve.o \
		opennurbs_polyline.o \
		opennurbs_polylinecurve.o \
		opennurbs_revsurface.o \
		opennurbs_sphere.o \
		opennurbs_string.o \
		opennurbs_sum.o \
		opennurbs_sumsurface.o \
		opennurbs_surface.o \
		opennurbs_surfaceproxy.o \
		opennurbs_textlog.o \
		opennurbs_torus.o \
		opennurbs_userdata.o \
		opennurbs_uuid.o \
		opennurbs_viewport.o \
		opennurbs_workspace.o \
		opennurbs_wstring.o \
		opennurbs_x.o \
		opennurbs_xform.o \
		opennurbs_zlib.o \
		opennurbs_zlib_memory.o

ZLIB_INC= ../opennurbs/zlib/crc32.h \
	../opennurbs/zlib/deflate.h \
	../opennurbs/zlib/inffast.h \
	../opennurbs/zlib/inffixed.h \
	../opennurbs/zlib/inflate.h \
	../opennurbs/zlib/inftrees.h \
	../opennurbs/zlib/trees.h \
	../opennurbs/zlib/zconf.h \
	../opennurbs/zlib/zlib.h \
	../opennurbs/zlib/zutil.h \
		
ZLIB_SRC= ../opennurbs/zlib/adler32.c \
	../opennurbs/zlib/compress.c \
	../opennurbs/zlib/crc32.c \
	../opennurbs/zlib/deflate.c \
	../opennurbs/zlib/infback.c \
	../opennurbs/zlib/inffast.c \
	../opennurbs/zlib/inflate.c \
	../opennurbs/zlib/inftrees.c \
	../opennurbs/zlib/trees.c \
	../opennurbs/zlib/uncompr.c \
	../opennurbs/zlib/zutil.c


ZLIB_OBJ=adler32.o \
	compress.o \
	crc32.o \
	deflate.o \
	infback.o \
	inffast.o \
	inflate.o \
	inftrees.o \
	trees.o \
	uncompr.o \
	zutil.o

KERNEL_INC=include/dynamicarray.h \
	include/mathglobals.h \
	include/matrix2x2.h \
	include/matrix3x3.h \
	include/matrix4x4.h \
	include/quaternion.h \
	include/vector2.h \
	include/vector3.h \
	include/vector4.h \
	include/aabb2.h \
	include/aabbtree.h \
	include/obb2.h \
	include/bctree.h \
	include/sbbtree.h \
	include/sbb.h \
	include/approxcurve.h \
	include/basiccurve.h \
	include/binstree.h \
	include/feat2ops.h \
	include/paramline.h \
	include/bspline.h \
	include/distops.h \
	include/exception.h \
	include/fileParser.h \
	include/geomops.h \
	include/graham.h \
	include/macros.h \
	include/nurbs.h \
	include/primops.h \
	include/rootfinder.h \
	include/interface.h

KERNEL_SRC= src/dynamicarray.cpp \
	src/matrix2x2.cpp \
	src/quaternion.cpp \
	src/vector2.cpp \
	src/vector3.cpp \
	src/aabb2.cpp \
	src/aabbtree.cpp \
	src/obb2.cpp \
	src/bctree.cpp \
	src/sbbtree.cpp \
	src/sbb.cpp \
	src/approxcurve.cpp \
	src/basiccurve.cpp \
	src/feat2ops.cpp \
	src/paramline.cpp \
	src/bspline.cpp \
	src/distops.cpp \
	src/exception.cpp \
	src/fileParser.cpp \
	src/geomops.cpp \
	src/graham.cpp \
	src/nurbs.cpp \
	src/primops.cpp \
	src/rootfinder.cpp \
	src/interface.cpp

KERNEL_OBJ=dynamicarray.o \
	matrix2x2.o \
	quaternion.o \
	vector2.o \
	vector3.o \
	aabbtree.o \
	aabb2.o \
	approxcurve.o \
	obb2.o \
	bctree.o \
	sbbtree.o \
	sbb.o \
	basiccurve.o \
	feat2ops.o \
	paramline.o \
	bspline.o \
	distops.o \
	exception.o \
	fileParser.o \
	geomops.o \
	graham.o \
	nurbs.o \
	primops.o \
	rootfinder.o \
	interface.o

$(ZLIB_OBJ) : $(ZLIB_INC)
	$(CC) $(CFLAGS) $(ZLIB_SRC)

zlib: $(ZLIB_OBJ)
	

opennurbs_memory.o : ../opennurbs/opennurbs_memory_default.c $(ON_INC)

$(ON_OBJ) : $(ON_INC)
	@echo "opennurbs" 
	$(CCC) $(CCFLAGS) $(LDFLAGS) $(ON_SRC)

$(KERNEL_OBJ): $(KERNEL_INC)
	@echo "kernel Inshape " 
	$(CCC) $(CCFLAGS) $(LDFLAGS) $(KERNEL_SRC) 

$(OPENNURBS_LIB_FILE) : $(KERNEL_OBJ) $(ON_OBJ) $(ZLIB_OBJ)
	@echo "Building <$(OPENNURBS_LIB_FILE)>" 
	-$(RM) $@
	$(AR) $@ $(KERNEL_OBJ) $(ON_OBJ) $(ZLIB_OBJ)  
	$(RANLIB) $@

kernel: $(KERNEL_OBJ)  
	@echo "only kernel Inshape " 
#	$(CCC) $(CCFLAGS) $(KERNEL_SRC) 

	
all: $(OPENNURBS_LIB_FILE)
	

#Inshape2d: $(KERNEL_OBJ) $(MATH_OBJ) $(UTILITY_OBJ) $(ON_OBJ) $(ZLIB_OBJ)
#	$(CCC) $(CCFLAGS) $(LDFLAGS) $(KERNEL_OBJ) $(MATH_OBJ) $(UTILITY_OBJ) $(ON_OBJ) $(ZLIB_OBJ) -o $@
	
	

