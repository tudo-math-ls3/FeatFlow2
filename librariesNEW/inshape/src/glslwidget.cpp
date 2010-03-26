#include "glslwidget.h"
#include "shader.h"
#include "visualdebug.h"
#include "quaternion.h"

CVisualDebug g_Renderer;


GlslWidget::GlslWidget(QWidget *parent) : GLWidget(parent)
{

}//end constructor

GlslWidget::~GlslWidget(void)
{

}//end deconstructor

void GlslWidget::initializeGL()
{

 	/* set OpenGL clear color to black */
	qglClearColor( Qt::white );

	glewInit();
	if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader)
			printf("Ready for GLSL\n");
	else
	{
			printf("Not totally ready :( \n");
			exit(1);
	}

	MakeTextures();
	
	/* enable smooth shading */
	glShadeModel( GL_SMOOTH );

	float shine =  20.0f;

	//enable lighting
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT , light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE , light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glMaterialfv(GL_FRONT, GL_SPECULAR, light_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS,&shine);
	glEnable(GL_COLOR_MATERIAL);

	//enable z-Buffer
    glEnable(GL_DEPTH_TEST);

	//initialize the shader
	m_GLSLShader.InitShaders("lighting.vert", "lighting.frag");
	m_GLSLShader.DisableShader();

	//disable the standard lighting
	glDisable(GL_LIGHTING);

}//end initializeGL

void GlslWidget::resizeGL(int w, int h)
{

	//establish viewport settings
 	glViewport(0,0,w,h);
	m_nWidth  = w;
	m_nHeight = h;
	
	//set a perspective projection
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluPerspective(60,(float)w/h,0.1,100);
	
	//reset the modelview matrix
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();

}//resizeGL

void GlslWidget::paintGL() 
{
	GLWidget::paintGL();

 //   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 //   glLoadIdentity();

	//glTranslatef(m_TransX, m_TransY, m_TransZ);

	//CQuaternionf qRotX;

	//qRotX.AxisAngleToQuat(0,1,0, m_RotX);

	//CMatrix4x4f matrix;

	//qRotX.CreateMatrix(matrix);

	//CQuaternionf qRotY;

	//qRotY.AxisAngleToQuat(1,0,0, m_RotY);

	//glMultMatrixf(matrix.Get());

	//qRotY.CreateMatrix(matrix);

	//glMultMatrixf(matrix.Get());

	////draw all objects
	//for(int i = 0; i < m_nNumObjects; i++)
	//{
	//	g_Renderer.DrawApproxCurve(ApproxCurves[i], 0,0,0);
	//}//end for

	//if(m_bDistance && !ApproxCurves.empty())
	//{

	//	GLWidget::GetFictitiousBoundaryValues();
	//	//GLWidget::RegLubSingle();

	//}//end if

	//if(m_bDistanceField)
	//{
	//	g_Grid.m_nMaxDist = 10.0;
	//	glEnable(GL_TEXTURE_2D);
	//	glBindTexture(GL_TEXTURE_2D, g_Textures[0]);
	//	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	//	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	//	glTexCoordPointer(2, GL_FLOAT, 0, g_Grid.m_pTCoords);
	//	glEnableClientState(GL_VERTEX_ARRAY);
	//	glVertexPointer(2, GL_FLOAT, 0, g_Grid.m_pVertices);
	//	g_Renderer.DrawtGrid(g_Grid, GRID_POINTS_X, GRID_POINTS_Y);
	//	glDisable(GL_TEXTURE_2D);
	//}

}//end paintGL
