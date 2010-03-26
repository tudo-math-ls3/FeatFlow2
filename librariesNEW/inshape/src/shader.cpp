/***************************************************************************
 *   Copyright (C) 2006 by Raphael Münster   *
 *   raphael@Cortez   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include<fstream>
#include "shader.h"


CShader::CShader()
{
	m_VertexShader   = 0;
	m_FragmentShader = 0;
	m_ProgramObject  = 0;
}//end CShader

////////////////////////////////////////////////////
//
//
//	loads the shader source code from the file
//
//	specified in strFile
//
//
///////////////////////////////////////////////////

string CShader::LoadShaderSource(std::string strFile)
{

	ifstream inFile;
	inFile.open(strFile.c_str());

	if(inFile.fail())
		return "";

	string strLine = "";
	string strText = "";

	/* read in file line by line */
	while(getline(inFile, strLine))
	{
		/* concatenate every line in the file */
		strText = strText + "\n" + strLine;

	}//end while

	inFile.close();

	return strText;

}//end loadShaderSource

void printInfoLog(GLhandleARB obj)
	{
	    int infologLength = 0;
	    int charsWritten  = 0;
	    char *infoLog;
	
	    glGetObjectParameterivARB(obj, GL_OBJECT_INFO_LOG_LENGTH_ARB,
						 &infologLength);
	
	    if (infologLength > 0)
	    {
		infoLog = (char *)malloc(infologLength);
		glGetInfoLogARB(obj, infologLength, &charsWritten, infoLog);
			printf("%s\n",infoLog);
		free(infoLog);
	    }
	}

void CShader::InitShaders(std::string strVertex, std::string strFragment)
{

	/* stores the shader source */
	string sVertexShader, sFragmentShader;


	/* check for valid filenames */
	if( !strVertex.length() || !strFragment.length() )
	{
		cerr << "Invalid string parameter in initShaders" <<endl;
		return;
	}//end if

	if(m_VertexShader || m_FragmentShader || m_ProgramObject)
	{
		
		Destroy();
	}//end if
	
	cout<<"creating shader objects..."<<endl;

	/* Create a shader object for the vertex shader */
	m_VertexShader   = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
	
	/* Create a shader object for the fragment shader */
	m_FragmentShader = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);

	cout<<"loading shader source..."<<endl;

	/* load the shader source */
	sVertexShader    = LoadShaderSource(strVertex);

	/* load the shader source */
	sFragmentShader  = LoadShaderSource(strFragment);
	
	const char *szVertexShader   = sVertexShader.c_str();
	const char *szFragmentShader = sFragmentShader.c_str();

	/* Set shader source */
	glShaderSourceARB(m_VertexShader, 1, &szVertexShader, NULL);

	/* Set shader source */
	glShaderSourceARB(m_FragmentShader, 1, &szFragmentShader, NULL);

	/* Compile shader */
	glCompileShaderARB(m_VertexShader);
	printInfoLog(m_VertexShader);

	/* Compile shader */
	glCompileShaderARB(m_FragmentShader);
	printInfoLog(m_FragmentShader);

	/* Create a program object */
	m_ProgramObject  = glCreateProgramObjectARB();

	/* Attach shader to program */
	glAttachObjectARB(m_ProgramObject, m_VertexShader);

	/* Attach shader to program */
	glAttachObjectARB(m_ProgramObject, m_FragmentShader);

	/* Link the program */
	glLinkProgramARB(m_ProgramObject);

	/* Enable the program object */
	glUseProgramObjectARB(m_ProgramObject);

}//end initShaders

void CShader::Destroy()
{

	if(m_VertexShader)
	{
		glDetachObjectARB(m_ProgramObject, m_VertexShader);
		glDeleteObjectARB(m_VertexShader);
		m_VertexShader = 0;
	}

	if(m_FragmentShader)
	{
		glDetachObjectARB(m_ProgramObject, m_FragmentShader);
		glDeleteObjectARB(m_FragmentShader);
		m_FragmentShader = 0;
	}

	if(m_ProgramObject)
	{
		glDeleteObjectARB(m_ProgramObject);
		m_ProgramObject = 0;
	}

}//end deleteShaders

GLint CShader::GetVariable(string strVariable)
{
	/* error checking */
	if(!m_ProgramObject)
		return -1;

	// This returns the variable ID for a variable that is used to find
	// the address of that variable in memory.
	return glGetUniformLocationARB(m_ProgramObject, strVariable.c_str());
}
