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

#if !defined _CSHADER_H_
#define _CSHADER_H_

#include<iostream>
#include<string>
#include <GL/glew.h>



using namespace std;

class CShader
{

public:

	/* constructor */
	CShader();
	~CShader(){};

	/* member functions */

	/* initialize vertex and fragment shaders */
	void InitShaders(string strVertex, string strFragment);

	// This returns an ID for a variable in our shader
	GLint GetVariable(string strVariable);

	/* loads a source file */
	string LoadShaderSource(string strFile);

	/* shorthands to functions */
	inline void SetInt(GLint variable, int newValue) const							{ glUniform1iARB(variable, newValue);		};
	inline void SetFloat(GLint variable, float newValue) const						{ glUniform1fARB(variable, newValue);		};
	inline void SetFloat2(GLint variable, float v0, float v1) const					{ glUniform2fARB(variable, v0, v1);			};
	inline void SetFloat3(GLint variable, float v0, float v1, float v2) const		{ glUniform3fARB(variable, v0, v1, v2);		};

	inline void SetFloat4(GLint variable, float v0, float v1, float v2, float v3) const
	{ glUniform4fARB(variable, v0, v1, v2, v3);	};

	inline GLhandleARB& GetPObject() {return m_ProgramObject;};
	inline GLhandleARB& GetVShader() {return m_VertexShader;};
	inline GLhandleARB& GetFShader() {return m_FragmentShader;};

	/* functions to enable and disable the shader */
	inline void EnableShader()   {glUseProgramObjectARB(m_ProgramObject);};
	inline void DisableShader()  {glUseProgramObjectARB(0);};

	/* cleanUP */
	void Destroy();

private:

	/* handle to vertex shader */
	GLhandleARB m_VertexShader;

	/* handle to fragment shader */
	GLhandleARB m_FragmentShader;

	/* handle to shader program */
	GLhandleARB m_ProgramObject;


};

#endif