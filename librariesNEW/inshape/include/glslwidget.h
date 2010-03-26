/****************************************************************************
**
** Copyright (C) 2005-2007 Trolltech ASA. All rights reserved.
**
** This file is part of the example classes of the Qt Toolkit.
**
** This file may be used under the terms of the GNU General Public
** License version 2.0 as published by the Free Software Foundation
** and appearing in the file LICENSE.GPL included in the packaging of
** this file.  Please review the following information to ensure GNU
** General Public Licensing requirements will be met:
** http://www.trolltech.com/products/qt/opensource.html
**
** If you are unsure which license is appropriate for your use, please
** review the following information:
** http://www.trolltech.com/products/qt/licensing.html or contact the
** sales department at sales@trolltech.com.
**
** This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
** WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
**
****************************************************************************/

#ifndef GLSLWIDGET_H
#define GLSLWIDGET_H

#ifdef WIN32
#pragma once
#endif

#include "shader.h"
#include "glwidget.h"

class CShader;

class GlslWidget : public GLWidget
{
public:

    GlslWidget(QWidget *parent = 0);
    virtual ~GlslWidget();

    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int w, int h);

	CShader m_GLSLShader;

};

#endif
