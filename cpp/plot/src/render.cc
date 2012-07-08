/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012  Daniel Robinson
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#include "plot.h"
#include "render.h"

#include <vector>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>


namespace srvf
{


////////////////////////////////////////////////////////////////////////////
// OpenGlRenderer
////////////////////////////////////////////////////////////////////////////


OpenGlRenderer::OpenGlRenderer() 
  : mode_(LINES) 
{
}

void OpenGlRenderer::viewport(double x, double y, double w, double h)
{
  glViewport(x, y, w, h);
}

void OpenGlRenderer::ortho(double left, double right, 
           double bottom, double top, 
           double near, double far)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(left, right, bottom, top, near, far);
  glMatrixMode(GL_MODELVIEW);
}

void OpenGlRenderer::frustum(double left, double right, 
             double bottom, double top, 
             double near, double far)
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glFrustum(left, right, bottom, top, near, far);
  glMatrixMode(GL_MODELVIEW);
}

void 
OpenGlRenderer::begin(DrawingMode mode)
{
  mode_ = mode;

  if (mode == srvf::LINES)
  {
    glBegin(GL_LINE_STRIP);
  }
  else if (mode == srvf::POINTS)
  {
    glBegin(GL_POINTS);
  }
}

void 
OpenGlRenderer::vertex(double x, double y)
{
  if (mode_ == srvf::POINTS || mode_ == srvf::LINES)
  {
    glVertex2d(x, y);
  }
}

void 
OpenGlRenderer::vertex(double x, double y, double z)
{
  if (mode_ == srvf::POINTS || mode_ == srvf::LINES)
  {
    glVertex3d(x, y, z);
  }
}

void 
OpenGlRenderer::end()
{
  if (mode_ == srvf::POINTS || mode_ == srvf::LINES)
  {
    glEnd();
  }
}

void 
OpenGlRenderer::set_color(Color c)
{
  glColor3f(c.red, c.green, c.blue);
}

void 
OpenGlRenderer::set_thickness(double t)
{
  thickness_ = t;

  if (mode_ == srvf::LINES)
  {
    glLineWidth((GLfloat)t);
  }
}

void OpenGlRenderer::translate(double x, double y)
{
  glTranslated(x, y, 0.0);
}

void OpenGlRenderer::translate(double x, double y, double z)
{
  glTranslated(x, y, z);
}

void 
OpenGlRenderer::scale(double sfx, double sfy)
{
  glScaled(sfx, sfy, 1.0);
}

void 
OpenGlRenderer::scale(double sfx, double sfy, double sfz)
{
  glScaled(sfx, sfy, sfz);
}

void 
OpenGlRenderer::rotate(double angle)
{
  glRotated(0.0, 0.0, 1.0, (GLfloat)angle);
}

void 
OpenGlRenderer::rotate(double ax, double ay, double az, double angle)
{
  glRotated((GLfloat)ax, (GLfloat)ay, (GLfloat)az, (GLfloat)angle);
}



} // namespace srvf
