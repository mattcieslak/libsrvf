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
#include <vector>
#include <iostream>

#include <GL/gl.h>
#include <GL/glu.h>

#include "plotwin.h"
#include "render.h"
#include "plot.h"


namespace srvf
{


FltkGlPlotWindow::FltkGlPlotWindow(int w, int h, const char *l)
 : fltk::GlWindow(w, h, l)
{ }

FltkGlPlotWindow::FltkGlPlotWindow(int w, int h, int x, int y, const char *l)
 : fltk::GlWindow(w, h, x, y, l)
{ }

void 
FltkGlPlotWindow::draw()
{
  if (!valid())
  {
    glMatrixMode(GL_PROJECTION);
    glViewport(0, 0, (GLsizei)w(), (GLsizei)h());
    double r = ((double)h()) / ((double)w());
    glOrtho(-1.0, 1.0, -r, r, -0.1, 10.0);
    glMatrixMode(GL_MODELVIEW);
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(0.0, 0.0, -1.0);

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  for (size_t i=0; i<plots_.size(); ++i)
  {
    plots_[i]->render(renderer_);
  }

  //glColor3f(1.0, 1.0, 1.0);
  //glBegin(GL_QUADS);
  //glVertex2d(0.0, 0.0);
  //glVertex2d(1.0, 0.0);
  //glVertex2d(1.0, 1.0);
  //glVertex2d(0.0, 1.0);
  //glEnd();
}

void 
FltkGlPlotWindow::add_plot(Plot *p)
{
  plots_.push_back(p);
}

} // namespace srvf
