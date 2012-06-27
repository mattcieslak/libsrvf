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

#include <fltk/events.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "plotwin.h"
#include "render.h"
#include "plot.h"


namespace srvf
{


FltkGlPlotWindow::FltkGlPlotWindow(int w, int h, const char *l)
 : fltk::GlWindow(w, h, l), 
   prev_x_(0), prev_y_(0), button_state_(0), 
   camera_x_(0.0), camera_y_(0.0), camera_z_(1.0)
{ }

FltkGlPlotWindow::FltkGlPlotWindow(int w, int h, int x, int y, const char *l)
 : fltk::GlWindow(w, h, x, y, l),
   prev_x_(0), prev_y_(0), button_state_(0), 
   camera_x_(0.0), camera_y_(0.0), camera_z_(1.0)
{ }

void 
FltkGlPlotWindow::draw()
{
  //if (!valid())
  //{
  //  glMatrixMode(GL_PROJECTION);
  //  glLoadIdentity();
  //  glViewport(0, 0, (GLsizei)w(), (GLsizei)h());
  //  double r = ((double)h()) / ((double)w());
  //  glOrtho(-1.0, 1.0, -r, r, 0.1, 10.0);
  //  glMatrixMode(GL_MODELVIEW);
  //}

  renderer_.device_width(w());
  renderer_.device_height(h());
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  //glTranslatef(camera_x_, -camera_y_, -1.0);
  //glScalef(camera_z_, camera_z_, camera_z_);

  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  for (size_t i=0; i<plots_.size(); ++i)
  {
    glPushMatrix();
    plots_[i]->render(renderer_);
    glPopMatrix();
  }
}

void 
FltkGlPlotWindow::add_plot(Plot *p)
{
  plots_.push_back(p);
}


#define DX  0.02f
#define DY  0.02f
#define DZ  0.05f
int
FltkGlPlotWindow::handle(int event)
{
  switch(event) {
    case fltk::PUSH:
      prev_x_ = fltk::event_x();
      prev_y_ = fltk::event_y();
      if      ( fltk::event_button() == 1 ) button_state_ |= 1;
      else if ( fltk::event_button() == 2 ) button_state_ |= 2;
      else if ( fltk::event_button() == 3 ) button_state_ |= 4;
      return 1;
    case fltk::DRAG:
      if ( button_state_ == 1 ){
        camera_x_ += DX * (float)(fltk::event_x() - prev_x_);
        camera_y_ += DY * (float)(fltk::event_y() - prev_y_);
        prev_x_ = fltk::event_x();
        prev_y_ = fltk::event_y();
        redraw();
      } else if ( button_state_ == 2 ){
        camera_z_ += DZ * (float)(fltk::event_y() - prev_y_);
        prev_x_ = fltk::event_x();
        prev_y_ = fltk::event_y();
        redraw();
      }
      return 1;
    case fltk::RELEASE:   
      if      ( fltk::event_button() == 1 ) button_state_ &= 6;
      else if ( fltk::event_button() == 2 ) button_state_ &= 5;
      else if ( fltk::event_button() == 3 ) button_state_ &= 3;
      return 1;
    case fltk::FOCUS :
    case fltk::UNFOCUS :
      return 1;  // we want keyboard events
    case fltk::KEY:
      switch(fltk::event_key())
      {
        case fltk::EscapeKey:
        case 'q':
          this->destroy();
          break;
      }
      return 1;
    case fltk::SHORTCUT:
      return 0;
    default:
      return GlWindow::handle(event);
  }
}

} // namespace srvf
