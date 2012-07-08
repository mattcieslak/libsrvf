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
#include "plotwin.h"
#include "render.h"
#include "plot.h"

#include <vector>
#include <iostream>

#include <FL/Fl.H>
#include <FL/Enumerations.H>
#include <GL/gl.h>
#include <GL/glu.h>


namespace srvf
{


FltkGlPlotWindow::FltkGlPlotWindow(int w, int h, const char *l)
 : Fl_Gl_Window(w, h, l), 
   prev_x_(0), prev_y_(0), button_state_(0), 
   camera_x_(0.0), camera_y_(0.0), camera_z_(1.0)
{ }

FltkGlPlotWindow::FltkGlPlotWindow(int w, int h, int x, int y, const char *l)
 : Fl_Gl_Window(w, h, x, y, l),
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
  double sf = 1.0 / camera_z_;
  glTranslatef(camera_x_, -camera_y_, -1.0);
  glScalef(sf, sf, sf);

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
#define DZ  0.02f
int
FltkGlPlotWindow::handle(int event)
{
  switch(event) {
    case FL_PUSH:
      prev_x_ = Fl::event_x();
      prev_y_ = Fl::event_y();
      if      ( Fl::event_button() == 1 ) button_state_ |= 1;
      else if ( Fl::event_button() == 2 ) button_state_ |= 2;
      else if ( Fl::event_button() == 3 ) button_state_ |= 4;
      return 1;
    case FL_DRAG:
      if ( button_state_ == 1 ){
        camera_x_ += DX * (float)(Fl::event_x() - prev_x_);
        camera_y_ += DY * (float)(Fl::event_y() - prev_y_);
        prev_x_ = Fl::event_x();
        prev_y_ = Fl::event_y();
        redraw();
      } else if ( button_state_ == 2 ){
        camera_z_ += DZ * (float)(Fl::event_y() - prev_y_);
        camera_z_ = std::max(camera_z_, 0.05f);
        prev_x_ = Fl::event_x();
        prev_y_ = Fl::event_y();
        redraw();
      }
      return 1;
    case FL_RELEASE:   
      if      ( Fl::event_button() == 1 ) button_state_ &= 6;
      else if ( Fl::event_button() == 2 ) button_state_ &= 5;
      else if ( Fl::event_button() == 3 ) button_state_ &= 3;
      return 1;
    case FL_FOCUS :
    case FL_UNFOCUS :
      return 1;  // we want keyboard events
    case FL_KEYDOWN:
      switch(Fl::event_key())
      {
        case FL_Escape:
        case 'q':
          this->hide();
          break;
      }
      return 1;
    case FL_SHORTCUT:
      return 0;
    default:
      return Fl_Gl_Window::handle(event);
  }
}

} // namespace srvf
