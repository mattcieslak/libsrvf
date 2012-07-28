/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012  FSU Statistical Shape Analysis and Modeling Group
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
   prev_x_(0), prev_y_(0), button_state_(0)
{ }

FltkGlPlotWindow::FltkGlPlotWindow(int w, int h, int x, int y, const char *l)
 : Fl_Gl_Window(w, h, x, y, l),
   prev_x_(0), prev_y_(0), button_state_(0)
{ }

void 
FltkGlPlotWindow::draw()
{
  renderer_.device_width(w());
  renderer_.device_height(h());
  renderer_.viewport(2, 2, w()-2, h()-2);
  
  renderer_.clear_color(Color(1.0,1.0,1.0));
  renderer_.clear();

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


int
FltkGlPlotWindow::handle(int event)
{
  double dx, dy;

  switch(event) {
    case FL_PUSH:
      prev_x_ = Fl::event_x();
      prev_y_ = Fl::event_y();
      if      ( Fl::event_button() == 1 ) button_state_ |= 1; // left down
      else if ( Fl::event_button() == 2 ) button_state_ |= 2; // middle down
      else if ( Fl::event_button() == 3 ) button_state_ |= 4; // right down
      return 1;
    case FL_DRAG:
      dx = (double)(Fl::event_x() - prev_x_) / (double)w();
      dy = (double)(Fl::event_y() - prev_y_) / (double)h();
      prev_x_ = Fl::event_x();
      prev_y_ = Fl::event_y();

      for (size_t i=0; i<plots_.size(); ++i)
      {
        switch(button_state_)
        {
          std::cout << button_state_ << std::endl;
          case 1: // left button down only
            plots_[i]->pan_view(dx, dy);
            break;
          case 2: // middle button down only
            break;
          case 4: // right button down only
            plots_[i]->rotate_view(dx, dy);
            break;
          case 5: // left and right buttons down
            plots_[i]->scale_view(dx, dy);
            break;
        }
      }
      redraw();
      return 1;
    case FL_RELEASE:   
      if      ( Fl::event_button() == 1 ) button_state_ &= 6; // left up
      else if ( Fl::event_button() == 2 ) button_state_ &= 5; // middle up
      else if ( Fl::event_button() == 3 ) button_state_ &= 3; // right up
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
