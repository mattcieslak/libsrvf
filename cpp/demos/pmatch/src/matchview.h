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
#ifndef MATCHVIEW_H
#define MATCHVIEW_H 1

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <GL/gl.h>
#include <GL/glu.h>

#include <plf.h>
#include <srvf.h>
#include <qmap.h>
#include <plot.h>
#include <rotate.h>
#include <paretoset.h>


class MatchView : public Fl_Gl_Window
{
public:

  MatchView(int x, int y, int w, int h)
   : Fl_Gl_Window(x, y, w, h),
     prev_x_(0), prev_y_(0), button_state_(0), 
     camera_x_(0.0), camera_y_(0.0), camera_z_(1.0)
  { }

  MatchView(int x, int y, int w, int h, const char *l)
   : Fl_Gl_Window(x, y, w, h, l),
     prev_x_(0), prev_y_(0), button_state_(0), 
     camera_x_(0.0), camera_y_(0.0), camera_z_(1.0)
  { }

  void set_plot(const srvf::SuperimposedPlot &p)
  { 
    plot_ = p; 
    redraw();
  }

  void set_matches(const srvf::pmatch::ParetoSet &S)
  {
    matches_ = S;
  }

  void set_interval(size_t idx, double a, double b)
  {
    plot_.set_plf_interval(idx, std::pair<double,double>(a,b));
    redraw();
  }

  void set_match(size_t bucket_idx, size_t match_idx)
  {
    if ( (bucket_idx < matches_.nbuckets()) && 
         (match_idx < matches_[bucket_idx].size()) )
    {
      double a = matches_[bucket_idx][match_idx].a;
      double b = matches_[bucket_idx][match_idx].b;
      double c = matches_[bucket_idx][match_idx].c;
      double d = matches_[bucket_idx][match_idx].d;

      plot_.set_plf_interval (0, std::pair<double,double>(a,b));
      plot_.set_plf_interval (1, std::pair<double,double>(c,d));

      srvf::Srvf Q1 = srvf::plf_to_srvf(plot_.get_plf(0));
      srvf::Srvf Q2 = srvf::plf_to_srvf(plot_.get_plf(1));
      srvf::Matrix R = srvf::optimal_rotation(Q1, Q2, a, b, c, d);
      plot_.get_plf(1).rotate(R);
      redraw();
    }
  }

  size_t nbuckets()
  { return matches_.nbuckets(); }

  size_t bucket_size(size_t bucket_idx)
  { return matches_[bucket_idx].size(); }

  double match_dist(size_t bucket_idx, size_t match_idx)
  { return matches_[bucket_idx][match_idx].dist; }

  double match_length(size_t bucket_idx, size_t match_idx)
  { return matches_[bucket_idx][match_idx].length(); }

  virtual void draw()
  {
    renderer_.device_width(w());
    renderer_.device_height(h());
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    double sf = 1.0 / camera_z_;
    glTranslatef(camera_x_, -camera_y_, -1.0);
    glScalef(sf, sf, sf);

    glPushMatrix();
    plot_.render(renderer_);
    glPopMatrix();
  }

#define DX  0.02f
#define DY  0.02f
#define DZ  0.02f
  int handle(int event)
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
        return 0;
      case FL_SHORTCUT:
        return 0;
      default:
        return Fl_Gl_Window::handle(event);
    }
  }

private:  
  srvf::SuperimposedPlot plot_;
  srvf::OpenGlRenderer renderer_;
  srvf::pmatch::ParetoSet matches_;

  // Mouse state
  int prev_x_;
  int prev_y_;
  int button_state_;

  // Camera position
  float camera_x_;
  float camera_y_;
  float camera_z_;
};

#endif // MATCHVIEW_H
