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
#include <srvf/render.h>
#include <srvf/plot.h>

#include <iostream>

namespace srvf
{

namespace plot
{


static void render_plfs_ (
  const std::vector<Plf> &plfs, 
  const std::vector<bool> &visibilities,
  const std::vector<std::pair<double,double> > &intervals,
  const std::vector<DrawingMode> &modes,
  const std::vector<Color> &colors,
  const std::vector<double> &thicknesses,
  Renderer &r)
{
  for (size_t i=0; i<plfs.size(); ++i)
  {
    if (!visibilities[i]) continue;

    r.set_thickness(thicknesses[i]);
    r.begin(modes[i]);
    r.set_color(colors[i]);

    double dom_lb = intervals[i].first;
    double dom_ub = intervals[i].second;

    for (size_t j=0; j<plfs[i].ncp(); ++j)
    {
      double t1 = plfs[i].params()[j];

      if (j+1 == plfs[i].ncp())
      {
        // Last point.  Output it if it lies in the restricted domain.
        if (t1 < dom_ub + 1e-5)
        {
          double x = (plfs[i].dim() > 0 ? plfs[i].samps()[j][0] : 0.0);
          double y = (plfs[i].dim() > 1 ? plfs[i].samps()[j][1] : 0.0);
          double z = (plfs[i].dim() > 2 ? plfs[i].samps()[j][2] : 0.0);
          r.vertex(x, y, z);
        }
      }
      else
      {
        // Not the last point.  Several cases to handle here.
        double t2 = plfs[i].params()[j+1];

        if (t1 < dom_lb && t2 > dom_lb + 1e-5)
        {
          // ith point is outside (left of) the restricted domain, but the 
          // following point is inside.  Output a linear combination 
          // of the two.
          double dt = t2-t1;
          double w1 = (t2 - dom_lb) / dt;
          double w2 = (dom_lb - t1) / dt;

          Point p = (plfs[i].samps()[j] * w1) + (plfs[i].samps()[j+1] * w2);
          double x = (p.dim() > 0 ? p[0] : 0.0);
          double y = (p.dim() > 1 ? p[1] : 0.0);
          double z = (p.dim() > 2 ? p[2] : 0.0);
          r.vertex(x, y, z);
        }
        else if (t1 > dom_lb - 1e-5 && t1 < dom_ub + 1e-5)
        {
          // ith point is in the restricted domain, so output it.
          double x = (plfs[i].dim() > 0 ? plfs[i].samps()[j][0] : 0.0);
          double y = (plfs[i].dim() > 1 ? plfs[i].samps()[j][1] : 0.0);
          double z = (plfs[i].dim() > 2 ? plfs[i].samps()[j][2] : 0.0);
          r.vertex(x, y, z);

          if (t2 > dom_ub)
          {
            // ith point is in the restricted domain, but the following 
            // point is not.  Output linear combination of the two.
            double dt = t2-t1;
            double w1 = (t2 - dom_ub) / dt;
            double w2 = (dom_ub - t1) / dt;

            Point p = (plfs[i].samps()[j] * w1) + (plfs[i].samps()[j+1] * w2);
            x = (p.dim() > 0 ? p[0] : 0.0);
            y = (p.dim() > 1 ? p[1] : 0.0);
            z = (p.dim() > 2 ? p[2] : 0.0);
            r.vertex(x, y, z);
          }
        }
      }
    }

    r.end();
  }
}


////////////////////////////////////////////////////////////////////////////
// Plot2D
////////////////////////////////////////////////////////////////////////////

void Plot2D::render(Renderer &r)
{
  double asp_rat = ((double)r.device_height()) / ((double)r.device_width());
  r.ortho(-plot_radius_*1.5, plot_radius_*1.5, 
          -plot_radius_*asp_rat*1.5, plot_radius_*asp_rat*1.5, 
          0.01, 1.5);

  r.translate(trans_x_, trans_y_, -0.1);
  r.rotate(rot_theta_);
  r.scale(scale_, scale_);
  render_plfs_(plfs_, plf_visibilities_, 
    plf_intervals_, plf_modes_, plf_colors_, 
    plf_thicknesses_, r);
}

void Plot2D::pan_view(double dx, double dy)
{
  trans_x_ += 3.0*dx;
  trans_y_ -= 3.0*dy;  // invert the mouse!
}

void Plot2D::scale_view(double dx, double dy)
{
  scale_ = std::max(scale_-3*dy, 0.01);
}

void Plot2D::rotate_view(double dx, double dy)
{
  rot_theta_ -= dx * 360.0;
}

////////////////////////////////////////////////////////////////////////////
// Plot3D
////////////////////////////////////////////////////////////////////////////

void Plot3D::render(Renderer &r)
{
  double asp_rat = ((double)r.device_height()) / ((double)r.device_width());
  r.ortho(-plot_radius_*1.5, plot_radius_*1.5, 
          -plot_radius_*asp_rat*1.5, plot_radius_*asp_rat*1.5, 
          -plot_radius_*1.5, plot_radius_*1.5);

  r.rotate(hrot_, 0.0, 1.0, 0.0);
  r.rotate(vrot_, 1.0, 0.0, 0.0);
  r.scale(scale_, scale_, scale_);
  render_plfs_(plfs_, plf_visibilities_,
    plf_intervals_, plf_modes_, plf_colors_, 
    plf_thicknesses_, r);
}

void Plot3D::scale_view(double dx, double dy)
{
  scale_ = std::max(scale_ - 3.0*dy, 0.05);
}

void Plot3D::rotate_view(double dx, double dy)
{
  hrot_ += dx * 360.0;
  vrot_ -= dy * 360.0;
}

////////////////////////////////////////////////////////////////////////////
// FunctionPlot
////////////////////////////////////////////////////////////////////////////

void FunctionPlot::render(Renderer &r)
{
  double x_min=1e9, x_max=-1e9;
  double y_min=1e9, y_max=-1e9;

  for (size_t i=0; i<plfs_.size(); ++i)
  {
    x_min = std::min(x_min, plfs_[i].domain_lb());
    x_max = std::max(x_max, plfs_[i].domain_ub());

    y_min = std::min(y_min, bounding_boxes_[i][0][0]);
    y_max = std::max(y_max, bounding_boxes_[i][1][0]);
  }

  // Set up the rendering context
  double dev_width = (double)r.device_width();
  double dev_height = (double)r.device_height();
  r.viewport(2, 2, dev_width-2, dev_height-2);
  r.ortho(x_min, x_max, y_min, y_max, -1.0, 1.0);

  // Draw the Plf's
  for (size_t i=0; i<plfs_.size(); ++i)
  {
    r.set_color(plf_colors_[i]);
    r.begin(LINES);
    for (size_t j=0; j<plfs_[i].ncp(); ++j)
    {
      r.vertex(plfs_[i].params()[j], plfs_[i].samps()[j][0]);
    }
    r.end();
  }

  // Draw the Srvf's
  for (size_t i=0; i<srvfs_.size(); ++i)
  {
    r.set_color(srvf_colors_[i]);
    for (size_t j=1; j<srvfs_[i].ncp(); ++j)
    {
      r.begin(LINES);
      r.vertex(srvfs_[i].params()[j-1], srvfs_[i].samps()[j-1][0]);
      r.vertex(srvfs_[i].params()[j], srvfs_[i].samps()[j-1][0]);
      r.end();
    }
  }

  // Draw some lame axes
  r.set_color(Color(1.0,1.0,1.0));
  r.set_thickness(1.0);
  r.begin(LINES);
  r.vertex(x_min, 0.0);
  r.vertex(x_max, 0.0);
  r.end();
  r.begin(LINES);
  r.vertex(0.0, y_min);
  r.vertex(0.0, y_max);
  r.end();
}


} //namespace plot

} // namespace srvf
