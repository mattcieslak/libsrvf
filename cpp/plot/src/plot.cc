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
#include "render.h"
#include "plot.h"

namespace srvf
{


////////////////////////////////////////////////////////////////////////////
// MontagePlot
////////////////////////////////////////////////////////////////////////////

void MontagePlot::insert(
  const Plf &F, 
  Color c, 
  double thickness, 
  DrawingMode mode)
{
  plfs_.push_back(F);
  plf_colors_.push_back(c);
  plf_thicknesses_.push_back(thickness);
  plf_modes_.push_back(mode);
}

void MontagePlot::insert(
  const Srvf &Q, 
  Color c, 
  double thickness, 
  DrawingMode mode)
{
  srvfs_.push_back(Q);
  srvf_colors_.push_back(c);
  srvf_thicknesses_.push_back(thickness);
  srvf_modes_.push_back(mode);
}

void MontagePlot::render(Renderer &r)
{
  
}

////////////////////////////////////////////////////////////////////////////
// SuperimposedPlot
////////////////////////////////////////////////////////////////////////////
void SuperimposedPlot::insert(
  const Plf &F, 
  Color c, 
  double thickness, 
  DrawingMode mode)
{
  plfs_.push_back(F);
  plf_colors_.push_back(c);
  plf_thicknesses_.push_back(thickness);
  plf_modes_.push_back(mode);
  plf_intervals_.push_back(std::pair<double,double>( 
    F.domain_lb(), F.domain_ub()) );

  std::vector<Point> bbox = F.bounding_box();
  for (size_t i=0; i<F.dim(); ++i)
  {
    bball_rad_ = std::max(bball_rad_, fabs(bbox[0][i]));
    bball_rad_ = std::max(bball_rad_, fabs(bbox[1][i]));
  }
}

void SuperimposedPlot::insert(
  const Srvf &Q, 
  Color c, 
  double thickness, 
  DrawingMode mode)
{
  srvfs_.push_back(Q);
  srvf_colors_.push_back(c);
  srvf_thicknesses_.push_back(thickness);
  srvf_modes_.push_back(mode);
}

void SuperimposedPlot::render(Renderer &r)
{
  double dev_width = (double)r.device_width();
  double dev_height = (double)r.device_height();
  double asp_rat = dev_height / dev_width;
  r.viewport(2, 2, dev_width-2, dev_height-2);
  r.ortho(-bball_rad_, bball_rad_, 
          -bball_rad_*asp_rat, bball_rad_*asp_rat, 
          -bball_rad_, bball_rad_);

  for (size_t i=0; i<plfs_.size(); ++i)
  {
    r.begin(plf_modes_[i]);
    r.set_color(plf_colors_[i]);

    double dom_lb = plf_intervals_[i].first;
    double dom_ub = plf_intervals_[i].second;

    for (size_t j=0; j<plfs_[i].ncp(); ++j)
    {
      double t1 = plfs_[i].params()[j];

      if (j+1 == plfs_[i].ncp())
      {
        // Last point.  Output it if it lies in the restricted domain.
        if (t1 < dom_ub + 1e-5)
        {
          double x = (plfs_[i].dim() > 0 ? plfs_[i].samps()[j][0] : 0.0);
          double y = (plfs_[i].dim() > 1 ? plfs_[i].samps()[j][1] : 0.0);
          double z = (plfs_[i].dim() > 2 ? plfs_[i].samps()[j][2] : 0.0);
          r.vertex(x, y, z);
        }
      }
      else
      {
        // Not the last point.  Several cases to handle here.
        double t2 = plfs_[i].params()[j+1];

        if (t1 < dom_lb && t2 > dom_lb + 1e-5)
        {
          // ith point is outside (left of) the restricted domain, but the 
          // following point is inside.  Output a linear combination 
          // of the two.
          double dt = t2-t1;
          double w1 = (t2 - dom_lb) / dt;
          double w2 = (dom_lb - t1) / dt;

          Point p = (plfs_[i].samps()[j] * w1) + (plfs_[i].samps()[j+1] * w2);
          double x = (p.dim() > 0 ? p[0] : 0.0);
          double y = (p.dim() > 1 ? p[1] : 0.0);
          double z = (p.dim() > 2 ? p[2] : 0.0);
          r.vertex(x, y, z);
        }
        else if (t1 > dom_lb - 1e-5 && t1 < dom_ub + 1e-5)
        {
          // ith point is in the restricted domain, so output it.
          double x = (plfs_[i].dim() > 0 ? plfs_[i].samps()[j][0] : 0.0);
          double y = (plfs_[i].dim() > 1 ? plfs_[i].samps()[j][1] : 0.0);
          double z = (plfs_[i].dim() > 2 ? plfs_[i].samps()[j][2] : 0.0);
          r.vertex(x, y, z);

          if (t2 > dom_ub)
          {
            // ith point is in the restricted domain, but the following 
            // point is not.  Output linear combination of the two.
            double dt = t2-t1;
            double w1 = (t2 - dom_ub) / dt;
            double w2 = (dom_ub - t1) / dt;

            Point p = (plfs_[i].samps()[j] * w1) + (plfs_[i].samps()[j+1] * w2);
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
// GeodesicPlot
////////////////////////////////////////////////////////////////////////////

void GeodesicPlot::render(Renderer &r)
{
  
}


////////////////////////////////////////////////////////////////////////////
// FunctionPlot
////////////////////////////////////////////////////////////////////////////

void FunctionPlot::render(Renderer &r)
{
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

} // namespace srvf
