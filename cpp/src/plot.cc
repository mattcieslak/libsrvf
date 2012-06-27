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
  colors_.push_back(c);
  thicknesses_.push_back(thickness);
  modes_.push_back(mode);
}

void MontagePlot::insert(
  const Srvf &Q, 
  Color c, 
  double thickness, 
  DrawingMode mode)
{
  srvfs_.push_back(Q);
  colors_.push_back(c);
  thicknesses_.push_back(thickness);
  modes_.push_back(mode);
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
  colors_.push_back(c);
  thicknesses_.push_back(thickness);
  modes_.push_back(mode);
}

void SuperimposedPlot::insert(
  const Srvf &Q, 
  Color c, 
  double thickness, 
  DrawingMode mode)
{
  srvfs_.push_back(Q);
  colors_.push_back(c);
  thicknesses_.push_back(thickness);
  modes_.push_back(mode);
}

void SuperimposedPlot::render(Renderer &r)
{
  for (size_t i=0; i<plfs_.size(); ++i)
  {
    r.begin(modes_[i]);
    r.set_color(colors_[i]);

    if (plfs_[i].dim()==2)
    {
      for (size_t j=0; j<plfs_[i].ncp(); ++j)
      {
        r.vertex(plfs_[i].samps()[j][0], 
                 plfs_[i].samps()[j][1]);
      }
    }
    else if (plfs_[i].dim()>=3)
    {
      for (size_t j=0; j<plfs_[i].ncp(); ++j)
      {
        r.vertex(plfs_[i].samps()[j][0], 
                 plfs_[i].samps()[j][1], 
                 plfs_[i].samps()[j][2]);
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

  for (size_t i=0; i<plfs_.size(); ++i)
  {
    r.set_color(colors_[i]);
    r.begin(LINES);
    for (size_t j=0; j<plfs_[i].ncp(); ++j)
    {
      r.vertex(plfs_[i].params()[j], plfs_[i].samps()[j][0]);
    }
    r.end();
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
