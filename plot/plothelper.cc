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
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/plot.h>
#include <srvf/plotwin.h>
#include <srvf/plothelper.h>

#include <FL/Fl.H>


namespace srvf
{

namespace plot
{

/**
 * High-level convenience routine for plotting a vector of 1-D Plf's.
 */
void plot_1d_plfs(const std::vector<Plf> &v, 
  const std::vector<Color> &colors,
  size_t x, size_t y, size_t w, size_t h, const char *title)
{
  FunctionPlot plot;
  for (size_t i=0; i<v.size(); ++i){
    if (!colors.empty())
      plot.insert( v[i], colors[(i % colors.size())] );
    else
      plot.insert( v[i], Color(0.0, 0.0, 1.0) );
  }

  FltkGlPlotWindow plotwin(x, y, w, h, title);
  plotwin.add_plot(&plot);
  plotwin.show();
  Fl::run();
}

/**
 * High-level convenience routine for plotting a vector of 1-D Srvf's.
 */
void plot_1d_srvfs(const std::vector<Srvf> &v, 
  const std::vector<Color> &colors,
  size_t x, size_t y, size_t w, size_t h, const char *title)
{
  FunctionPlot plot;
  for (size_t i=0; i<v.size(); ++i){
    if (!colors.empty())
      plot.insert( v[i], colors[(i % colors.size())] );
    else
      plot.insert( v[i], Color(0.0, 0.0, 1.0) );
  }

  FltkGlPlotWindow plotwin(x, y, w, h, title);
  plotwin.add_plot(&plot);
  plotwin.show();
  Fl::run();
}

} // namespace plot
} // namespace srvf
