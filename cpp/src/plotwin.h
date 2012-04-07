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
#ifndef SRVF_PLOTWIN_H
#define SRVF_PLOTWIN_H 1

#include <vector>
#include <fltk/GlWindow.h>

#include "render.h"
#include "plot.h"


namespace srvf
{

/**
 * A basic plotting window based on FLTK GlWindow.
 */
class FltkGlPlotWindow : public fltk::GlWindow
{
public:
  
  FltkGlPlotWindow(int w, int h, const char *l);
  FltkGlPlotWindow(int w, int h, int x, int y, const char *l);

  virtual void draw();

  void add_plot(Plot *p);

private:
  
  OpenGlRenderer renderer_;
  std::vector<Plot*> plots_;
};

} // namespace srvf

#endif // SRVF_PLOTWIN_H
