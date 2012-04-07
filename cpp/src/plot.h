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
#ifndef SRVF_PLOT_H
#define SRVF_PLOT_H 1

#include "plf.h"
#include "srvf.h"
#include "render.h"

namespace srvf
{

/**
 * Abstract base for plot classes.
 */
class Plot
{ 
public:
  
  virtual void 
  render(Renderer &r) = 0;

  virtual void 
  insert(const Plf &F, Color c, double thickness=1.0, DrawingMode mode=LINES)
  {
    plfs_.push_back(F);
    colors_.push_back(c);
    thicknesses_.push_back(thickness);
    modes_.push_back(mode);
  }

  virtual void 
  insert(const Srvf &Q, Color c, double thickness=1.0, DrawingMode mode=LINES)
  {
    srvfs_.push_back(Q);
    colors_.push_back(c);
    thicknesses_.push_back(thickness);
    modes_.push_back(mode);
  }

protected:
  
  std::vector<Plf> plfs_;
  std::vector<Srvf> srvfs_;
  std::vector<Color> colors_;
  std::vector<double> thicknesses_;
  std::vector<DrawingMode> modes_;
};


class MontagePlot : public Plot
{ 
public:

  MontagePlot() 
   : nrows_(1), next_x_(0.0), next_y_(0.0)
  { }

  MontagePlot(int nrows) 
   : nrows_(nrows), next_x_(0.0), next_y_(0.0)
  { }

  virtual void 
  insert(const Plf &F, Color c, double thickness=1.0, DrawingMode mode=LINES);

  virtual void 
  insert(const Srvf &Q, Color c, double thickness=1.0, DrawingMode mode=LINES);
  
  virtual void render(Renderer &r);

private:
  
  int nrows_;
  double next_x_;
  double next_y_;
};


class SuperimposedPlot : public Plot
{ 
public:
  
  virtual void 
  insert(const Plf &F, Color c, double thickness=1.0, DrawingMode mode=LINES);

  virtual void 
  insert(const Srvf &Q, Color c, double thickness=1.0, DrawingMode mode=LINES);

  virtual void render(Renderer &r);
};


class GeodesicPlot : public Plot
{ 
public:

  virtual void render(Renderer &r);
};

} // namespace srvf

#endif // SRVF_PLOT_H
