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
 *
 * A \c Plot object is responsible for managing a collection of \c Plf 
 * and \c Srvf instances, and for defining the geometry of all of 
 * these functions.  \c Plot objects do not handle viewing.
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
    plf_colors_.push_back(c);
    plf_thicknesses_.push_back(thickness);
    plf_modes_.push_back(mode);
  }

  virtual void 
  insert(const Srvf &Q, Color c, double thickness=1.0, DrawingMode mode=LINES)
  {
    srvfs_.push_back(Q);
    srvf_colors_.push_back(c);
    srvf_thicknesses_.push_back(thickness);
    srvf_modes_.push_back(mode);
  }

protected:
  
  std::vector<Plf> plfs_;
  std::vector<Srvf> srvfs_;
  std::vector<Color> plf_colors_;
  std::vector<Color> srvf_colors_;
  std::vector<double> plf_thicknesses_;
  std::vector<double> srvf_thicknesses_;
  std::vector<DrawingMode> plf_modes_;
  std::vector<DrawingMode> srvf_modes_;
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


/**
 * A specialized plot for 1-D functions.
 */
class FunctionPlot : public Plot
{ 
public:

  virtual void 
  insert(const Plf &F, Color c, double thickness=1.0, DrawingMode mode=LINES)
  {
    plfs_.push_back(F);
    plf_colors_.push_back(c);
    plf_thicknesses_.push_back(thickness);
    plf_modes_.push_back(mode);

    std::vector<Point> bbox = F.bounding_box();
    x_min = std::min(x_min, F.params().front());
    x_max = std::max(x_max, F.params().back());
    y_min = std::min(y_min, bbox[0][0]);
    y_max = std::max(y_max, bbox[1][0]);
  }

  virtual void 
  insert(const Srvf &Q, Color c, double thickness=1.0, DrawingMode mode=LINES)
  {
    srvfs_.push_back(Q);
    srvf_colors_.push_back(c);
    srvf_thicknesses_.push_back(thickness);
    srvf_modes_.push_back(mode);
  }

  virtual void render(Renderer &r);

private:
  
  double x_min, x_max;
  double y_min, y_max;
};


} // namespace srvf

#endif // SRVF_PLOT_H
