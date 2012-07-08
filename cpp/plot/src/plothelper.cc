#include "plf.h"
#include "srvf.h"
#include "plot.h"
#include "plotwin.h"

#include <FL/Fl.H>

#include "plothelper.h"

namespace srvf
{

/**
 * High-level convenience routine for plotting a vector of 1-D Plf's.
 */
void plot_1d_plfs(const std::vector<Plf> &v, 
  size_t x, size_t y, size_t w, size_t h, const char *title)
{
  srvf::FunctionPlot plot;
  for (size_t i=0; i<v.size(); ++i){
    plot.insert( v[i], srvf::Color(0.0, 0.0, 1.0) );
  }

  srvf::FltkGlPlotWindow plotwin(x, y, w, h, title);
  plotwin.add_plot(&plot);
  plotwin.show();
  Fl::run();
}

/**
 * High-level convenience routine for plotting a vector of 1-D Srvf's.
 */
void plot_1d_srvfs(const std::vector<Srvf> &v, 
  size_t x, size_t y, size_t w, size_t h, const char *title)
{
  srvf::FunctionPlot plot;
  for (size_t i=0; i<v.size(); ++i){
    plot.insert( v[i], srvf::Color(0.0, 0.0, 1.0) );
  }

  srvf::FltkGlPlotWindow plotwin(x, y, w, h, title);
  plotwin.add_plot(&plot);
  plotwin.show();
  Fl::run();
}

} // namespace srvf
