#include <boost/test/unit_test.hpp>

#include <srvf/plot.h>
#include <srvf/plotwin.h>
#include <srvf/plf.h>
#include <srvf/srvf.h>

#include <FL/Fl.H>
#include <vector>


BOOST_AUTO_TEST_SUITE(plot_tests)

BOOST_AUTO_TEST_CASE(function_plot_test1)
{
  double samps_data[] = { 0.0, 1.0, 0.0, 1.0, 0.0 };
  size_t npts = sizeof(samps_data) / sizeof(double);
  srvf::Pointset samps(1, npts, samps_data);
  std::vector<double> params = srvf::util::linspace(0.0, 1.0, npts);
  srvf::Plf F(samps, params);

  srvf::plot::FunctionPlot plot;
  plot.insert(F, srvf::plot::Color(0.0, 0.0, 1.0));
  srvf::plot::FltkGlPlotWindow plotwin(800, 400, "plot_tests/function_plot_test1");
  plotwin.add_plot(static_cast<srvf::plot::Plot*>(&plot));
  plotwin.show();
  Fl::run();
}

BOOST_AUTO_TEST_SUITE_END()
