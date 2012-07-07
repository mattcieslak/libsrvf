#include <boost/test/unit_test.hpp>
#include <vector>
#include <fstream>

#include <fltk/run.h>
#include "plotwin.h"
#include "plf.h"
#include "srvf.h"
#include "qmap.h"
#include "rotate.h"
#include "reparam.h"
#include "fileio.h"

BOOST_AUTO_TEST_SUITE(plotwin_tests)

BOOST_AUTO_TEST_CASE(superimposed_test1)
{
  std::ifstream is1("data/horse-1.csv");
  std::ifstream is2("data/horse-2.csv");
  std::vector<srvf::Matrix> samps1_data = srvf::io::load_csv(is1);
  std::vector<srvf::Matrix> samps2_data = srvf::io::load_csv(is2);
  srvf::Pointset samps1(samps1_data[0],srvf::Pointset::POINT_PER_COLUMN);
  srvf::Pointset samps2(samps2_data[0],srvf::Pointset::POINT_PER_COLUMN);

  srvf::Plf F1(samps1);
  F1.scale_to_unit_arc_length();
  F1.translate_to_origin();

  srvf::Plf F2(samps2);
  F2.scale_to_unit_arc_length();
  F2.translate_to_origin();

  srvf::Srvf Q1=srvf::plf_to_srvf(F1);
  srvf::Srvf Q2=srvf::plf_to_srvf(F2);
  srvf::Matrix R=srvf::optimal_rotation(Q1, Q2);
  Q2.rotate(R);
  F2.rotate(R);
  srvf::Plf G=srvf::opencurves::optimal_reparam(Q1, Q2);
  srvf::Plf F2r=composition(F2, G);
  srvf::Srvf Q2r=gamma_action(Q2, G);

  srvf::SuperimposedPlot plot1;
  plot1.insert(F1, srvf::Color(0.0f,0.0f,1.0f));
  plot1.insert(F2r, srvf::Color(1.0f,0.0f,0.0f));

  srvf::FltkGlPlotWindow win(400, 400, "SuperimposedPlot Test 1");
  win.add_plot(&plot1);
  win.show();
  fltk::run();
}

BOOST_AUTO_TEST_SUITE_END()
