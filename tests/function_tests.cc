#include <cmath>
#include <boost/test/unit_test.hpp>

#include <srvf/matrix.h>
#include <srvf/pointset.h>
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/functions.h>
using srvf::functions::match_vertex_t;

#include <srvf/util.h>
#include <srvf/render.h>
#include <srvf/plotwin.h>
#include <FL/Fl.H>

#include <vector>
#include <deque>
#include <iostream>
#include <iterator>
#include <time.h>
#include <cstdlib>


#define MY_CHECK_CLOSE(a,b) \
 BOOST_CHECK_EQUAL(srvf::numeric::almost_equal((a),(b)), true)

#define MY_REQUIRE_CLOSE(a,b) \
 BOOST_REQUIRE_EQUAL(srvf::numeric::almost_equal((a),(b)), true)


BOOST_AUTO_TEST_SUITE(functions_tests)

BOOST_AUTO_TEST_CASE(edge_variation_test1)
{
  double samps_data[] = { 1.0, -1.0, 1.0, -1.0 };
  std::vector<double> params = srvf::util::linspace(0.0, 1.0, 5);
  size_t i1 = 1;  // params[1] is a peak
  size_t i2 = 4;  // params[4] is a valley
  double expected = -0.5;

  srvf::Pointset samps(1, 4, samps_data);
  srvf::Srvf Q(samps, params);

  double actual = srvf::functions::TestAccess::edge_variation(Q, i1, i2);
  BOOST_CHECK_SMALL(actual-expected,1e-6);
}

BOOST_AUTO_TEST_CASE(edge_score_test1)
{
  double samps1_data[] = { 1.0, -1.0, 1.0, -1.0 };
  double samps2_data[] = { -1.0, 1.0 };
  std::vector<double> params1 = srvf::util::linspace(0.0, 1.0, 5);
  std::vector<double> params2 = srvf::util::linspace(0.0, 1.0, 3);
  size_t sc=1, sr=0;
  size_t tc=4, tr=1;
  double expected = 0.5;

  srvf::Pointset samps1(1, 4, samps1_data);
  srvf::Pointset samps2(1, 2, samps2_data);
  srvf::Srvf Q1(samps1, params1);
  srvf::Srvf Q2(samps2, params2);

  double actual = srvf::functions::TestAccess::edge_score(
    Q1, Q2, sc, sr, tc, tr);
  BOOST_CHECK_SMALL(actual-expected,1e-6);
}

BOOST_AUTO_TEST_CASE(build_gamma_segment_test1)
{
  double samps1_data[] = { 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0 };
  double samps2_data[] = { 1.0, -1.0, 1.0 };
  double G1exp[] = { 0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875 };
  double G2exp[] = { 0.0, 0.166667, 0.166667, 0.333333, 0.666667, 
                     0.833333, 0.833333, 1.0 };
  size_t Glen = sizeof(G1exp) / sizeof(double);
  size_t ncp1 = sizeof(samps1_data) / sizeof(double) + 1;
  size_t ncp2 = sizeof(samps2_data) / sizeof(double) + 1;
  std::vector<double> params1 = srvf::util::linspace(0.0, 1.0, ncp1);
  std::vector<double> params2 = srvf::util::linspace(0.0, 1.0, ncp2);
  size_t sc=0, sr=0;
  size_t tc=7, tr=3;

  srvf::Pointset samps1(1, ncp1-1, samps1_data);
  srvf::Pointset samps2(1, ncp2-1, samps2_data);
  srvf::Srvf Q1(samps1, params1);
  srvf::Srvf Q2(samps2, params2);
  std::vector<double> G1samps(1,0.0);
  std::vector<double> G2samps(1,0.0);

  srvf::functions::TestAccess::build_gamma_segment(
    Q1, Q2, sc, sr, tc, tr, G1samps, G2samps);

  BOOST_REQUIRE_EQUAL(G1samps.size(), Glen);
  BOOST_REQUIRE_EQUAL(G2samps.size(), Glen);

  for (size_t i=0; i<Glen; ++i)
  {
    BOOST_CHECK_SMALL(G1exp[i] - G1samps[i], 1e-3);
    BOOST_CHECK_SMALL(G2exp[i] - G2samps[i], 1e-3);
  }
}


BOOST_AUTO_TEST_CASE(calculate_scores_test1)
{
  double samps1_data[] = { 0.0, 0.2, 0.0, 0.2, 0.0, 0.2 };
  double samps2_data[] = { 1.0/3.0, 0.0, 1.0/3.0, 0.0 };
  size_t nsamps1 = sizeof(samps1_data) / sizeof(double);
  size_t nsamps2 = sizeof(samps2_data) / sizeof(double);
  srvf::Plf F1(srvf::Pointset(1, nsamps1, samps1_data), 
               srvf::util::linspace(0.0, 1.0, nsamps1) );
  srvf::Plf F2(srvf::Pointset(1, nsamps2, samps2_data), 
               srvf::util::linspace(0.0, 1.0, nsamps2) );
  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);
  std::map<match_vertex_t,double> score;
  std::map<match_vertex_t,match_vertex_t> preds;
  match_vertex_t start_vertex(nsamps1,nsamps2);
  score[start_vertex] = 0.0;
  score[match_vertex_t(0,1)] = 0.0;
  score[match_vertex_t(1,0)] = 0.0;
  preds[match_vertex_t(0,1)] = start_vertex;
  preds[match_vertex_t(1,0)] = start_vertex;
  srvf::functions::TestAccess::calculate_scores(Q1,Q2,score,preds);

}


BOOST_AUTO_TEST_CASE(segment_translate_test1)
{
  double samps1_data[]={ 1.0, -1.0, 1.0, -1.0 };
  double samps2_data[]={ 1.0, -1.0, 1.0, -1.0, 1.0 };
  size_t nsamps1 = sizeof(samps1_data) / sizeof(double);
  size_t nsamps2 = sizeof(samps2_data) / sizeof(double);
  srvf::Srvf Q1(srvf::Pointset(1,nsamps1,samps1_data));
  srvf::Srvf Q2(srvf::Pointset(1,nsamps2,samps2_data));

  std::deque<match_vertex_t> path;
  path.push_back(match_vertex_t(0,0));
  path.push_back(match_vertex_t(3,1));
  path.push_back(match_vertex_t(4,4));

  // Cases for translation from Q1 to Q2
  size_t segs12[] = { 0, 0, 0, 0, 0, 0,
                      1, 1, 1, 1 };
  double ts12[]= { 0.0, 0.125, 0.25, 0.5, 0.625, 0.75, 
                   0.75, 10.0/12.0, 11.0/12.0, 1.0 };
  double exps12[] = { 0.0, 0.05, 0.1, 0.1, 0.15, 0.2, 
                      0.2, 1.0/3.0, 2.0/3.0, 0.8 };
  size_t ncases12 = sizeof(segs12) / sizeof(size_t);

  for (size_t i=0; i<ncases12; ++i)
  {
    double s = srvf::functions::TestAccess::segment_translate(
      Q1, Q2, path[segs12[i]].first, path[segs12[i]].second, 
      path[segs12[i]+1].first, path[segs12[i]+1].second, ts12[i]);
    MY_CHECK_CLOSE(s, exps12[i]);
  }
}


BOOST_AUTO_TEST_CASE(match_test1)
{
  unsigned int seed = (unsigned int)time(NULL);
  srandom(seed);
  std::vector<double> F1samps_data=srvf::util::random_vector(120, 0.5);
  std::vector<double> F2samps_data=srvf::util::random_vector(130, 0.5);
  srvf::Pointset F1samps(1,F1samps_data.size(),F1samps_data);
  srvf::Pointset F2samps(1,F2samps_data.size(),F2samps_data);
  std::vector<double> F1params = 
    srvf::util::linspace(0.0,1.0,F1samps_data.size());
  std::vector<double> F2params = 
    srvf::util::linspace(0.0,1.0,F2samps_data.size());
  srvf::Plf F1(F1samps,F1params);
  srvf::Plf F2(F2samps,F2params);
  F1.scale_to_unit_arc_length();
  F2.scale_to_unit_arc_length();

  srvf::plot::FltkGlPlotWindow plotwin_orig(600, 300, 
    "match_test1: original functions");
  srvf::plot::FunctionPlot plot_orig;
  plot_orig.insert(F1,srvf::plot::Color(0.0,0.0,1.0));
  plot_orig.insert(F2,srvf::plot::Color(1.0,0.0,0.0));
  plotwin_orig.add_plot(&plot_orig);
  plotwin_orig.show();

  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);
  Q1.scale_to_unit_norm();
  Q2.scale_to_unit_norm();

  std::cout << "match_test1: distance between original parametrizations = "
            << srvf::sphere_distance(Q1, Q2) << std::endl;

  Q1 = srvf::constant_speed_param(Q1);
  Q2 = srvf::constant_speed_param(Q2);
  F1 = srvf::srvf_to_plf(Q1);
  F2 = srvf::srvf_to_plf(Q2);

  std::cout << "match_test1: normalized F1 has " 
            << F1.ncp() << " changepoints"<< std::endl;
  std::cout << "match_test1: normalized F2 has " 
            << F2.ncp() << " changepoints"<< std::endl;

  std::vector<srvf::Plf> Gv = srvf::functions::optimal_reparam(Q1, Q2);
  srvf::Plf F1r = srvf::composition(F1, Gv[0]);
  srvf::Plf F2r = srvf::composition(F2, Gv[1]);
  srvf::Srvf Q1r = srvf::plf_to_srvf(F1r);
  srvf::Srvf Q2r = srvf::plf_to_srvf(F2r);
  std::cout << "match_test1:  shape distance = "
            << srvf::sphere_distance(Q1r, Q2r) << std::endl;

  //std::cout << "G1:  ";
  //for (size_t i=0; i<Gv[0].samps().npts(); std::cout << Gv[0].samps()[i++][0] << " ");
  //std::cout << std::endl;
  //for (size_t i=0; i<Gv[0].params().size(); std::cout << Gv[0].params()[i++] << " ");
  //std::cout << std::endl;
  //std::cout << "G2:  ";
  //for (size_t i=0; i<Gv[1].samps().npts(); std::cout << Gv[1].samps()[i++][0] << " ");
  //std::cout << std::endl;
  //for (size_t i=0; i<Gv[1].params().size(); std::cout << Gv[1].params()[i++] << " ");
  //std::cout << std::endl;

  srvf::plot::FltkGlPlotWindow plotwin(600, 300, 
    "match_test1: aligned functions");
  srvf::plot::FunctionPlot plot;
  plot.insert(F1r,srvf::plot::Color(0.0,0.0,1.0));
  plot.insert(F2r,srvf::plot::Color(1.0,0.0,0.0));
  plotwin.add_plot(&plot);
  plotwin.show();
  Fl::run();
}

BOOST_AUTO_TEST_CASE(groupwise_optimal_reparam_test1)
{
  double Msamps_data[] = { 1.0, -1.0, 1.0 };
  double Q1samps_data[] = { 1.0, -1.0 };
  double Q2samps_data[] = { -1.0, 1.0 };
  double Q3samps_data[] = { -1.0, 1.0 };
  double Mparams_data[] = { 0.0, 1.0/3.0, 2.0/3.0, 1.0 };
  double Q1params_data[] = { 0.0, 5.0/6.0, 1.0 };
  double Q2params_data[] = { 0.0, 1.0/6.0, 1.0 };
  double Q3params_data[] = { 0.0, 3.0/4.0, 1.0 };
  size_t Mncp = sizeof(Mparams_data) / sizeof(double);
  size_t Q1ncp = sizeof(Q1params_data) / sizeof(double);
  size_t Q2ncp = sizeof(Q2params_data) / sizeof(double);
  size_t Q3ncp = sizeof(Q3params_data) / sizeof(double);
  double Eexp[] = { 0.0, 0.0, 1.0/3.0, 2.0/3.0, 1.0, 1.0 };
  double G1exp[] = { 0.0, 0.0, 5.0/12.0, 5.0/12.0, 5.0/6.0, 1.0 };
  double G2exp[] = { 0.0, 1.0/6.0, 7/12.0, 7/12.0, 1.0, 1.0 };
  double G3exp[] = { 0.0, 0.0, 0.0, 3.0/4.0, 1.0, 1.0 };
  size_t gamma_ncp = sizeof(Eexp) / sizeof(double);

  std::vector<std::deque<srvf::functions::match_vertex_t> > paths(3);
  paths[0].push_back(srvf::functions::match_vertex_t(0,0));
  paths[0].push_back(srvf::functions::match_vertex_t(3,1));
  paths[1].push_back(srvf::functions::match_vertex_t(Mncp,Q2ncp));
  paths[1].push_back(srvf::functions::match_vertex_t(0,1));
  paths[1].push_back(srvf::functions::match_vertex_t(3,2));
  paths[2].push_back(srvf::functions::match_vertex_t(Mncp,Q3ncp));
  paths[2].push_back(srvf::functions::match_vertex_t(1,0));
  paths[2].push_back(srvf::functions::match_vertex_t(2,1));
  paths[2].push_back(srvf::functions::match_vertex_t(3,2));

  srvf::Srvf Mu(
    srvf::Pointset(1,Mncp-1,Msamps_data), 
    std::vector<double>(&(Mparams_data[0]), &(Mparams_data[Mncp])) );

  std::vector<srvf::Srvf> Qs;
  Qs.push_back( srvf::Srvf(
    srvf::Pointset(1,Q1ncp-1,Q1samps_data), 
    std::vector<double>(&(Q1params_data[0]), &(Q1params_data[Q1ncp]))) );
  Qs.push_back( srvf::Srvf(
    srvf::Pointset(1,Q2ncp-1,Q2samps_data), 
    std::vector<double>(&(Q2params_data[0]), &(Q2params_data[Q2ncp]))) );
  Qs.push_back( srvf::Srvf(
    srvf::Pointset(1,Q3ncp-1,Q3samps_data), 
    std::vector<double>(&(Q3params_data[0]), &(Q3params_data[Q3ncp]))) );

  std::vector<srvf::Plf> GE = 
    srvf::functions::groupwise_build_gammas(Mu, Qs, paths);

  for (size_t i=0; i<4; ++i)
  { BOOST_REQUIRE_EQUAL(GE[i].ncp(), gamma_ncp); }

  for (size_t i=0; i<gamma_ncp; ++i)
  {
    MY_CHECK_CLOSE(GE[0].samps()[i][0], G1exp[i]);
    MY_CHECK_CLOSE(GE[1].samps()[i][0], G2exp[i]);
    MY_CHECK_CLOSE(GE[2].samps()[i][0], G3exp[i]);
    MY_CHECK_CLOSE(GE[3].samps()[i][0], Eexp[i]);
  }
}

BOOST_AUTO_TEST_CASE(karcher_mean_test1)
{
  double Q1samps_data[] = { 1.0, -1.0 };
  double Q2samps_data[] = { -1.0, 1.0 };
  double Q3samps_data[] = { -1.0, 1.0 };
  double Q1params_data[] = { 0.0, 5.0/6.0, 1.0 };
  double Q2params_data[] = { 0.0, 1.0/6.0, 1.0 };
  double Q3params_data[] = { 0.0, 3.0/4.0, 1.0 };
  size_t Q1ncp = sizeof(Q1params_data) / sizeof(double);
  size_t Q2ncp = sizeof(Q2params_data) / sizeof(double);
  size_t Q3ncp = sizeof(Q3params_data) / sizeof(double);

  std::vector<srvf::Srvf> Qs;
  Qs.push_back( srvf::Srvf(
    srvf::Pointset(1,Q1ncp-1,Q1samps_data), 
    std::vector<double>(&(Q1params_data[0]), &(Q1params_data[Q1ncp]))) );
  Qs.push_back( srvf::Srvf(
    srvf::Pointset(1,Q2ncp-1,Q2samps_data), 
    std::vector<double>(&(Q2params_data[0]), &(Q2params_data[Q2ncp]))) );
  Qs.push_back( srvf::Srvf(
    srvf::Pointset(1,Q3ncp-1,Q3samps_data), 
    std::vector<double>(&(Q3params_data[0]), &(Q3params_data[Q3ncp]))) );

  srvf::Srvf Mu = srvf::functions::karcher_mean(Qs, 0.00001, 100);
  std::vector<srvf::Plf> GE=srvf::functions::groupwise_optimal_reparam(Mu,Qs);

  Mu = gamma_action(Mu, GE.back());
  for (size_t i=0; i<Qs.size(); ++i)
  {
    Qs[i] = gamma_action(Qs[i], GE[i]);
  }

  srvf::plot::FltkGlPlotWindow plotwin(600, 300, "karcher_mean_test1");
  srvf::plot::FunctionPlot plot;

  for (size_t i=0; i<Qs.size(); ++i)
  {
    plot.insert(srvf::srvf_to_plf(Qs[i]),srvf::plot::Color(0.0,0.0,1.0),2.0);
  }
  plot.insert(srvf::srvf_to_plf(Mu),srvf::plot::Color(1.0,0.0,0.0),2.0);

  plotwin.add_plot(&plot);
  plotwin.show();
  Fl::run();
}

BOOST_AUTO_TEST_SUITE_END()
