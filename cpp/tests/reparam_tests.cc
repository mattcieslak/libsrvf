#include <boost/test/unit_test.hpp>

#include "reparam.h"
#include "srvf.h"
#include "plf.h"
#include "pointset.h"

BOOST_AUTO_TEST_SUITE(reparam_tests)

BOOST_AUTO_TEST_CASE(match_cost_test1)
{
  double samps1_data[]={1.0, -1.0};
  double samps2_data[]={1.0,  1.0};
  srvf::Pointset samps1(1, 2, samps1_data);
  srvf::Pointset samps2(1, 2, samps2_data);
  srvf::Srvf Q1(samps1);
  srvf::Srvf Q2(samps2);
  double exp_costs[]=
  {0.08578643763, 1.5, 2.0, 1.5, 2.914213562};
  size_t idx[5][4]={
    {0, 1, 0, 2},
    {0, 2, 1, 2},
    {0, 2, 0, 2},
    {0, 2, 0, 1},
    {1, 2, 0, 2}
  };

  for (size_t i=0; i<5; ++i)
  {
    double c = match_cost(Q1, idx[i][0], idx[i][1], Q2, idx[i][2], idx[i][3]);
    BOOST_CHECK_CLOSE(c, exp_costs[i], 1e-4);
  }
}

BOOST_AUTO_TEST_CASE(match_cost_test2)
{
  double samps1_data[] =
  {
    1.0/3.0, -2.0, 1.0,
    2.0,     -1.0, -1.0
  };
  double samps2_data[] =
  {
    0.5, -0.5, 2.0,  0.9,
    0.2,  0.1, 0.3, -1.1
  };
  srvf::Pointset samps1(2, 3, samps1_data, srvf::Pointset::POINT_PER_COLUMN);
  srvf::Pointset samps2(2, 4, samps2_data, srvf::Pointset::POINT_PER_COLUMN);
  srvf::Srvf Q1(samps1);
  srvf::Srvf Q2(samps2);
  
  double c = srvf::match_cost(Q1, 1, 3, Q2, 1, 4);
  BOOST_CHECK_CLOSE(c, 3.1716, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
