#include <boost/test/unit_test.hpp>

#include "interp.h"
#include "util.h"

BOOST_AUTO_TEST_SUITE(interp_tests)

BOOST_AUTO_TEST_CASE(lookup_test1)
{
  double table[]={ 0.0, 0.00001, 0.5, 0.55, 1.0 };
  double tv[]={ -5.0, -0.0001, 0.0, 0.00001, 0.49999, 0.5, 0.9999, 1.0, 23.0 };
  int expv[]={ -1, -1, 0, 1, 1, 2, 3, 4, 4 };
  int ntable=sizeof(table)/sizeof(double);
  int ncases=sizeof(expv)/sizeof(int);
  for (int i=0; i<ncases; ++i)
  {
    int idx=srvf::interp::lookup(table,ntable,tv[i]);
    BOOST_CHECK_EQUAL(idx,expv[i]);
  }
}

BOOST_AUTO_TEST_CASE(interp_const_test1)
{
  double samps_data[] = {0.25, -0.5, 1.0, -1.5};
  double tv_data[] =
  {
    0.0, 0.1,   0.24999, 0.25, 0.3, 0.499, 
    0.5, 0.749, 0.75,    0.8,  0.99, 1.0
  };
  double exp_data[] =
  {
    0.25, 0.25, 0.25, -0.5, -0.5, -0.5, 
    1.0,  1.0,  -1.5, -1.5, -1.5, -1.5
  };
  int ntv=sizeof(tv_data)/sizeof(double);

  srvf::Matrix samps(1,4,samps_data);
  srvf::Matrix params=srvf::util::linspace(0.0,1.0,5);
  srvf::Matrix tv(1,ntv,tv_data);
  srvf::Matrix result(1,ntv);

  srvf::interp::interp_const(samps,params,tv,result);
  for (int i=0; i<ntv; ++i)
  {
    BOOST_CHECK_CLOSE(result(i),exp_data[i],1e-9);
  }
}

BOOST_AUTO_TEST_SUITE_END()
