#include <boost/test/unit_test.hpp>

#include "interp.h"

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

BOOST_AUTO_TEST_SUITE_END()
