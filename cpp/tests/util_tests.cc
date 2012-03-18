#include <boost/test/unit_test.hpp>

#include "util.h"
#include "matrix.h"

BOOST_AUTO_TEST_SUITE(util_tests)

BOOST_AUTO_TEST_CASE(linspace_test1)
{
  double expdata[]={ 0.0, 0.25, 0.5, 0.75, 1.0 };
  srvf::Matrix A=srvf::linspace(0.0,1.0,5);
  BOOST_CHECK_EQUAL(A.rows(),1);
  BOOST_CHECK_EQUAL(A.cols(),5);
  for (int i=0; i<5; ++i)
  {
    BOOST_CHECK_EQUAL(A(i),expdata[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
