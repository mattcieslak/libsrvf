#include <boost/test/unit_test.hpp>

#include "util.h"
#include "matrix.h"

BOOST_AUTO_TEST_SUITE(util_tests)

BOOST_AUTO_TEST_CASE(linspace_test1)
{
  double expdata[]={ 0.0, 0.25, 0.5, 0.75, 1.0 };
  srvf::Matrix A=srvf::util::linspace(0.0,1.0,5);
  BOOST_CHECK_EQUAL(A.rows(),1);
  BOOST_CHECK_EQUAL(A.cols(),5);
  for (int i=0; i<5; ++i)
  {
    BOOST_CHECK_EQUAL(A(i),expdata[i]);
  }
}

BOOST_AUTO_TEST_CASE(unique_test1)
{
  double v1_data[]={ 0.67, 0.23, 0.9, 0.42 };
  double v2_data[]={ 1.0, 0.0, 0.23, 0.89999, 0.4200001, 0.67  };
  double exp_data[]={ 0.0, 0.23, 0.42, 0.4200001, 0.67, 0.89999, 0.9, 1.0 };
  int nv1=sizeof(v1_data)/sizeof(double);
  int nv2=sizeof(v2_data)/sizeof(double);
  int nexp=sizeof(exp_data)/sizeof(double);
  srvf::Matrix v1(1,nv1,v1_data);
  srvf::Matrix v2(1,nv2,v2_data);
  srvf::Matrix v=srvf::util::unique(v1,v2);
  BOOST_REQUIRE_EQUAL(v.size(),nexp);
  for (int i=0; i<v.size(); ++i)
  {
    BOOST_CHECK_CLOSE(v(i),exp_data[i],1e-9);
  }
}

BOOST_AUTO_TEST_CASE(diff_test1)
{
  srvf::Matrix X(2,10,1.523);
  srvf::Matrix dX=srvf::util::diff(X);
  BOOST_REQUIRE_EQUAL(dX.rows(),2);
  BOOST_REQUIRE_EQUAL(dX.cols(),9);
  for (int i=0; i<dX.size(); ++i)
  {
    BOOST_CHECK_CLOSE(dX(i),0.0,1e-9);
  }
}

BOOST_AUTO_TEST_CASE(diff_test2)
{
  srvf::Matrix X(3,11);
  srvf::Matrix tv=srvf::util::linspace(0.0,1.0,11);
  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<11; ++j)
    {
      X(i,j)=(double)j;
    }
  }
  srvf::Matrix dX=srvf::util::diff(X,tv);
  BOOST_REQUIRE_EQUAL(dX.rows(),3);
  BOOST_REQUIRE_EQUAL(dX.cols(),10);
  for (int i=0; i<dX.size(); ++i)
  {
    BOOST_CHECK_CLOSE(dX(i),10.0,1e-6);
  }
}

BOOST_AUTO_TEST_SUITE_END()
