#include <boost/test/unit_test.hpp>

#include "qmap.h"
#include "matrix.h"
#include "plf.h"
#include "srvf.h"
#include "util.h"

BOOST_AUTO_TEST_SUITE(qmap_tests)

BOOST_AUTO_TEST_CASE(plf_to_srvf_test1)
{
  double samps_data[]={0.0, 1.0/3.0, 0.0, 1.0/3.0};
  double exp_vals[]={1.0, -1.0, 1.0};
  srvf::Matrix params=srvf::util::linspace(0.0,1.0,4);
  srvf::Matrix samps(1,4,samps_data);
  srvf::Plf F(samps,params);
  srvf::Srvf Q=srvf::plf_to_srvf(F);

  BOOST_REQUIRE_EQUAL(Q.dim(),1);
  BOOST_REQUIRE_EQUAL(Q.nsamps(),4);

  // The Srvf should have same parameters as the Plf
  for (int i=0; i<4; ++i)
  {
    BOOST_CHECK_CLOSE(Q.params()(i),params(i),1e-9);
  }
  // Check samples
  for (int i=0; i<3; ++i)
  {
    BOOST_CHECK_CLOSE(Q.samps()(i),exp_vals[i],1e-9);
  }
}

BOOST_AUTO_TEST_CASE(plf_to_srvf_test2)
{
  double samps_data[]={
    0.0, 4.0, 4.0, 3.0, 2.0,
    0.0, 0.0, 2.0, 2.0, 1.0
  };
  double exp_vals[]={
    4.0, 0.0,      -2.0, -1.6817928,
    0.0, 2.828427,  0.0, -1.6817928
  };
  srvf::Matrix params=srvf::util::linspace(0.0,1.0,5);
  srvf::Matrix samps(2,5,samps_data);
  srvf::Plf F(samps,params);
  srvf::Srvf Q=srvf::plf_to_srvf(F);

  BOOST_REQUIRE_EQUAL(Q.dim(),2);
  BOOST_REQUIRE_EQUAL(Q.nsamps(),5);

  // The Srvf should have same parameters as the Plf
  for (int i=0; i<5; ++i)
  {
    BOOST_CHECK_CLOSE(Q.params()(i),params(i),1e-9);
  }
  // Check samples
  for (int i=0; i<8; ++i)
  {
    BOOST_CHECK_CLOSE(Q.samps()(i),exp_vals[i],1e-4);
  }
}

BOOST_AUTO_TEST_SUITE_END()
