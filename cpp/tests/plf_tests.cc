#include <boost/test/unit_test.hpp>

#include "plf.h"
#include "matrix.h"
#include "util.h"

BOOST_AUTO_TEST_SUITE(plf_tests)

BOOST_AUTO_TEST_CASE(evaluate_test1)
{
  srvf::Matrix A=srvf::util::linspace(0.0,1.0,2);
  srvf::Plf F(A,A);
  srvf::Matrix tv=srvf::util::linspace(0.0,1.0,5);
  srvf::Matrix Ftv1(1,5);
  srvf::Matrix Ftv2(1,1);
  double expdata[]={0.0, 0.25, 0.5, 0.75, 1.0};
  F.evaluate(tv,Ftv1);
  for (int i=0; i<5; ++i)
  {
    F.evaluate(tv(i),Ftv2);
    BOOST_CHECK_EQUAL(Ftv1(i),expdata[i]);
    BOOST_CHECK_EQUAL(Ftv2(0),expdata[i]);
  }
}

BOOST_AUTO_TEST_CASE(evaluate_test2)
{
  srvf::Matrix uv=srvf::util::linspace(0.0,1.0,500);
  srvf::Matrix uvl=srvf::util::linspace(0.0,1.0,999);
  srvf::Matrix X(3,500);
  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<500; ++j)
    {
      X(i,j)=(double)j;
    }
  }
  srvf::Plf F(X,uv);
  srvf::Matrix res(3,999);
  F.evaluate(uvl,res);
  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<999; ++j)
    {
      double ev=(double)(j) * 499.0 / 998.0;
      BOOST_CHECK_CLOSE(ev,res(i,j),1e-5);
    }
  }
}

BOOST_AUTO_TEST_CASE(preimages_test1)
{
  srvf::Matrix samps=srvf::util::linspace(0.0,1.0,2);
  srvf::Matrix uv=srvf::util::linspace(0.0,1.0,13);
  srvf::Matrix Fiuv(1,uv.size());
  srvf::Plf F(samps);
  F.preimages(uv,Fiuv);
  for (int i=0; i<uv.size(); ++i)
  {
    BOOST_CHECK_CLOSE(uv(i),Fiuv(i),1e-5);
  }
}

BOOST_AUTO_TEST_CASE(preimages_test2)
{
  double samps_data[]={0.0, 0.5, 0.5, 1.0};
  srvf::Matrix samps(1,4,samps_data);
  srvf::Matrix tv=srvf::util::linspace(0.0,1.0,5);
  srvf::Matrix Fitv(1,tv.size());
  double exp_data[]={0.0, 1.0/6.0, 2.0/3.0, 5.0/6.0, 1.0};
  srvf::Plf F(samps);
  F.preimages(tv,Fitv);
  for (int i=0; i<tv.size(); ++i)
  {
    BOOST_CHECK_CLOSE(exp_data[i],Fitv(i),1e-5);
  }
}

BOOST_AUTO_TEST_CASE(arc_length_test1)
{
  double samps_data[]=
  {
    0.0, 1.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 1.0
  };
  srvf::Plf F(srvf::Matrix(2,4,samps_data));
  double ev=3.0;
  double av=F.arc_length();
  BOOST_CHECK_CLOSE(ev,av,1e-5);
}

BOOST_AUTO_TEST_CASE(translate_test1)
{
  srvf::Matrix samps(3,20,0.0);
  srvf::Plf F(samps);
  double v_data[]={1.0,2.0,3.0};
  srvf::Matrix v(3,1,v_data);
  F.translate(v);
  for (int i=0; i<F.samps().rows(); ++i)
  {
    for (int j=0; j<F.samps().cols(); ++j)
    {
      BOOST_CHECK_EQUAL(F.samps()(i,j),v(i));
    }
  }
}

BOOST_AUTO_TEST_CASE(rotate_test1)
{
  // Rotation by 90 degrees counter-clockwise
  double R_data[]={
    0.0, -1.0,
    1.0, 0.0
  };
  srvf::Matrix R(2,2,R_data);
  srvf::Matrix samps(2,3,1.0);
  srvf::Plf F(samps);
  F.rotate(R);
  for (int i=0; i<F.samps().cols(); ++i)
  {
    BOOST_CHECK_EQUAL(F.samps()(0,i),-1.0);
    BOOST_CHECK_EQUAL(F.samps()(1,i),1.0);
  }
}

BOOST_AUTO_TEST_CASE(scale_test1)
{
  srvf::Matrix samps(3,57,1.0);
  srvf::Plf F(samps);
  F.scale(2.0);
  for (int i=0; i<F.samps().size(); ++i)
  {
    BOOST_CHECK_EQUAL(F.samps()(i),2.0);
  }
}

BOOST_AUTO_TEST_CASE(linear_combination_test1)
{
  double params1_data[]= { 0.0, 0.1, 0.2, 0.9, 1.0 };
  double params2_data[]= { 0.0, 0.3, 0.4, 0.6, 1.0 };
  srvf::Matrix samps1(1,5,1.0);
  srvf::Matrix samps2(1,5,-1.0);
  srvf::Matrix params1(1,5,params1_data);
  srvf::Matrix params2(1,5,params2_data);
  double exp_params_data[]={ 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.9, 1.0 };
  srvf::Plf F1(samps1,params1);
  srvf::Plf F2(samps2,params2);
  srvf::Plf Fr=srvf::linear_combination(F1,F2,0.75,0.25);

  BOOST_REQUIRE_EQUAL(Fr.samps().size(),8);
  BOOST_REQUIRE_EQUAL(Fr.nsamps(),8);
  for (int i=0; i<8; ++i)
  {
    BOOST_CHECK_EQUAL(Fr.samps()(i),0.5);
    BOOST_CHECK_EQUAL(Fr.params()(i),exp_params_data[i]);
  }
}

BOOST_AUTO_TEST_CASE(composition_test1)
{
  double params1_data[]={0.0, 0.25, 0.5, 1.0};
  double params2_data[]={0.0, 1.0/3.0, 2.0/3.0, 1.0};
  double exp_params[]={0.0, 1.0/6.0, 1.0/3.0, 2.0/3.0, 1.0};
  double samps1_data[]={0.0, 0.5, 0.0, 0.5};
  double samps2_data[]={0.0, 0.5, 0.5, 1.0};
  double exp_samps[]={0.0, 0.5, 0.0, 0.0, 0.5};
  int n1=sizeof(params1_data)/sizeof(double);
  int n2=sizeof(params2_data)/sizeof(double);
  int n3=sizeof(exp_params)/sizeof(double);
  srvf::Matrix params1(1,n1,params1_data);
  srvf::Matrix params2(1,n2,params2_data);
  srvf::Matrix samps1(1,n1,samps1_data);
  srvf::Matrix samps2(1,n2,samps2_data);
  srvf::Plf F1(samps1,params1);
  srvf::Plf F2(samps2,params2);
  srvf::Plf F12=srvf::composition(F1,F2);
  
  BOOST_REQUIRE_EQUAL(F12.nsamps(),n3);
  for (int i=0; i<n3; ++i)
  {
    BOOST_CHECK_EQUAL(F12.params()(i),exp_params[i]);
    BOOST_CHECK_EQUAL(F12.samps()(i),exp_samps[i]);
  }
}

BOOST_AUTO_TEST_CASE(inverse_test1)
{
  double params_data[]={0.0, 0.25, 0.5, 0.75, 1.0};
  double samps_data[]={0.0, 1.0/3.0, 1.0/3.0, 2.0/3.0, 1.0};
  int n1=sizeof(params_data)/sizeof(double);
  int nexp=n1;
  srvf::Matrix params(1,n1,params_data);
  srvf::Matrix samps(1,n1,samps_data);
  srvf::Plf F(samps,params);
  srvf::Plf Fi=srvf::inverse(F);
  
  BOOST_REQUIRE_EQUAL(Fi.nsamps(),nexp);
  for (int i=0; i<nexp; ++i)
  {
    BOOST_CHECK_EQUAL(Fi.params()(i),samps_data[i]);
    BOOST_CHECK_EQUAL(Fi.samps()(i),params_data[i]);
  }
}


BOOST_AUTO_TEST_SUITE_END()
