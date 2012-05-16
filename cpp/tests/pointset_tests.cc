#include <boost/test/unit_test.hpp>

#include "pointset.h"
#include "matrix.h"

BOOST_AUTO_TEST_SUITE(pointset_tests)

BOOST_AUTO_TEST_CASE(ctor_test1)
{
  double points_data[] = {
    0.0, 1.0, 2.0, 
    3.0, 4.0, 5.0
  };
  srvf::Matrix points(2,3,points_data);
  srvf::Pointset X1(points,srvf::Pointset::POINT_PER_ROW);
  srvf::Pointset X2(points,srvf::Pointset::POINT_PER_COLUMN);

  BOOST_REQUIRE_EQUAL(X1.dim(),3);
  BOOST_REQUIRE_EQUAL(X1.npts(),2);
  BOOST_REQUIRE_EQUAL(X2.dim(),2);
  BOOST_REQUIRE_EQUAL(X2.npts(),3);

  for (int i=0; i<2; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      BOOST_CHECK_CLOSE(X1(i,j),points_data[3*i+j],1e-9);
      BOOST_CHECK_CLOSE(X2(j,i),points_data[3*i+j],1e-9);
    }
  }
}

BOOST_AUTO_TEST_CASE(push_test1)
{
  std::vector<double> v(3,0.2);
  srvf::Pointset S;

  for (int i=1; i<=5; ++i)
  {
    S.push(v);
    BOOST_CHECK_EQUAL(S.dim(),3);
    BOOST_CHECK_EQUAL(S.npts(),i);
    for (size_t j=0; j<v.size(); ++j)
    {
      BOOST_CHECK_EQUAL(S(S.npts()-1,j), v[j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(distance_and_dot_product_test1)
{
  double points_data[] = 
  {
    0.0, 1.0, 2.0, 
    3.0, 4.0, 5.0
  };
  double dists1[2][2] = 
  {
    {0.0, 5.196152423},
    {5.196152423, 0.0}
  };
  double dists2[3][3] = 
  {
    {0.0, 1.414213562, 2.828427125},
    {1.414213562, 0.0, 1.414213562},
    {2.828427125, 1.414213562, 0.0}
  };
  double dots1[2][2] = 
  {
    {5.0, 14.0},
    {14.0, 50.0}
  };
  double dots2[3][3] = 
  {
    {9.0, 12.0, 15.0},
    {12.0, 17.0, 22.0},
    {15.0, 22.0, 29.0}
  };
  srvf::Matrix data(2,3,points_data);
  srvf::Pointset X1(data,srvf::Pointset::POINT_PER_ROW);
  srvf::Pointset X2(data,srvf::Pointset::POINT_PER_COLUMN);

  for (int i=0; i<2; ++i)
  {
    for (int j=0; j<2; ++j)
    {
      BOOST_CHECK_CLOSE(X1.distance(i,j),dists1[i][j],1e-4);
      BOOST_CHECK_CLOSE(X1.dot_product(i,j),dots1[i][j],1e-4);
    }
  }

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      BOOST_CHECK_CLOSE(X2.distance(i,j),dists2[i][j],1e-4);
      BOOST_CHECK_CLOSE(X2.dot_product(i,j),dots2[i][j],1e-4);
    }
  }
}


BOOST_AUTO_TEST_SUITE_END()
