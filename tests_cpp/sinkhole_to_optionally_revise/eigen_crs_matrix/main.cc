
#include <iostream>
#include <iomanip>
#include <vector>
#include "Eigen/Sparse"

template<class T>
typename T::InnerIterator row(T & M, int k)
{
  return typename T::InnerIterator(M,k);
}

int main(int argc, char *argv[])
{
  using mat_type = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;
  using Tr = Eigen::Triplet<double>;

  std::vector<Tr> v;
  v.push_back( Tr(0, 0, 1.) );
  v.push_back( Tr(0, 2, 2.) );
  v.push_back( Tr(2, 1, 3.) );
  v.push_back( Tr(2, 2, 4.) );
  v.push_back( Tr(2, 3, 5.) );
  v.push_back( Tr(3, 3, 6.) );
  v.push_back( Tr(3, 4, 7.) );
  v.push_back( Tr(3, 5, 8.) );
  v.push_back( Tr(5, 2, 9.) );
  v.push_back( Tr(5, 3, 10.) );
  v.push_back( Tr(5, 4, 11.) );
  v.push_back( Tr(8, 1, 12.) );
  v.push_back( Tr(8, 2, 13.) );
  v.push_back( Tr(8, 3, 14.) );
  v.push_back( Tr(8, 4, 15.) );
  v.push_back( Tr(10, 5, 16.) );
  v.push_back( Tr(10, 6, 17.) );

  mat_type M(11,7);
  M.setFromTriplets(v.begin(), v.end());
  M.makeCompressed();

  // for (int k=0; k<M.outerSize(); ++k)
  //   for (mat_type::InnerIterator it(M,k); it; ++it)
  //     {
  // 	std::cout << it.value() <<  " "
  // 		  << it.row() << " "
  // 		  << it.col() << std::endl;
  // 	// it.value();
  // 	// it.row();   // row index
  // 	// it.col();   // col index (here it is equal to k)
  // 	// it.index(); // inner index, here it is equal to it.row()
  //     }

  const auto nnz = M.nonZeros();
  std::cout << "innersize = " <<M.innerSize() << '\n';
  std::cout << "outersize = " <<M.outerSize() << '\n';

  // column indices of nonzeros
  std::cout << "\n inner \n";
  auto ipt = M.innerIndexPtr();
  for (int i=0; i<M.innerSize(); ++i){
    std::cout << ipt[i] << '\n';
  }

  std::cout << "\n outer \n";
  auto opt = M.outerIndexPtr();
  for (int i=0; i<M.outerSize(); ++i){
    std::cout << opt[i] << '\n';
  }

  std::cout << "\n values \n";
  auto values = M.valuePtr();
  for (int i=0; i<nnz; ++i){
    std::cout << values[i] << '\n';
  }

  auto it0 = row(M, 2);
  auto it1 = row(M, 3);
  for (; it0,it1; ++it0, ++it1){
    //it0.valueRef() *= 2;
    std::cout << it0.value() <<  " "
	      << it0.row() << " "
	      << it0.col() << std::endl;
    std::cout << it1.value() <<  " "
	      << it1.row() << " "
	      << it1.col() << std::endl;
  }

  // using mt = Eigen::Map< Eigen::Matrix<double, -1, 1> >;
  // int rowindex = 1;
  // mt row(&values[opt[rowindex]], 3);
  // std::cout << " row : " << row[0] << std::endl;
  // std::cout << " row : " << row[1] << std::endl;

  std::cout << M << '\n';

  std::puts("PASS");
  return 0;
}
