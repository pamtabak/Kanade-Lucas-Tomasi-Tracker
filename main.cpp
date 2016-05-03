#include "libs/CImg.h"
#include <iostream>
#include <string>
#include <vector>
#include "lucasKanade.hpp"
#include "matrix.hpp"

#define cimg_use_magick

using namespace cimg_library;

// g++ main.cpp -o main.out -L/opt/X11/lib -lX11 -pthread

int main(int argc, char * argv[]) 
{
	std::vector<std::vector<double> > matrix;
	std::vector<double> row1;
	std::vector<double> row2;
	std::vector<double> row3;

	row1.push_back(1);
	row1.push_back(2);
	row1.push_back(3);
	row2.push_back(4);
	row2.push_back(5);
	row2.push_back(6);
	row3.push_back(7);
	row3.push_back(8);
	row3.push_back(9);


	matrix.push_back(row1);
	matrix.push_back(row2);
	matrix.push_back(row3);

	Matrix m;
	std::vector<std::vector<double> > transposedMatrix = m.getTransposedMatrix(matrix);
	for (int i = 0; i < matrix.size(); i++)
	{
		for (int j = 0; j < matrix.size(); j++)
		{
			std::cout << transposedMatrix[i][j] << std::endl;
		}
	}
	// CImg<double> image1("images/input1.png");
	// CImg<double> image2("images/input2.png");

    // boost::numeric::ublas::matrix<double> m (3, 3);
    // for (unsigned i = 0; i < m.size1 (); ++ i)
    //     for (unsigned j = 0; j < m.size2 (); ++ j)
    //         m (i, j) = 3 * i + j;
    // std::cout << m << std::endl;

	return 0;
}