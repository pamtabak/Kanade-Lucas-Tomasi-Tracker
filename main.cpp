#include "libs/CImg.h"
#include <iostream>
#include <string>
#include <vector>
#include "lucasKanade.hpp"

#define cimg_use_magick

using namespace cimg_library;

// g++ main.cpp -o main.out -L/opt/X11/lib -lX11 -pthread

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

int main(int argc, char * argv[]) 
{
	// CImg<double> image1("images/input1.png");
	// CImg<double> image2("images/input2.png");

    // boost::numeric::ublas::matrix<double> m (3, 3);
    // for (unsigned i = 0; i < m.size1 (); ++ i)
    //     for (unsigned j = 0; j < m.size2 (); ++ j)
    //         m (i, j) = 3 * i + j;
    // std::cout << m << std::endl;

	return 0;
}