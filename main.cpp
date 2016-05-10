#include "libs/CImg.h"
#include <iostream>
#include <string>
#include <vector>
#include "lucasKanade.hpp"
#include "matrix.hpp"
#include "harrisdetector.hpp"

#define cimg_use_magick

using namespace cimg_library;

// g++ main.cpp -o main.out -L/opt/X11/lib -lX11 -pthread -std=c++11

int main(int argc, char * argv[]) 
{
	std::vector<CImg<double> > images;
	CImg<double> image1("images/input1.png");
	CImg<double> image2("images/input2.png");
	images.push_back(image1);
	images.push_back(image2);
	
	// CImg<double> image1(3,3,1,1,0);
	// image1(0,0) = 1.0;
	// image1(0,1) = 2.0;
	// image1(0,2) = 3.0;
	// image1(1,0) = 4.0;
	// image1(1,1) = 5.0;
	// image1(1,2) = 6.0;
	// image1(2,0) = 7.0;
	// image1(2,1) = 8.0;
	// image1(2,2) = 9.0;
	LucasKanade lucasKanade;
	lucasKanade.algorithm(images);
	// HarrisDetector harris;
	
	// harris.algorithm(image1);

    // boost::numeric::ublas::matrix<double> m (3, 3);
    // for (unsigned i = 0; i < m.size1 (); ++ i)
    //     for (unsigned j = 0; j < m.size2 (); ++ j)
    //         m (i, j) = 3 * i + j;
    // std::cout << m << std::endl;

	return 0;
}