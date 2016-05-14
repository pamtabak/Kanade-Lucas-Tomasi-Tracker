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
	
	LucasKanade lucasKanade;
	lucasKanade.algorithm(images);
	// HarrisDetector harris;
	
	// harris.algorithm(image1);

	return 0;
}