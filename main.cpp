#include "libs/CImg.h"
#include <iostream>
#include <string>
#include <vector>
// #include "lucasKanade.hpp"
// #include "harrisdetector.hpp"
#include "lib.hpp"

#define cimg_use_magick

using namespace cimg_library;

// g++ main.cpp -o main.out -I/usr/local/Cellar/eigen/3.2.4/include/eigen3/ -L/opt/X11/lib -lX11 -std=c++11 -w

int main(int argc, char * argv[]) 
{
	std::vector<CImg<double> > images;

	// for (int i = 1; i <= 2; i++)
	for (int i = 1; i <= 252; i++)
	{
		std::string fileName = "images/CarScale/";
		std::string numberOfFile = std::to_string(i);
		for (int l = numberOfFile.size(); l < 4; l++)
		{
			numberOfFile = "0" + numberOfFile;
		}
		fileName += numberOfFile + ".jpg";

		const char * c = fileName.c_str();
		CImg<double> image(c);
		images.push_back(image);
	}
	
	Tracking tracking;
	tracking.lucasKanade(images);

	// LucasKanade lucasKanade;
	// lucasKanade.pyramidAlgorithm(images);

	return 0;
}