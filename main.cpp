#include "libs/CImg.h"
#include <iostream>
#include <string>
#include <vector>
#include "lucasKanade.hpp"
#include "harrisdetector.hpp"

#define cimg_use_magick

using namespace cimg_library;

// g++ main.cpp -o main.out -L/opt/X11/lib -lX11 -pthread -std=c++11

int main(int argc, char * argv[]) 
{
	std::vector<CImg<double> > images;

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
	
	LucasKanade lucasKanade;
	
	lucasKanade.pyramidAlgorithm(images);

	return 0;
}