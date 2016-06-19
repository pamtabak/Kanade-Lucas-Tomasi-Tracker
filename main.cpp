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
	// CImg<double> image1("images/Dumptruck/frame07.png");
	// CImg<double> image2("images/Dumptruck/frame08.png");
	// CImg<double> image3("images/Dumptruck/frame09.png");
	// CImg<double> image4("images/Dumptruck/frame10.png");
	// CImg<double> image5("images/Dumptruck/frame11.png");
	// CImg<double> image6("images/Dumptruck/frame12.png");
	// CImg<double> image7("images/Dumptruck/frame13.png");
	// CImg<double> image8("images/Dumptruck/frame14.png");
	// CImg<double> image1("images/input1.png");
	// CImg<double> image2("images/input2.png");
	// images.push_back(image1);
	// images.push_back(image2);
	// images.push_back(image3);
	// images.push_back(image4);
	// images.push_back(image5);
	// images.push_back(image6);
	// images.push_back(image7);
	// images.push_back(image8);

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