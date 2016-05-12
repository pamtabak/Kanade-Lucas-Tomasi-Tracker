#include "libs/CImg.h"
#include <iostream>
#include <vector>

using namespace cimg_library;

class LucasKanade 
{
public:
	LucasKanade()
	{

	}

	~LucasKanade()
	{

	}

	void algorithm (std::vector<CImg<double> > images)
	{
		// I(x+∆x,y+∆y,t+∆t)                 = I(x,y,t) + (∂I/∂x)∆x + (∂I/∂y)∆y + (∂I/∂t)∆t

		// Image 1                           = Image 2, assuming intensity of the pixel doesnt change
		// I (x + ∆x, y + ∆y, t + ∆t )       = I (x, y, t)

		// (∂I/∂x)∆x + (∂I/∂y)∆y + (∂I/∂t)∆t = 0     ->     (Ix,Iy)·(vx,vy)=−It
		// Ix = ∂I/∂x ,Iy = ∂I/∂y ,It = ∂I/∂t


		// First, we are only working with 2 images
		CImg<double> image1 = images[0];
		CImg<double> image2 = images[1];

		// Calculating Ix and Iy for image1
		std::vector<CImg<double> > derived = derive(image1);
		CImg<double> ix = derived[0];
		CImg<double> iy = derived[1];

		// Calculating It, difference between image 1 and image 2
	}

	CImg<double> calculateIt (CImg<double> image1, CImg<double> image2)
	{
		// Setting image`s attributes
		int height    = image1.height();
		int width     = image1.width();
		int depth     = 1;
		int channel   = 1;
		int initValue = 0;

		// Initializing return object
		CImg<double> it(width, height,depth,channel,initValue);

		// Iterating over each pixel
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				it(x,y) = image2(x,y) - image1(x,y);
			}
		}

		return it;
	}
	
	// Calculate image`s derived
	std::vector<CImg<double> > derive (CImg<double> image)
	{
		// Initilazing return object		
		std::vector<CImg<double> > derived;

		// Setting image`s attributes
		int height    = image.height();
		int width     = image.width();
		int depth     = 1;
		int channel   = 1;
		int initValue = 0;

		CImg<double> ix(width, height,depth,channel,initValue);
		CImg<double> iy(width, height,depth,channel,initValue);
		
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				if ((y != 0) && (y != (height - 1)))
				{
					ix(x,y) = image(x,y+1) - image(x,y-1);	
				}
				else
				{
					// Avoiding dealing with border
					ix(x,y) = image(x,y);
				}
				
				if ((x != 0) && (x != (width - 1)))
				{
					iy(x,y) = image(x+1,y) - image(x-1,y);
				}
				else
				{
					// Avoiding dealing with border
					iy(x,y) = image(x,y);
				}
			}
		}

		derived.push_back(ix);
		derived.push_back(iy);

		return derived;
	}
	
// private:
};