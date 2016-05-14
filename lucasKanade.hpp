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

		int height    = image1.height();
		int width     = image1.width();

		// Calculating Ix and Iy for image1
		std::vector<CImg<double> > derived = derive(image1);
		CImg<double> ix = derived[0];
		CImg<double> iy = derived[1];

		// Calculating It, difference between image 1 and image 2
		CImg<double> it = calculateIt(image1, image2);

		// (Ix,Iy)·(vx,vy)= −It
		// To solve this equation, we need to use neighboor pixels

		// Creating Matrix A (Ix,Iy) and Matrix b (It) for each pixel, and solving linear system
		// Neighboors: we are gonna get all pixels around the select pixel. And ignore borders
		CImg<double> minEigenValues(width,height,depth,channel,0);
		double maximumMinEigenValue = 0.0;
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				// Initializing A matrix	
				CImg<double> a(2, 9,depth,channel,initValue);
				CImg<double> b(1, 9,depth,channel,initValue);

				// We are going to apply Gaussian weights
				a(0,0) = 1/16 * ix(x-1,y-1);
				a(0,1) = 1/16 * iy(x-1,y-1);
				b(0,0) = 1/16 * it(x-1,y-1);

				a(1,0) = 2/16 * ix(x,y-1);
				a(1,1) = 2/16 * iy(x,y-1);
				b(0,1) = 2/16 * it(x,y-1);

				a(2,0) = 1/16 * ix(x+1,y-1);
				a(2,1) = 1/16 * iy(x+1,y-1);
				b(0,2) = 1/16 * it(x+1,y-1);

				a(3,0) = 2/16 * ix(x-1,y);
				a(3,1) = 2/16 * iy(x-1,y);
				b(0,3) = 2/16 * it(x-1,y);

				a(4,0) = 4/16 * ix(x,y); // current pixel
				a(4,1) = 4/16 * iy(x,y); // current pixel
				b(4,0) = 4/16 * it(x,y); // current pixel

				a(5,0) = 2/16 * ix(x+1,y);
				a(5,1) = 2/16 * iy(x+1,y);
				b(0,5) = 2/16 * it(x+1,y);

				a(6,0) = 1/16 * ix(x-1,y+1);
				a(6,1) = 1/16 * iy(x-1,y+1);
				b(0,6) = 1/16 * it(x-1,y+1);

				a(7,0) = 2/16 * ix(x,y+1);
				a(7,1) = 2/16 * iy(x,y+1);
				b(0,7) = 2/16 * it(x,y+1);

				a(8,0) = 1/16 * ix(x+1,y+1);
				a(8,1) = 1/16 * iy(x+1,y+1);
				b(0,8) = 1/16 * it(x+1,y+1);

				// Solving linear system Av = b

				CImg<double> transposedATimesA = a.get_transpose() * a;
				// Calculating eigen values
				CImgList<double> eigen      = transposedATimesA.get_eigen();
				CImg<double> eigenValuesImg = eigen(0);
				double lambda0              = eigenValuesImg(0,0);
				double lambda1              = eigenValuesImg(0,1);

				if (lambda0 > 0 && lambda1 > 0)
				{
					// We are only choosing pixels where both eigen values are positive
					if (lambda0 > lambda1)
					{
						minEigenValues(x,y) = lambda1;
						if (lambda1 > maximumMinEigenValue)
						{
							maximumMinEigenValue = lambda1;
						}
					}
					else
					{
						minEigenValues(x,y) = lambda0;
						if (lambda0 > maximumMinEigenValue)
						{
							maximumMinEigenValue = lambda0;
						}
					}
				}


				CImg<double> v(1, 9, depth, channel, initValue);
				v = ((a.get_transpose() * a).get_invert())*a.get_transpose()*b;

			}
		}
	}






	CImg<double> calculateIt (CImg<double> image1, CImg<double> image2)
	{
		// Setting image`s attributes
		int height    = image1.height();
		int width     = image1.width();

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

		CImg<double> ix(width, height,depth,channel,initValue);
		CImg<double> iy(width, height,depth,channel,initValue);
		
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				if ((y != 0) && (y != (height - 1)))
				{
					ix(x,y) = 0.5*(image(x,y+1) - image(x,y-1));	
				}
				else
				{
					// Avoiding dealing with border
					ix(x,y) = image(x,y);
				}
				
				if ((x != 0) && (x != (width - 1)))
				{
					iy(x,y) = 0.5*(image(x+1,y) - image(x-1,y));
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
	
private:
	int depth     = 1;
	int channel   = 1;
	int initValue = 0;
};