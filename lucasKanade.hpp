#include "libs/CImg.h"
#include <iostream>
#include <vector>
#include <algorithm>    // std::min_element, std::max_element

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
		CImg<double> it = getIt(image1, image2);

		// (Ix,Iy)·(vx,vy)= −It
		// To solve this equation, we need to use neighboor pixels

		// Creating Matrix A (Ix,Iy) and Matrix b (It) for each pixel, and solving linear system
		// Neighboors: we are gonna get all pixels around the select pixel. And ignore borders
		CImg<double> minEigenValues(width,height,depth,channel,0);

		std::vector<std::vector<std::vector<CImg<double> >>> matrixes = getMatrixes(width, height, minEigenValues, ix, iy, it);
		std::vector<std::vector<CImg<double> >> allA = matrixes[0];
		std::vector<std::vector<CImg<double> >> allB = matrixes[1];
		minEigenValues = matrixes[2][0][0];

		const unsigned char white[] = { 255,255,255 };
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				if (minEigenValues(x,y) > 0.0)
				{
					// CHOSEN POINT! Calculating vector
					CImg<double> v = ((allA[x][y].get_transpose() * allA[x][y]).get_invert())*allA[x][y].get_transpose()*allB[x][y];
					image1.draw_line(x, y , x + (int) v(0,0), y + (int) v(0,1), white);
				}
			}
		}

		image1.display();
		// image1.save("NOME.png");
	}



	// Returns all A matrixes (for all pixels), all B matrixes and minEigenValue matrix (for each pixel,
	// the minimum eigen value, in order to choose which pixels we are going to track)
	std::vector<std::vector<std::vector<CImg<double> >>> getMatrixes(int width, int height, CImg<double> minEigenValues, CImg<double> ix, CImg<double> iy, CImg<double> it)
	{
		// Initializing return object
		std::vector<std::vector<std::vector<CImg<double> >>> matrixes;	

		std::vector<std::vector<CImg<double> >> allA; // all A matrixes
		std::vector<std::vector<CImg<double> >> allB; // all B matrixes
		std::vector<std::vector<CImg<double> >> modifiedMinEigenValue; // creating a vector just so we can return minEigenValue object

		std::vector<CImg<double> > fixingBorderIssue;
		fixingBorderIssue.push_back(minEigenValues);
		allA.push_back(fixingBorderIssue);
		allB.push_back(fixingBorderIssue);

		for (int x = 1; x < width - 1; x++)
		{
			std::vector<CImg<double> > rowA;
			std::vector<CImg<double> > rowB;

			rowA.push_back(minEigenValues);
			rowB.push_back(minEigenValues);

			for (int y = 1; y < height - 1; y++)
			{
				// Initializing A matrix	
				CImg<double> a(2, 9,depth,channel,initValue);
				CImg<double> b(1, 9,depth,channel,initValue);

				// We are going to apply Gaussian weights
				
				a(0,0) = ix(x-1,y-1)/16;
				a(0,1) = iy(x-1,y-1)/16;
				b(0,0) = it(x-1,y-1)/16;

				a(1,0) = 2 * ix(x,y-1)/16;
				a(1,1) = 2 * iy(x,y-1)/16;
				b(0,1) = 2 * it(x,y-1)/16;

				a(2,0) = ix(x+1,y-1)/16;
				a(2,1) = iy(x+1,y-1)/16;
				b(0,2) = it(x+1,y-1)/16;

				a(3,0) = 2 * ix(x-1,y)/16;
				a(3,1) = 2 * iy(x-1,y)/16;
				b(0,3) = 2 * it(x-1,y)/16;

				a(4,0) = 4 * ix(x,y)/16; // current pixel
				a(4,1) = 4 * iy(x,y)/16; // current pixel
				b(4,0) = 4 * it(x,y)/16; // current pixel

				a(5,0) = 2 * ix(x+1,y)/16;
				a(5,1) = 2 * iy(x+1,y)/16;
				b(0,5) = 2 * it(x+1,y)/16;

				a(6,0) = ix(x-1,y+1)/16;
				a(6,1) = iy(x-1,y+1)/16;
				b(0,6) = it(x-1,y+1)/16;

				a(7,0) = 2 * ix(x,y+1)/16;
				a(7,1) = 2 * iy(x,y+1)/16;
				b(0,7) = 2 * it(x,y+1)/16;

				a(8,0) = ix(x+1,y+1)/16;
				a(8,1) = iy(x+1,y+1)/16;
				b(0,8) = it(x+1,y+1)/16;

				rowA.push_back(a);
				rowB.push_back(b);

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
			}

			allA.push_back(rowA);
			allB.push_back(rowB);
		}

		double threshold = 0.1*maximumMinEigenValue;
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				// Choosing possible points to be tracked
				if (minEigenValues(x,y) <= threshold)
				{
					// Dont choose this point
					minEigenValues(x,y) = 0.0;
				}
			}
		}

		// Only staying with one chosen pixel per window
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				double myEigen[] = {minEigenValues(x-1,y-1), minEigenValues(x-1,y), minEigenValues(x-1,y+1), minEigenValues(x,y-1), minEigenValues(x,y), minEigenValues(x,y+1)
					, minEigenValues(x+1,y-1), minEigenValues(x+1,y), minEigenValues(x+1,y+1)};

				double * maxValue = std::max_element(myEigen, myEigen+9);
				
				if (*maxValue > minEigenValues(x-1,y-1))
					minEigenValues(x-1,y-1) = 0.0;
				if (*maxValue > minEigenValues(x-1,y))
					minEigenValues(x-1,y) = 0.0;
				if (*maxValue > minEigenValues(x-1,y+1))
					minEigenValues(x-1,y+1) = 0.0;
				if (*maxValue > minEigenValues(x,y-1))
					minEigenValues(x,y-1) = 0.0;
				if (*maxValue > minEigenValues(x,y))
					minEigenValues(x,y) = 0.0;
				if (*maxValue > minEigenValues(x,y+1))
					minEigenValues(x,y+1) = 0.0;
				if (*maxValue > minEigenValues(x+1,y-1))
					minEigenValues(x+1,y-1) = 0.0;
				if (*maxValue > minEigenValues(x+1,y))
					minEigenValues(x+1,y) = 0.0;
				if (*maxValue > minEigenValues(x+1,y+1))
					minEigenValues(x+1,y+1) = 0.0;
			}
		}

		std::vector<CImg<double> > minEigenValueVector;
		minEigenValueVector.push_back(minEigenValues);
		modifiedMinEigenValue.push_back(minEigenValueVector);

		matrixes.push_back(allA);
		matrixes.push_back(allB);
		matrixes.push_back(modifiedMinEigenValue);

		return matrixes;
	}



	// Calculate It between two images
	CImg<double> getIt (CImg<double> image1, CImg<double> image2)
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
				
				if ((x != 0) && (x != (width - 1)))
				{
					iy(x,y) = 0.5*(image(x+1,y) - image(x-1,y));
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
	double maximumMinEigenValue = 0.0;
};