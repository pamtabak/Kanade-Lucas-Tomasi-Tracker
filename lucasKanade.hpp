#include "libs/CImg.h"
#include <iostream>
#include <vector>
#include <algorithm>    // std::min_element, std::max_element
#include "gaussianpyramid.hpp"

using namespace cimg_library;

struct chosenPoint
{
	int x;
	int y;
	CImg<double> *v;
};

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
					image1.draw_line(x, y ,x + (int) v(0,0),y + (int) v(0,1), white);
				}
			}
		}

		image1.display();
		image1.save("NOME.png");
	}

	void pyramidAlgorithm (std::vector<CImg<double> > images)
	{
		std::vector<std::vector<CImg<double> >> pyramids = getGaussianPyramids(images);
		std::vector<CImg<double> > flowPyramid;

		// Calculating flow for highest level of the pyramid
		CImg<double> image1 = pyramids[0][pyramidSize - 1];
		CImg<double> image2 = pyramids[1][pyramidSize - 1];

		int height    = image1.height();
		int width     = image1.width();

		// Calculating Ix and Iy for image1
		std::vector<CImg<double> > derived = derive(image1);
		CImg<double> ix                    = derived[0];
		CImg<double> iy                    = derived[1];

		// Calculating It, difference between image 1 and image 2
		CImg<double> it = getIt(image1, image2);

		CImg<double> minEigenValues(width,height,depth,channel,0);

		std::vector<std::vector<std::vector<CImg<double> >>> matrixes = getMatrixes(width, height, minEigenValues, ix, iy, it);
		std::vector<std::vector<CImg<double> >> allA                  = matrixes[0];
		std::vector<std::vector<CImg<double> >> allB                  = matrixes[1];
		minEigenValues                                                = matrixes[2][0][0];

		const unsigned char white[] = { 255,255,255 };
		CImg<double> flowImage(width, height, depth, 3, 0); // rgb image.
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				if (minEigenValues(x,y) > 0.0)
				{
					CImg<double> v = ((allA[x][y].get_transpose() * allA[x][y]).get_invert())*allA[x][y].get_transpose()*allB[x][y];
					flowImage(x,y,0,0) = v(0,0); // r
					flowImage(x,y,0,2) = v(0,1); // b
					image1.draw_line(x, y ,x + (int) v(0,0),y + (int) v(0,1), white);
				}
			}
		}
		flowPyramid.push_back(flowImage);

		int iteration = 0;
		for (int p = pyramidSize - 2; p >= 0; p--)
		{
			// std::cout << p << std::endl;

			image1 = pyramids[0][p];
			image2 = pyramids[1][p];

			height    = image1.height();
			width     = image1.width();

			minEigenValues(width,height,depth,channel,0);
			flowImage.assign(width,height,depth,3);
			maximumMinEigenValue = 0.0;

			// Calculating Ix and Iy for image1
			derived = derive(image1);
			ix      = derived[0];
			iy      = derived[1];

			std::vector<std::vector<CImg<double> >> allA;
			std::vector<std::vector<CImg<double> >> allB;

			// only getting even x and y, so this point exists on lower level image
			int kkk = 0;
			for (int x = 0; x < width - 1; x++)
			{
				std::vector<CImg<double> > rowA;
				std::vector<CImg<double> > rowB;
				for (int y = 0; y < height - 1; y++)
				{
					CImg<double> a(2,9,depth,channel,initValue);
					CImg<double> b(1,9,depth,channel,initValue);
					if (x > 0 && y > 0 && (x % 2 == 0) && (y % 2 == 0))
					{
						// checking if flow is different than 0 in this position
						if (flowPyramid[iteration](x/2, y/2, 0,0) > 0.0 || flowPyramid[iteration](x/2, y/2, 0,2) > 0.0)
						{
							kkk++;
							// Calculating it
							it = getIt(ix, iy, flowPyramid[iteration](x/2, y/2, 0,0), flowPyramid[iteration](x/2, y/2, 0,2));
							a  = applyGaussianWeightsA(ix,iy,x,y);
							b  = applyGaussianWeightsB(it,x,y); 

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
						else
						{
							minEigenValues(x,y) = 0.0;
						}
					}

					rowA.push_back(a);
					rowB.push_back(b);
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
			// for (int x = 1; x < width - 1; x++)
			// {
			// 	for (int y = 1; y < height - 1; y++)
			// 	{
			// 		minEigenValues = reducingAmountOfPointsToTrack (minEigenValues, x, y);
			// 	}
			// }

			
			for (int x = 0; x < width - 1; x++)
			{
				for (int y = 0; y < height - 1; y++)
				{
					if (x > 0 && y > 0 && (x % 2 == 0) && (y % 2 == 0))
					{
						if (minEigenValues(x,y) > 0.0)
						{
							CImg<double> v    = ((allA[x][y].get_transpose() * allA[x][y]).get_invert())*allA[x][y].get_transpose()*allB[x][y];	
							if (!std::isnan(v(0,0)) && !std::isnan(v(0,1)))
							{
								std::cout << x << "," << y << std::endl;
								flowImage(x,y,0,0) = 2*(flowPyramid[iteration](x/2, y/2, 0,0)) + v(0,0); // r
								flowImage(x,y,0,1) = 0.0; // g
								// std::cout << 2*(flowPyramid[iteration](x/2, y/2, 0,2)) + v(0,1) << std::endl;
								flowImage(x,y,0,2) = 2*(flowPyramid[iteration](x/2, y/2, 0,2)) + v(0,1); // b
								// flowImage(x,y,0,2) = 0.0; 
								image1.draw_line(x, y ,x + (int) flowImage(x,y,0,0),y + (int) flowImage(x,y,0,2), white);
							}
						}
						else if (flowPyramid[iteration](x/2, y/2, 0,0) > 0.0 || flowPyramid[iteration](x/2, y/2, 0,2) > 0.0)
						{
							flowImage(x,y,0,0) = 2*(flowPyramid[iteration](x/2, y/2, 0,0)); // r
							flowImage(x,y,0,1) = 0.0; // g
							flowImage(x,y,0,2) = 2*(flowPyramid[iteration](x/2, y/2, 0,2)); // b
							image1.draw_line(x, y ,x + (int) flowImage(x,y,0,0),y + (int) flowImage(x,y,0,2), white);
						}
					}
					else if (minEigenValues(x,y) > 0.0)
					{
						CImg<double> v    = ((allA[x][y].get_transpose() * allA[x][y]).get_invert())*allA[x][y].get_transpose()*allB[x][y];	
						if (!std::isnan(v(0,0)) && !std::isnan(v(0,1)))
						{
							image1.draw_line(x, y ,x + (int) v(0,0),y + (int) v(0,1), white);
						}
					}
				}
			}

			flowPyramid.push_back(flowImage);
			iteration++;
			// std::cout << iteration << std::endl;
			// std::cout << kkk << std::endl;
			image1.display();
		}
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
				// CImg<double> a(2, 9,depth,channel,initValue);
				// CImg<double> b(1, 9,depth,channel,initValue);

				CImg<double> a = applyGaussianWeightsA(ix,iy,x,y);
				CImg<double> b = applyGaussianWeightsB(it,x,y);

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
				minEigenValues = reducingAmountOfPointsToTrack (minEigenValues, x, y);
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

	CImg<double> reducingAmountOfPointsToTrack (CImg<double> minEigenValues, int x, int y)
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

		return minEigenValues;
	}

	CImg<double> applyGaussianWeightsA (CImg<double> ix, CImg<double> iy, int x, int y)
	{
		CImg<double> a(2, 9,depth,channel,initValue);

		a(0,0) = ix(x-1,y-1)/16;
		a(0,1) = iy(x-1,y-1)/16;

		a(1,0) = 2 * ix(x,y-1)/16;
		a(1,1) = 2 * iy(x,y-1)/16;

		a(2,0) = ix(x+1,y-1)/16;
		a(2,1) = iy(x+1,y-1)/16;

		a(3,0) = 2 * ix(x-1,y)/16;
		a(3,1) = 2 * iy(x-1,y)/16;

		a(4,0) = 4 * ix(x,y)/16; // current pixel
		a(4,1) = 4 * iy(x,y)/16; // current pixel

		a(5,0) = 2 * ix(x+1,y)/16;
		a(5,1) = 2 * iy(x+1,y)/16;

		a(6,0) = ix(x-1,y+1)/16;
		a(6,1) = iy(x-1,y+1)/16;

		a(7,0) = 2 * ix(x,y+1)/16;
		a(7,1) = 2 * iy(x,y+1)/16;

		a(8,0) = ix(x+1,y+1)/16;
		a(8,1) = iy(x+1,y+1)/16;

		return a;
	}

	CImg<double> applyGaussianWeightsB (CImg<double> it, int x, int y)
	{
		CImg<double> b(1, 9,depth,channel,initValue);

		b(0,0) = it(x-1,y-1)/16;

		b(0,1) = 2 * it(x,y-1)/16;

		b(0,2) = it(x+1,y-1)/16;

		b(0,3) = 2 * it(x-1,y)/16;

		b(4,0) = 4 * it(x,y)/16; // current pixel

		b(0,5) = 2 * it(x+1,y)/16;

		b(0,6) = it(x-1,y+1)/16;

		b(0,7) = 2 * it(x,y+1)/16;

		b(0,8) = it(x+1,y+1)/16;

		return b;
	}

	std::vector<std::vector<CImg<double> >> getGaussianPyramids(std::vector<CImg<double>> images)
	{
		// Generating pyramid for each image.
		GaussianPyramid gPyramid;
		double filter[5] = {1.0/16, 4.0/16, 6.0/16, 4.0/16, 1.0/16};
		gPyramid.generateFilter(filter);

		// Generating pyramids
		std::vector<std::vector<CImg<double> >> pyramids;
		for (int i = 0; i < images.size(); i++)
		{
			std::vector<CImg<double> > gaussianPyramid;
			gaussianPyramid.push_back(images[i]);
			
			CImg<double> reducedImage  = images[i]; // original image
			for (int p = 1; p < pyramidSize; p++)
			{
				CImg<double> imageGenerated = gPyramid.reduce(reducedImage);
				gaussianPyramid.push_back(imageGenerated);
				reducedImage = imageGenerated;
			}

			pyramids.push_back(gaussianPyramid);
		}

		return pyramids;
	}

	// Calculate It between two images
	CImg<double> getIt (CImg<double> image1, CImg<double> image2)
	{
		// Setting image`s attributes
		int height    = image1.height();
		int width     = image1.width();

		// Initializing return object
		CImg<double> it(width, height, depth, channel, initValue);

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
	
	// Calculate It between two images, when working with pyramids
	CImg<double> getIt (CImg<double> image1, CImg<double> image2, double xFlow, double yFlow)
	{
		// Setting image`s attributes
		int height    = image1.height();
		int width     = image1.width();

		// Initializing return object
		CImg<double> it(width, height, depth, channel, initValue);

		// Iterating over each pixel
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				it(x,y) = image2((int)(x + xFlow), (int)(y + yFlow)) - image1(x,y);
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
	int depth                   = 1;
	int channel                 = 1;
	int initValue               = 0;
	int pyramidSize             = 2;
	double maximumMinEigenValue = 0.0;
};