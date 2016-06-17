#include "libs/CImg.h"
#include <iostream>
#include <vector>
#include <algorithm>    // std::min_element, std::max_element
#include "gaussianpyramid.hpp"
#include "point.hpp"

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
		CImg<double> image1 = images[0];
		CImg<double> image2 = images[1];

		int height    = image1.height();
		int width     = image1.width();

		std::vector<CImg<double> > derived = derive(image1);
		CImg<double> ix = derived[0];
		CImg<double> iy = derived[1];

		CImg<double> it = getIt(image1, image2);

		CImg<double> minEigenValues(width,height,depth,channel,0);
		
		std::vector<std::vector<std::vector<CImg<double> >>> matrixes = getMatrixes(width, height, minEigenValues, ix, iy, it);
		std::vector<std::vector<CImg<double> >> allA = matrixes[0];
		std::vector<std::vector<CImg<double> >> allB = matrixes[1];
		minEigenValues = matrixes[2][0][0];

		const unsigned char white[] = { 255,255,255 };
		double xf = 0.0;
		double yf = 0.0;
		int times = 0;
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
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
	}

	void pyramidAlgorithm (std::vector<CImg<double> > images)
	{
		int numberOfFrames = images.size() - 1;
		int frame          = 0; // for iteration purposes
		// std::vector<ChosenPoint*> points;
		std::vector<ChosenPoint> points;

		CImg<double> image1 = images[0];
		CImg<double> image2 = images[1];

		std::vector<CImg<double> > imagesBeingUsed;
		imagesBeingUsed.push_back(image1);
		imagesBeingUsed.push_back(image2);
		std::vector<std::vector<CImg<double> >> pyramids = getGaussianPyramids(imagesBeingUsed);

		const int height    = image1.height();
		const int width     = image1.width();

		// Calculating Ix and Iy for image1
		std::vector<CImg<double> > derived = derive(image1);
		CImg<double> ix                    = derived[0];
		CImg<double> iy                    = derived[1];

		// Calculating It, difference between image 1 and image 2
		CImg<double> it = getIt(image1, image2); // We are not actually using this information at this point

		CImg<double> minEigenValues(width,height,depth,channel,0); // Matrix that helps us decide which points should be chosen


		std::cout << "oie" << std::endl;
		getMatrixes(width, height, minEigenValues, ix, iy, it);
		// std::vector<std::vector<std::vector<CImg<double> >>> matrixes = getMatrixes(width, height, &minEigenValues, ix, iy, it, &allA, &allB);
		// std::vector<std::vector<CImg<double> >> allA                  = matrixes[0];
		// std::vector<std::vector<CImg<double> >> allB                  = matrixes[1];
		// minEigenValues                                                = matrixes[2][0][0];
		std::cout << "oie2" << std::endl;

		const unsigned char white[] = { 255,255,255 };
		for (int y = 1; y < height - 1; y++)
		{
			for (int x = 1; x < width - 1; x++)
			{
				if (minEigenValues(x,y) > 0.0)
				{
					// choose this point
					ChosenPoint chosenP;
					chosenP.setNumberOfFrames(numberOfFrames);
					chosenP.setPoint(x, y);
					points.push_back(chosenP);

					// CImg<double> v = ((allA[x][y].get_transpose() * allA[x][y]).get_invert())*allA[x][y].get_transpose()*allB[x][y];
					// image1.draw_line(x, y ,x + (int) v(0,0),y + (int) v(0,1), white);
				}
			}
		}

		// Once the points are choosen from the original image, we build the pyramid
		for (int level = pyramidSize - 1; level >= 0; level--)
		{
			std::cout << level << std::endl;
			// STILL NEED TO DO FOR MORE THAN 2 IMAGES
			// STILL NEED TO DO SOMETHING TO CHANGE CHOSEN POINTS AFTER X FRAMES
			for (int p = 0; p < points.size(); p++)
			{
				int xOnLevel = points[p].getPoint().x / pow(2,level);
				int YOnLevel = points[p].getPoint().y / pow(2,level);

				if (xOnLevel > 0 && xOnLevel < pyramids[0][level].width() && YOnLevel > 0 && YOnLevel < pyramids[0][level].height())
				{
					if (level == pyramidSize - 1)
					{
						std::vector<CImg<double> > derived = derive(pyramids[0][level]);
						ix                    = derived[0];
						iy                    = derived[1];

						it = getIt(pyramids[0][level], pyramids[1][level], 2*points[p].getFlow()[frame].x, 2*points[p].getFlow()[frame].y);

						// Calculating matrix A and B, at this point
						CImg<double> a = applyGaussianWeightsA(ix, iy, xOnLevel, YOnLevel);
						CImg<double> b = applyGaussianWeightsB(it, xOnLevel ,YOnLevel);

						CImg<double> v = ((a.get_transpose() * a).get_invert())*a.get_transpose()*b;

						if (!std::isnan(v(0,0)) && !std::isnan(v(0,1)))
						{
							// just checking if everything went ok with all matrixes transformations
							points[p].setFlow(v(0,0), v(0,1), frame);
						}
					}
					else
					{
						points[p].updateFlow(frame);
					}
				}
				else
				{
					points[p].updateFlow(0.0, 0.0, frame);
				}
			}
		}

		double meanX = 0.0;
		double meanY = 0.0;
		for (int p = 0; p < points.size(); p++)
		{
			double finalX = points[p].getPoint().x + points[p].getFlow()[frame].x;
			double finalY = points[p].getPoint().y + points[p].getFlow()[frame].y;
			point finalPoint = bilinearInterpolation(finalX, finalY);
			// We need to  verify if pixel is outside the image

			image1.draw_line(points[p].getPoint().x, points[p].getPoint().y, (int) finalPoint.x, (int) finalPoint.y, white);

			meanX += points[p].getFlow()[frame].x;
			meanY += points[p].getFlow()[frame].y;
		}

		std::cout << meanX/points.size() << "," << meanY/points.size() << std::endl;

		// image1.save("testePiramide.png");
		image1.display();
	}

	CImg<double>* getMatrixes(const int width,const int height, CImg<double> minEigenValues, CImg<double> ix, CImg<double> iy, CImg<double> it)
	{
		CImg<double> matrixes[3];

		CImg<double> allA[width][height];
		CImg<double> allB[width][height];

		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				CImg<double> a = applyGaussianWeightsA(ix,iy,x,y);
				CImg<double> b = applyGaussianWeightsB(it,x,y);

				allA[x][y] = a;
				allB[x][y] = b;

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
		}

		// double threshold = 0.1*maximumMinEigenValue;
		// for (int x = 1; x < width - 1; x++)
		// {
		// 	for (int y = 1; y < height - 1; y++)
		// 	{
		// 		// Choosing possible points to be tracked
		// 		if (minEigenValues(x,y) <= threshold)
		// 		{
		// 			minEigenValues(x,y) = 0.0;
		// 		}
		// 	}
		// }

		// // // Only staying with one chosen pixel per window
		// for (int x = 1; x < width - 1; x++)
		// {
		// 	for (int y = 1; y < height - 1; y++)
		// 	{
		// 		minEigenValues = reducingAmountOfPointsToTrack (minEigenValues, x, y);
		// 	}
		// }

		matrixes[0] = allA;
		matrixes[1] = allB;
		matrixes[2] = minEigenValues;

		return matrixes;
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
				CImg<double> a = applyGaussianWeightsA(ix,iy,x,y);
				CImg<double> b = applyGaussianWeightsB(it,x,y);

				rowA.push_back(a);
				rowB.push_back(b);

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
					minEigenValues(x,y) = 0.0;
				}
			}
		}

		// // Only staying with one chosen pixel per window
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
		a(1,0) = iy(x-1,y-1)/16;

		a(0,1) = 2 * ix(x,y-1)/16;
		a(1,1) = 2 * iy(x,y-1)/16;

		a(0,2) = ix(x+1,y-1)/16;
		a(1,2) = iy(x+1,y-1)/16;

		a(0,3) = 2 * ix(x-1,y)/16;
		a(1,3) = 2 * iy(x-1,y)/16;

		a(0,4) = 4 * ix(x,y)/16; // current pixel
		a(1,4) = 4 * iy(x,y)/16; // current pixel

		a(0,5) = 2 * ix(x+1,y)/16;
		a(1,5) = 2 * iy(x+1,y)/16;

		a(0,6) = ix(x-1,y+1)/16;
		a(1,6) = iy(x-1,y+1)/16;

		a(0,7) = 2 * ix(x,y+1)/16;
		a(1,7) = 2 * iy(x,y+1)/16;

		a(0,8) = ix(x+1,y+1)/16;
		a(1,8) = iy(x+1,y+1)/16;

		return a;
	}

	CImg<double> applyGaussianWeightsB (CImg<double> it, int x, int y)
	{
		CImg<double> b(1, 9,depth,channel,initValue);

		b(0,0) = it(x-1,y-1)/16;

		b(0,1) = 2 * it(x,y-1)/16;

		b(0,2) = it(x+1,y-1)/16;

		b(0,3) = 2 * it(x-1,y)/16;

		b(0,4) = 4 * it(x,y)/16; // current pixel

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
		double newX, newY;
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				newX = x + xFlow;
				if (newX < 0)
					newX = x;
				if (newX >= width)
					newX = x;
				

				newY = y + yFlow;
				if (newY < 0)
					newY = y;
				if (newY >= height)
					newY = y;

				point newPoint = bilinearInterpolation(newX, newY);
				if (newPoint.y == height)
					newPoint.y--;
				if (newPoint.x == width)
					newPoint.x--;
				it(x,y) = image2(newPoint.x, newPoint.y) - image1(x,y);
			}
		}

		return it;
	}

	point bilinearInterpolation(double x, double y)
	{
		point a,b,c,d;
		a.x = floor(x);
		a.y = floor(y);

		b.x = floor(x);
		b.y = ceil(y);

		c.x = ceil(x);
		c.y = floor(y);

		d.x = ceil(x);
		d.y = ceil(y);

		double distanceA = sqrt(pow((x - a.x), 2) + pow ((y - a.y), 2));
		double distanceB = sqrt(pow((x - b.x), 2) + pow ((y - b.y), 2));
		double distanceC = sqrt(pow((x - c.x), 2) + pow ((y - c.y), 2));
		double distanceD = sqrt(pow((x - d.x), 2) + pow ((y - d.y), 2));

		if (distanceA <= distanceB && distanceA <= distanceC && distanceA <= distanceD)
			return a;
		if (distanceB <= distanceA && distanceB <= distanceC && distanceB <= distanceD)
			return b;
		if (distanceC <= distanceA && distanceC <= distanceB && distanceC <= distanceD)
			return c;
		if (distanceD <= distanceA && distanceD <= distanceB && distanceD <= distanceC)
			return d;
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
		
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				ix (x,y,0) = 0.5 * (image(x+1, y) - image(x-1, y));
				iy (x,y,0) = 0.5 * (image(x, y+1) - image(x, y-1));
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
	int pyramidSize             = 1;
	double maximumMinEigenValue = 0.0;
};