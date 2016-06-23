#include "libs/CImg.h"
#include <iostream>
#include <vector>
#include <algorithm>    // std::min_element, std::max_element
#include "gaussianpyramid.hpp"
#include "point.hpp"

using namespace cimg_library;

typedef struct matrix
{
	double **m;
} matrix;

typedef struct pointInArray
{
	bool inArray;
	int position;
} pointInArray;


class LucasKanade 
{
public:
	LucasKanade()
	{
		
	}

	~LucasKanade()
	{
		delete allA;
		delete allB;
	}

	void initArrays(int width, int height)
	{
		allA = new matrix*[width];
		allB = new matrix*[width];

		for (int i = 0; i < width; i++)
		{
			allA[i] = new matrix[height];
			allB[i] = new matrix[height];
		}
	}

	// void algorithm (std::vector<CImg<double> > images)
	// {
	// 	CImg<double> image1 = images[0];
	// 	CImg<double> image2 = images[1];

	// 	int height    = image1.height();
	// 	int width     = image1.width();

	// 	std::vector<CImg<double> > derived = derive(image1);
	// 	CImg<double> ix = derived[0];
	// 	CImg<double> iy = derived[1];

	// 	CImg<double> it = getIt(image1, image2);

	// 	minEigenValues.assign(width,height,depth,channel,0);
		
	// 	initArrays(width, height);

	// 	getMatrixes(width, height, ix, iy, it);

	// 	const unsigned char white[] = { 255,255,255 };
	// 	double xf = 0.0;
	// 	double yf = 0.0;
	// 	int times = 0;
	// 	for (int x = 1; x < width - 1; x++)
	// 	{
	// 		for (int y = 1; y < height - 1; y++)
	// 		{
	// 			if (minEigenValues(x,y) > 0.0)
	// 			{
	// 				// CHOSEN POINT! Calculating vector
	// 				CImg<double> v = ((allA[x][y].get_transpose() * allA[x][y]).get_invert())*allA[x][y].get_transpose()*allB[x][y];
	// 				image1.draw_line(x, y ,x + (int) v(0,0),y + (int) v(0,1), white);
	// 			}
	// 		}
	// 	}

	// 	image1.display();
	// }

	pointInArray pointInArrayFunction(std::vector<ChosenPoint> points, double pointX, double pointY)
	{
		pointInArray p;
		for (int i = 0; i < points.size(); i++)
		{
			if (points[i].getPoint().x == pointX && points[i].getPoint().y == pointY)
			{
				p.inArray = true;
				p.position = i;
				return p;
			}
		}

		p.inArray = false;
		return p;
	}

	void pyramidAlgorithm (std::vector<CImg<double> > images)
	{
		int numberOfFrames = images.size() - 1;
		// std::vector<ChosenPoint> points;
		for (int frame = 0; frame < numberOfFrames; frame++)
		{
			CImg<double> image1 = images[frame];
			CImg<double> image2 = images[frame + 1];

			std::vector<CImg<double> > imagesBeingUsed;
			imagesBeingUsed.push_back(image1);
			imagesBeingUsed.push_back(image2);
			std::vector<std::vector<CImg<double> >> pyramids = getGaussianPyramids(imagesBeingUsed);

			const unsigned char white[] = { 255,255,255 };
			const unsigned char pink [] = { 255,0,255 };

			const int height    = image1.height();
			const int width     = image1.width();

			matrix ix;
			matrix iy;
			matrix it;
			if (frame % 10 == 0)
			{
				// points.clear();

				// we are only recalculating points to track in this case
				// Calculating Ix and Iy for image1
				std::vector<matrix> derived = derive(image1);
				ix                    = derived[0];
				iy                    = derived[1];

				// Calculating It, difference between image 1 and image 2
				it = getIt(image1, image2); // We are not actually using this information at this point

				// minEigenValues.assign(width,height,depth,channel,0); // Matrix that helps us decide which points should be chosen
				minEigenValues.m = new double*[width];
				for (int i = 0; i < width; i++)
				{
					minEigenValues.m[i] = new double[height];
					for (int j = 0; j < height; j++)
					{
						minEigenValues.m[i][j] = 0.0;
					}
				}

				initArrays(width, height);

				getMatrixes(width, height, ix, iy, it);

				for (int y = 1; y < height - 1; y++)
				{
					for (int x = 1; x < width - 1; x++)
					{
						if (minEigenValues.m[x][y] > 0.0)
						{
							pointInArray pointIn = pointInArrayFunction(points, (double) x, (double) y);
							if (!pointIn.inArray)
							{
								// choose this point
								ChosenPoint chosenP;
								chosenP.setNumberOfFrames(numberOfFrames);
								chosenP.setPoint(x, y);
								points.push_back(chosenP);
							}
						}
					}
				}

				std::cout << points.size() << std::endl;
			}

			// Once the points are choosen from the original image, we build the pyramid
			for (int level = pyramidSize - 1; level >= 0; level--)
			{
				std::cout << level << std::endl;
				for (int p = 0; p < points.size(); p++)
				{
					int xOnLevel = points[p].getPoint().x / pow(2,level);
					int YOnLevel = points[p].getPoint().y / pow(2,level);

					if (xOnLevel > 0 && xOnLevel < pyramids[0][level].width() - 1 && YOnLevel > 0 && YOnLevel < pyramids[0][level].height() - 1)
					{
						if (level == pyramidSize - 1)
						{
							std::vector<matrix> derived = derive(pyramids[0][level]);
							
							ix                    = derived[0];
							iy                    = derived[1];

							it = getIt(pyramids[0][level], pyramids[1][level], 2*points[p].getFlow()[frame].x, 2*points[p].getFlow()[frame].y);

							// Calculating matrix A and B, at this point
							matrix a = applyGaussianWeightsA(ix, iy, xOnLevel, YOnLevel);
							matrix b = applyGaussianWeightsB(it, xOnLevel ,YOnLevel);

							matrix v = calculateFlow(a,b);
							if (!std::isnan(v.m[0][0]) && !std::isnan(v.m[1][0]))
							{
								// just checking if everything went ok with all matrixes transformations
								points[p].setFlow(v.m[0][0], v.m[1][0], frame);
							}
							else
							{
								points.erase(points.begin() + p - 1);
							}

							delete v.m;
							delete a.m;
							delete b.m;
						}
						else
						{
							points[p].updateFlow(frame);
						}
					}
					else
					{
						// points[p].updateFlow(0.0, 0.0, frame);
						points.erase(points.begin() + p - 1);
					}
				}
			}

			for (int p = 0; p < points.size(); p++)
			{
				if (points[p].getFlow()[frame].x > 0.0 || points[p].getFlow()[frame].y > 0.0)
				{
					double finalX = points[p].getPoint().x - points[p].getFlow()[frame].x;
					double finalY = points[p].getPoint().y - points[p].getFlow()[frame].y;
					point finalPoint = bilinearInterpolation(finalX, finalY);


					double initFlowX = points[p].getPoint().x;
					double initFlowY = points[p].getPoint().y;
					double lastFlowX = finalPoint.x;
					double lastFlowY = finalPoint.y;
					
					if (frame >= 10)
					{
						for (int f = frame - 10; f < frame; f++)
						{
							image1.draw_line((int) initFlowX, (int) initFlowY, (int) lastFlowX, (int) lastFlowY, pink);
							initFlowX = lastFlowX;
							initFlowY = lastFlowY;
							lastFlowX = lastFlowX - points[p].getFlow()[f].x;
							lastFlowY = lastFlowY - points[p].getFlow()[f].y;
						}
					}
					else
					{
						image1.draw_line((int) initFlowX, (int) initFlowY, (int) lastFlowX, (int) lastFlowY, pink);
					}
					
					points[p].setPoint(finalPoint.x, finalPoint.y);
				}
			}
			image1.save("images/output/Segments/piramide.png", frame);

			// delete minEigenValues.m;
			// delete ix.m;
			// delete iy.m;
			// delete it.m;
		}

		images[numberOfFrames].save("images/output/Segments/piramide.png", numberOfFrames);
	}

	void getMatrixes(const int width,const int height, matrix ix, matrix iy, matrix it)
	{
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				matrix a = applyGaussianWeightsA(ix,iy,x,y);
				matrix b = applyGaussianWeightsB(it,x,y);

				allA[x][y] = a;
				allB[x][y] = b;				

				matrix aT = getTranspose(a, 9, 2);
				matrix aTa = multiplyMatrix(aT,a,2,2,9);

				double lambda0 = getMinEigenValue2x2(aTa.m[0][0], aTa.m[0][1], aTa.m[1][0], aTa.m[1][1]);
				if (lambda0 > 0)
				{
					minEigenValues.m[x][y] = lambda0;	
					if (lambda0 > maximumMinEigenValue)
					{
						maximumMinEigenValue = lambda0;
					}	
				}

				delete aTa.m;
				delete aT.m;
			}
		}

		double threshold = 0.1*maximumMinEigenValue;
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				// Choosing possible points to be tracked
				if (minEigenValues.m[x][y] <= threshold)
				{
					minEigenValues.m[x][y] = 0.0;
				}
			}
		}

		// Only staying with one chosen pixel per window
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				reducingAmountOfPointsToTrack (x, y);
			}
		}
	}

	matrix getTranspose(matrix a, int width, int height)
	{
		matrix transpose;
		transpose.m = new double*[height];
		for (int i = 0; i < height; i++)
		{
			transpose.m[i] = new double[width];
		}

		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				transpose.m[y][x] = a.m[x][y];
			}
		}

		return transpose;
	}

	// knowing it`s a 2x2 matrix
	matrix getInvert (matrix a)
	{
		matrix invert;
		invert.m    = new double*[2];
		invert.m[0] = new double[2];
		invert.m[1] = new double[2];

		double d = 1 / (a.m[0][0]*a.m[1][1] - a.m[0][1]*a.m[1][0]);

		invert.m[0][0] =   d*a.m[1][1];
		invert.m[0][1] = - d*a.m[0][1];
		invert.m[1][0] = - d*a.m[1][0];
		invert.m[1][1] =   d*a.m[0][0];

		return invert;
	}

	matrix calculateFlow (matrix a, matrix b)
	{
		matrix aTranspose = getTranspose(a, 9, 2);
		matrix aTa = multiplyMatrix(aTranspose,a,2,2,9);

		matrix ataInvert = getInvert(aTa);

		matrix ataInvertATranpose = multiplyMatrix(ataInvert, aTranspose,2,9,2);
		matrix answer = multiplyMatrix(ataInvertATranpose, b,2,1,9);
		
		delete aTranspose.m;
		delete aTa.m;
		delete ataInvert.m;
		delete ataInvertATranpose.m;

		return answer;
	}

	matrix multiplyMatrix(matrix a, matrix b, int widthA, int heightB, int heightA)
	{
		matrix answer;
		answer.m = new double*[widthA];
		for (int i = 0; i < widthA; i++)
		{
			answer.m[i] = new double[heightB];
		}

		for (int x = 0; x < widthA; x++)
		{
			for (int y = 0; y < heightB; y++)
			{
				for (int z = 0; z < heightA; z++)
				{
					answer.m[x][y] += a.m[x][z]*b.m[z][y];
				}
			}
		}

		return answer;
	}

	double getMinEigenValue2x2(double& matA, double& matB, double& matC, double& matD) {
	    double b = matA+matD;
	    double c = matA*matD - matB*matC;

	    //b^2 - 4ac, where a=1.0
	    const double delta = b*b - 4*c;
	    if (delta < 0)
	        return -1.0;

	    double ev1 = 0.0;
	    double ev2 = 0.0;
	    double sqrtDelta = std::sqrt(delta);
	    
	    ev1 = b - sqrtDelta;
	    ev2 = b + sqrtDelta;

	    ev1 *= 0.5;
	    ev2 *= 0.5;

	    return std::min(ev1, ev2);
	}

	void reducingAmountOfPointsToTrack (int x, int y)
	{	
		double myEigen[] = {minEigenValues.m[x-1][y-1], minEigenValues.m[x-1][y], minEigenValues.m[x-1][y+1], 
			minEigenValues.m[x][y-1], minEigenValues.m[x][y], minEigenValues.m[x][y+1],
			minEigenValues.m[x+1][y-1], minEigenValues.m[x+1][y], minEigenValues.m[x+1][y+1]};

		double * maxValue = std::max_element(myEigen, myEigen+9);
		
		if (*maxValue > minEigenValues.m[x-1][y-1])
			minEigenValues.m[x-1][y-1] = 0.0;
		if (*maxValue > minEigenValues.m[x-1][y])
			minEigenValues.m[x-1][y] = 0.0;
		if (*maxValue > minEigenValues.m[x-1][y+1])
			minEigenValues.m[x-1][y+1] = 0.0;
		if (*maxValue > minEigenValues.m[x][y-1])
			minEigenValues.m[x][y-1] = 0.0;
		if (*maxValue > minEigenValues.m[x][y])
			minEigenValues.m[x][y] = 0.0;
		if (*maxValue > minEigenValues.m[x][y+1])
			minEigenValues.m[x][y+1] = 0.0;
		if (*maxValue > minEigenValues.m[x+1][y-1])
			minEigenValues.m[x+1][y-1] = 0.0;
		if (*maxValue > minEigenValues.m[x+1][y])
			minEigenValues.m[x+1][y] = 0.0;
		if (*maxValue > minEigenValues.m[x+1][y+1])
			minEigenValues.m[x+1][y+1] = 0.0;
	}

	matrix applyGaussianWeightsA (matrix ix, matrix iy, int x, int y)
	{	
		matrix a;
		a.m = new double*[9];
		for (int i = 0; i < 9; i++)
		{
			a.m[i] = new double[2];	
		}
		
		a.m[1] = new double[9];

		a.m[0][0] = ix.m[x-1][y-1]/16;
		a.m[0][1] = iy.m[x-1][y-1]/16;

		a.m[1][0] = 2 * ix.m[x][y-1]/16;
		a.m[1][1] = 2 * iy.m[x][y-1]/16;

		a.m[2][0] = ix.m[x+1][y-1]/16;
		a.m[2][1] = iy.m[x+1][y-1]/16;

		a.m[3][0] = 2 * ix.m[x-1][y]/16;
		a.m[3][1] = 2 * iy.m[x-1][y]/16;

		a.m[4][0] = 4 * ix.m[x][y]/16; // current pixel
		a.m[4][1] = 4 * iy.m[x][y]/16; // current pixel

		a.m[5][0] = 2 * ix.m[x+1][y]/16;
		a.m[5][1] = 2 * iy.m[x+1][y]/16;

		a.m[6][0] = ix.m[x-1][y+1]/16;
		a.m[6][1] = iy.m[x-1][y+1]/16;

		a.m[7][0] = 2 * ix.m[x][y+1]/16;
		a.m[7][1] = 2 * iy.m[x][y+1]/16;

		a.m[8][0] = ix.m[x+1][y+1]/16;
		a.m[8][1] = iy.m[x+1][y+1]/16;

		return a;
	}

	matrix applyGaussianWeightsB (matrix it, int x, int y)
	{
		// CImg<double> b(1, 9,depth,channel,initValue);
		matrix b;
		b.m = new double*[9];
		for (int i = 0; i < 9; i++)
		{
			b.m[i] = new double[1];	
		}
		

		b.m[0][0] = it.m[x-1][y-1]/16;

		b.m[1][0] = 2 * it.m[x][y-1]/16;

		b.m[2][0] = it.m[x+1][y-1]/16;

		b.m[3][0] = 2 * it.m[x-1][y]/16;

		b.m[4][0] = 4 * it.m[x][y]/16; // current pixel

		b.m[5][0] = 2 * it.m[x+1][y]/16;

		b.m[6][0] = it.m[x-1][y+1]/16;

		b.m[7][0] = 2 * it.m[x][y+1]/16;

		b.m[8][0] = it.m[x+1][y+1]/16;

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
	matrix getIt (CImg<double> image1, CImg<double> image2)
	{
		// Setting image`s attributes
		int height    = image1.height();
		int width     = image1.width();

		// Initializing return object
		matrix it;
		it.m = new double*[width];

		// Iterating over each pixel
		for (int x = 0; x < width; x++)
		{
			it.m[x] = new double[height];
			for (int y = 0; y < height; y++)
			{
				it.m[x][y] = image2(x,y) - image1(x,y);
			}
		}

		return it;
	}
	
	// Calculate It between two images, when working with pyramids
	matrix getIt (CImg<double> image1, CImg<double> image2, double xFlow, double yFlow)
	{
		// Setting image`s attributes
		int height = image1.height();
		int width  = image1.width();

		// Initializing return object
		matrix it;
		it.m = new double*[width];

		// Iterating over each pixel
		double newX, newY;
		for (int x = 0; x < width; x++)
		{
			it.m[x] = new double[height];
			for (int y = 0; y < height; y++)
			{
				newX = x + xFlow;
				newY = y + yFlow;
				if (newX < 0 || newX >= width || newY < 0 || newY >= height)
				{
					pointInArray p = pointInArrayFunction(points, (double) x, (double) y);
					if (p.inArray)
					{
						points.erase(points.begin() + p.position - 1);
					}
					it.m[x][y] = 0.0;
					continue;
				}

				point newPoint = bilinearInterpolation(newX, newY);
				if (newPoint.y == height || newPoint.x == width)
				{
					pointInArray p = pointInArrayFunction(points, (double) x, (double) y);
					if (p.inArray)
					{
						points.erase(points.begin() + p.position - 1);
					}
					it.m[x][y] = 0.0;
					continue;
				}

				it.m[x][y] = image2(newPoint.x, newPoint.y) - image1(x,y);
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
	std::vector<matrix> derive (CImg<double> image)
	{
		// Initilazing return object		
		std::vector<matrix > derived;

		// Setting image`s attributes
		int height    = image.height();
		int width     = image.width();

		matrix ix;
		matrix iy;

		ix.m = new double*[width];
		iy.m = new double*[width];

		for (int x = 0; x < width; x++)
		{
			ix.m[x] = new double[height];
			iy.m[x] = new double[height];
			for (int y = 0; y < height; y++)
			{
				if (x == 0 || x == (width - 1) || y == 0 || y == (height - 1))
				{
					ix.m[x][y] = 0.0;
					iy.m[x][y] = 0.0;
				}
				else
				{
					ix.m[x][y] = 0.5 * (image(x+1, y) - image(x-1, y));	
					iy.m[x][y] = 0.5 * (image(x, y+1) - image(x, y-1));
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
	int pyramidSize             = 4;
	double maximumMinEigenValue = 0.0;
	matrix ** allA;
	matrix ** allB;
	matrix minEigenValues;
	std::vector<ChosenPoint> points;
};