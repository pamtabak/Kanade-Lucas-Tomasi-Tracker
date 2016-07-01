#include "libs/CImg.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "gaussianpyramid.hpp"
#include "point.hpp"
#include <tuple>

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
		delete minEigenValues.m;
	}

	void initArrays(int width, int height)
	{
		allA = new matrix*[width];
		allB = new matrix*[width];

		minEigenValues.m = new double*[width];

		for (int i = 0; i < width; i++)
		{
			allA[i] = new matrix[height];
			allB[i] = new matrix[height];

			minEigenValues.m[i] = new double[height];
			for (int j = 0; j < height; j++)
			{
				minEigenValues.m[i][j] = 0.0;
			}
		}
	}

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

		images[0].save("images/output/Segments/piramide.png", 0);

		int height    = images[0].height();
		int width     = images[0].width();

		initArrays(width, height);

		for (int frame = 0; frame < numberOfFrames; frame++)
		{
			std::cout << "Frame: " << frame << std::endl;
			
			CImg<double> image1 = images[frame];
			CImg<double> image2 = images[frame + 1];

			std::vector<CImg<double> > imagesBeingUsed;
			imagesBeingUsed.push_back(image1);
			imagesBeingUsed.push_back(image2);
			std::vector<std::vector<CImg<double> >> pyramids = getGaussianPyramids(imagesBeingUsed);

			const unsigned char white[] = { 255,255,255 };
			const unsigned char pink [] = { 255,0,255 };

			matrix ix;
			matrix iy;
			matrix it;
			
			if (frame % 10 == 0)
			{
				std::vector<matrix> derived = derive(image1, 0, width - 1, 0, height - 1);
				ix                    = derived[0];
				iy                    = derived[1];

				it = getIt(image1, image2, 0, width - 1, 0, height - 1);

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
								chosenP.setInitialPoint(x,y);
								points.push_back(chosenP);
							}
						}
					}
				}
				std::cout << points.size() << std::endl;
			}

			for (int level = pyramidSize - 1; level >= 0; level--)
			{
				for (int p = 0; p < points.size(); p++)
				{
					double auxX = ((double) points[p].getPoint().x / pow(2,level));
					double auxY = ((double) points[p].getPoint().y / pow(2,level));

					point interpolatedPoint = bilinearInterpolation(auxX, auxY);
					int xOnLevel = (int) interpolatedPoint.x;
					int yOnLevel = (int) interpolatedPoint.y;

					if ((xOnLevel) <= 0 || (xOnLevel >= (pyramids[0][level].width() - 1)) || (yOnLevel <= 0) || (yOnLevel >= (pyramids[0][level].height() - 1)))
					{
						// points.erase(points.begin() + p - 1);
						points[p].updateFlow(0.0, 0.0, p);
						continue;
					}

					double initialFlowX = 2*points[p].getFlow()[frame].x;
					double initialFlowY = 2*points[p].getFlow()[frame].y;

					std::vector<matrix> derived = derive(pyramids[0][level], xOnLevel - 1, xOnLevel + 1, yOnLevel - 1, yOnLevel + 1);
					
					ix = derived[0];
					iy = derived[1];

					it = getIt(pyramids[0][level], pyramids[1][level], xOnLevel - 1, xOnLevel + 1, yOnLevel - 1, yOnLevel + 1 ,initialFlowX, initialFlowY);

					matrix a = applyGaussianWeightsA(ix, iy, xOnLevel, yOnLevel);
					matrix b = applyGaussianWeightsB(it, xOnLevel ,yOnLevel);

					matrix v = calculateFlow(a,b);
					if (!std::isnan(v.m[0][0]) && !std::isnan(v.m[1][0]))
					{
						// just checking if everything went ok with all matrixes transformations
						points[p].setFlow(v.m[0][0] + initialFlowX, v.m[1][0] + initialFlowY, frame);
					}
					else
					{
						points.erase(points.begin() + p - 1);
					}

					delete v.m;
					delete a.m;
					delete b.m;
				}
			}

			for (int p = 0; p < points.size(); p++)
			{
				// points[p].setPoint(points[p].getPoint().x + points[p].getFlow()[frame].x, points[p].getPoint().y + points[p].getFlow()[frame].y);
				points[p].setPoint(points[p].getPoint().x - points[p].getFlow()[frame].x, points[p].getPoint().y - points[p].getFlow()[frame].y);
				
				double finalX    = points[p].getPoint().x - points[p].getFlow()[frame].x;
				double finalY    = points[p].getPoint().y - points[p].getFlow()[frame].y;
				point finalPoint = bilinearInterpolation(finalX, finalY);

				double initFlowX = points[p].getPoint().x;
				double initFlowY = points[p].getPoint().y;
				double lastFlowX = finalPoint.x;
				double lastFlowY = finalPoint.y;
				
				if (frame >= 30)
				{
					for (int f = frame - 1; f >= frame - 30; f--)
					{
						image2.draw_line((int) initFlowX, (int) initFlowY, (int) lastFlowX, (int) lastFlowY, pink);
						initFlowX = lastFlowX;
						initFlowY = lastFlowY;
						lastFlowX = lastFlowX - points[p].getFlow()[f].x;
						lastFlowY = lastFlowY - points[p].getFlow()[f].y;
					}
				}
				else
				{
					for (int f = frame - 1; f >= 0; f--)
					{
						image2.draw_line((int) initFlowX, (int) initFlowY, (int) lastFlowX, (int) lastFlowY, pink);
						initFlowX = lastFlowX;
						initFlowY = lastFlowY;
						lastFlowX = lastFlowX - points[p].getFlow()[f].x;
						lastFlowY = lastFlowY - points[p].getFlow()[f].y;
					}					
				}
			}
			image2.save("images/output/Segments/piramide.png", frame + 1);

			delete ix.m;
			delete iy.m;
			delete it.m;
		}
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

				matrix aT  = getTranspose(a, 9, 2);
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
				if (minEigenValues.m[x][y] <= threshold)
				{
					minEigenValues.m[x][y] = 0.0;
					pointInArray p = pointInArrayFunction(points, (double) x, (double) y);
					if (p.inArray)
					{
						points.erase(points.begin() + p.position - 1);
					}
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

	matrix getInvert2x2 (matrix a)
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
		matrix aTa        = multiplyMatrix(aTranspose,a,2,2,9);

		matrix ataInvert  = getInvert2x2(aTa);

		matrix ataInvertATranpose = multiplyMatrix(ataInvert, aTranspose,2,9,2);
		matrix answer             = multiplyMatrix(ataInvertATranpose, b,2,1,9);
		
		delete aTranspose.m;
		delete aTa.m;
		delete ataInvert.m;
		delete ataInvertATranpose.m;

		return answer;
	}

	matrix multiplyMatrix(matrix a, matrix b, int widthA, int heightB, int inCommon)
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
				for (int z = 0; z < inCommon; z++)
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
		std::tuple<int, int> tuples[9] = { std::make_tuple(x-1, y-1), std::make_tuple(x-1, y), std::make_tuple(x-1, y+1),
								           std::make_tuple(x, y-1), std::make_tuple(x, y), std::make_tuple(x-1, y+1),
								           std::make_tuple(x+1, y-1), std::make_tuple(x+1, y), std::make_tuple(x+1, y+1)};

		for (int i = 0; i < 9; i++)
		{
			if (*maxValue > minEigenValues.m[std::get<0>(tuples[i])][std::get<1>(tuples[i])] &&  minEigenValues.m[std::get<0>(tuples[i])][std::get<1>(tuples[i])] > 0.0)
			{
				minEigenValues.m[std::get<0>(tuples[i])][std::get<1>(tuples[i])] = 0.0;
				pointInArray p = pointInArrayFunction(points, (double) x, (double) y);
				if (p.inArray)
				{
					points.erase(points.begin() + p.position - 1);
				}
			}	
		}						           
	}

	matrix applyGaussianWeightsA (matrix ix, matrix iy, int x, int y)
	{	
		matrix a;
		a.m = new double*[9];
		for (int i = 0; i < 9; i++)
		{
			a.m[i] = new double[2];	
		}

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
		GaussianPyramid gPyramid;
		double filter[5] = {1.0/16, 4.0/16, 6.0/16, 4.0/16, 1.0/16};
		gPyramid.generateFilter(filter);

		std::vector<std::vector<CImg<double> >> pyramids;
		for (int i = 0; i < images.size(); i++)
		{
			std::vector<CImg<double> > gaussianPyramid;
			gaussianPyramid.push_back(images[i]);
			
			CImg<double> reducedImage  = images[i];
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

	matrix getIt (CImg<double> image1, CImg<double> image2, int xFrom, int xTo, int yFrom, int yTo, double xFlow = 0.0, double yFlow = 0.0)
	{
		int height    = image1.height();
		int width     = image1.width();

		matrix it;
		it.m = new double*[width];

		for (int x = xFrom; x <= xTo; x++)
		{
			it.m[x] = new double[height];
			for (int y = yFrom; y <= yTo; y++)
			{
				double image2x = x + xFlow;
				double image2y = y + yFlow;
				if (image2x < 0 || image2x >= width || image2y < 0 || image2y >= height)
				{
					pointInArray p = pointInArrayFunction(points, (double) x, (double) y);
					if (p.inArray)
					{
						points.erase(points.begin() + p.position - 1);
					}
					it.m[x][y] = 0.0;
					continue;
				}
				
				point newPoint = bilinearInterpolation(image2x, image2y);
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

				it.m[x][y] = - (image2(newPoint.x, newPoint.y) - image1(x,y));
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

	std::vector<matrix> derive (CImg<double> image, int xFrom, int xTo, int yFrom, int yTo)
	{
		std::vector<matrix > derived;

		int height    = image.height();
		int width     = image.width();

		matrix ix;
		matrix iy;

		ix.m = new double*[width];
		iy.m = new double*[width];

		for (int x = xFrom; x <= xTo; x++)
		{
			ix.m[x] = new double[height];
			iy.m[x] = new double[height];
			for (int y = yFrom; y <= yTo; y++)
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