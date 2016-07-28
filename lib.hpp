#undef Success 

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "libs/CImg.h"
#include "point.hpp"
#include "gaussianpyramid.hpp"



using namespace cimg_library;
using namespace Eigen;

class Tracking
{
public:
	Tracking()
	{

	}

	~Tracking()
	{

	}

	void lucasKanade (std::vector<CImg<double> > images)
	{
		// Initializing variables
		width  = images[0].width();
		height = images[0].height();

		const unsigned char pink [] = { 255,0,255 };
		
		MatrixXd ix(width, height);
		MatrixXd iy(width, height);
		MatrixXd it(width, height);
		MatrixXd minEigenValue(width, height);

		int numberOfFrames = images.size() - 1;
		images[0].save("images/output/Segments/piramide.png", 0);

		ChosenPoint ** points = new ChosenPoint*[width];
		for (int i = 0; i < width; i++)
		{
			points[i] = new ChosenPoint[height];
			for (int j = 0; j < height; j++)
			{
				points[i][j].setPoint(i, j, false);
				points[i][j].setNumberOfFrames(numberOfFrames);
			}
		}

		for (int frame = 0; frame < numberOfFrames; frame++)
		{
			// DO THE MAGIC!

			std::cout << "Frame: " << frame << std::endl;
			
			CImg<double> image1 = images[frame];
			CImg<double> image2 = images[frame + 1];

			std::vector<CImg<double>> image1pyramids = getGaussianPyramids(image1);
			std::vector<CImg<double>> image2pyramids = getGaussianPyramids(image2);

			if (frame % 10 == 0)
			{
				// recalculate points to be tracked
				point from, to, flow;
				from.x = 0.0;
				from.y = 0.0;
				to.x   = width  - 1;
				to.y   = height - 1;
				flow.x = 0.0;
				flow.y = 0.0;

				calculateIxAndIy(image1, from, to, ix, iy);
				calculateIt(image1, image2, from, to, flow, it);

				calculateMinEigenValue(ix, it, iy, minEigenValue);

				for (int y = 1; y < height - 1; y++)
				{
					for (int x = 1; x < width - 1; x++)
					{
						// Verify if point is already inside array
						if (minEigenValue(x, y) > 0.0)
						{
							points[x][y].setPoint(points[x][y].getPoint().x, points[x][y].getPoint().y, true);
						}
						else
						{
							points[x][y].setPoint(points[x][y].getPoint().x, points[x][y].getPoint().y, false);	
						}
					}
				}
			}

			for (int level = pyramidSize - 1; level >= 0; level--)
			{
				// Calculate flow
				for (int y = 1; y < height - 1; y++)
				{
					for (int x = 1; x < width - 1; x++)
					{
						if (!points[x][y].getPoint().isValid)
						{
							continue;
						}

						double auxX = (double) points[x][y].getPoint().x / pow (2, level);
						double auxY = (double) points[x][y].getPoint().y / pow (2, level);
						
						point interpolatedPoint = bilinearInterpolation(auxX, auxY);
						int xOnLevel            = interpolatedPoint.x;
						int yOnLevel            = interpolatedPoint.y;

						// Verify if pixel is outside image
						if (xOnLevel == 0 || xOnLevel == (image1pyramids[level].width() - 1) || yOnLevel == 0 || yOnLevel == (image1pyramids[level].height() - 1))
						{
							points[x][y].updateFlow(0.0, 0.0, frame);
							points[x][y].setPoint(x, y, false);
							continue;
						}

						point from, to, flow;
						from.x = xOnLevel - 1;
						from.y = yOnLevel - 1;
						to.x   = xOnLevel + 1;
						to.y   = yOnLevel + 1;
						flow.x = 2 * points[x][y].getFlow()[frame].x;
						flow.y = 2 * points[x][y].getFlow()[frame].y;

						calculateIxAndIy(image1pyramids[level], from, to, ix, iy);
						calculateIt(image1pyramids[level], image2pyramids[level], from, to, flow, it);

						MatrixXd a = calculateA(ix, iy, xOnLevel, yOnLevel);
						MatrixXd b = calculateB(it, xOnLevel, yOnLevel);

						MatrixXd aT      = a.transpose();
						MatrixXd aTa     = aT * a; 
						MatrixXd inverse = aTa.inverse();

						MatrixXd v = inverse * aT * b;
						if (v.rows() == 2 && v.cols() == 1)
						{
							// std::cout << v(0,0) << "," << v(1,0) << std::endl;
							// points[x][y].setFlow(finalX, finalY, frame);
							points[x][y].setFlow(v(0,0), v(1,0), frame);
						}
						else
						{
							points[x][y].updateFlow(0.0, 0.0, frame);
							points[x][y].setPoint(points[x][y].getPoint().x, points[x][y].getPoint().y, false);
						}
					}
				}
			}

			// Draw
			for (int y = 1; y < height - 1; y++)
			{
				for (int x = 1; x < width - 1; x++)
				{
					if (points[x][y].getPoint().isValid)
					{
						// bilinear interpolation
						point finalPoint = bilinearInterpolation(points[x][y].getPoint().x + points[x][y].getFlow()[frame].x, points[x][y].getPoint().y + points[x][y].getFlow()[frame].y);
						image2.draw_line(points[x][y].getPoint().x, points[x][y].getPoint().y, finalPoint.x, finalPoint.y, pink);
						points[x][y].setPoint(finalPoint.x, finalPoint.y, true);
					}
				}
			}
			// image2.display();
			image2.save("images/output/Segments/piramide.png", frame + 1);
		}

		for (int i = 0; i < height; i++)
		{
			delete [] points[i];
		}
		delete [] points;
	}

	void calculateIxAndIy (CImg<double> image, point from, point to, MatrixXd &ix, MatrixXd &iy)
	{
		for (int x = from.x; x <= to.x; x++)
		{
			for (int y = from.y; y <= to.y; y++)
			{
				if ((x == 0) || (x == (width - 1)))
					ix(x,y) = 0.0;
				else
					ix(x,y) = 0.5 * (image(x+1, y) - image(x-1, y));


				if ((y == 0) || (y == (height - 1)))
					iy(x,y) = 0.0;
				else
					iy(x,y) = 0.5 * (image(x, y+1) - image(x, y-1));
			}
		}
	}

	void calculateIt (CImg<double> image1, CImg<double> image2, point from, point to, point flow, MatrixXd &it)
	{
		for (int x = from.x; x <= to.x; x++)
		{
			for (int y = from.y; y <= to.y; y++)
			{
				double finalX = x + flow.x;
				double finalY = y + flow.y;

				if ((finalX < 0) || (finalX >= width) || (finalY < 0) || (finalY >= height))
				{
					// points[x][y].setPoint(x, y, false); // ?
					it(x,y) = 0.0; // Is this correct?
					continue;
				}

				point p = bilinearInterpolation(finalX, finalY);

				it(x,y) = - (image2(p.x, p.y) - image1(x,y));
			}
		}
	}

	void calculateMinEigenValue (MatrixXd &ix, MatrixXd &iy, MatrixXd &it, MatrixXd &minEigenValues)
	{
		// maxMinEigenValue = 0.0;
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				MatrixXd a = calculateA(ix, iy, x, y);
				MatrixXd b = calculateB(it, x, y);

				MatrixXd aTa = a.transpose() * a;
				double lambda0 = getMinEigenValue2x2(aTa(0,0), aTa(0,1), aTa(1,0), aTa(1,1));
				if (lambda0 > 0.0)
				{
					minEigenValues(x,y) = lambda0;
					if (lambda0 > maxMinEigenValue)
					{
						maxMinEigenValue = lambda0;
					}
				}				
			}
		}

		// THRESHOLD
		double threshold = 0.1*maxMinEigenValue;
		for (int y = 1; y < height - 1; y++)
		{
			for (int x = 1; x < width - 1; x++)
			{
				if (minEigenValues(x,y) < threshold)
				{
					minEigenValues(x,y) = 0.0;
				}
			}
		}

		// CLEANING MINEIGENVALUE
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				reducingAmountOfPointsToTrack (x, y, minEigenValues);
			}
		}
	}

	double getMinEigenValue2x2(double& a00, double& a01, double& a10, double& a11) 
	{
	    double b = a00 + a11;
	    double c = (a00 * a11) - (a01 * a10);

	    const double delta = pow(b, 2) - (4 * c);
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

	void reducingAmountOfPointsToTrack (int x, int y, MatrixXd &minEigenValues)
	{
		double myEigen[] = {minEigenValues(x-1,y-1), minEigenValues(x-1,y), minEigenValues(x-1,y+1),
			                minEigenValues(x,y-1),   minEigenValues(x,y),   minEigenValues(x,y+1),
			                minEigenValues(x+1,y-1), minEigenValues(x+1,y), minEigenValues(x+1,y+1)};

		double * maxValue = std::max_element(myEigen, myEigen+9);
		std::tuple<int, int> tuples[9] = { std::make_tuple(x-1, y-1), std::make_tuple(x-1, y), std::make_tuple(x-1, y+1),
								           std::make_tuple(x, y-1), std::make_tuple(x, y), std::make_tuple(x-1, y+1),
								           std::make_tuple(x+1, y-1), std::make_tuple(x+1, y), std::make_tuple(x+1, y+1)};

		for (int i = 0; i < 9; i++)
		{
			if (*maxValue > minEigenValues(std::get<0>(tuples[i]),std::get<1>(tuples[i])) &&  minEigenValues(std::get<0>(tuples[i]),std::get<1>(tuples[i])) > 0.0)
			{
				minEigenValues(std::get<0>(tuples[i]),std::get<1>(tuples[i])) = 0.0;
			}	
		}		
	}

	MatrixXd calculateA (MatrixXd &ix, MatrixXd &iy, int x, int y)
	{
		MatrixXd a(9,2);

		a(0,0) = 1 * ix(x-1,y-1)/16;
		a(0,1) = 1 * iy(x-1,y-1)/16;
		a(1,0) = 2 * ix(x,y-1)/16;
		a(1,1) = 2 * iy(x,y-1)/16;
		a(2,0) = 1 * ix(x+1,y-1)/16;
		a(2,1) = 1 * iy(x+1,y-1)/16;
		a(3,0) = 2 * ix(x-1,y)/16;
		a(3,1) = 2 * iy(x-1,y)/16;
		a(4,0) = 4 * ix(x,y)/16; // current pixel
		a(4,1) = 4 * iy(x,y)/16; // current pixel
		a(5,0) = 2 * ix(x+1,y)/16;
		a(5,1) = 2 * iy(x+1,y)/16;
		a(6,0) = 1 * ix(x-1,y+1)/16;
		a(6,1) = 1 * iy(x-1,y+1)/16;
		a(7,0) = 2 * ix(x,y+1)/16;
		a(7,1) = 2 * iy(x,y+1)/16;
		a(8,0) = 1 * ix(x+1,y+1)/16;
		a(8,1) = 1 * iy(x+1,y+1)/16;

		return a;
	}

	MatrixXd calculateB (MatrixXd &it, int x, int y)
	{
		MatrixXd b(9,1);

		b(0,0) = 1 * it(x-1,y-1)/16;
		b(1,0) = 2 * it(x,y-1)/16;
		b(2,0) = 1 * it(x+1,y-1)/16;
		b(3,0) = 2 * it(x-1,y)/16;
		b(4,0) = 4 * it(x,y)/16; // current pixel
		b(5,0) = 2 * it(x+1,y)/16;
		b(6,0) = 1 * it(x-1,y+1)/16;
		b(7,0) = 2 * it(x,y+1)/16;
		b(8,0) = 1 * it(x+1,y+1)/16;

		return b;
	}

	point bilinearInterpolation(double x, double y)
	{
		point a,b,c,d;
		a.x = floor(x);
		a.y = floor(y);

		b.x = floor(x);
		b.y = ceil(y);
		if (b.y == height)
			b.y -= 1;

		c.x = ceil(x);
		if (c.x == width)
			c.x -= 1;
		c.y = floor(y);

		d.x = ceil(x);
		if (d.x == width)
			d.x -= 1;
		d.y = ceil(y);
		if (d.y == height)
			d.y -= 1;

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

	std::vector<CImg<double>> getGaussianPyramids (CImg<double> image)
	{
		GaussianPyramid gPyramid;
		double filter[5] = {1.0/16, 4.0/16, 6.0/16, 4.0/16, 1.0/16};
		gPyramid.generateFilter(filter);

		std::vector<CImg<double> > gaussianPyramid;
		gaussianPyramid.push_back(image);
		
		CImg<double> reducedImage  = image;
		for (int p = 1; p < pyramidSize; p++)
		{
			CImg<double> imageGenerated = gPyramid.reduce(reducedImage);
			gaussianPyramid.push_back(imageGenerated);
			reducedImage = imageGenerated;
		}

		return gaussianPyramid;
	}

private:
	int width;
	int height;
	int pyramidSize = 4;
	double maxMinEigenValue = 0.0;
};