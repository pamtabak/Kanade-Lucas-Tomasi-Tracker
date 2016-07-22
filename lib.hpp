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
		
		MatrixXd ix(width, height);
		MatrixXd iy(width, height);
		MatrixXd it(width, height);
		MatrixXd minEigenValue(width, height);

		int numberOfFrames = images.size() - 1;
		images[0].save("images/output/Segments/piramide.png", 0);

		for (int frame = 0; frame < 1; frame++)
		// for (int frame = 0; frame < numberOfFrames; frame++)
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

				// maxMinEigenValue = 0.0;
				calculateMinEigenValue(ix, it, iy, minEigenValue);
			}

			for (int level = pyramidSize - 1; level >= 0; level--)
			{
				// Calculate flow
			}
		}
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
				int finalX = x + flow.x;
				int finalY = y + flow.y;

				if ((finalX < 0) || (finalX >= width) || (finalY < 0) || (finalY >= height))
				{
					// point is outside image. 
					// Decide what to do
					it(x,y) = 0.0; // Is this correct?
					continue;
				}

				point p = bilinearInterpolation(finalX, finalY);

				it(x,y) = - (image2(p.x, p.y) - image1(x,y));
			}
		}
	}

	void calculateMinEigenValue (MatrixXd &ix, MatrixXd &iy, MatrixXd &it, MatrixXd &minEigenValue)
	{
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				MatrixXd a = calculateA(ix, iy, x, y);
				MatrixXd b = calculateB(it, x, y);

				MatrixXd aTa = a.transpose() * a;
				// EigenSolver<MatrixXd> es(aTa);
				double lambda0 = getMinEigenValue2x2(aTa(0,0), aTa(0,1), aTa(1,0), aTa(1,1));
				if (lambda0 > 0)
				{
					minEigenValues(x,y) = lambda0;
					if (lambda0 > maxMinEigenValue)
					{
						maxMinEigenValue = lambda0;
					}	
				}

				// std::cout << es.eigenvalues() << std::endl;
				// double lambda0 = es.eigenvalues();
			}
		}

		// THRESHOLD

		// CLEANING MINEIGENVALUE
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