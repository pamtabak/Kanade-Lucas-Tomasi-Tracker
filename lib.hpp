#undef Success // in order to make Eigen work

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

		points = new ChosenPoint*[width];
		for (int i = 0; i < width; i++)
		{
			points[i] = new ChosenPoint[height];
			for (int j = 0; j < height; j++)
			{
				points[i][j].setPoint(i, j, false);
				points[i][j].setNumberOfFrames(numberOfFrames);
			}
		}

		// for (int frame = 0; frame < 1; frame++)
		for (int frame = 0; frame < numberOfFrames; frame++)
		{
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
						if (minEigenValue(x, y) > 0.0)
						{
							points[x][y].setPoint(points[x][y].getPoint().x, points[x][y].getPoint().y, true);
						}
						// else
						// {
						// 	points[x][y].setPoint(points[x][y].getPoint().x, points[x][y].getPoint().y, false);	
						// }
					}
				}
			}

			// for (int level = 0; level >= 0; level--)
			for (int level = pyramidSize - 1; level >= 0; level--)
			{
				// Calculate flow
				for (int y = 1; y < height - 1; y++)
				{
					for (int x = 1; x < width - 1; x++)
					{
						if (points[x][y].getPoint().isValid)
						{
							double auxX = (double) points[x][y].getPoint().x / pow (2, level);
							double auxY = (double) points[x][y].getPoint().y / pow (2, level);
							// double auxX = (double) x / pow (2, level);
							// double auxY = (double) y / pow (2, level);
							
							point interpolatedPoint = bilinearInterpolation(auxX, auxY);
							int xOnLevel            = interpolatedPoint.x;
							int yOnLevel            = interpolatedPoint.y;

							// Verify if pixel is outside image
							if (xOnLevel <= 0 || xOnLevel >= (image1pyramids[level].width() - 1) || yOnLevel <= 0 || yOnLevel >= (image1pyramids[level].height() - 1))
							{
								points[x][y].updateFlow(0.0, 0.0, frame);
								points[x][y].setPoint(points[x][y].getPoint().x, points[x][y].getPoint().y, false);
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

							MatrixXd aTa     = calculateAta(ix, iy, xOnLevel, yOnLevel);
							MatrixXd aTb     = calculateAtB(ix, iy, it, xOnLevel, yOnLevel);
							MatrixXd inverse = aTa.inverse();

							MatrixXd v = inverse * aTb;
							if (v.rows() == 2 && v.cols() == 1)
							{
								flow.x += v(0,0);
								flow.y += v(1.0);
								points[x][y].setFlow(flow.x, flow.y, frame);

								if (level == 0)
								{
									// DRAW
									double finalX    = points[x][y].getPoint().x - points[x][y].getFlow()[frame].x;
									double finalY    = points[x][y].getPoint().y - points[x][y].getFlow()[frame].y;
									point finalPoint = bilinearInterpolation(finalX, finalY);

									double initFlowX = points[x][y].getPoint().x;
									double initFlowY = points[x][y].getPoint().y;
									double lastFlowX = finalPoint.x;
									double lastFlowY = finalPoint.y;

									points[x][y].setPoint(finalPoint.x, finalPoint.y, true);

									int finalFrame;
									if (frame >= 30)
									{
										finalFrame = frame - 30;
									}
									else
									{
										finalFrame = 0;
									}

									image2.draw_line((int) initFlowX, (int) initFlowY, (int) lastFlowX, (int) lastFlowY, pink);
									for (int f = frame - 1; f >= finalFrame; f--)
									{
										image2.draw_line((int) initFlowX, (int) initFlowY, (int) lastFlowX, (int) lastFlowY, pink);
										initFlowX = lastFlowX;
										initFlowY = lastFlowY;
										lastFlowX = lastFlowX - points[x][y].getFlow()[f].x;
										lastFlowY = lastFlowY - points[x][y].getFlow()[f].y;
									}
									image2.save("images/output/Segments/piramide.png", frame + 1);
								}
							}
							else
							{
								points[x][y].updateFlow(0.0, 0.0, frame);
								points[x][y].setPoint(points[x][y].getPoint().x, points[x][y].getPoint().y, false);
					     	}
						}
					}
				}
			}
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
				
				point p = bilinearInterpolation(finalX, finalY);

				if ((p.x < 0) || (p.x >= width) || (p.y < 0) || (p.y >= height))
				{
					points[x][y].setPoint(points[x][y].getPoint().x, points[x][y].getPoint().y, false);
					it(x,y) = 0.0;
				}
				else
				{
					it(x,y) = - (image2(p.x, p.y) - image1(x,y));
				}
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
				MatrixXd aTa = calculateAta(ix, iy, x, y);
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
		double myEigen[]  = { minEigenValues(x-1,y-1), minEigenValues(x-1,y), minEigenValues(x-1,y+1),
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

	MatrixXd calculateAta (MatrixXd &ix, MatrixXd &iy, int x, int y)
	{
		MatrixXd aTa(2,2);
		aTa(0,0) = 0.0;
		aTa(0,1) = 0.0;
		aTa(1,0) = 0.0;
		aTa(1,1) = 0.0;

		double weight[3][3] = { { 1.0/16, 2.0/16, 1.0/16 }, { 2.0/16, 4.0/16, 2.0/16 }, { 1.0/16, 2.0/16, 1.0/16 } };

		int iCount = 0;
		for (int i = x - 1; i <= x + 1; i++)
		{
			int jCount = 0;
			for (int j = y - 1; j <= y + 1; j++)
			{
				aTa(0,0) += ix(i,j) * ix(i,j) * weight[iCount][jCount];
				aTa(0,1) += ix(i,j) * iy(i,j) * weight[iCount][jCount];
				aTa(1,0) += ix(i,j) * iy(i,j) * weight[iCount][jCount];
				aTa(1,1) += iy(i,j) * iy(i,j) * weight[iCount][jCount];
				jCount++;
			}
			iCount++;
		}

		return aTa;
	}

	MatrixXd calculateAtB (MatrixXd &ix, MatrixXd &iy, MatrixXd &it, int x, int y)
	{
		MatrixXd aTb(2,1);
		aTb(0,0) = 0.0;
		aTb(1,0) = 0.0;

		double weight[3][3] = { { 1.0/16, 2.0/16, 1.0/16 }, { 2.0/16, 4.0/16, 2.0/16 }, { 1.0/16, 2.0/16, 1.0/16 } };

		int iCount = 0;
		for (int i = x - 1; i <= x + 1; i++)
		{
			int jCount = 0;
			for (int j = y - 1; j <= y + 1; j++)
			{
				aTb(0,0) += ix(i,j) * it(i,j) * weight[iCount][jCount];
				aTb(1,0) += iy(i,j) * it(i,j) * weight[iCount][jCount];
				jCount++;
			}

			aTb(0,0) = -1.0 * aTb(0,0);
			aTb(1,0) = -1.0 * aTb(1,0);

			iCount++;
		}

		return aTb;
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
	ChosenPoint ** points;
};