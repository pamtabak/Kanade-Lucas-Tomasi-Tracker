#include "libs/CImg.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace cimg_library;

class HarrisDetector
{
public:
	HarrisDetector()
	{
		generateFilter();
	}
	~HarrisDetector()
	{

	}

	void algorithm (CImg<double> image)
	{
		int height    = image.height();
		int width     = image.width();
		int depth     = 1;
		int channel   = 1;
		int initValue = 0;

		// 1. Calculate Ix, Iy
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

		// 2. Calculating Ix2 = IxIx, Iy2 = IyIy e Ixy = IxIy
		CImg<double> ixx = ix*ix;
		CImg<double> iyy = iy*iy;
		CImg<double> ixy = ix*iy;

		// 3. Filter images with Gaussian Filter
		CImg<double> filterIxx = filterImage(ixx);
		CImg<double> filterIyy = filterImage(iyy);
		CImg<double> filterIxy = filterImage(ixy);

		// 4. para cada pixel: encontrar os autovalores e utilizar uma das medidas
		std::vector<std::vector<double> > eigenValues;
		std::vector<double> sortedValues;
		double maxValue = 0.0;
		for (int x = 0; x < width; x++)
		{
			std::vector<double> row;
			for (int y = 0; y < height; y++)
			{
				// [[Ixx Ixy],[Ixy Iyy]]
				CImg<double> matrix(2,2,depth,channel,initValue);
				matrix(0,0) = filterIxx(x,y);
				matrix(0,1) = filterIxy(x,y);
				matrix(1,0) = filterIxy(x,y);
				matrix(1,1) = filterIyy(x,y);

				CImgList<double> eigen      = matrix.get_eigen();
				CImg<double> eigenValuesImg = eigen(0);
				double lambda0              = eigenValuesImg(0,0);
				double lambda1              = eigenValuesImg(0,1);

				double measure;
				// measure = abs(lambda0*lambda1 - 0.06*(pow(lambda0+lambda1, 2)));
				if (lambda0 < lambda1)
				{
					measure = abs(lambda0 - 0.05*lambda1);
				}
				else
				{
					measure = abs(lambda1 - 0.05*lambda0);	
				}
				sortedValues.push_back(measure);
				row.push_back(measure);

				// Using this variable to define a threshold
				if (measure > maxValue)
				{
					maxValue = measure;
				}
			}
			eigenValues.push_back(row);	
		}

		//  5. limiarizar para encontrar maximos
		std::sort(sortedValues.begin(), sortedValues.end());
		std::unique(sortedValues.begin(),sortedValues.end());
		int thresholdIndex = (int) sortedValues.size() * 0.9;
		double threshold = sortedValues[thresholdIndex];
		// CImg<double> points(width, height, depth, channel, initValue);
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				if (threshold <= eigenValues[x][y])
				{
					image(x,y,0,0) = 255.0;
					image(x,y,0,1) = 0.0;
					image(x,y,0,2) = 255.0;
				}
			}
		}

		// 6. limitar numero de maximos por regiao
		// for (int x = 0; x < width - 5; x++)
		// {
		// 	for (int y = 0; y < height - 5; y++)
		// 	{
		// 		double maxValue = 0.0;
		// 		for (int i = 0; i < 5; i++)
		// 		{
		// 			for (int j = 0; j < 5; j++)
		// 			{
		// 				if (points(x + i, y + j) == 255.0)
		// 				{
		// 					if (eigenValues[x+i][i+j] > maxValue)
		// 					{
		// 						maxValue = eigenValues[x+i][i+j];
		// 					}
		// 				}
		// 			}
		// 		}
		// 		for (int i = 0; i < 5; i++)
		// 		{
		// 			for (int j = 0; j < 5; j++)
		// 			{
		// 				if (eigenValues[x+i][y+j] < maxValue)
		// 				{
		// 					points(x+i,y+j) = 0.0;
		// 				}
		// 			}
		// 		}
		// 	}
		// }

		image.display();
		image.save("test.png");
	}

	CImg<double> filterImage (CImg<double> image)
	{
		int height = image.height();
		int width  = image.width();

		CImg<double> filteredImage(width, height, 1, 1, 0);

		for (int y = 0; y < height; y++)
		{
			for (int x = 0; x < width; x++)
			{
				filteredImage(x,y) = 0.0;
				for (int i = -2; i <= 2; i++)
				{
					for (int j = -2; j <= 2; j++) 
					{
						// we need to verify the edges
						int newValueOfX = ((x + i) < 0) ? (-x - i - 1) : x + i;
						newValueOfX = (newValueOfX >= (width)) ? ((2 * (width)) - x - i - 1) : newValueOfX;
						int newValueOfY = ((y + j) < 0) ? (-y - j - 1) : y + j;
						newValueOfY = (newValueOfY >= (height)) ? ((2 * (height)) - y - j - 1) : newValueOfY;
						filteredImage(x,y) += image(newValueOfX, newValueOfY) * filter2d[i+2][j+2];
					}
				}
			}
		}
		return filteredImage;	
	}

	void generateFilter () {
		for (int i = 0; i < 5; i++) {
			for (int j = 0; j < 5; j++) {
				filter2d[i][j] = filter[i]*filter[j];
			}
		}
	}

private:	
	double filter2d[5][5];
	double filter[5] = {1.0/16, 4.0/16, 6.0/16, 4.0/16, 1.0/16};
};