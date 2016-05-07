#include "libs/CImg.h"
#include <iostream>
#include <vector>
#include <math.h>

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
		int height = image.height();
		int width  = image.width();

		/* 1. calcular Ix e Iy : gradiente na direcao x e y pode usar 
		por exemplo: convolucao da derivada de uma Gaussiana */

		// Calculating Ix, Iy
		CImg<double> ix(width, height,1,1,0);
		CImg<double> iy(width, height,1,1,0);
		// Avoiding dealing with border
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				ix(x,y) = image(x+1,y) - image(x-1,y);
				iy(x,y) = image(x,y+1) - image(x,y-1);
			}
		}
		ix(0,0)                  = image(0,0);
		iy(0,0)                  = image(0,0);
		ix(0,height - 1)         = image(0,height - 1);
		iy(0,height - 1)         = image(0,height - 1);
		ix(width - 1,0)          = image(width - 1,0);
		iy(width - 1,0)          = image(width - 1,0);
		ix(width - 1,height - 1) = image(width - 1,height - 1);
		iy(width - 1,height - 1) = image(width - 1,height - 1);

		// 2. calcular os produtos Ix2 = IxIx, Iy2 = IyIy e Ixy = IxIy
		CImg<double> ixx = ix*ix;
		CImg<double> iyy = iy*iy;
		CImg<double> ixy = ix*iy;

		// 3. convoluir as 3 imagens com uma Gaussiana
		CImg<double> filterIxx = filterImage(ixx);
		CImg<double> filterIyy = filterImage(iyy);
		CImg<double> filterIxy = filterImage(ixy);

		// 4. para cada pixel: encontrar os autovalores e utilizar uma das medidas
		std::vector<std::vector<double> > eigenValues;
		double maxValue = 0.0;
		for (int x = 0; x < width; x++)
		{
			std::vector<double> row;
			for (int y = 0; y < height; y++)
			{
				// [[Ixx Ixy],[Ixy Iyy]]
				CImg<double> matrix(2,2,1,1,0);
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

				std::cout << measure << std::endl;
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
		double threshold = 0.5*maxValue;
		CImg<double> points(width, height, 1, 1, 0);
		for (int x = 0; x < width; x++)
		{
			for (int y = 0; y < height; y++)
			{
				if (threshold <= eigenValues[x][y])
				{
					image(x,y,0,0) = 255.0;
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