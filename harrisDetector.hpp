#include "libs/CImg.h"
#include <iostream>

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
		CImg<double> ix(width, height, 1, 1, 0);
		CImg<double> iy(width, height, 1, 1, 0);
		// Avoiding dealing with border
		for (int x = 1; x < width - 1; x++)
		{
			for (int y = 1; y < height - 1; y++)
			{
				ix(x,y,0,1) = image(x+1,y,0,1) - image(x-1,y,0,1);
				iy(x,y,0,1) = image(x,y+1,0,1) - image(x,y-1,0,1);
			}
		}
		ix(0,0,0,1)          = image(0,0,0,1);
		iy(0,0,0,1)          = image(0,0,0,1);
		ix(0,height,0,1)     = image(0,height,0,1);
		iy(0,height,0,1)     = image(0,height,0,1);
		ix(width,0,0,1)      = image(width,0,0,1);
		iy(width,0,0,1)      = image(width,0,0,1);
		ix(width,height,0,1) = image(width,height,0,1);
		iy(width,height,0,1) = image(width,height,0,1);


		// 2. calcular os produtos Ix2 = IxIx, Iy2 = IyIy e Ixy = IxIy
		CImg<double> ixx = ix*ix;
		CImg<double> iyy = iy*iy;
		CImg<double> ixy = ix*iy;

		// 3. convoluir as 3 imagens com uma Gaussiana
		CImg<double> filterIxx = filterImage(ixx);
		CImg<double> filterIyy = filterImage(iyy);
		CImg<double> filterIxy = filterImage(ixy);

		// 4. para cada pixel: encontrar os autovalores e utilizar uma das medidas

		//  5. limiarizar para encontrar maximos

		// 6. limitar numero de maximos por regiao
	}

	CImg<double> filterImage (CImg<double> image)
	{
		int height = image.height();
		int width  = image.width();

		CImg<double> filteredImage(width, height, 1, 1,0);

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