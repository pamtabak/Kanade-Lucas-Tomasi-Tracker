
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
					flowImage(x,y,0,1) = 0.0;    // g
					flowImage(x,y,0,2) = v(0,1); // b
					image1.draw_line(x, y ,x + (int) v(0,0),y + (int) v(0,1), white);
				}
			}
		}
		flowPyramid.push_back(flowImage);

		int iteration = 0;
		for (int p = pyramidSize - 2; p >= 0; p--)
		{
			image1 = pyramids[0][p];
			image2 = pyramids[1][p];

			height    = image1.height();
			width     = image1.width();

			minEigenValues.assign(width,height,depth,channel,0);
			flowImage.assign(width,height,depth,3,0);
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
					if (x > 0 && y > 0)
					{
						if ((x % 2 == 0) && (y % 2 == 0))
						{
							// checking if flow is different than 0 in this position
							if (flowPyramid[iteration](x/2, y/2, 0,0) > 0.0 || flowPyramid[iteration](x/2, y/2, 0,2) > 0.0)
							{
								// Calculating it
								it = getIt(image1, image2, flowPyramid[iteration](x/2, y/2, 0,0), flowPyramid[iteration](x/2, y/2, 0,2));
							}
						}
						else
						{
							it = getIt(image1, image2);	
						}

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
			for (int x = 1; x < width - 1; x++)
			{
				for (int y = 1; y < height - 1; y++)
				{
					minEigenValues = reducingAmountOfPointsToTrack (minEigenValues, x, y);
				}
			}
			
			for (int x = 0; x < width - 1; x++)
			{
				for (int y = 0; y < height - 1; y++)
				{
					if (x > 0 && y > 0 && (x % 2 == 0) && (y % 2 == 0))
					{
						// std::cout << minEigenValues(x,y) << std::endl;
						if (minEigenValues(x,y) > 0.0)
						{
							CImg<double> v    = ((allA[x][y].get_transpose() * allA[x][y]).get_invert())*allA[x][y].get_transpose()*allB[x][y];	
							if (!std::isnan(v(0,0)) && !std::isnan(v(0,1)))
							{
								flowImage(x,y,0,0) = 2*(flowPyramid[iteration](x/2, y/2, 0,0)) + v(0,0); // r
								flowImage(x,y,0,1) = 0.0; // g
								flowImage(x,y,0,2) = 2*(flowPyramid[iteration](x/2, y/2, 0,2)) + v(0,1); // b
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
							flowImage(x,y,0,0) = v(0,0); // r
							flowImage(x,y,0,1) = 0.0;    // g
							flowImage(x,y,0,2) = v(0,1); // b
							image1.draw_line(x, y ,x + (int) flowImage(x,y,0,0),y + (int) flowImage(x,y,0,2), white);
						}
					}
				}
			}

			flowPyramid.push_back(flowImage);
			iteration++;
			// std::cout << iteration << std::endl;
			image1.display();
		}
	}