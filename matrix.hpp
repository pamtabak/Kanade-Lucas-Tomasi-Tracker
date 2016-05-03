#include <iostream>
#include <vector>

using namespace std;

class Matrix
{
public:
	Matrix()
	{

	}
	~Matrix()
	{

	}

	vector<vector<double> > getTransposedMatrix (vector<vector<double> > matrix)
	{
		// Initializing return object
		vector<vector<double> > transposedMatrix;
		
		// Iterating over each column
		for (int j = 0; j < matrix[0].size(); j++)
		{
			// Iterating over each row
			vector<double> row;
			for (int i = 0; i < matrix.size(); i++)
			{
				row.push_back(matrix[i][j]);
			}	
			transposedMatrix.push_back(row);
		}

		return transposedMatrix;
	}

};