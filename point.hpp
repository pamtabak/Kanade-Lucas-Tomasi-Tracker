typedef struct point
{
	double x;
	double y;
	bool isValid;
} point;

class ChosenPoint 
{
public:	
	ChosenPoint()
	{

	}

	~ChosenPoint()
	{
	}

	void setNumberOfFrames(int numberOfFrames)
	{
		this->flow = new point[numberOfFrames];

		for (int i = 0; i < numberOfFrames; i++)
		{
			point f;
			f.x = 0.0;
			f.y = 0.0;
			this->flow[i] = f;
		}
	}

	void setPoint(int x, int y, bool isValid)
	{
		this->pt.x = (double) x;
		this->pt.y = (double) y;
		this->pt.isValid = isValid;
	}

	point getPoint()
	{
		return this->pt;
	}

	point* getFlow()
	{
		return this->flow;
	}

	void setFlow(double &xFlow, double &yFlow, int &position)
	{
		this->flow[position].x = xFlow;
		this->flow[position].y = yFlow;
	}

	void updateFlow(int position)
	{
		this->flow[position].x = 2*this->flow[position].x;
		this->flow[position].y = 2*this->flow[position].y;
	}

	void updateFlow(double xFlow, double yFlow, int position)
	{
		this->flow[position].x = xFlow;
		this->flow[position].y = yFlow;
	}

	point pt;
	point *flow;
};