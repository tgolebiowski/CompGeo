#include "QuasiCrystal.h"


QuasiCrystal::QuasiCrystal(int waveCount, float period)
{
	this->waveCount = waveCount;
	this->period = period;

	for(int i = 0; i < waveCount; i++)
	{
		waveAngles[i] = ((double)i) * ((3.1415926f * 2.0f) / waveCount) + (3.1415926f/2.0f);
	}

	offset = (((float)(rand() % 360)) / 360.0f) * (3.1415926f * 2.0f);
}


QuasiCrystal::~QuasiCrystal(void)
{
}

double QuasiCrystal::get(double x, double y, bool contrastOn)
{
	double total = 0.0f;
	for(int i = 0; i < waveCount; i++)
	{
		double theta = (2.0 * 3.1415926) * (x / (period / std::sin(waveAngles[i])));
		double yOffset = (2.0 * 3.1415926) * (y / (period / std::cos(waveAngles[i])));

		double value = std::cos(theta + yOffset + offset);
		value = (value + 1.0)/2.0;

		//if(value < 0.9f && value >= 0.5f) value = 1.0;
		//if(value > 0.1 && value <= 0.5f) value = 0.0;

		total += value;
	}
	total /= waveCount;

	if(contrastOn)
	{
		if(total <= 0.7) total = 0.0;
		else total = 1.0;
	}

	return total;
}

void QuasiCrystal::generatePointSet(double gridSize, double cellSize, float badPointBias, std::vector<Vector2>* set, bool addJitter, Vector2 offset)
{
	float badPointDistBias = badPointBias;
	int setDoubleCheck = 36;
	int count = gridSize / cellSize;

	struct Bounds
	{
		int start;
		int end;
	};
	std::vector<Bounds> lastRow;
	std::vector<Bounds> currentRow;

	for(int y = 0; y < count; y++)
	{
		bool regionStartedFlag = false;
		for(int x = 0; x < count; x++)
		{
			double val = get((double)x + offset.x, (double)y + offset.y);
			if(val > 0.9 && !regionStartedFlag)
			{
				currentRow.push_back(Bounds());
				currentRow.back().start = x;
				regionStartedFlag = true;
			}
			else if(val < 0.9 && regionStartedFlag)
			{
				currentRow.back().end = x;
				regionStartedFlag = false;
			}
		}

		if(!lastRow.empty() && !currentRow.empty())
		{
			int lastRowIndex = 0;
			int currentRowIndex = 0;

			while(lastRowIndex < lastRow.size() && currentRowIndex < currentRow.size())
			{
				Bounds lastBounds = lastRow.at(lastRowIndex);
				Bounds currentBounds = currentRow.at(currentRowIndex);

				if((currentBounds.end <= lastBounds.end && currentBounds.end >= lastBounds.start) ||
					(currentBounds.start <= lastBounds.end && currentBounds.start >= lastBounds.start))
				{
					//overlap between regions confirmed

					int lastBoundsMiddle = (lastBounds.start + lastBounds.end) / 2;
					int currentBoundsMiddle = (currentBounds.start + currentBounds.end) / 2;

					int lastBoundsLength = lastBounds.end - lastBounds.start;
					int currentBoundsLength = currentBounds.end - currentBounds.start;

					if(currentBoundsLength < lastBoundsLength)
					{
						Vector2 position = Vector2(lastBoundsMiddle, y);

						int setIndex = set->size() - 1;
						bool badPointFlag = false;
						while(setIndex > 0 && setIndex > set->size() - setDoubleCheck)
						{
							float dist = (position - set->at(setIndex)).length();
							if(dist < badPointDistBias)
							{
								badPointFlag = true;
								break;
							}

							setIndex--;
						}

						if(!badPointFlag)
						{
							if(addJitter)
							{
								Vector2 jitter = Vector2(((float)(rand()%1000))/1000.0f,((float)(rand()%1000))/1000.0f);
								jitter -= Vector2(0.5f, 0.5f);
								position += jitter;
							}
							set->push_back(position);
						}
					}

					lastRowIndex++;
					currentRowIndex++;
				}
				else
				{
					if(currentBounds.end <= lastBounds.end && currentRowIndex <= currentRow.size() - 1) currentRowIndex++;
					else lastRowIndex++;
				}
			}
		}
		
		lastRow.clear();
		for(int i = 0; i < currentRow.size(); i++)
		{
			lastRow.push_back(currentRow.at(i));
		}
		currentRow.clear();
	}
}
