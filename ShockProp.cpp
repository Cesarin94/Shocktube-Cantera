/**
* @file ShockProp.cpp
*/

#include "cantera/shocktubeProp/ShockProp.h"

using namespace std;
namespace Cantera
{

void ShockProp::initializeValues(double area)
{
//Initialize the gasSpeed and labTime to default values of zero
//The cross-sectional area of the tube is set as given the input
gasSpeed = 0.0;
labTime = 0.0;
Area = area;
}

void ShockProp::setValues(double speed, double time)
{
labTime = time;
gasSpeed = speed;
}

double* ShockProp::getCurrentValues()
{
	static double r[3];
	r[0]=labTime;
	r[1]=gasSpeed;
	r[2]=Area;
	return r;
}

};

