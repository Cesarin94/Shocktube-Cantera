/**
*@file ShockProp.h
*/


#ifndef SHOCKPROP_H
#define SHOCKPROP_H

namespace Cantera
{
class ShockProp 
{
	double gasSpeed,labTime,Area;
public:
	void initializeValues(double area);
	void setValues(double speed, double time);
	double* getCurrentValues();

};

}
#endif
