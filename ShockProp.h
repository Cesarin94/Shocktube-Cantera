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
	void initializeShockProp();
	void setShockProp();
	int* getCurrentShockProp();

};

}
#endif