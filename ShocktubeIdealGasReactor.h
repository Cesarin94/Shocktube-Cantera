/**
 *  @file ShocktubeIdealGasReactor.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_SHOCKTUBEIDEALGASREACTOR_H
#define CT_SHOCKTUBEIDEALGASREACTOR_H

#include "Reactor.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera
{

/**
 * Class ShocktubeIdealGasReactor is a class to simulate a batch reactor
 * optimized for ideal gases. More specifically, a gas mixture is excited by 
 * an incident shockwave to raise its temperature for reactions to be observed.
 */
class ShocktubeIdealGasReactor : public Reactor
{
public:
    ShocktubeIdealGasReactor() {}
	ShockProp shock;
	ShockProp *m_shock;
	m_shock = &shock;

    virtual int type() const {
        return ShocktubeIdealGasReactorType;
    }

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getInitialConditions(doublereal t0, size_t leny,
                                      doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);

    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    virtual void updateState(doublereal* y);

    virtual size_t componentIndex(const std::string& nm) const;


protected:
    vector_fp m_uk; //!< Species molar internal energies
	vector_fp m_hk; //!< Species molar enthalpies
};

}

#endif