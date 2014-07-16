/**
 *  @file ShocktubeIdealGasReactor.cpp A zero-dimensional reactor
 */

#include "cantera/shocktubeProp/ShockProp.h"
#include "cantera/zeroD/ShocktubeIdealGasReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/zeroD/ReactorNet.h"

#include <cfloat>
#include <cmath>

using namespace std;

namespace Cantera
{

void ShocktubeIdealGasReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.eosType() != cIdealGas) {
        throw CanteraError("ShocktubeIdealGasReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void ShocktubeIdealGasReactor::getInitialConditions(double t0, size_t leny, double* y)
{
    m_init = true;
    if (m_thermo == 0) {
        cout << "Error: reactor is empty." << endl;
        return;
    }
    m_thermo->restoreState(m_state);
	//ADDED BY ME:
	m_rho0 = m_thermo->density();
	//
    // set the first component to the total mass
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_mass;

    // set the second component to the total volume
    y[1] = m_vol;

    // Set the third component to the temperature
    y[2] = m_thermo->temperature();

	// Set the fourth component to the density
	y[3] = m_thermo->density();

	// Set the fifth component to the gas velocity, by default 0.0
	y[4] = 0.0;

	// Set the sixth component to the lab time, by default 0.0
	y[5] = 0.0;
    // set components y+6 ... y+K+5 to the mass fractions of each species
    m_thermo->getMassFractions(y+6);

    // set the remaining components to the surface species
    // coverages on the walls
    size_t loc = m_nsp + 6;
    SurfPhase* surf;
    for (size_t m = 0; m < m_nwalls; m++) {
        surf = m_wall[m]->surface(m_lr[m]);
        if (surf) {
            m_wall[m]->getCoverages(m_lr[m], y + loc);
            loc += surf->nSpecies();
        }
    }
}

void ShocktubeIdealGasReactor::initialize(doublereal t0)
{
    m_thermo->restoreState(m_state);
    m_sdot.resize(m_nsp, 0.0);
    m_wdot.resize(m_nsp, 0.0);
    m_uk.resize(m_nsp, 0.0);
    m_nv = m_nsp + 6;
    for (size_t w = 0; w < m_nwalls; w++)
        if (m_wall[w]->surface(m_lr[w])) {
            m_nv += m_wall[w]->surface(m_lr[w])->nSpecies();
        }

    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();

    size_t nt = 0, maxnt = 0;
    for (size_t m = 0; m < m_nwalls; m++) {
        m_wall[m]->initialize();
        if (m_wall[m]->kinetics(m_lr[m])) {
            nt = m_wall[m]->kinetics(m_lr[m])->nTotalSpecies();
            if (nt > maxnt) {
                maxnt = nt;
            }
            if (m_wall[m]->kinetics(m_lr[m])) {
                if (&m_kin->thermo(0) !=
                        &m_wall[m]->kinetics(m_lr[m])->thermo(0)) {
                    throw CanteraError("ShocktubeIdealGasReactor::initialize",
                                       "First phase of all kinetics managers must be"
                                       " the gas.");
                }
            }
        }
    }
    m_work.resize(maxnt);
    std::sort(m_pnum.begin(), m_pnum.end());
    m_init = true;
}

void ShocktubeIdealGasReactor::updateState(doublereal* y)
{
    for (size_t i = 0; i < m_nv; i++) {
        AssertFinite(y[i], "ShocktubeIdealGasReactor::updateState",
                     "y[" + int2str(i) + "] is not finite");
    }

    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, y[3] the density, y[4] the gas velocity, y[5] the lab time, [6...K+6] are the mass fractions of each species,
    // and [K+6...] are the coverages of surface species on each wall.
    m_mass = y[0];
    m_vol = y[1];
	//ADDED BY ME:
	m_shock->gasSpeed = y[4];
	m_shock->labTime = y[5];
	//
    m_thermo->setMassFractions_NoNorm(y+6);
    m_thermo->setState_TR(y[2], y[3]);

    size_t loc = m_nsp + 6;
    SurfPhase* surf;
    for (size_t m = 0; m < m_nwalls; m++) {
        surf = m_wall[m]->surface(m_lr[m]);
        if (surf) {
            m_wall[m]->setCoverages(m_lr[m], y+loc);
            loc += surf->nSpecies();
        }
    }

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}



void ShocktubeIdealGasReactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    m_thermo->restoreState(m_state);
//---------------------------------------------------------------------------------
    // process sensitivity parameters                                             |
    if (params) {                                                                 
        size_t npar = m_pnum.size();
        for (size_t n = 0; n < npar; n++) {
            double mult = m_kin->multiplier(m_pnum[n]);
            m_kin->setMultiplier(m_pnum[n], mult*params[n]);
        }
        size_t ploc = npar;
        for (size_t m = 0; m < m_nwalls; m++) {
            if (m_nsens_wall[m] > 0) {
                m_wall[m]->setSensitivityParameters(m_lr[m], params + ploc);
                ploc += m_nsens_wall[m];
            }
        }
    } //                                                                           |
//----------------------------------------------------------------------------------
    m_vdot = 0.0; // change in volume. (Array Index = 1) 
    m_Q    = 0.0; // Heat exchange term.
    double dTdt = 0.0; // m * c_v * dT/dt (Array Index = 2)
    double dmdt = 0.0; // dm/dt (gas phase)  (Array index = 0)
    double* dYdt = ydot + 6; // changed index from 3 to 6. (Pointer Because of it being array.)
	double dR/dt = 0.0; // R here represents the density. (Array index = 3)
	double dvdt = 0.0; //v here represents the gas velocity (m/s) (Array index = 4)
	double dtdt = 0.0; // here represents the lab time (secs) (Array index = 5)

    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
	// ADDED BY ME:
	m_thermo->getPartialMolarIntEnergies(&m_hk[0]); // J/kmol
//----------------------------------------------------------------------------------
    // compute wall terms                                                          |
    size_t loc = m_nsp+6; // Changed Index from 3 to 6.
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    for (size_t i = 0; i < m_nwalls; i++) {
        int lr = 1 - 2*m_lr[i];
        double vdot = lr*m_wall[i]->vdot(time);
        m_vdot += vdot;
        m_Q += lr*m_wall[i]->Q(time);
        Kinetics* kin = m_wall[i]->kinetics(m_lr[i]);
        SurfPhase* surf = m_wall[i]->surface(m_lr[i]);
        if (surf && kin) {
            double rs0 = 1.0/surf->siteDensity();
            size_t nk = surf->nSpecies();
            double sum = 0.0;
            surf->setTemperature(m_state[0]);
            m_wall[i]->syncCoverages(m_lr[i]);
            kin->getNetProductionRates(DATA_PTR(m_work));
            size_t ns = kin->surfacePhaseIndex();
            size_t surfloc = kin->kineticsSpeciesIndex(0,ns);
            for (size_t k = 1; k < nk; k++) {
                ydot[loc + k] = m_work[surfloc+k]*rs0*surf->size(k);
                sum -= ydot[loc + k];
            }
            ydot[loc] = sum;
            loc += nk;

            double wallarea = m_wall[i]->area();
            for (size_t k = 0; k < m_nsp; k++) {
                m_sdot[k] += m_work[k]*wallarea;
            }
        }
    } //                                                                          |
// --------------------------------------------------------------------------------
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    double mdot_surf = 0.0; // net mass flux from surfaces
    for (size_t k = 0; k < m_nsp; k++) {
        // production in gas phase and from surfaces
        dYdt[k] = (m_wdot[k] * m_vol + m_sdot[k]) * mw[k] / m_mass;
        mdot_surf += m_sdot[k] * mw[k];
    }
    dmdt += mdot_surf;

	// MAJOR CHANGES START HERE:
	// Add the computation for the time derivative of density (kg/m^3)
	for (size_t n = 0; n < m_nsp; n++) {
		dRdt += m_wdot[n]*( m_hk[n] - m_thermo->cp_mole() * m_thermo->temperature() )  /( m_thermo->cp_mass() * m_thermo->temperature() - m_thermo->cv_mole() * pow((m_shock->gasSpeed),2) / 8314.4621);
    // compression work and external heat transfer
    // mcvdTdt += - m_pressure * m_vdot - m_Q; 

	//Temperature contribution from non-summation terms:
	dTdt += pow(m_shock->gasSpeed,2)* dRdt / (m_thermo->density() * m_thermo->cp_mass());
	
    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        // mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol;
        // mcvdTdt -= m_sdot[n] * m_uk[n];
		dvdt += (-1 / m_thermo->density() ) * dRdt;
		dTdt -= m_wdot[n] * m_hk[n] / ( m_thermo->density() * m_thermo->cp_mass() ) ;
        // dilution by net surface mass flux
        dYdt[n] -= Y[n] * mdot_surf / m_mass;
    }
	//-----------------------------------------------------------------------------------------------
    // add terms for open system                                                                    |
    if (m_open) {
        // outlets
        for (size_t i = 0; i < m_nOutlets; i++) {
            double mdot_out = m_outlet[i]->massFlowRate(time);
            dmdt -= mdot_out; // mass flow out of system
            //mcvdTdt -= mdot_out * m_pressure * m_vol / m_mass; // flow work
        }

        // inlets
        for (size_t i = 0; i < m_nInlets; i++) {
            double mdot_in = m_inlet[i]->massFlowRate(time);
            dmdt += mdot_in; // mass flow into system
            //mcvdTdt += m_inlet[i]->enthalpy_mass() * mdot_in;
            for (size_t n = 0; n < m_nsp; n++) {
                double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
                // flow of species into system and dilution by other species
                dYdt[n] += (mdot_spec - mdot_in * Y[n]) / m_mass;

				// INCLUDE THE COMPUTATION FOR THE LAB TIME HERE:
				// Note that the Actual Formula is :
				// dt/dt = (rho0 * Area0) / (rho * Area)
				// But here the area is constant.
				dtdt += (m_rho0)/m_thermo->density();

                // In combintion with h_in*mdot_in, flow work plus thermal
                // energy carried with the species
                //mcvdTdt -= m_uk[n] / mw[n] * mdot_spec;
            }
        }
    }//                                                                                               |
	//-------------------------------------------------------------------------------------------------
    ydot[0] = dmdt;
    ydot[1] = m_vdot;
    if (m_energy) {
        ydot[2] = dTdt;
			//mcvdTdt / (m_mass * m_thermo->cv_mass());
    } else {
        ydot[2] = 0;
    }
	ydot[3] = dRdt;
	ydot[4] = dvdt;
	ydot[5] = dtdt;

    for (size_t i = 0; i < m_nv; i++) {
        AssertFinite(ydot[i], "ShocktubeIdealGasReactor::evalEqs",
                     "ydot[" + int2str(i) + "] is not finite");
    }

    // reset sensitivity parameters
    if (params) {
        size_t npar = m_pnum.size();
        for (size_t n = 0; n < npar; n++) {
            double mult = m_kin->multiplier(m_pnum[n]);
            m_kin->setMultiplier(m_pnum[n], mult/params[n]);
        }
        size_t ploc = npar;
        for (size_t m = 0; m < m_nwalls; m++) {
            if (m_nsens_wall[m] > 0) {
                m_wall[m]->resetSensitivityParameters(m_lr[m]);
                ploc += m_nsens_wall[m];
            }
        }
    }
}

size_t ShocktubeIdealGasReactor::componentIndex(const string& nm) const
{
    if (nm == "T") {
        return 2;
    } else {
        return Reactor::componentIndex(nm);
    }
}

}
