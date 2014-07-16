Shocktube-Cantera
=================

Cantera Version 2.1.1 modifications for Shocktube Ignition
This repository is for additional classes that can be added to the already existing Cantera-2.1.1 library. (For more information about Cantera please see the attached linked to the citation at the end of this README file.)
The two classes that are included in this repository are the ShocktubeIdealGasReactor and ShockProp as well as the necessary changes to the parent class ReactorBase(only minor changes were made to this already existing file in the Cantera 2.1.1 library to allow for the addition of the new derived class ShocktubeIdealGasReactor).

Descriptions of each file:

--ShockProp:
    This is a basic class that initializes, sets and gets the current values of the Shock Property variables that are necessary for the solution of the differential equations evaluated in the Shocktube IdealGasPhaseReactor Class. An instance of this class is called directly in the aforementioned reactor class.
    

--ShocktubeIdealGasReactor:
  This is a modified reactor class modeled after the already exisiting IdealGasPhaseReactor class but with changes made to include the differential equations necessary to describe the mass,volume,density, temperature, gas speed and lab time variations with time (zero dimensional).

Where to put these files inside Cantera:

1) ReactorBase.h and ShocktubeIdealGasPhase.h should be placed inside the "zeroD" include folder directory. Note: the ReactorBase class is still fully compatible with other ReactorTypes. It has simply been modified to include this new Reactortype.
2) ReactorBase.cpp and ShokctubeIdealGasPhase.cpp should be placed inside the src directory under the "zeroD" file.
3) Inside the src directory a folder named "ShocktubeProp" should be created to contain the ShockProp.cpp file.  Lastly, a "ShocktubeProp" folder should be included in the include directory of Cantera to contain the ShockProp.h file .



CITATIONS:

Dave Goodwin, Nicholas Malaya, Harry Moffat, and Raymond Speth. Cantera: An object-oriented software toolkit for chemical kinetics, thermodynamics, and transport processes. Version 2.1.1, available at https://code.google.com/p/cantera/
