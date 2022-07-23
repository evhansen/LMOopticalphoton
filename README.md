
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
(Including the licences and such in individual files drives me mad so here you go. Based on some example.)



EM is encapsulatory material, e.g. aluminium reflectors.
LD is light detector.
ARC is antireflective coating.

ROOT files saved to data/
Flags set by files in options/ (I know, I know, I know, I know)


Main issue: muons getting stuck in crystal/ LD boundaries (when G4OpticalPhysics is registered -- fixed? OOM still kills with certain time segments.).



TODO: fix all of the derp flags .....

TODO: change "lmo" to "crys" everywhere ... teo2 is sad.

    char prng; // if 'c', C PRNG, otherwise G4 PRNG

    char lmos; 	// if e, LMO exists
   		// if s, LMO exists and is a sensitive detector 
		// else, not

TODO: float crysd; // cubic crystal dimension (in mm)

    float ldg;	// value of LD gap! (in mm)
    float emg;	// value of EM gap! (in mm)

    char arc; // whether the antireflective coating is being used
    float arct; // thickness of antireflective coating (in mm)

    char face1; char face2; char face3; char face4; char face5; char face6;
    // l is an LD, and e is an EM, else nothing

    char crysmat; // t is teo2, else lmo
    char ldmat; // s is Si, else Ge
    char emat; // c is Cu, m is mysterious (ooOoOOo),  else Al
    char arcmat; // switch sio2, sio, and nb2o5 





note, for troubleshooting:
	options/misc
		b: beta
		g: gamma
		else: muon

	options/gamma
		1: 1.17 MeV
		2: 1.33 MeV
		else: both



note, iMom is always primary's (seconds get 0,0,0)

TODO: fix scintillation, cerenkov
