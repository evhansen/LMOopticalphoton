


Note: hardware troubles, so the repo is _not_ in a polished state currently.

In order to swap between projects, need to manually move the DetectorConstructions and PrimaryGeneratorActions into src/ (and move the other ones out) for right now.

Compiled with geant4.10.7.p03 (not yet tested with "geant4.10.7.03").

	TODO: add compiler options!



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



Some terminology:

	EM is encapsulatory material, e.g. aluminium reflectors.
	
	LD is light detector.
	
	ARC is antireflective coating.


Some data locations:
	
	ROOT files saved to data/.
	
	Flags set by files in options/ (to be improved at a later date!).
	
	Material properties are potentially changed by files in options/mpt/ (note BOMs and strange behaviour generally -- better to hardcode your data).



TODOs / fixes:
	
	PbWO4 birefringence
	
	Change "lmo" to "crys" everywhere ... 
	
	Add float crysd; // cubic crystal dimension (in mm)
	
	Fix rootfile data (e.g. iMom is always primary's (seconds get 0,0,0).).
	
	Better modularization 
	
	Inaccuracy with OP capture in TES module?
	
	TES module fibre geometry (snake horizontal fibre)
	
	Add more photos to this README. :)


/////////////////
// **************
// Option files for the LMO runs
// **************
/////////////////


// DetectorConstruction //

    in /options/lmos
		
		e: LMO exists
		
		s: LMO exists and is a sensitive detector 
		
		else: no crystal


    in /options/ldg (float)
		
		value of LD gap! (in mm)


    in /options/emg (float)
		
		value of EM gap! (in mm)


    in /options/arct (float)
		
		thickness of antireflective coating (in mm)


	in /options/face1,

	/options/face2,

	/options/face3,

	/options/face4,

	/options/face5,

	/options/face6
		
		l: LD
		
		e: EM
		
		else: nothing on that side


	in /options/crysmat
		
		t: TeO2
		
		else: LMO


	in /options/ldmat 
		
		s: Si
		
		else: Ge


	in /options/emat 
		
		c: Cu
		
		m: mysterious (ooOoOOo, meant for variable rindices)
		
		else: Al
	
	Note that EMs are automatically marked sensitive.


	in /options/arcmat
		
		2: sio2
		
		1: sio
		
		n: nb2o5
		
		else: none
	
	Note that ARC are automatically marked sensitive.
	
	
// PrimaryGeneratorAction //

	in /options/prng
		
		c: C PRNG
		
		else: G4 PRNG

	in /options/misc:
		
		b: beta (0.31 MeV)
		
		1: gamma (1.17 MeV)
		
		2: gamma (1.33 MeV)
		
		t: thorium (2615 keV gammas)
		
		else: muons (5 GeV)



/////////////////
// **************
// Option files for the Cocalib runs
// **************
/////////////////

None -- geometry is fixed.



/////////////////
// **************
// Option files for the TES module runs
// **************
/////////////////

// DetectorConstruction //

	in /options/TES_topgap (float)
		
		0: the default 10 mm gap 
		
		all else: 10+ that value gap in mm. (Negative values allowed :))


	in /options/TES_arc
		
		1: ARC included
		
		else: none


	in /options/TES_ref
		
		1: reflector on the underside of the top slab included
		
		else: none


	in /options/TES_can
		
		1: Cu can encapsulating the entire system included
		
		else: none


	in /options/TES_canref
		
		1: 70% can reflectivity included
		
		else: none


// PrimaryGeneratorAction //

	in /options/TES_escape
		
		a: OPs escape from the wires
		
		else: no OPs escape


	in /options/TES_opticalfibre
		
		h: horizontal fibres (at the moment, separate horizontal ones)
![Image of horizontal fibre.](https://github.com/SporadicAnyonsInhabitGregariousAsylums/LMOopticalphoton/blob/ModularProj/imgs/f2.JPG)
		
		v: vertical fibres
![Image of vertical fibre.](https://github.com/SporadicAnyonsInhabitGregariousAsylums/LMOopticalphoton/blob/ModularProj/imgs/f1.JPG)
		
		else: no fibre at all, i.e. no generation!
	
	Note that geometry details are in the documentation.

	Note that fibres all have diffusion implemented.
	
	![Diffusion length and power percentage image.](https://github.com/SporadicAnyonsInhabitGregariousAsylums/LMOopticalphoton/blob/ModularProj/imgs/diffusionfibreplot.JPG)


/////////////////
// **************
// Description of relevant homemade functions and structures
// **************
/////////////////


There are three main includes in this vein, that aren't G4 classes:
	
	DefMats.hh -- defines the materials for the runs
	
	DefGeos.hh -- defines the geometries for the runs
	
	DefSrcs.hh -- defines the sources for the runs
	
	
DefMats isn't used in a special way really, just define the functions you want in the appropriate header (probably in DetectorConstruction.hh) and you can throw your function into this file.
	
	
DefGeos isn't used in a special way either. Similar to DefMats.
	
	
* DefSrcs _is_ used in a special way. You do not have to use it or its functions though. :)
	
	
First, we define the "process" (the particles generated) then we define the "source" (or we use one of the predefined functions available) and that's it! Let's start with the "source" function.
		
		
The "source" function will look (something) like:
				
		void PrimaryGeneratorAction::SOURCE_NAME(
					
			G4int tot, G4int num,
					
			G4double PSHalfSizes[],G4double PSPos[], 
					
			void(*Process)(float[], float[], G4Event*, G4ParticleGun*, char),
					
			G4Event* anEvent, char swch)
	
	
SOURCE_NAME is taken from a rectangle source but with a little bit of modification, we can create lots of other geometries (see DefSrcs.hh for an example of a circular source). 


For such, we pass PSHalfSizes to delineate the halfsizes for the source (so half the lengths of the sides of the abstract rectangle in which we will generate particles). 


We pass PSPos to set the coordinate of the centre of our rectangle.


Note that void(*Process)(float[], float[], G4Event*, G4ParticleGun*, char) looks a lot like a function -- it is! We pass a function pointer to our "process" function (PROCESS_NAME below) to be used in SOURCE_NAME. So we can swap out PROCESS_NAME easily to change the type of particles generated, how they're generated, etc. :)


We pass anEvent from our main PrimaryActionGenerator function.


We can optionally pass a character: swch. That's for if we have a option file read somewhere in PrimaryActionGenerator that influences either SOURCE_NAME or PROCESS_NAME. It's absolutely not necessary.


Now most importantly we have the two G4ints: tot and num. tot is the "total number of decays" for that event and num is the "total number of particles [per decay]". We use this to compute the initial positions and momenta for each particle that we generate -- we have num * 3 positions / momenta (so num particles, each with an x, y, z position / momenta) and we iterate over the computation of said positons and momenta tot times. Effectively we generate tot * num number of particles, just allowing for num to have different proportions of different particles!
			
		
The "process" function will look like:
				
		void PrimaryGeneratorAction::PROCESS_NAME(
					
			float ConstituentPos[], float ConstituentDir[],
					
			G4Event* anEvent, G4ParticleGun*PG, char swch)
					
					
Here ConstituentPos[],  ConstituentDir[], anEvent, PG, and swch are all passed to PROCESS_NAME by you and the source function. 
	
	
More importantly here are ConstituentPos and ConstituentDir -- these are the starting positions and momenta of the particle(s) generated. They are computed in SOURCE_NAME then passed to PROCESS_NAME to use when generating. 

			
Mind that you do _not_ have to use these files at all! They're only there to make the files shorter, more readable, and make it easier to build simulations!
	


See the images below for graphical representations of the functions:

![Visual representation of some classes.](https://github.com/SporadicAnyonsInhabitGregariousAsylums/LMOopticalphoton/blob/ModularProj/imgs/classes.jpeg)
