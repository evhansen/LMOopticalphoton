
#include "PrimaryGeneratorAction.hh"


static char readExtFile(const char* filename);
void appendP(std::string filename, std::string line);
static std::string readExtFileLine(const char* filename);





//Rejection sampling

void PrimaryGeneratorAction::Cobalt60(float ConstituentPos[],float ConstituentDir[],G4Event* anEvent,G4ParticleGun*PG,char swch){
	//99.88% of the radiation will be a 0.31 MeV e-, 1.17 MeV g, and 1.33 MeV g.
	
	static std::string fname = "./data/" + readExtFileLine("./options/name");

	if(swch == '1'){	
		PG->SetParticleEnergy(1.17*MeV);
		PG->SetParticleDefinition(G4Gamma::Definition());
		PG->SetParticlePosition(G4ThreeVector(ConstituentPos[3]*m, ConstituentPos[4]*m, ConstituentPos[5]*m));
		PG->SetParticleMomentumDirection(G4ThreeVector(ConstituentDir[3], ConstituentDir[4], ConstituentDir[5]));
			
		PG->GeneratePrimaryVertex(anEvent);


		appendP(fname, 
	std::to_string(anEvent->GetEventID())  
	+ "," + std::to_string(ConstituentDir[3]) 
	+ "," + std::to_string(ConstituentDir[4]) 
	+ "," + std::to_string(ConstituentDir[5])
	+ "," + std::to_string(ConstituentPos[3]) 
	+ "," + std::to_string(ConstituentPos[4]) 
	+ "," + std::to_string(ConstituentPos[5]));

		
	} else if (swch == '2'){
		PG->SetParticleEnergy(1.33*MeV);
		PG->SetParticleDefinition(G4Gamma::Definition());
		PG->SetParticlePosition(G4ThreeVector(ConstituentPos[6]*m, ConstituentPos[7]*m, ConstituentPos[8]*m));
		PG->SetParticleMomentumDirection(G4ThreeVector(ConstituentDir[6], ConstituentDir[7], ConstituentDir[8]));
			
		PG->GeneratePrimaryVertex(anEvent);


		appendP(fname, 
	std::to_string(anEvent->GetEventID()) 
	+ "," + std::to_string(ConstituentDir[6]) 
	+ "," + std::to_string(ConstituentDir[7]) 
	+ "," + std::to_string(ConstituentDir[8])
	+ "," + std::to_string(ConstituentPos[6]) 
	+ "," + std::to_string(ConstituentPos[7]) 
	+ "," + std::to_string(ConstituentPos[8]));


	} else { // beta

		PG->SetParticleEnergy(0.31*MeV);
		PG->SetParticleDefinition(G4Electron::Definition());
		PG->SetParticlePosition(G4ThreeVector(ConstituentPos[0]*m, ConstituentPos[1]*m, ConstituentPos[2]*m));
		PG->SetParticleMomentumDirection(G4ThreeVector(ConstituentDir[0], ConstituentDir[1], ConstituentDir[2]));
		
		PG->GeneratePrimaryVertex(anEvent);


		appendP(fname, 
	std::to_string(anEvent->GetEventID())
	+ "," + std::to_string(ConstituentDir[0]) 
	+ "," + std::to_string(ConstituentDir[1]) 
	+ "," + std::to_string(ConstituentDir[2])
	+ "," + std::to_string(ConstituentPos[0]) 
	+ "," + std::to_string(ConstituentPos[1]) 
	+ "," + std::to_string(ConstituentPos[2]));


	}

}




void PrimaryGeneratorAction::CosmogenicMuons(float ConstituentPos[],float ConstituentDir[],G4Event* anEvent,G4ParticleGun*PG,char swch){
	//The muons will be 0.2-1000 GeV
	//ratio of m+/m- ~ 1.268 +/- (0.008 + 0.0002 p) where p is momentum.
	// https://doi.org/10.1016/0370-2693(71)90741-6


	static std::default_random_engine eng; //meh
	std::discrete_distribution<> d {125,121,360,259,96,6030,224,270,524,537,685,471,337,292,7670,6470,5931,3754,1767,613,198,49,11}; //23 vals
	
	int val = d(eng);
	float a,b,energy;

	switch(d(eng)){ // U(a,b) <= a + x(b-a) :U(0,1)
		case 0:
			a=0.2;b=0.5;break;
		case 1:
			a=0.5;b=0.8;break;
		case 2:
			a=0.8;b=2.0;break;
		case 3:
			a=2.0;b=4.0;break;
		case 4:
			a=4.0;b=6.0;break;
		case 5:
			a=0.98;b=1.24;break;
		case 6:
			a=0.2;b=0.4;break;
		case 7:
			a=0.4;b=0.6;break;
		case 8:
			a=0.6;b=1.0;break;
		case 9:
			a=1.0;b=1.5;break;
		case 10:
			a=1.5;b=2.5;break;
		case 11:
			a=2.5;b=4.0;break;
		case 12:
			a=4.0;b=6.0;break;
		case 13:
			a=6.0;b=10.0;break;
		case 14:
			a=10;b=13;break;
		case 15:
			a=13;b=17;break;
		case 16:
			a=17;b=25;break;
		case 17:
			a=25;b=40;break;
		case 18:
			a=40;b=70;break;
		case 19:
			a=70;b=128;break;
		case 20:
			a=128;b=250;break;
		case 21:
			a=250;b=450;break;
		case 22:
			a=450;b=1000;break;
		default:
			a=0.2;b=1000;
	}

	energy = (a + (((float)rand()/RAND_MAX)*(b - a)));
	PG->SetParticleEnergy( energy*GeV );

	
	double ratio,c; 
	
	//m+:m- - 1.268 +/- (0.008 + 0.0002 p)  
	c = 0.008 + 0.0002 * energy;
	ratio = 1.268 + (-c + (((float)rand()/RAND_MAX)*(c - (-c)))); // 1.268 + U(-c,c) 

	std::discrete_distribution<> r {ratio,1};
	
	if(r(eng) == 0){
		PG->SetParticleDefinition(G4MuonPlus::Definition());
	} else {
		PG->SetParticleDefinition(G4MuonMinus::Definition());
	}
	
	PG->SetParticlePosition(G4ThreeVector(
			ConstituentPos[0]*m, 
			ConstituentPos[1]*m, 
			ConstituentPos[2]*m));
	PG->SetParticleMomentumDirection(G4ThreeVector(
			ConstituentDir[0], 
			ConstituentDir[1], 
			ConstituentDir[2]));
	
	
	PG->GeneratePrimaryVertex(anEvent);


	static std::string fname = "./data/" + readExtFileLine("./options/name");

	appendP(fname,	
	std::to_string(anEvent->GetEventID()) 
	+ "," + std::to_string(ConstituentDir[0]) 
	+ "," + std::to_string(ConstituentDir[1]) 
	+ "," + std::to_string(ConstituentDir[2])
	+ "," + std::to_string(ConstituentPos[0]) 
	+ "," + std::to_string(ConstituentPos[1]) 
	+ "," + std::to_string(ConstituentPos[2]));

}






// tot: number of decays
// num: number of constituent particles .. YES I KNOW, DUMB REEEEEEEEEE

void PrimaryGeneratorAction::RectanglePlateSource(G4int tot,G4int num,G4double PSHalfSizes[],G4double PSPos[],void (*Process)(float[],float[],G4Event*,G4ParticleGun*,char),G4Event* anEvent,char swch){
	
	float ConstituentPos[num*3];
	float ConstituentDir[num*3];
	
	for(int i = 0; i < tot; i++){
		
		// rand in (0,1) or [0,1]
		// need affine ~U(0,1) => ~U(-halfsize,halfsize)
		// for x\in (0,1) 
		// => y = -halfsize + x (halfsize - (-halfsize)) 
		if(prng == 'c'){ 
		
			for(int j = 0; j < num*3; j += 3){
				
				ConstituentPos[j] = PSPos[0] + -PSHalfSizes[0] 
					+ ((float)rand()/RAND_MAX) * (PSHalfSizes[0] - (-PSHalfSizes[0]));
				
				ConstituentPos[j+1] = PSPos[1] + -PSHalfSizes[1] 
					+ ((float)rand()/RAND_MAX) * (PSHalfSizes[1] - (-PSHalfSizes[1]));
				
				ConstituentPos[j+2] =  PSPos[2];
				
				ConstituentDir[j] = 0;
				ConstituentDir[j+1] = 0;				
				ConstituentDir[j+2] = -1;

			}
			
		} else {
			
			for(int j = 0; j < num*3; j += 3){
				
				ConstituentPos[j] = PSPos[0] + -PSHalfSizes[0] 
					+ G4UniformRand() * (PSHalfSizes[0] - (-PSHalfSizes[0]));
				
				ConstituentPos[j+1] = PSPos[1] + -PSHalfSizes[1] 
					+ G4UniformRand() * (PSHalfSizes[1] - (-PSHalfSizes[1]));
				
				ConstituentPos[j+2] = PSPos[2];
				
				ConstituentDir[j] = 0;
				ConstituentDir[j+1] = 0;
				ConstituentDir[j+2] = -1;
				
			}
		}

		Process(ConstituentPos,ConstituentDir,anEvent,fParticleSource,swch); // num us encoded in the Process
	}
	
}





void PrimaryGeneratorAction::CircularPlateSource(G4int tot,G4int num,G4double Radius,G4double PSPos[],void (*Process)(float[],float[],G4Event*,G4ParticleGun*,char),G4Event* anEvent,char swch){
	
	float ConstituentPos[num*3];
	float ConstituentDir[num*3];
	float vx,vy,vz;
	
	for(int i = 0; i < tot; i++){
		
		// rand in (0,1) or [0,1]
		if(prng == 'c'){ 
		
			for(int j = 0; j < num*3; j += 3){
				
				ConstituentPos[j] = PSPos[0] + (Radius * std::sqrt(((float)rand()/RAND_MAX))) * std::cos(((float)rand()/RAND_MAX) * CLHEP::twopi);
				
				ConstituentPos[j+1] = PSPos[1] + (Radius * std::sqrt(((float)rand()/RAND_MAX))) * std::sin(((float)rand()/RAND_MAX) * CLHEP::twopi);
				
				ConstituentPos[j+2] = PSPos[2]; 
				 
				
				//want U(-1,1) x,y:   -1 + x (1 - (-1))   :   -1 + 2x
				vx = -1 + 2*((float)rand()/RAND_MAX);
				vy = -1 + 2*((float)rand()/RAND_MAX);

				//U(0,1) ~ vz = ((float)rand()/RAND_MAX);
				vz = ((float)rand()/RAND_MAX);
				
				ConstituentDir[j] = vx / std::sqrt((vx*vx)+(vy*vy)+(vz*vz));
				
				ConstituentDir[j+1] = vy / std::sqrt((vx*vx)+(vy*vy)+(vz*vz));
				
				ConstituentDir[j+2] = vz / std::sqrt((vx*vx)+(vy*vy)+(vz*vz));
			}
			
		} else {
			
			for(int j = 0; j < num*3; j += 3){
				
				ConstituentPos[j] = PSPos[0] + (Radius * std::sqrt(G4UniformRand())) * std::cos(G4UniformRand() * CLHEP::twopi);
				
				ConstituentPos[j+1] = PSPos[1] + (Radius * std::sqrt(G4UniformRand())) * std::sin(G4UniformRand() * CLHEP::twopi);
				
				ConstituentPos[j+2] = PSPos[2]; 
				
				vx = -1 + 2*G4UniformRand();
				vy = -1 + 2*G4UniformRand();
			
				vz = G4UniformRand();
	
				ConstituentDir[j] = vx / std::sqrt((vx*vx)+(vy*vy)+(vz*vz));
				
				ConstituentDir[j+1] = vy / std::sqrt((vx*vx)+(vy*vy)+(vz*vz));
				
				ConstituentDir[j+2] = vz / std::sqrt((vx*vx)+(vy*vy)+(vz*vz));
				
			}
		}
		
		Process(ConstituentPos,ConstituentDir,anEvent,fParticleSource,swch); // num is encoded in the Process
	}
}






char readExtFile(const char* filename){
	FILE *p = fopen(filename,"r");
	if(p == NULL){return '0';}
	return fgetc(p);
}

void appendP(std::string filename, std::string line){
	std::ofstream p;
	p.open(filename.c_str(),std::ios_base::app);
	p << line << "\n";
	p.close();
}

std::string readExtFileLine(const char* filename){

	std::string out;
	std::ifstream p;

	p.open(filename);

	std::getline(p,out);

	return out;
}



