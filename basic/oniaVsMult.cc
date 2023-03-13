/**
 *
 * SIMPLIFIED VERSION OF MACRO
 *
 * This macro generates events with PYTHIA, looks for quarkonia in the events,
 * and stores the information about the events and probes in TTree's that it writes into a root file.
 * The trees can be later used to determine e.g. the multiplicity dependence of J/psi production.
 *
 *
 *
 * General structure:
 *
 *
 *
 *
 *
 * Details:
 *
 * The code was written to be run on a computing farm (like GSI kronos)
 *
 *
 *
 *
 * To simulate only J/psi in the EE decay channel use the settings file settings_EEdecay.cmnd
 *
 *
 *
 *
 *
 * It takes arguments from the command line:
 *
 * argument 1:  number of events to produce (default: 10000)
 * argument 2:  whether or not to scale down the MB events (0=no, 1=yes) (default: 0)
 * argument 3:  whether to print out more information while running (default: 0)
 * argument 4:  the settings file to use (default: "settings.cmd")
 * argument 5:  the output filename (without the final .root) (default: "out")
 * argument 6:  an initial seed for the random number generator (default:1)
 *
 */


#include <time.h>
#include "TFile.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "Pythia8/Pythia.h"
// You need to include this to get access to the HIInfo object for HeavyIons.
#include "Pythia8/HeavyIons.h"

using namespace Pythia8;


/**
 *
 * This is the information stored for each Quarkonium in the tree:
 *  - the pdg code (see http://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf)
 *  - the pT and rapidity
 *  - the multiplicities in the towards ("Region1"), transverse ("Region2") and away ("Region3") region
 *    these multiplicities are counted either
 *    - for all etas
 *    - for |eta| < 1
 *    - in the V0A acceptance at positive or negative rapidity
 *    - in the V0C acceptance at positive or negative rapidity
 *
 */
struct Quarkonium {
    int pdg;
    float pt,y;
    int initialPdg;
    int decayChannel;
    unsigned short multRegion1;
    unsigned short multRegion2;
    unsigned short multRegion3;
    unsigned short multEta09Region1;
    unsigned short multEta09Region2;
    unsigned short multEta09Region3;
    unsigned short multEta1Region1;
    unsigned short multEta1Region2;
    unsigned short multEta1Region3;
    unsigned short multV0APosRegion1;
    unsigned short multV0APosRegion2;
    unsigned short multV0APosRegion3;
    unsigned short multV0ANegRegion1;
    unsigned short multV0ANegRegion2;
    unsigned short multV0ANegRegion3;
    unsigned short multV0CPosRegion1;
    unsigned short multV0CPosRegion2;
    unsigned short multV0CPosRegion3;
    unsigned short multV0CNegRegion1;
    unsigned short multV0CNegRegion2;
    unsigned short multV0CNegRegion3;
    Quarkonium():
        pdg(0),
        pt(0.),
        y(0.),
        initialPdg(0),
        decayChannel(0),
        multEta09Region1(0),
        multEta09Region2(0),
        multEta09Region3(0),
        multEta1Region1(0),
        multEta1Region2(0),
        multEta1Region3(0),
        multV0APosRegion1(0),
        multV0APosRegion2(0),
        multV0APosRegion3(0),
        multV0ANegRegion1(0),
        multV0ANegRegion2(0),
        multV0ANegRegion3(0),
        multV0CPosRegion1(0),
        multV0CPosRegion2(0),
        multV0CPosRegion3(0),
        multV0CNegRegion1(0),
        multV0CNegRegion2(0),
        multV0CNegRegion3(0)
    {};
};

struct Track {
    double eta;
    Track():
        eta(0.)
    {};
};

struct ChEvent {
    unsigned short multFull;
    TClonesArray* tracks;
    ChEvent():
        multFull(0),
        tracks(0x0)
    {
      tracks = new TClonesArray("Track", 100000);
    };
};


void traceBackQuarkonium(Quarkonium &found, unsigned short particle, Pythia &pythia);
bool isOnium(int pdg);

void fillRegions(unsigned short &multRegion1,
                 unsigned short &multEta09Region1,  unsigned short &multEta1Region1,
                 unsigned short &multV0APosRegion1, unsigned short &multV0ANegRegion1,
                 unsigned short &multV0CPosRegion1, unsigned short &multV0CNegRegion1,
                 unsigned short &multRegion2,
                 unsigned short &multEta09Region2,  unsigned short &multEta1Region2,
                 unsigned short &multV0APosRegion2, unsigned short &multV0ANegRegion2,
                 unsigned short &multV0CPosRegion2, unsigned short &multV0CNegRegion2,
                 unsigned short &multRegion3,
                 unsigned short &multEta09Region3,  unsigned short &multEta1Region3,
                 unsigned short &multV0APosRegion3, unsigned short &multV0ANegRegion3,
                 unsigned short &multV0CPosRegion3, unsigned short &multV0CNegRegion3,
                 double phi, double eta);
void fillRegionRandom(unsigned short &multRegion,
                      unsigned short &multEta09Region,  unsigned short &multEta1Region,
                      unsigned short &multV0APosRegion, unsigned short &multV0ANegRegion,
                      unsigned short &multV0CPosRegion, unsigned short &multV0CNegRegion,
                      double phi, double eta);
void fillRegionsWrtQuarkonium(Quarkonium &found, Pythia& pythia, double phi);

bool isPrimaryChargedALICE(unsigned short idx, Pythia &pythia);
bool isLongLived(unsigned int pdg);

bool inV0APosAcceptance(float eta);
bool inV0ANegAcceptance(float eta);
bool inV0CPosAcceptance(float eta);
bool inV0CNegAcceptance(float eta);

/**
 *
 * The main function.
 *
 * - processes the command line arguments
 * - starts parallel threads using openMP
 * - created the output trees
 * - initialized pythia
 * - runs the event loop
 *   - here it calculates the multplicities, looks for the hard probes, and calls helper function
 *     to get the multiplicities in the different regions, and information about the hard probes.
 *
 *
 */
int main(int argc, char** argv) {

    int nev            = argc > 1 ? std::stoi(argv[1]) : 10000;
    string settings    = argc > 2 ? argv[2]            : "settings.cmnd";
    string outfilename = argc > 3 ? argv[3]            : "out";
    int taskid         = argc > 4 ? std::stoi(argv[4]) : 1;

    // prepare output file
    outfilename += ".root";
    TFile* fout;
    fout = TFile::Open(outfilename.c_str(), "RECREATE");
    fout->cd();


    // prepare trees
    unsigned short multFull, multEta09, multEta1, multV0APos, multV0ANeg, multV0CPos, multV0CNeg, nMPI;
    unsigned char  type;
    unsigned short multRegionRnd;
    unsigned short multEta09RegionRnd;
    unsigned short multEta1RegionRnd;
    unsigned short multV0APosRegionRnd;
    unsigned short multV0ANegRegionRnd;
    unsigned short multV0CPosRegionRnd;
    unsigned short multV0CNegRegionRnd;
    double weight;
    unsigned short nPartTarg;
    unsigned short nAbsTarg;
    unsigned short nDiffTarg;

    // define the tree for minimum bias events
    TTree* eventTree = new TTree("eventTree", "event information");
    eventTree->Branch("multFull",            &multFull);
    eventTree->Branch("multEta09",           &multEta09);
    eventTree->Branch("multEta1",            &multEta1);
    eventTree->Branch("multV0APos",          &multV0APos);
    eventTree->Branch("multV0ANeg",          &multV0ANeg);
    eventTree->Branch("multV0CPos",          &multV0CPos);
    eventTree->Branch("multV0CNeg",          &multV0CNeg);
    eventTree->Branch("type",                &type);
    eventTree->Branch("nMPI",                &nMPI);
    eventTree->Branch("multRegionRnd",       &multRegionRnd);
    eventTree->Branch("multEta09RegionRnd",  &multEta09RegionRnd);
    eventTree->Branch("multEta1RegionRnd",   &multEta1RegionRnd);
    eventTree->Branch("multV0APosRegionRnd", &multV0APosRegionRnd);
    eventTree->Branch("multV0ANegRegionRnd", &multV0ANegRegionRnd);
    eventTree->Branch("multV0CPosRegionRnd", &multV0CPosRegionRnd);
    eventTree->Branch("multV0CNegRegionRnd", &multV0CNegRegionRnd);
    eventTree->Branch("weight",              &weight);     // weight of the current event.
    eventTree->Branch("nPartTarg",           &nPartTarg);  // nof interacting target nucleons
    eventTree->Branch("nAbsTarg",            &nAbsTarg);   // nof absorptively wounded target nucleons
    eventTree->Branch("nDiffTarg",           &nDiffTarg);  // nof diffrectively wounded target nucleons

    TTree* eventTestTree = new TTree("eventTestTree", "test event information");
    ChEvent* testEvent = new ChEvent();
    eventTestTree->Branch("testEvent", &testEvent, 16000, 99);

    TTree* oniumTree = new TTree("oniumTree", "Onia information");
    Quarkonium onium;
    oniumTree->Branch("multFull",           &multFull);
    oniumTree->Branch("multEta09",          &multEta09);
    oniumTree->Branch("multEta1",           &multEta1);
    oniumTree->Branch("multV0APos",         &multV0APos);
    oniumTree->Branch("multV0ANeg",         &multV0ANeg);
    oniumTree->Branch("multV0CPos",         &multV0CPos);
    oniumTree->Branch("multV0CNeg",         &multV0CNeg);
    oniumTree->Branch("nMPI",               &nMPI);
    oniumTree->Branch("type",               &type);                                                                 
    oniumTree->Branch("onium.pdg",          &onium.pdg);
    oniumTree->Branch("onium.decayChannel", &onium.decayChannel);
    oniumTree->Branch("onium.pt",           &onium.pt);
    oniumTree->Branch("onium.y",            &onium.y);
    oniumTree->Branch("onium.initialPdg",   &onium.initialPdg);
    oniumTree->Branch("multRegion1",        &onium.multRegion1 );
    oniumTree->Branch("multRegion2",        &onium.multRegion2);
    oniumTree->Branch("multRegion3",        &onium.multRegion3);
    oniumTree->Branch("multEta09Region1",   &onium.multEta09Region1);
    oniumTree->Branch("multEta09Region2",   &onium.multEta09Region2);
    oniumTree->Branch("multEta09Region3",   &onium.multEta09Region3);
    oniumTree->Branch("multEta1Region1",    &onium.multEta1Region1);
    oniumTree->Branch("multEta1Region2",    &onium.multEta1Region2);
    oniumTree->Branch("multEta1Region3",    &onium.multEta1Region3);
    oniumTree->Branch("multV0APosRegion1",  &onium.multV0APosRegion1);
    oniumTree->Branch("multV0APosRegion2",  &onium.multV0APosRegion2);
    oniumTree->Branch("multV0APosRegion3",  &onium.multV0APosRegion3);
    oniumTree->Branch("multV0ANegRegion1",  &onium.multV0ANegRegion1);
    oniumTree->Branch("multV0ANegRegion2",  &onium.multV0ANegRegion2);
    oniumTree->Branch("multV0ANegRegion3",  &onium.multV0ANegRegion3);
    oniumTree->Branch("multV0CPosRegion1",  &onium.multV0CPosRegion1);
    oniumTree->Branch("multV0CPosRegion2",  &onium.multV0CPosRegion2);
    oniumTree->Branch("multV0CPosRegion3",  &onium.multV0CPosRegion3);
    oniumTree->Branch("multV0CNegRegion1",  &onium.multV0CNegRegion1);
    oniumTree->Branch("multV0CNegRegion2",  &onium.multV0CNegRegion2);
    oniumTree->Branch("multV0CNegRegion3",  &onium.multV0CNegRegion3);


    // start pythia
    gRandom = new TRandom3();
    Pythia pythia;
    pythia.readFile(settings);

    pythia.readString("Random:setSeed = on");
    int seed = time(NULL) * taskid % 900000000;
    std::stringstream sstm;
    sstm << "Random:seed =" << seed;
    std::string seedString = sstm.str();
    pythia.readString(seedString);

    pythia.init();


    // Start the event loop
    for(int iev=0; iev<nev; iev++) {

        if(!pythia.next()) continue;
        multFull   = 0; // mult in eta +- inf.
        multEta09  = 0; // mult in eta +-0.9
        multEta1   = 0; // mult in eta +-1.0
        multV0APos = 0; // mult in V0APos acceptance
        multV0ANeg = 0; // mult in V0ANeg acceptance
        multV0CPos = 0; // mult in V0CPos acceptance
        multV0CNeg = 0; // mult in V0CNeg acceptance

        // multiplicity in a randomly chosen 120 degree phi region
        multRegionRnd       = 0;
        multEta09RegionRnd  = 0;
        multEta1RegionRnd   = 0;
        multV0APosRegionRnd = 0;
        multV0ANegRegionRnd = 0;
        multV0CPosRegionRnd = 0;
        multV0CNegRegionRnd = 0;
        
        // collision type: single, double, central or non-diffractive,
        //   see http://home.thep.lu.se/~torbjorn/pythia82html/QCDProcesses.html
        type = pythia.info.code();
        // number of hard interactions (0 is for elastic and diffractive events)
        nMPI = pythia.info.nMPI();
        
        weight    = pythia.info.hiInfo->weight();
        nPartTarg = pythia.info.hiInfo->nPartTarg();
        nAbsTarg  = pythia.info.hiInfo->nAbsTarg();
        nDiffTarg = pythia.info.hiInfo->nDiffTarg();

        vector <Quarkonium> foundQuarkoniumPerEvent;
        // here the loop over the particles starts
        for(int iPart=0; iPart<pythia.event.size(); iPart++) {

            Particle* part = &pythia.event[iPart];
            if(isPrimaryChargedALICE(iPart, pythia)) {

                TClonesArray& chTracks = *(testEvent->tracks);
                Track* chTrack = NULL;
                chTrack = new(chTracks[testEvent->multFull]) Track();
                chTrack->eta = part->eta();
                testEvent->multFull += 1;

                ++multFull;
                if(abs(part->eta())<0.9)                 ++multEta09;
                if(abs(part->eta())<1.0)                 ++multEta1;
                else if(inV0APosAcceptance(part->eta())) ++multV0APos;
                else if(inV0ANegAcceptance(part->eta())) ++multV0ANeg;
                if(inV0CPosAcceptance(part->eta()))      ++multV0CPos;
                else if(inV0CNegAcceptance(part->eta())) ++multV0CNeg;
                fillRegionRandom(multRegionRnd, multEta09RegionRnd, multEta1RegionRnd,
                                 multV0APosRegionRnd, multV0ANegRegionRnd,
                                 multV0CPosRegionRnd, multV0CNegRegionRnd,
                                 abs(part->phi())/M_PI, part->eta());
            }

            // we found a quarkonium
            else if(isOnium(part->id())) {
                int pdg = part->id();
                int pdgDau1 = pythia.event[part->daughter1()].id();
                int pdgDau2 = pythia.event[part->daughter2()].id();
                // avoid double counting ("J/psi -> J/psi")
                if(!(pdgDau1==pdg) && !(pdgDau2==pdg)) {
                    Quarkonium found;
                    found.pt    = part->pT();
                    found.y     = part->y();
                    found.pdg   = pdg;
                    // look for electrons(11) and muons(13)
                    if     (abs(pdgDau1)==11 && abs(pdgDau2)==11) found.decayChannel = 1;
                    else if(abs(pdgDau1)==13 && abs(pdgDau2)==13) found.decayChannel = -1;
                    // get information on where the J/psi came from, i.e. prompt or non-prompt
                    traceBackQuarkonium(found, iPart, pythia);

                    // fill in the multiplicity values in the 3 regions
                    fillRegionsWrtQuarkonium(found, pythia, part->phi()/M_PI);

                    // add the found quarkonium to the list of found quarkonia in this event
                    foundQuarkoniumPerEvent.push_back(found);
                }
            }

        }  // end of particle loop

        eventTree->Fill();
        // we found some quarkonia in this event -> fill them in the tree
        if(foundQuarkoniumPerEvent.size()) {
            for(std::vector<Quarkonium>::iterator iOnium=foundQuarkoniumPerEvent.begin(); iOnium!=foundQuarkoniumPerEvent.end(); iOnium++) {
                onium = *iOnium;
                oniumTree->Fill();
            }
            foundQuarkoniumPerEvent.clear();
        }
        
        // pythia.event.list();
    }  // end event loop

    pythia.stat();
    
    eventTree->Write();
    oniumTree->Write();
    fout->Close();

}


bool inV0APosAcceptance(float eta) {
    return eta > 2.8 && eta < 5.1;
}


bool inV0ANegAcceptance(float eta) {
    return eta < -2.8 && eta > -5.1;
}


bool inV0CPosAcceptance(float eta) {
    return eta < 3.7 && eta > 1.7;
}


bool inV0CNegAcceptance(float eta) {
    return eta > -3.7 && eta < -1.7;
}


/**
 *
 * Function that finds out, where the quarkonia comes from. In this simplified version it only checks
 * if it is a prompt or non-prompt J/psi (i.e. if one of it's ancestor's was a b-meson or baryon)
 * The function gets the mother of the quarkonium, checks it's id and status code, and decides if it has to
 * go further back to find the initially produced particle.
 * Check http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html for the meaning of the status
 * codes.
 *
 */
void traceBackQuarkonium(Quarkonium &ret, unsigned short nextParticle, Pythia & pythia) {

    while(nextParticle) {
        Particle* part      = &pythia.event[nextParticle];
        unsigned int status = part->statusAbs();
        unsigned int m1     = part->mother1();
        unsigned int m2     = part->mother2();

        // particle is carbon-copy of mother -> keep going back
        if(m1 == m2) {
            nextParticle = m1;
        }
        // particle was produced in hard scattering -> we found the initially produced particle!
        else if(status==33 || status==23) {
            ret.initialPdg = part->id();
            nextParticle   = 0;
        }
        // particle was produced in string fragmentation step -> we found the initially produced particle!
        else if(status>80 && status<90) {
            ret.initialPdg = part->id();
            nextParticle   = 0;
        }
        // we are at some intermediate state (particle is rescattering or decaying) -> keep going back
        else {
            nextParticle = m1;
        }
    }

    return;
}


// Find out the quark content of a particle from it's PDG identifier
// see http://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf


int firstQuark(int pdg) {
    return (abs(pdg)%100 - abs(pdg)%10) / 10;
}

int secondQuark(int pdg) {
    return (abs(pdg)%1000 - abs(pdg)%100) / 100;
}

int thirdQuark(int pdg) {
    return (abs(pdg)%10000 - abs(pdg)%1000) / 1000;
}

bool isCharmonium(int pdg) {
    return (firstQuark(pdg)==4 && secondQuark(pdg)==4 && thirdQuark(pdg)==0);
}

bool isBottomonium(int pdg) {
    return (firstQuark(pdg)==5 && secondQuark(pdg)==5 && thirdQuark(pdg)==0);
}

bool isOnium(int pdg) {
    return (isCharmonium(pdg) || isBottomonium(pdg));
}


/**
 *
 * Find out if given particle is a primary charged particle according to ALICE definition,
 * see https://cds.cern.ch/record/2270008
 *
 */
bool isPrimaryChargedALICE(unsigned short idx, Pythia &pythia) {

    if(!pythia.event[idx].isCharged()) return false;

    int status = pythia.event[idx].statusHepMC();
    if(status==0 || status ==4 || (status>=11 && status<=200)) return false;
    unsigned int pdg = pythia.event[idx].idAbs();
    if(!isLongLived(pdg)) return false;
    while(idx = pythia.event[idx].mother1()) {
        status = pythia.event[idx].statusHepMC();
        if(status == 4 || status==13) return true;
        pdg = pythia.event[idx].idAbs();
        if(isLongLived(pdg)) return false;
    }
    return true;
}


bool isLongLived(unsigned int pdg) {
    if(pdg > 1000000000) return true;

    switch (pdg) {
    case 13:
    case 11:
    case 22:
    case 211:
    case 321:
    case 130:
    case 310:
    case 2212:
    case 2112:
    case 3122:
    case 3112:
    case 3222:
    case 3312:
    case 3322:
    case 3334:
    case 12:
    case 14:
    case 16:
        return true;
    }
    return false;
}


void fillRegions(unsigned short &multRegion1,
                 unsigned short &multEta09Region1,  unsigned short &multEta1Region1,
                 unsigned short &multV0APosRegion1, unsigned short &multV0ANegRegion1,
                 unsigned short &multV0CPosRegion1, unsigned short &multV0CNegRegion1,
                 unsigned short &multRegion2,
                 unsigned short &multEta09Region2,  unsigned short &multEta1Region2,
                 unsigned short &multV0APosRegion2, unsigned short &multV0ANegRegion2,
                 unsigned short &multV0CPosRegion2, unsigned short &multV0CNegRegion2,
                 unsigned short &multRegion3,
                 unsigned short &multEta09Region3,  unsigned short &multEta1Region3,
                 unsigned short &multV0APosRegion3, unsigned short &multV0ANegRegion3,
                 unsigned short &multV0CPosRegion3, unsigned short &multV0CNegRegion3,
                 double phi, double eta) {
    if(phi<1./3.) {
        ++multRegion1;
        if(abs(eta)<0.9)                 ++multEta09Region1;
        if(abs(eta)<1.0)                 ++multEta1Region1;
        else if(inV0APosAcceptance(eta)) ++multV0APosRegion1;
        else if(inV0ANegAcceptance(eta)) ++multV0ANegRegion1;
        if(inV0CPosAcceptance(eta))      ++multV0CPosRegion1;
        else if(inV0CNegAcceptance(eta)) ++multV0CNegRegion1;
    }
    else if(phi<2./3) {
        ++multRegion2;
        if(abs(eta)<0.9)                 ++multEta09Region2;
        if(abs(eta)<1.0)                 ++multEta1Region2;
        else if(inV0APosAcceptance(eta)) ++multV0APosRegion2;
        else if(inV0ANegAcceptance(eta)) ++multV0ANegRegion2;
        if(inV0CPosAcceptance(eta))      ++multV0CPosRegion2;
        else if(inV0CNegAcceptance(eta)) ++multV0CNegRegion2;
    }
    else {
        ++multRegion3;
        if(abs(eta)<0.9)                 ++multEta09Region3;
        if(abs(eta)<1.0)                 ++multEta1Region3;
        else if(inV0APosAcceptance(eta)) ++multV0APosRegion3;
        else if(inV0ANegAcceptance(eta)) ++multV0ANegRegion3;
        if(inV0CPosAcceptance(eta))      ++multV0CPosRegion3;
        else if(inV0CNegAcceptance(eta)) ++multV0CNegRegion3;
    }
}

void fillRegionRandom(unsigned short &multRegion,
                      unsigned short &multEta09Region,  unsigned short &multEta1Region,
                      unsigned short &multV0APosRegion, unsigned short &multV0ANegRegion,
                      unsigned short &multV0CPosRegion, unsigned short &multV0CNegRegion,
                      double phi, double eta) {
    if(phi<1./3.) {
        ++multRegion;
        if(abs(eta)<0.9)                 ++multEta09Region;
        if(abs(eta)<1.0)                 ++multEta1Region;
        else if(inV0APosAcceptance(eta)) ++multV0APosRegion;
        else if(inV0ANegAcceptance(eta)) ++multV0ANegRegion;
        if(inV0CPosAcceptance(eta))      ++multV0CPosRegion;
        else if(inV0CNegAcceptance(eta)) ++multV0CNegRegion;
    }
}

void fillRegionsWrtQuarkonium(Quarkonium &found, Pythia& pythia, double phi_quarkonium) {
    for(int iPart2=0; iPart2<pythia.event.size(); iPart2++) {
        Particle* part2 = &pythia.event[iPart2];
        if(isPrimaryChargedALICE(iPart2, pythia)) {
            float phi_assoc = part2->phi() / M_PI;
            float deltaPhi  = abs(phi_assoc - phi_quarkonium);
            if(deltaPhi>1) deltaPhi = 2 - deltaPhi;
            fillRegions(found.multRegion1,
                        found.multEta09Region1,  found.multEta1Region1,
                        found.multV0APosRegion1, found.multV0ANegRegion1,
                        found.multV0CPosRegion1, found.multV0CNegRegion1,
                        found.multRegion2,
                        found.multEta09Region2,  found.multEta1Region2,
                        found.multV0APosRegion2, found.multV0ANegRegion2,
                        found.multV0CPosRegion2, found.multV0CNegRegion2,
                        found.multRegion3,
                        found.multEta09Region3,  found.multEta1Region3,
                        found.multV0APosRegion3, found.multV0ANegRegion3,
                        found.multV0CPosRegion3, found.multV0CNegRegion3,
                        deltaPhi, part2->eta());
        }
    }
}
