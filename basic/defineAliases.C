/*
 * This macro defines the aliases in the event and quarkonium trees that are used
 * in the plotRelOrigin macro.
 * 
 */
void defineAliases(TTree *eventTree, TTree *oniumTree) {
    
    
    // V0M = V0A + V0C. 
    eventTree->SetAlias("multV0MpPb", "(multV0APos+multV0CNeg)");
    oniumTree->SetAlias("multV0MpPb", "(multV0APos+multV0CNeg)");
    eventTree->SetAlias("multV0MPbp", "(multV0ANeg+multV0CPos)");
    oniumTree->SetAlias("multV0MPbp", "(multV0ANeg+multV0CPos)");
    
    eventTree->SetAlias("multV0MpPbRegion1", "(multV0APosRegion1+multV0CNegRegion1)");
    oniumTree->SetAlias("multV0MpPbRegion1", "(multV0APosRegion1+multV0CNegRegion1)");
    eventTree->SetAlias("multV0MPbpRegion1", "(multV0ANegRegion1+multV0CPosRegion1)");
    oniumTree->SetAlias("multV0MPbpRegion1", "(multV0ANegRegion1+multV0CPosRegion1)");
    
    eventTree->SetAlias("multV0MpPbRegion2", "(multV0APosRegion2+multV0CNegRegion2)");
    oniumTree->SetAlias("multV0MpPbRegion2", "(multV0APosRegion2+multV0CNegRegion2)");
    eventTree->SetAlias("multV0MPbpRegion2", "(multV0ANegRegion2+multV0CPosRegion2)");
    oniumTree->SetAlias("multV0MPbpRegion2", "(multV0ANegRegion2+multV0CPosRegion2)");
    
    eventTree->SetAlias("multV0MpPbRegion3", "(multV0APosRegion3+multV0CNegRegion3)");
    oniumTree->SetAlias("multV0MpPbRegion3", "(multV0APosRegion3+multV0CNegRegion3)");
    eventTree->SetAlias("multV0MPbpRegion3", "(multV0ANegRegion3+multV0CPosRegion3)");
    oniumTree->SetAlias("multV0MPbpRegion3", "(multV0ANegRegion3+multV0CPosRegion3)");
    
    
    // for the event tree, the 3 regions are chosen randomly. 
    // in the trees only one value of the multiplicity in a randomy chosen region is stored
    // and used for all 3 regions
    
    eventTree->SetAlias("multRegion1", "multRegionRnd");
    eventTree->SetAlias("multRegion2", "multRegionRnd");
    eventTree->SetAlias("multRegion3", "multRegionRnd");
    
    eventTree->SetAlias("multEta1Region1", "multEta1RegionRnd");
    eventTree->SetAlias("multEta1Region2", "multEta1RegionRnd");
    eventTree->SetAlias("multEta1Region3", "multEta1RegionRnd");
    
    eventTree->SetAlias("multV0APosRegion1", "multV0APosRegionRnd");
    eventTree->SetAlias("multV0APosRegion2", "multV0APosRegionRnd");
    eventTree->SetAlias("multV0APosRegion3", "multV0APosRegionRnd");
    eventTree->SetAlias("multV0ANegRegion1", "multV0ANegRegionRnd");
    eventTree->SetAlias("multV0ANegRegion2", "multV0ANegRegionRnd");
    eventTree->SetAlias("multV0ANegRegion3", "multV0ANegRegionRnd");
    
    eventTree->SetAlias("multV0CPosRegion1", "multV0CPosRegionRnd");
    eventTree->SetAlias("multV0CPosRegion2", "multV0CPosRegionRnd");
    eventTree->SetAlias("multV0CPosRegion3", "multV0CPosRegionRnd");
    eventTree->SetAlias("multV0CNegRegion1", "multV0CNegRegionRnd");
    eventTree->SetAlias("multV0CNegRegion2", "multV0CNegRegionRnd");
    eventTree->SetAlias("multV0CNegRegion3", "multV0CNegRegionRnd");
    
//     // V0AND cut
//     eventTree->SetAlias("V0ANDpPb", "(multV0APos>0 && multV0CNeg>0)");
//     oniumTree->SetAlias("V0ANDpPb", "(multV0APos>0 && multV0CNeg>0)");
//     eventTree->SetAlias("V0ANDPbp", "(multV0ANeg>0 && multV0CPos>0)");
//     oniumTree->SetAlias("V0ANDPbp", "(multV0ANeg>0 && multV0CPos>0)");

    // INEL>0 cut
    eventTree->SetAlias("inel0", "(multEta1>0)");
    oniumTree->SetAlias("inel0", "(multEta1>0)");

    // Inelastic events
    eventTree->SetAlias("inelastic", "(type!=102)");
    oniumTree->SetAlias("inelastic", "(type!=102)");

    // NSD events (double diffractive + non-diffractive)
    eventTree->SetAlias("NSD", "(type==101 || type==105)");
    oniumTree->SetAlias("NSD", "(type==101 || type==105)");
    
    // Diffractive events
    eventTree->SetAlias("diffractive", "(type==103 || type==104 || type==105 || type==106)");  // for pPb?
    oniumTree->SetAlias("diffractive", "(type==103 || type==104 || type==105 || type==106)");  // for pPb?
//     eventTree->SetAlias("diffractive", "(type!=101)");  // for pp with SoftQCD:inelatsic = on
//     oniumTree->SetAlias("diffractive", "(type!=101)");  // for pp with SoftQCD:inelatsic = on
    
    
    // number of MPI (nMPI is one for events with one hard scattering, 2 for events with 1 additional
    // scattering)
    // so subtract 1 to get number of MPI
    eventTree->SetAlias("mpi", "(nMPI-1)");
    oniumTree->SetAlias("mpi", "(nMPI-1)");
    
    
    oniumTree->SetAlias("midRapidity",     "TMath::Abs(onium.y)<0.9");
    oniumTree->SetAlias("muonRapidity",    "TMath::Abs(onium.y)>2.5 && TMath::Abs(onium.y)<4.");
    oniumTree->SetAlias("muonRapidityPos", "onium.y>2.5 && onium.y<4.");
    oniumTree->SetAlias("muonRapidityNeg", "onium.y<-2.5 && onium.y>-4.");
    
    
    oniumTree->SetAlias("dilepton",   "abs(onium.decayChannel)==1");
    oniumTree->SetAlias("dielectron", "onium.decayChannel==1");
    
    
    // absolute value of the PDG code of the initially produced particle
    oniumTree->SetAlias("pdg", "(onium.initialPdg<0?-1*onium.initialPdg:onium.initialPdg)");
    
    // absolute value of the PDG code of the final particle
    oniumTree->SetAlias("PDG", "(onium.pdg<0?-1*onium.pdg:onium.pdg)");
    
    oniumTree->SetAlias("jpsi", "PDG==443");
    
    oniumTree->SetAlias("charmonia",
        "((PDG%1000 - PDG%100)/100 == 4 && (PDG%100 - PDG%10)/10 == 4 && (PDG%10000 - PDG%1000)/1000 == 0)");
    
    oniumTree->SetAlias("bottomonia",
        "((PDG%1000 - PDG%100)/100 == 5 && (PDG%100 - PDG%10)/10 == 5 && (PDG%10000 - PDG%1000)/1000 == 0)");
    
    
    // directly produced J/psi, i.e. no non-prompt and no feeddown of higher mass charmonia
    // but J/psi can be initially produced in a color octet state, these are the codes starting with 99...
    // see http://home.thep.lu.se/~torbjorn/pythia82html/OniaProcesses.html
    oniumTree->SetAlias("direct", "(pdg==443 || pdg == 9940003 || pdg == 9941003 || pdg == 9942003)");
    oniumTree->SetAlias("nonprompt",
        "(pdg < 9940000) && ((pdg%10000 - pdg%1000)/1000 == 5 || (pdg%1000 - pdg%100)/100 == 5)");
    oniumTree->SetAlias("cfeeddown", "(!direct && !nonprompt)");
    oniumTree->SetAlias("octet",     "(pdg > 9900000)");
    
}
