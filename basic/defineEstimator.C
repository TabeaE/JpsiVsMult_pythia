/**
 * 
 * This macro defines the multiplicity estimator according to the code given to it.
 * The codes are the following:
 * 
 * 
 * 0: nch in full rapidity
 * 1: nch in |eta| < 0.9
 * 2: nch in |eta| < 1.0
 * 3: nch in V0A (pPb)
 * 4: nch in V0A (Pbp)
 * 5: nch in V0C (pPb)
 * 6: nch in V0C (Pbp)
 * 7: nch in V0M = V0A+V0C (pPb)
 * 8: nch in V0M = V0A+V0C (Pbp)
 * 9: number of MPI
 * 
 * 10: nch in full rapidity in towards region
 * 11: nch in |eta| < 0.9 in towards region
 * 12: nch in |eta| < 1.0 in towards region
 * 13: nch in V0A in towards region (pPb)
 * 14: nch in V0A in towards region (Pbp)
 * 15: nch in V0C in towards region (pPb)
 * 16: nch in V0C in towards region (Pbp)
 * 17: nch in V0M in towards region (pPb)
 * 18: nch in V0M in towards region (Pbp)
 * 
 * 20: nch in full rapidity in transverse region
 * 21: nch in |eta| < 0.9 in transverse region
 * 22: nch in |eta| < 1.0 in transverse region
 * 23: nch in V0A in transverse region (pPb)
 * 24: nch in V0A in transverse region (Pbp)
 * 25: nch in V0C in transverse region (pPb)
 * 26: nch in V0C in transverse region (Pbp)
 * 27: nch in V0M in transverse region (pPb)
 * 28: nch in V0M in transverse region (Pbp)
 * 
 * 30: nch in full rapidity in away region
 * 31: nch in |eta| < 0.9 in away region
 * 32: nch in |eta| < 1.0 in away region
 * 33: nch in V0A in away region (pPb)
 * 34: nch in V0A in away region (Pbp)
 * 35: nch in V0C in away region (pPb)
 * 36: nch in V0C in away region (Pbp)
 * 37: nch in V0M in away region (pPb)
 * 38: nch in V0M in away region (Pbp)
 * 
 * 
 * 
 * The macro returns 3 strings:
 *  - estimatorString:    to be given to the tree for drawing
 *  - stringForAxisTitle: used for the x axis title
 *  - regionString:       also used in the x axis title, describes the phi region
 * 
 * In addition it returns a factor by which to scale the x axis (because multiplicity in V0M is higher than in
 * |eta| < 1 and so on).
 * 
 * 
 * 
 * 
 */

void defineEstimator(int estimator, TString &estimatorString, TString &stringForAxisTitle, TString &regionString, 
                     double &factor) {
  
  factor = 1.;
  switch(estimator) {
    case 9:
      stringForAxisTitle = "#it{N}_{MPI}";
      estimatorString    = "nMPI";
      return;
      break;
    default:
      stringForAxisTitle = "#it{N}_{ch}";
      estimatorString    = "mult";
      break;
  }

  int last_digit = estimator %10;
  switch(last_digit ) {
    case 0:
      regionString = "full rapidity";
      factor *= 2.;
      break;
    case 1:
      regionString     = "|#it{#eta}|<0.9";
      estimatorString += "Eta09";
      break;
    case 2:
      regionString     = "|#it{#eta}|<1";
      estimatorString += "Eta1";
      break;
    case 3:
      regionString     = "V0ApPb";
      estimatorString += "V0ApPb";
      break;
    case 4:
      regionString     = "V0APbp";
      estimatorString += "V0APbp";
      break;
    case 5:
      regionString     = "V0CpPb";
      estimatorString += "V0CpPb";
      break;
    case 6:
      regionString     = "V0CPbp";
      estimatorString += "V0CPbp";
      break;
    case 7:
      regionString     = "V0MpPb";
      estimatorString += "V0MpPb";
      factor *= 1.5;
      break;
    case 8:
      regionString     = "V0MPbp";
      estimatorString += "V0MPbp";
      factor *= 1.5;
      break;
  }
  
  int second_digit = (estimator%100 - estimator%10) / 10;
  switch(second_digit) {
    case 1:
      regionString    += ", toward region";
      estimatorString += "Region1";
      factor /= 2.;
      break;
    case 2:
      regionString    += ", transverse region";
      estimatorString += "Region2";
      factor /= 2.;
      break;
    case 3:
      regionString    += ", away region";
      estimatorString += "Region3";
      factor /= 2.;
      break;
  }
  
}
