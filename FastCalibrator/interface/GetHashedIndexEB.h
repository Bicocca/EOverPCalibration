#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <vector>
#include "hChain.h"
#include "h2Chain.h"
#include <TGraphErrors.h>

#include <TLorentzVector.h>
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"



int GetHashedIndexEB(int iEta, int iPhi, int Zside);

int GetIetaFromHashedIndex(int Index);

int GetIphiFromHashedIndex(int Index);
