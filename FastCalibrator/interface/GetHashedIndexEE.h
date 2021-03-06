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

int GetHashedIndexEE(int, int, int);

int GetIxFromHashedIndex(int);

int GetIyFromHashedIndex(int);

int GetZsideFromHashedIndex(int);

// Essential values to get EE geometry
static const int IX_MIN = 1;
static const int IY_MIN = 1;
static const int IX_MAX = 100;
static const int IY_MAX = 100;
static const int kEEhalf = 7324;
