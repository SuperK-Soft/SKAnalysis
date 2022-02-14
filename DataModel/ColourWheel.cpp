#ifndef COLOUR_WHEEL_CLASS
#define COLOUR_WHEEL_CLASS

#include "ColourWheel.h"

const std::vector<EColor> ColourWheel::colours{kBlack, kBlue, (EColor)(kGreen+1), kRed, kViolet, kOrange, kMagenta,(EColor)(kAzure+2),(EColor)(kOrange+4),(EColor)(kViolet-6),(EColor)(kTeal-6)};
const std::vector<std::string> ColourWheel::colournames{"kBlack", "kBlue", "kGreen", "kRed", "kViolet", "kOrange", "kMagenta","kAzure","kOrange","kViolet","kTeal"};
// these are the raw colours, for ref.
//enum EColor { kWhite =0, kBlack =1, kGray=920, kRed =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432, kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };

#endif
