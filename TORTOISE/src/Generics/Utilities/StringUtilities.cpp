// Copyright © 2019, Marco Battiato <marco.battiato@ntu.edu.sg; battiato.marco@gmail.com>, All rights reserved.
//
// Licensed under the GNU GENERAL PUBLIC LICENSE Version 3 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   https://www.gnu.org/licenses/
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.
//
//
// DISCLAIMER: This is a version under active development and testing.
// Not all features have been sufficiently tested, and several features
// are only partially implemented. Moreover both the core and the interface
// may change. You are discouraged from using this version to publish
// results without the supervision of the developer.
//
// If you want, you are welcome to act as a beta tester. In that case
// please contact Marco Battiato at:
// marco.battiato@ntu.edu.sg or battiato.marco@gmail.com
//
// Check if a newer, stable and tested version has, in the meanwhile,
// been made available at:
// https://github.com/MarcoBattiato/TORTOISE
//
//
//  StringUtilities.cpp
//  TORTOISE
//
//  Created by Marco Battiato on 15/10/20.
//

#include <Generics/Utilities/StringUtilities.hpp>

// Taken from https://stackoverflow.com/questions/23457305/compare-strings-with-wildcard
bool matchStringRec(const std::string &pattern, const std::string &candidate, int p, int c) {
  if (pattern[p] == '\0') {
    return candidate[c] == '\0';
  } else if (pattern[p] == '*') {
    for (; candidate[c] != '\0'; c++) {
      if (matchStringRec(pattern, candidate, p+1, c))
        return true;
    }
    return matchStringRec(pattern, candidate, p+1, c);
  } else if (pattern[p] != '?' && pattern[p] != candidate[c]) {
    return false;
  }  else {
    return matchStringRec(pattern, candidate, p+1, c+1);
  }
}
bool matchString(const std::string &pattern, const std::string &candidate) {
    return matchStringRec(pattern,candidate,0,0);
}


std::string toString(int value,int digitsCount)
{
    std::ostringstream os;
    os << std::setfill('0') << std::setw(digitsCount)<<value;
    return os.str();
}
