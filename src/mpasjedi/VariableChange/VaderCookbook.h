/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <map>
#include <string>
#include <vector>

namespace mpas {

/// This static method is used to define the default value for the "vader custom cookbook" Parameter
/// for the VariableChange and LinearVariableChange classes.
// -------------------------------------------------------------------------------------------------

  static std::map<std::string, std::vector<std::string>> vaderMPASCustomCookbook() {
    return {
      // pm25 from fine smoke and dust in RRFS
      {"mass_density_of_particulate_matter_2p5_in_air", {"ParticulateMatter2p5_B"}}
    };
  }

}  // namespace mpas
