/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include <ostream>
#include <string>
#include <vector>

#include "atlas/util/Config.h"

#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/ModelData/ModelData.h"
// (TODO) #include "mpasjedi/Utilities/Constants.h"

// -------------------------------------------------------------------------------------------------

namespace mpas {

// -------------------------------------------------------------------------------------------------

ModelData::ModelData(const Geometry & geometry) {}

// -------------------------------------------------------------------------------------------------

ModelData::~ModelData() {}

// -------------------------------------------------------------------------------------------------

const eckit::LocalConfiguration ModelData::modelData() const {
  eckit::LocalConfiguration modelData;

//  (TODO) BJJ: We need this if we actively use vader.
//  // Add all constants to modelData config
//  std::vector<std::string> allConstantsNames = getAllConstantsNames();
//  for (std::string allConstantsName : allConstantsNames) {
//    modelData.set(allConstantsName, getConstant(allConstantsName));
//  }

  // BJJ: Hard-coded for testing
  modelData.set("epsilon", 0.60836236933);

  return modelData;
}

// -------------------------------------------------------------------------------------------------

void ModelData::print(std::ostream & os) const {
  os << "mpas::ModelData::modelData(): " << modelData();
}

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
