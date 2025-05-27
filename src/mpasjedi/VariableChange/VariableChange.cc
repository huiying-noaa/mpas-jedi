/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/State/State.h"
#include "mpasjedi/Traits.h"
#include "mpasjedi/VariableChange/VariableChange.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------
  typedef VariableChangeParameters Parameters_;

VariableChange::VariableChange(const eckit::Configuration & config, const Geometry & geometry)
  : vader_(), run_vader_(), run_mpasjedi_()
{
  // Create the variable change
  VariableChangeParametersWrapper params;
  params.deserialize(config);
    run_vader_ = params.run_vader.value();
  run_mpasjedi_ = params.run_mpasjedi.value();
  eckit::LocalConfiguration variableChangeConfig = params.toConfiguration();
  ModelData modelData{geometry};
  eckit::LocalConfiguration vaderConfig;
  vaderConfig.set(vader::configCookbookKey,
                  variableChangeConfig.getSubConfiguration("vader custom cookbook"));
  vaderConfig.set(vader::configModelVarsKey, modelData.modelData());

  // Create vader with mpas-jedi custom cookbook
  vader_.reset(new vader::Vader(params.variableChangeParameters.value().vader,
                                vaderConfig));
  variableChange_.reset(VariableChangeFactory::create(geometry,
                        params.variableChangeParameters.value()));
}

// -------------------------------------------------------------------------------------------------

VariableChange::~VariableChange() {}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVar(State & x, const oops::Variables & vars) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVar starting" << std::endl;

//clthtink vars should be long names ,s in mpasjedi , howow
  // If all variables already in incoming state just remove the no longer needed fields
  if (vars == x.variables()) {
    oops::Log::info() << "VariableChange::changeVar done (identity)" << std::endl;
    return;
  }
  oops::Log::trace() << "VariableChange::changeVar, vars" << vars << std::endl;
  // Call Vader to perform first set of variable transforms
  // ------------------------------------------------------

  // Record start variables
  oops::Variables varsFilled = x.variables();

  oops::Variables varsVader = vars;
  varsVader -= varsFilled;  // Pass only the needed variables


  // Call Vader. On entry, varsVader holds the vars requested from Vader; on exit,
  // it holds the vars NOT fulfilled by Vader, i.e., the vars still to be requested elsewhere.
  // vader_->changeVar also returns the variables fulfilled by Vader. These variables are
  // allocated and populated and added to the FieldSet (xfs).
  if (run_vader_) {
    atlas::FieldSet xfs;
    x.toFieldSet(xfs);
    const oops::Variables varsVaderPopulated = vader_->changeVar(xfs, varsVader);
    if (varsVaderPopulated.size() > 0) {
      varsFilled += varsVaderPopulated;
//cltthink      x.updateFields(varsFilled);
      x.fromFieldSet(xfs);
    }
  }


  // Create output state
  State xout(x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVar(x, xout);

  // Copy data from temporary state
  x = xout;

  // Trace
  oops::Log::trace() << "VariableChange::changeVar done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::changeVarInverse(State & x, const oops::Variables & vars) const {
  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse starting" << std::endl;

  // If all variables already in incoming state just remove the no longer needed fields
  if (vars == x.variables()) {
    oops::Log::info() << "VariableChange::changeVarInverse done (identity)" << std::endl;
    return;
  }

  // Call Vader to perform first set of variable transforms
  // ------------------------------------------------------

  // Record start variables
  oops::Variables varsFilled = x.variables();

  oops::Variables varsVader = vars;
  varsVader -= varsFilled;  // Pass only the needed variables

  // Call Vader. On entry, varsVader holds the vars requested from Vader; on exit,
  // it holds the vars NOT fulfilled by Vader, i.e., the vars still to be requested elsewhere.
  // vader_->changeVar also returns the variables fulfilled by Vader. These variables are
  // allocated and populated and added to the FieldSet (xfs).
  atlas::FieldSet xfs;
  x.toFieldSet(xfs);
  const oops::Variables varsVaderPopulated = vader_->changeVar(xfs, varsVader);
  if (varsVaderPopulated.size() > 0) {
    varsFilled += varsVaderPopulated;
//cltthink    x.updateFields(varsFilled);
    x.fromFieldSet(xfs);
  }

  // Create output state
  State xout(x.geometry(), vars, x.time());

  // Call variable change
  variableChange_->changeVarInverse(x, xout);

  // Copy data from temporary state
  x = xout;

  // Trace
  oops::Log::trace() << "VariableChange::changeVarInverse done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void VariableChange::print(std::ostream & os) const {
  os << *variableChange_;
}

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
