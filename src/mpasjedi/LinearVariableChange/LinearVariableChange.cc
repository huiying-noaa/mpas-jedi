/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"

#include "oops/util/Logger.h"

#include "mpasjedi/Geometry/Geometry.h"
#include "mpasjedi/Increment/Increment.h"
#include "mpasjedi/LinearVariableChange/LinearVariableChange.h"
#include "mpasjedi/State/State.h"

namespace mpas {

// -------------------------------------------------------------------------------------------------

LinearVariableChange::LinearVariableChange(const Geometry & geom,
                                           const eckit::Configuration & config)
  : geom_(geom), linearVariableChange_(),vader_() {
  params_.deserialize(config);
  
  eckit::LocalConfiguration variableChangeConfig = params_.toConfiguration();
  ModelData modelData{geom};
  eckit::LocalConfiguration vaderConfig;
  vaderConfig.set(vader::configCookbookKey,
                  variableChangeConfig.getSubConfiguration("vader custom cookbook"));
  vaderConfig.set(vader::configModelVarsKey, modelData.modelData());

  // Create vader with mpas-jedi custom cookbook
  vader_.reset(new vader::Vader(params_.linearVariableChangeParameters.value().vader,
                                vaderConfig));
}

// -------------------------------------------------------------------------------------------------

LinearVariableChange::~LinearVariableChange() {}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarTraj(const State & xfg, const oops::Variables & vars) {
  oops::Log::trace() << "LinearVariableChange::changeVarTraj starting" << std::endl;
  // Call Vader's changeVarTraj to populate its initial trajectory FieldSet

  oops::Variables varsVader = vars; //cltthink it should be long names
  
  atlas::FieldSet xfgfs;
  xfg.toFieldSet(xfgfs);
  vader_->changeVarTraj(xfgfs, varsVader);
  // If input and output variables are specified in the yaml, we use those variables to finish
  // initializing vader's linear variable change now. Otherwise we have to wait until changeVarTL or
  // changeVarAD is called to find out the ingredient/increment vars.

  const auto &lvc_params = params_.linearVariableChangeParameters.value();
  if (lvc_params.inputVariables.value() != boost::none &&
      lvc_params.outputVariables.value() != boost::none) {
    oops::Variables inputVars = *lvc_params.inputVariables.value();
    oops::Variables outputVars = *lvc_params.outputVariables.value();
    ASSERT_MSG(outputVars == vars, "outputVariables in config file must match output "
          "variables passed to changeVarTraj");
    oops::Variables ingredientVars =inputVars;
    initVaderTLAD(ingredientVars);
  }

  // Create the variable change
  linearVariableChange_.reset(LinearVariableChangeFactory::create(xfg, xfg, geom_,
             params_.linearVariableChangeParameters.value()));
  oops::Log::trace() << "LinearVariableChange::changeVarTraj done" << std::endl;
}

// -------------------------------------------------------------------------------------------------
void LinearVariableChange::initVaderTLAD(oops::Variables & ingredientVars) const {
  oops::Log::trace() << "LinearVariableChange::initVaderTLAD starting" << std::endl;
  oops::Variables originalIngredientVars = ingredientVars;
  varsVaderPopulates_ = vader_->initTLAD(ingredientVars);
  varsVaderPopulates_ -= originalIngredientVars;
  oops::Log::trace() << "LinearVariableChange::initVaderTLAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarTL(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarTL starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::changeVarTL done (identity)" << std::endl;
    return;
  }
  // Make sure this object is fully initialized
  if (vader_->needsTLADInit()) {
    oops::Variables ingredientVars = dx.variables();
    initVaderTLAD(ingredientVars);
  }
  // If Vader is doing anything, call Vader
  if (varsVaderPopulates_.size() > 0) {
    atlas::FieldSet dxfs;
    dx.toFieldSet(dxfs);
    vader_->changeVarTL(dxfs);

    // Set intermediate state for the Increment containing original fields plus the ones
    // Vader has done
// now updateFields is not available in mpasjedi::Increment
// assume vadar only works on a part (like tv ) pf the variables
    oops::Variables varsVader = dx.variables();
    varsVader += varsVaderPopulates_;
    dx.updateFields(varsVader);
    dx.fromFieldSet(dxfs);
  }

  // Create output increment
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->changeVarTL(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarTL done" << dx << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseTL(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarInverseTL starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::changeVarInverseTL done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->changeVarInverseTL(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarInverseTL done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarAD(Increment & dx, const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarAD starting" << std::endl;

 // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::changeVarAD done (identity)" << std::endl;
    return;
  }

  // Make sure this object is fully initialized
  if (vader_->needsTLADInit()) {
    oops::Variables ingredientVars(vars);
    initVaderTLAD(ingredientVars);
  }
  // Create dxin as a copy of dx, minus the variables created by Vader (in the forward direction)
  // This way we ensure the model code will not be able to do the adjoint for these vars
  Increment dxin(dx, true);  // true => full copy
  oops::Variables varsVaderDidntPopulate = dx.variables();
  varsVaderDidntPopulate -= varsVaderPopulates_;
  dxin.updateFields(varsVaderDidntPopulate);

  dx.updateFields(varsVaderPopulates_);



  // Create output increment
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->changeVarAD(dxin, dxout);

  // dxout needs to temporarily have the variables that Vader populated put into it before
  // being passed into vader_.changeVarAD, so Vader can do its adjoints.
  atlas::FieldSet dxout_fs;
  dxout.toFieldSet(dxout_fs);
  oops::Variables varsVaderWillAdjoint = varsVaderPopulates_;
  if (varsVaderWillAdjoint.size() > 0) {
    atlas::FieldSet dx_fs;
    dx.toFieldSet(dx_fs);
    for (const auto field : dx_fs) {
      dxout_fs.add(field);
    }

    oops::Variables varsAdjointed = vader_->changeVarAD(dxout_fs);
    varsVaderWillAdjoint -= varsAdjointed;
    // After changeVarAD, vader should have removed everything from varsVaderWillAdjoint,
    // indicating it did all the adjoints we expected it to.
    ASSERT(varsVaderWillAdjoint.size() == 0);
  }


  // Copy data from temporary increment
  dx.updateFields(vars);
  dx.fromFieldSet(dxout_fs);

  oops::Log::trace() << "LinearVariableChange::changeVarAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::changeVarInverseAD(Increment & dx,
                                             const oops::Variables & vars) const {
  oops::Log::trace() << "LinearVariableChange::changeVarInverseAD starting" << std::endl;

  // If all variables already in incoming increment just remove the no longer needed fields
  if (vars == dx.variables()) {
    oops::Log::trace() << "LinearVariableChange::changeVarInverseAD done (identity)" << std::endl;
    return;
  }

  // Create output increment
  Increment dxout(dx.geometry(), vars, dx.time());

  // Call variable change
  linearVariableChange_->changeVarInverseAD(dx, dxout);

  // Copy data from temporary increment
  dx = dxout;

  oops::Log::trace() << "LinearVariableChange::changeVarInverseAD done" << std::endl;
}

// -------------------------------------------------------------------------------------------------

void LinearVariableChange::print(std::ostream & os) const {
  os << classname() << " linear variable change";
}

// -------------------------------------------------------------------------------------------------

}  // namespace mpas
