#include <TMB.hpp>
#include "sparse_joint_model.h"
#include "sparse_joint_cohort_model.h"
#include "sparse_joint_series_model.h"
#include "sparse_joint_pairwise_model.h"
#include "sparse_joint_state_space_model.h"

// One DLL, several models: the data list names the model in `model`.
template<class Type>
Type objective_function<Type>::operator() () {
  SEXP model_name = getListElement(this->data, "model");
  if (Rf_isNull(model_name)) Rf_error("The TMB data list must name its model in `model`.");
  std::string model = CHAR(STRING_ELT(model_name, 0));
  if (model == "sparse_joint") return sparse_joint_model(this);
  if (model == "sparse_joint_cohort") return sparse_joint_cohort_model(this);
  if (model == "sparse_joint_series") return sparse_joint_series_model(this);
  if (model == "sparse_joint_pairwise") return sparse_joint_pairwise_model(this);
  if (model == "sparse_joint_state_space") return sparse_joint_state_space_model(this);
  Rf_error("Unknown TMB model: %s", model.c_str());
  return Type(0);
}
