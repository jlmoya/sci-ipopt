//////////////////////////////////////////////////////////////////////////////////////////
// sciipopt: a scilab interface to the ipopt non linear constrained optimization solver //
//////////////////////////////////////////////////////////////////////////////////////////

//  Copyright (C) 2008-2010 Yann Collette
//  Copyright (C) 2020 - UTC - Stéphane MOTTELET
//
//  SCIIPOPT is free software; you can redistribute it and/or modify it
//  under the terms of the GNU General Public License as published by the
//  Free Software Foundation; either version 2, or (at your option) any
//  later version.
//
//  This part of code is distributed with the FURTHER condition that it 
//  can be compiled and linked with the Matlab libraries and it can be 
//  used within the Matlab environment.
//
//  SCIIPOPT is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with SciCoinOR; see the file COPYING.  If not, write to the Free
//  Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "IpIpoptApplication.hpp"
#include "scilabjournal.hpp"
#include <string>
#include <new>
#include <limits>

#include "struct.hxx"
#include "function.hxx"
#include "parameters.hxx"

extern "C"
{
#include <sciprint.h>
#include <Scierror.h>
}

int manage_ipopt_params(types::Struct *pStructOptions, Ipopt::SmartPtr<Ipopt::OptionsList> options, int iLog)
{
  // Get the parameters stored in the options struct
  int     iValue;
  char *  pStrValue;
  double  dblValue;
  bool bOk;

#ifdef DEBUG
  DBGPRINTF("sciipopt: processing options\n");
#endif


  ////////////
  // Output //
  ////////////
  // verbosity level

  bOk = getIntInPList(pStructOptions, L"print_level", &iValue, 5, iLog, CHECK_BOTH, 0, 12);
  if (bOk) options->SetIntegerValue("print_level", iValue);  
  // File name of options file (default: ipopt.opt)
  bOk = getStringInPList(pStructOptions, L"option_file_name", &pStrValue, "ipopt.opt", iLog, CHECK_NONE);
  if (bOk) options->SetStringValue("option_file_name", pStrValue);
  FREE(pStrValue);

  /////////////////
  // Termination //
  /////////////////
  // Desired convergence tolerance (relative)
  bOk = getDoubleInPList(pStructOptions, L"tol", &dblValue, 1e-8, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("tol", dblValue);
  // Maximum number of iterations
  bOk = getIntInPList(pStructOptions, L"max_iter", &iValue, 3000, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("max_iter", iValue);
  // Desired threshold for the dual infeasibility
  // Absolute tolerance on the dual infeasibility. Successful termination requires that the max-norm of the (unscaled) dual infeasibility is less than this threshold.
  bOk = getDoubleInPList(pStructOptions, L"dual_inf_tol", &dblValue, 1.0, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("dual_inf_tol", dblValue);
  // Desired threshold for the constraint violation
  // Absolute tolerance on the constraint violation. Successful termination requires that the max-norm of the (unscaled) constraint violation is less than this threshold.
  bOk = getDoubleInPList(pStructOptions, L"constr_viol_tol", &dblValue, 1e-4, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("constr_viol_tol", dblValue);
  // Acceptance threshold for the complementarity conditions
  // Absolute tolerance on the complementarity. Successful termination requires that the max-norm of the (unscaled) complementarity is less than this threshold.
  bOk = getDoubleInPList(pStructOptions, L"compl_inf_tol", &dblValue, 1e-4, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("compl_inf_tol", dblValue);
  // Acceptable convergence tolerance (relative)
  // Determines which (scaled) overall optimality error is considered to be "acceptable." 
  // There are two levels of termination criteria. If the usual "desired" tolerances (see tol, dual_inf_tol etc) are satisfied at an iteration, 
  // the algorithm immediately terminates with a success message. On the other hand, if the algorithm encounters "acceptable_iter" many iterations 
  // in a row that are considered "acceptable", it will terminate before the desired convergence tolerance is met. This is useful in cases where
  // the algorithm might not be able to achieve the "desired" level of accuracy.
  bOk = getDoubleInPList(pStructOptions, L"acceptable_tol", &dblValue, 1e-6, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("acceptable_tol", dblValue);
  // Acceptance threshold for the constraint violation
  // Absolute tolerance on the constraint violation. "Acceptable" termination requires that the max-norm of the (unscaled) constraint violation is
  // less than this threshold; see also acceptable_tol.
  bOk = getDoubleInPList(pStructOptions, L"acceptable_constr_viol_tol", &dblValue, 1e-2, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("acceptable_constr_viol_tol", dblValue);
  // Acceptance threshold for the dual infeasibility
  // Absolute tolerance on the dual infeasibility. "Acceptable" termination requires that the (max-norm of the unscaled) dual infeasibility is less than this 
  // threshold; see also acceptable_tol.
  bOk = getDoubleInPList(pStructOptions, L"acceptable_dual_inf_tol", &dblValue, 1e10, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("acceptable_dual_inf_tol", dblValue);
  // Acceptance threshold for the complementarity conditions
  // Absolute tolerance on the complementarity. "Acceptable" termination requires that the max-norm of the (unscaled) complementarity is less than this 
  // threshold; see also acceptable_tol.
  bOk = getDoubleInPList(pStructOptions, L"acceptable_compl_inf_tol", &dblValue, 1e-2, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("acceptable_compl_inf_tol", dblValue);
  // "Acceptance" stopping criterion based on objective function change.
  // If the relative change of the objective function (scaled by Max(1,|f(x)|)) is less than this value, this part of the acceptable tolerance termination 
  // is satisfied; see also acceptable_tol. This is useful for the quasi-Newton option, which has trouble to bring down the dual infeasibility.
  bOk = getDoubleInPList(pStructOptions, L"acceptable_obj_change_tol", &dblValue, 1e20, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("acceptable_obj_change_tol", dblValue);
  // Threshold for maximal value of primal iterates
  // If any component of the primal iterates exceeded this value (in absolute terms), the optimization is aborted with the exit message that the iterates
  // seem to be diverging.
  bOk = getDoubleInPList(pStructOptions, L"diverging_iterates_tol", &dblValue, 1e20, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("diverging_iterates_tol", dblValue);
  // Number of "acceptable" iterates before triggering termination.
  // If the algorithm encounters this many successive "acceptable" iterates (see "acceptable_tol"), it terminates, assuming that the problem 
  // has been solved to best possible accuracy given round-off. If it is set to zero, this heuristic is disabled.
  bOk = getIntInPList(pStructOptions, L"acceptable_iter", &iValue, 15, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("acceptable_iter", iValue);
  // Indicates if all variable bounds should be replaced by inequality
  // constraints
  //   This option must be set for the inexact algorithm
  // Possible values:
  //  - no                      [leave bounds on variables]
  //  - yes                     [replace variable bounds by inequality constraints]
  bOk = getStringInPList(pStructOptions, L"replace_bounds", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("replace_bounds", pStrValue);
  FREE(pStrValue);
  // Scaling threshold for the NLP error.
  //   (See paragraph after Eqn. (6) in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"s_max", &dblValue, 100, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("s_max", dblValue);
  // Maximum number of CPU seconds.
  //   A limit on CPU seconds that Ipopt can use to solve one problem.  If
  //   during the convergence check this limit is exceeded, Ipopt will terminate
  //   with a corresponding error message.
  bOk = getDoubleInPList(pStructOptions, L"max_cpu_time", &dblValue, 1e6, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("max_cpu_time", dblValue);

  /////////////////
  // NLP Scaling //
  /////////////////
  // Scaling factor for the objective function
  bOk = getDoubleInPList(pStructOptions, L"obj_scaling_factor", &dblValue, 1, iLog, CHECK_NONE);
  if (bOk) options->SetNumericValue("obj_scaling_factor", dblValue);
  // Select the technique used for scaling the NLP
  //  none: no problem scaling will be performed
  //  user-scaling: scaling parameters will come from the user
  //  gradient-based: scale the problem so the maximum gradient at the starting point is scaling_max_gradient
  //  equilibration-based: scale the problem so that first derivatives are of order 1 at random points (only available with MC19)
  // Selects the technique used for scaling the problem internally before it is solved. For user-scaling, the parameters come from the NLP.
  //  If you are using AMPL, they can be specified through suffixes ("scaling_factor")
  bOk = getStringInPList(pStructOptions, L"nlp_scaling_method", &pStrValue, "gradient-based", iLog, 
		       CHECK_VALUES, 4, 
		       "none",
		       "user_scaling",
		       "gradient-based",
		       "equilibration-based");
  if (bOk) options->SetStringValue("nlp_scaling_method", pStrValue);
  FREE(pStrValue);
  // Maximum gradient after scaling
  // This is the gradient scaling cut-off. If the maximum gradient is above this value, then gradient based scaling will be performed.
  // Scaling parameters are calculated to scale the maximum gradient back to this value. (This is g_max in Section 3.8 of the implementation paper.)
  //  Note: This option is only used if "nlp_scaling_method" is chosen as "gradient-based".
  bOk = getDoubleInPList(pStructOptions, L"nlp_scaling_max_gradient", &dblValue, 100, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("nlp_scaling_max_gradient", dblValue);
  // Target value for objective function gradient size.
  //   If a positive number is chosen, the scaling factor the objective function
  //   is computed so that the gradient has the max norm of the given size at
  //   the starting point.  This overrides nlp_scaling_max_gradient for the
  //   objective function.
  bOk = getDoubleInPList(pStructOptions, L"nlp_scaling_obj_target_gradient", &dblValue, 0, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("nlp_scaling_obj_target_gradient", dblValue);
  // Target value for constraint function gradient size.
  //   If a positive number is chosen, the scaling factor the constraint
  //   functions is computed so that the gradient has the max norm of the given
  //   size at the starting point.  This overrides nlp_scaling_max_gradient for
  //   the constraint functions.
  bOk = getDoubleInPList(pStructOptions, L"nlp_scaling_constr_target_gradient", &dblValue, 0, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("nlp_scaling_constr_target_gradient", dblValue);

  /////////
  // NLP //
  /////////
  // Factor for initial relaxation of the bounds
  // Before start of the optimization, the bounds given by the user are relaxed.  This option sets the factor for this relaxation. If it 
  // is set to zero, then then bounds relaxation is disabled. (See Eqn.(35) in implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"bound_relax_factor", &dblValue, 1e-8, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("bound_relax_factor", dblValue);
  // Indicates whether final points should be projected into original bounds.
  // * no: Leave final point unchanged
  // * yes: Project final point back into original bounds
  // Ipopt might relax the bounds during the optimization (see, e.g., option "bound_relax_factor"). This option determines whether the final 
  // point should be projected back into the user-provide original bounds after the optimization.
  bOk = getStringInPList(pStructOptions, L"honor_original_bounds", &pStrValue, "yes", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("honor_original_bounds", pStrValue);
  FREE(pStrValue);
  // Indicates whether it is desired to check for Nan/Inf in derivative matrices
  // * no: Don't check (faster).
  // * yes: Check Jacobians and Hessian for Nan and Inf.
  // Activating this option will cause an error if an invalid number is detected in the constraint Jacobians or the Lagrangian Hessian. 
  // If this is not activated, the test is skipped, and the algorithm might proceed with invalid numbers and fail.
  bOk = getStringInPList(pStructOptions, L"check_derivatives_for_naninf", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("check_derivatives_for_naninf", pStrValue);
  FREE(pStrValue);
  // any bound less or equal this value will be considered -inf (i.e. not lower bounded).
  bOk = getDoubleInPList(pStructOptions, L"nlp_lower_bound_inf", &dblValue, -1e19, iLog, CHECK_NONE);
  if (bOk) options->SetNumericValue("nlp_lower_bound_inf", dblValue);
  // any bound greater or this value will be considered +inf (i.e. not upper bounded).
  bOk = getDoubleInPList(pStructOptions, L"nlp_upper_bound_inf", &dblValue, 1e19, iLog, CHECK_NONE);
  if (bOk) options->SetNumericValue("nlp_upper_bound_inf", dblValue);
  // Determines how fixed variables should be handled.
  // * make_parameter: Remove fixed variable from optimization variables
  // * make_constraint: Add equality constraints fixing variables
  // * relax_bounds: Relax fixing bound constraints
  // The main difference between those options is that the starting point in the "make_constraint" case still has the fixed variables at 
  // their given values, whereas in the case "make_parameter" the functions are always evaluated with the fixed values for those variables.  
  // Also, for "relax_bounds", the fixing bound constraints are relaxed (according to "bound_relax_factor"). For both "make_constraints"
  //  and "relax_bounds", bound multipliers are computed for the fixed variables.
  bOk = getStringInPList(pStructOptions, L"fixed_variable_treatment", &pStrValue, "make_parameter", iLog, 
		       CHECK_VALUES, 3, 
		       "make_parameter",
		       "make_constraint",
		       "relax_bounds");
  if (bOk) options->SetStringValue("fixed_variable_treatment", pStrValue);
  FREE(pStrValue);
  //Indicates which linear solver should be used to detect linearly dependent equality constraints.
  // none: don't check; no extra work at beginning
  // mumps: use MUMPS
  // wsmp: use WSMP
  // ma28: use MA28
  // The default and available choices depend on how Ipopt has been compiled. This is experimental and does not work well.
  bOk = getStringInPList(pStructOptions, L"dependency_detector", &pStrValue, "none", iLog, 
		       CHECK_VALUES, 4, 
		       "none",
		       "mumps",
		       "wsmp",
		       "ma28");
  if (bOk) options->SetStringValue("dependency_detector", pStrValue);
  FREE(pStrValue);
  // Indicates if the right hand sides of the constraints should be considered during dependency detection
  // no: only look at gradients
  // yes: also consider right hand side
  bOk = getStringInPList(pStructOptions, L"dependency_detection_with_rhs", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("dependency_detection_with_rhs", pStrValue);
  FREE(pStrValue);
  // Number of linear variables
  // When the Hessian is approximated, it is assumed that the first num_linear_variables variables are linear. The Hessian is then not 
  // approximated in this space. If the get_number_of_nonlinear_variables method in the TNLP is implemented, this option is ignored.
  bOk = getIntInPList(pStructOptions, L"num_linear_variables", &iValue, 0, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("num_linear_variables", iValue);
  // Weight for linear damping term (to handle one-sided bounds).
  //   (see Section 3.7 in implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"kappa_d", &dblValue, 1e-5, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("kappa_d", dblValue);
  // Indicates whether all equality constraints are linear
  // * no: Don't assume that all equality constraints are linear
  // * yes: Assume that equality constraints Jacobian are constant
  // Activating this option will cause Ipopt to ask for the Jacobian of the equality constraints only once from the NLP and reuse this information later.
  bOk = getStringInPList(pStructOptions, L"jac_c_constant", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("jac_c_constant", pStrValue);
  FREE(pStrValue);
  // Indicates whether all inequality constraints are linear
  // * no: Don't assume that all inequality constraints are linear
  // * yes: Assume that equality constraints Jacobian are constant
  // Activating this option will cause Ipopt to ask for the Jacobian of the inequality constraints only once from the NLP and reuse this information later.
  bOk = getStringInPList(pStructOptions, L"jac_d_constant", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("jac_d_constant", pStrValue);
  FREE(pStrValue);
  // Indicates whether the problem is a quadratic problem
  // * no: Assume that Hessian changes
  // * yes: Assume that Hessian is constant
  // Activating this option will cause Ipopt to ask for the Hessian of the Lagrangian function only once from the NLP and reuse this information later.
  bOk = getStringInPList(pStructOptions, L"hessian_constant", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("hessian_constant", pStrValue);
  FREE(pStrValue);

  ////////////////////
  // Initialization //
  ////////////////////
  // Desired minimal relative distance of initial point to bound
  // Determines how much the initial point might have to be modified in order to be sufficiently inside the bounds (together with "bound_push").
  // (This is kappa_2 in Section 3.6 of implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"bound_frac", &dblValue, 1e-2, iLog, CHECK_BOTH, 0., 0.5);
  if (bOk) options->SetNumericValue("bound_frac", dblValue);
  // Desired minimal absolute distance of initial point to bound
  // Determines how much the initial point might have to be modified in order to be sufficiently inside the bounds (together with "bound_frac").
  // (This is kappa_1 in "Section 3.6 of implementation paper.
  bOk = getDoubleInPList(pStructOptions, L"bound_push", &dblValue, 1e-2, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("bound_push", dblValue);
  // Desired minimal relative distance of initial slack to bound
  // Determines how much the initial slack variables might have to be modified in order to be sufficiently inside the inequality bounds
  // (together with "slack_bound_push").  (This is kappa_2 in Section 3.6 of implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"slack_bound_frac", &dblValue, 1e-2, iLog, CHECK_BOTH, 0., 0.5);
  if (bOk) options->SetNumericValue("slack_bound_frac", dblValue);
  // Desired minimal absolute distance of initial slack to bound
  // Determines how much the initial slack variables might have to be modified in order to be sufficiently inside the inequality bounds
  //  (together with "slack_bound_frac").  (This is kappa_1 in Section 3.6 of implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"slack_bound_push", &dblValue, 1e-2, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("slack_bound_push", dblValue);
  // Initial value for the bound multipliers
  // All dual variables corresponding to bound constraints are initialized to this value.
  bOk = getDoubleInPList(pStructOptions, L"bound_mult_init_val", &dblValue, 1, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("bound_mult_init_val", dblValue);
  // Maximal allowed least-square guess of constraint multipliers
  // Determines how large the initial least-square guesses of the constraint multipliers are allowed to be (in max-norm). If the guess is larger
  // than this value, it is discarded and all constraint multipliers are set to zero. This options is also used when initializing the 
  // restoration phase. By default, "resto.constr_mult_init_max" (the one used in RestoIterateInitializer) is set to zero.
  bOk = getDoubleInPList(pStructOptions, L"constr_mult_init_max", &dblValue, 1000, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("constr_mult_init_max", dblValue);
  // Initialization method for bound multipliers
  // * constant: set all bound multipliers to the value of bound_mult_init_val
  // * mu-based: initialize to mu_init/x_slack
  // This option defines how the iterates for the bound multipliers are initialized.  If "constant" is chosen, then all bound multipliers 
  // are initialized to the value of "bound_mult_init_val".  If "mu-based" is chosen, the each value is initialized to the the value 
  // of "mu_init" divided by the corresponding slack variable. This latter option might be useful if the starting point is close to the
  // optimal solution.
  bOk = getStringInPList(pStructOptions, L"bound_mult_init_method", &pStrValue, "constant", iLog, 
		       CHECK_VALUES, 2, 
		       "constant",
		       "mu-based");
  if (bOk) options->SetStringValue("bound_mult_init_method", pStrValue);
  FREE(pStrValue);
  // Least square initialization of the primal variables
  // no: take user-provided point
  // yes: overwrite user-provided point with least-square estimates
  // If set to yes, Ipopt ignores the user provided point and solves a least square problem for the primal variables (x and s), to fit the 
  // linearized equality and inequality constraints.  This might be useful if the user doesn't know anything about the starting point, or for 
  // solving an LP or QP.
  bOk = getStringInPList(pStructOptions, L"least_square_init_primal", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("least_square_init_primal", pStrValue);
  FREE(pStrValue);
  // Least square initialization of all dual variables
  // no: use bound_mult_init_val and least-square equality constraint multipliers
  // yes: overwrite user-provided point with least-square estimates
  // If set to yes, Ipopt tries to compute least-square multipliers (considering ALL dual variables). If successful, the bound 
  // multipliers are possibly corrected to be at least bound_mult_init_val. This might be useful if the user doesn't know anything 
  // about the starting point, or for solving an LP or QP. This overwrites option "bound_mult_init_method".
  bOk = getStringInPList(pStructOptions, L"least_square_init_duals", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("least_square_init_duals", pStrValue);
  FREE(pStrValue);

  ///////////////////////
  // Barrier Parameter //
  ///////////////////////
  // Indicates if we want to do Mehrotra's algorithm.
  // * no: Do the usual Ipopt algorithm.
  // * yes: Do Mehrotra's predictor-corrector algorithm.
  // If set to yes, Ipopt runs as Mehrotra's predictor-corrector algorithm. This works usually very well for LPs and convex QPs. This
  // automatically disables the line search, and chooses the (unglobalized) adaptive mu strategy with the "probing" oracle, and uses 
  // "corrector_type=affine" without any safeguards; you should not set any of those options explicitly in addition. Also, unless 
  // otherwise specified, the values of "bound_push", "bound_frac", and "bound_mult_init_val" are set more aggressive, and sets 
  // "alpha_for_y=bound_mult".
  bOk = getStringInPList(pStructOptions, L"mehrotra_algorithm", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("mehrotra_algorithm", pStrValue);
  FREE(pStrValue);
  // Indicates if the linear system should be solved quickly.
  //   If set to yes, the algorithm assumes that the linear system that is
  //   solved to obtain the search direction, is solved sufficiently well. In
  //   that case, no residuals are computed, and the computation of the search
  //   direction is a little faster.
  // Possible values:
  //  - no                      [Verify solution of linear system by computing residuals.]
  //  - yes                     [Trust that linear systems are solved well.]
  bOk = getStringInPList(pStructOptions, L"fast_step_computation", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("fast_step_computation", pStrValue);
  FREE(pStrValue);
  // Update strategy for barrier parameter
  // * monotone: use the monotone (Fiacco-McCormick) strategy
  // * adaptive: use the adaptive update strategy
  // Determines which barrier parameter update strategy is to be used.
  bOk = getStringInPList(pStructOptions, L"mu_strategy", &pStrValue, "monotone", iLog, 
		       CHECK_VALUES, 2,
		       "monotone",
		       "adaptive");
  if (bOk) options->SetStringValue("mu_strategy", pStrValue);
  FREE(pStrValue);
  // Oracle for a new barrier parameter in the adaptive strategy
  // * probing: Mehrotra's probing heuristic
  // * loqo: LOQO's centrality rule
  // * quality-function: minimize a quality function
  // Determines how a new barrier parameter is computed in each "free-mode" iteration of the adaptive barrier parameter strategy.
  // (Only considered if "adaptive" is selected for option "mu_strategy").
  bOk = getStringInPList(pStructOptions, L"mu_oracle", &pStrValue, "quality-function", iLog, 
		       CHECK_VALUES, 3,
		       "probing",
		       "loqo",
		       "quality-function");
  if (bOk) options->SetStringValue("mu_oracle", pStrValue);
  FREE(pStrValue);
  bOk = getIntInPList(pStructOptions, L"quality_function_max_section_steps", &iValue, 8, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("quality_function_max_section_steps", iValue);
  // Oracle for the barrier parameter when switching to fixed mode.
  // * probing: Mehrotra's probing heuristic
  // * loqo: LOQO's centrality rule
  // * quality-function: minimize a quality function
  // * average_compl: base on current average complementarity
  // Determines how the first value of the barrier parameter should be computed when switching to the "monotone mode" in the adaptive strategy. 
  // (Only considered if "adaptive" is selected for option "mu_strategy")
  bOk = getStringInPList(pStructOptions, L"fixed_mu_oracle", &pStrValue, "average_compl", iLog, 
		       CHECK_VALUES, 4,
		       "probing",
		       "loqo",
		       "quality-function",
		       "average_compl");
  if (bOk) options->SetStringValue("fixed_mu_oracle", pStrValue);
  FREE(pStrValue);
  // Initial value for the barrier parameter
  // This option determines the initial value for the barrier parameter (mu).It is only relevant in the monotone, Fiacco-McCormick 
  // version of the algorithm. (i.e., if "mu_strategy" is chosen as "monotone")
  bOk = getDoubleInPList(pStructOptions, L"mu_init", &dblValue, 0.1, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("mu_init", dblValue);
  // Factor for initialization of maximum value for barrier parameter.
  // This option determines the upper bound on the barrier parameter. This upper bound is computed as the average complementarity at the initial 
  // point times the value of this option. (Only used if option "mu_strategy" is chosen as "adaptive".)
  bOk = getDoubleInPList(pStructOptions, L"mu_max_fact", &dblValue, 1000, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("mu_max_fact", dblValue);
  // Maximal value for barrier parameter for adaptive strategy
  // This option specifies an upper bound on the barrier parameter in the adaptive mu selection mode. If this option is set, it overwrites the 
  // effect of mu_max_fact. (Only used if option "mu_strategy" is chosen as "adaptive".)
  bOk = getDoubleInPList(pStructOptions, L"mu_max", &dblValue, 100000, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("mu_max", dblValue);
  // Minimum value for barrier parameter.
  // This option specifies the lower bound on the barrier parameter in the adaptive mu selection mode. By default, it is set to the minimum of 1e-11 and 
  // min("tol","compl_inf_tol")/("barrier_tol_factor"+1), which should be a reasonable value. (Only used if option "mu_strategy" is chosen as "adaptive".)
  bOk = getDoubleInPList(pStructOptions, L"mu_min", &dblValue, 1e-11, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("mu_min", dblValue);
  // Factor for mu in barrier stop test.
  // The convergence tolerance for each barrier problem in the monotone mode is the value of the barrier parameter times "barrier_tol_factor".
  // This option is also used in the adaptive mu strategy during the monotone mode. (This is kappa_epsilon in implementation paper).
  bOk = getDoubleInPList(pStructOptions, L"barrier_tol_factor", &dblValue, 10, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("barrier_tol_factor", dblValue);
  // Determines linear decrease rate of barrier parameter.
  // For the Fiacco-McCormick update procedure the new barrier parameter mu is obtained by taking the minimum of mu*"mu_linear_decrease_factor"
  // and mu^"superlinear_decrease_power". (This is kappa_mu in implementation paper.) This option is also used in the adaptive mu 
  // strategy during the monotone mode.
  bOk = getDoubleInPList(pStructOptions, L"mu_linear_decrease_factor", &dblValue, 0.2, iLog, CHECK_BOTH, 0., 1.);
  if (bOk) options->SetNumericValue("mu_linear_decrease_factor", dblValue);
  // Determines superlinear decrease rate of barrier parameter.
  // For the Fiacco-McCormick update procedure the new barrier parameter mu is obtained by taking the minimum of mu*"mu_linear_decrease_factor"
  // and mu^"superlinear_decrease_power". (This is theta_mu in implementation paper.) This option is also used in the adaptive mu 
  // strategy during the monotone mode.
  bOk = getDoubleInPList(pStructOptions, L"mu_superlinear_decrease_power", &dblValue, 1.5, iLog, CHECK_BOTH, 1., 2.);
  if (bOk) options->SetNumericValue("mu_superlinear_decrease_power", dblValue);
  // Allow skipping of barrier problem if barrier test is already met.
  // no: Take at least one iteration per barrier problem
  // yes: Allow fast decrease of mu if barrier test it met
  // If set to "no", the algorithm enforces at least one iteration per barrier problem, even if the barrier test is already met for the updated barrier parameter.
  bOk = getStringInPList(pStructOptions, L"mu_allow_fast_monotone_decrease", &pStrValue, "yes", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("mu_allow_fast_monotone_decrease", pStrValue);
  FREE(pStrValue);
  // Lower bound on fraction-to-the-boundary parameter tau
  // (This is tau_min in the implementation paper.)  This option is also used in the adaptive mu strategy during the monotone mode.
  bOk = getDoubleInPList(pStructOptions, L"tau_min", &dblValue, 0.99, iLog, CHECK_BOTH, 0., 1.);
  if (bOk) options->SetNumericValue("tau_min", dblValue);
  // Maximum value of the centering parameter.
  //   This is the upper bound for the centering parameter chosen by the quality
  //   function based barrier parameter update. (Only used if option "mu_oracle"
  //   is set to "quality-function".)
  bOk = getDoubleInPList(pStructOptions, L"sigma_max", &dblValue, 100, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("sigma_max", dblValue);
  // Minimum value of the centering parameter.
  //   This is the lower bound for the centering parameter chosen by the quality
  //   function based barrier parameter update. (Only used if option "mu_oracle"
  //   is set to "quality-function".)
  bOk = getDoubleInPList(pStructOptions, L"sigma_min", &dblValue, 100, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("sigma_min", dblValue);
  // Norm used for components of the quality function.
  //   (Only used if option "mu_oracle" is set to "quality-function".)
  // Possible values:
  //  - 1-norm                  [use the 1-norm (abs sum)]
  //  - 2-norm-squared          [use the 2-norm squared (sum of squares)]
  //  - max-norm                [use the infinity norm (max)]
  //  - 2-norm                  [use 2-norm]
  bOk = getStringInPList(pStructOptions, L"quality_function_norm_type", &pStrValue, "2-norm-squared", iLog, 
		       CHECK_VALUES, 4, 
		       "1-norm", 
		       "2-norm-squared",
		       "max-norm",
		       "2-norm");
  if (bOk) options->SetStringValue("quality_function_norm_type", pStrValue);
  FREE(pStrValue);
  // The penalty term for centrality that is included in quality function.
  //   This determines whether a term is added to the quality function to
  //   penalize deviation from centrality with respect to complementarity.  The
  //   complementarity measure here is the xi in the Loqo update rule. (Only
  //   used if option "mu_oracle" is set to "quality-function".)
  // Possible values:
  //  - none                    [no penalty term is added]
  //  - log                     [complementarity * the log of the centrality measure]
  //  - reciprocal              [complementarity * the reciprocal of the centrality measure]
  //  - cubed-reciprocal        [complementarity * the reciprocal of the centrality measure cubed]
  bOk = getStringInPList(pStructOptions, L"quality_function_centrality", &pStrValue, "none", iLog, 
		       CHECK_VALUES, 4, 
		       "none", 
		       "log",
		       "reciprocal",
		       "cubed-reciprocal");
  if (bOk) options->SetStringValue("quality_function_centrality", pStrValue);
  FREE(pStrValue);
  // The balancing term included in the quality function for centrality.
  //   This determines whether a term is added to the quality function that
  //   penalizes situations where the complementarity is much smaller than dual
  //   and primal infeasibilities. (Only used if option "mu_oracle" is set to
  //   "quality-function".)
  // Possible values:
  //  - none                    [no balancing term is added]
  //  - cubic                   [Max(0,Max(dual_inf,primal_inf)-compl)^3]
  bOk = getStringInPList(pStructOptions, L"quality_function_balancing_term", &pStrValue, "none", iLog, 
		       CHECK_VALUES, 2, 
		       "none", 
		       "cubic");
  if (bOk) options->SetStringValue("quality_function_balancing_term", pStrValue);
  FREE(pStrValue);
  // Tolerance for the section search procedure determining the optimal
  // centering parameter (in sigma space).
  //   The golden section search is performed for the quality function based mu
  //   oracle. (Only used if option "mu_oracle" is set to "quality-function".)
  bOk = getDoubleInPList(pStructOptions, L"quality_function_section_sigma_tol", &dblValue, 0.01, iLog, CHECK_BOTH, 0., 1.);
  if (bOk) options->SetNumericValue("quality_function_section_sigma_tol", dblValue);
  // Tolerance for the golden section search procedure determining the optimal
  // centering parameter (in the function value space).
  //   The golden section search is performed for the quality function based mu
  //   oracle. (Only used if option "mu_oracle" is set to "quality-function".)
  bOk = getDoubleInPList(pStructOptions, L"quality_function_section_qf_tol", &dblValue, 0, iLog, CHECK_BOTH, 0., 1.);
  if (bOk) options->SetNumericValue("quality_function_section_qf_tol", dblValue);
  // Factor limiting the deviation of dual variables from primal estimates.
  // If the dual variables deviate from their primal estimates, a correction is performed. (See Eqn. (16) in the implementation paper.) 
  // Setting the value to less than 1 disables the correction.
  bOk = getDoubleInPList(pStructOptions, L"kappa_sigma", &dblValue, 1e10, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("kappa_sigma", dblValue);
  // Globalization method used in backtracking line search
  // filter: Filter method
  // cg-penalty: Chen-Goldfarb penalty function
  // penalty: Standard penalty function
  bOk = getStringInPList(pStructOptions, L"line_search_method", &pStrValue, "cg-penalty", iLog, 
		       CHECK_VALUES, 3, 
		       "filter",
		       "cg-penalty",
		       "penalty");
  if (bOk) options->SetStringValue("line_search_method", pStrValue);
  FREE(pStrValue);
  // Globalization strategy for the adaptive mu selection mode.
  // kkt-error: nonmonotone decrease of kkt-error
  // obj-constr-filter: 2-dim filter for objective and constraint violation
  // never-monotone-mode: disables globalization
  // To achieve global convergence of the adaptive version, the algorithm has to switch to the monotone mode (Fiacco-McCormick approach) when 
  // convergence does not seem to appear.  This option sets the criterion used to decide when to do this switch. (Only used if option 
  // "mu_strategy" is chosen as "adaptive".)
  bOk = getStringInPList(pStructOptions, L"adaptive_mu_globalization", &pStrValue, "obj-constr-filter", iLog, 
		       CHECK_VALUES, 3, 
		       "kkt-error",
		       "obj-constr-filter",
		       "never-monotone-mode");
  if (bOk) options->SetStringValue("adaptive_mu_globalization", pStrValue);
  FREE(pStrValue);
  // Maximum number of iterations requiring sufficient progress.
  // For the "kkt-error" based globalization strategy, sufficient progress must be made for "adaptive_mu_kkterror_red_iters
  // iterations. If this number of iterations is exceeded, the globalization strategy switches to the monotone mode.
  bOk = getIntInPList(pStructOptions, L"adaptive_mu_kkterror_red_iters", &iValue, 4, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("adaptive_mu_kkterror_red_iters", iValue);
  // Sufficient decrease factor for "kkt-error" globalization strategy.
  // For the "kkt-error" based globalization strategy, the error must decrease by this factor to be deemed sufficient decrease.
  bOk = getDoubleInPList(pStructOptions, L"adaptive_mu_kkterror_red_fact", &dblValue, 0.9999, iLog, CHECK_BOTH, 0., 1.);
  if (bOk) options->SetNumericValue("adaptive_mu_kkterror_red_fact", dblValue);
  // Factor determining width of margin for obj-constr-filter adaptive globalization strategy.
  // When using the adaptive globalization strategy, "obj-constr-filter" sufficient progress for a filter entry is defined as 
  // follows: (new obj) < (filter obj) - filter_margin_fact*(new constr-viol) OR (new constr-viol) < (filter constr-viol) - 
  // filter_margin_fact*(new constr-viol).  For the description of the "kkt-error-filter" option see "filter_max_margin.
  bOk = getDoubleInPList(pStructOptions, L"filter_margin_fact", &dblValue, 1e-5, iLog, CHECK_BOTH, 0., 1.);
  if (bOk) options->SetNumericValue("filter_margin_fact", dblValue);
  // Maximum width of margin in obj-constr-filter adaptive globalization strategy.
  bOk = getDoubleInPList(pStructOptions, L"filter_max_margin", &dblValue, 1.0, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("filter_max_margin", dblValue);
  // Indicates if the previous iterate should be restored if the monotone mode is entered.
  // no: don't restore accepted iterate
  // yes: restore accepted iterate
  // When the globalization strategy for the adaptive barrier algorithm switches to the monotone mode, it can either start 
  // from the most recent iterate (no), or from the last iterate that was accepted (yes).
  bOk = getStringInPList(pStructOptions, L"adaptive_mu_restore_previous_iterate", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("adaptive_mu_restore_previous_iterate", pStrValue);
  FREE(pStrValue);
  // Determines the initial value of the barrier parameter when switching to the monotone mode.
  // When the globalization strategy for the adaptive barrier algorithm switches to the monotone mode and fixed_mu_oracle is chosen as 
  // "average_compl", the barrier parameter is set to the current average complementarity times the value of "adaptive_mu_monotone_init_factor"
  bOk = getDoubleInPList(pStructOptions, L"adaptive_mu_monotone_init_factor", &dblValue, 0.8, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("adaptive_mu_monotone_init_factor", dblValue);
  // Norm used for the KKT error in the adaptive mu globalization strategies.
  // 1-norm: use the 1-norm (abs sum)
  // 2-norm-squared: use the 2-norm squared (sum of squares)
  // max-norm: use the infinity norm (max)
  // 2-norm: use 2-norm
  // When computing the KKT error for the globalization strategies, the norm to be used is specified with this option. Note, this options is also used 
  // in the QualityFunctionMuOracle.
  bOk = getStringInPList(pStructOptions, L"adaptive_mu_kkt_norm_type", &pStrValue, "2-norm-squared", iLog, 
		       CHECK_VALUES, 4, 
		       "1-norm",
		       "2-norm-squared",
		       "max-norm",
		       "2-norm");
  if (bOk) options->SetStringValue("adaptive_mu_kkt_norm_type", pStrValue);
  FREE(pStrValue);

  ////////////////////////
  // Multiplier Updates //
  ////////////////////////

  // Fractional reduction of the trial step size in the backtracking line search.
  //   At every step of the backtracking line search, the trial step size is
  //   reduced by this factor.
  bOk = getDoubleInPList(pStructOptions, L"alpha_red_factor", &dblValue, 0.5, iLog, CHECK_BOTH, 0., 1.);
  if (bOk) options->SetNumericValue("alpha_red_factor", dblValue);
  // Always accept the first trial step.
  //   Setting this option to "yes" essentially disables the line search and
  //   makes the algorithm take aggressive steps, without global convergence guarantees.
  // Possible values:
  //  - no                      [don't arbitrarily accept the full step]
  //  - yes                     [always accept the full step]
  bOk = getStringInPList(pStructOptions, L"accept_every_trial_step", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("accept_every_trial_step", pStrValue);
  FREE(pStrValue);
  // Accept a trial point after maximal this number of steps.
  //   Even if it does not satisfy line search conditions.
  bOk = getDoubleInPList(pStructOptions, L"accept_after_max_steps", &dblValue, -1, iLog, CHECK_MIN, -1.0);
  if (bOk) options->SetNumericValue("accept_after_max_steps", dblValue);

  // Method to determine the step size for constraint multipliers.
  // * primal: use primal step size
  // * bound_mult: use step size for the bound multipliers (good for LPs)
  // * min: use the min of primal and bound multipliers
  // * max: use the max of primal and bound multipliers
  // * full: take a full step of size one
  // * min_dual_infeas: choose step size minimizing new dual infeasibility
  // * safe_min_dual_infeas: like "min_dual_infeas", but safeguarded by "min" and "max"
  // * primal-and-full: use the primal step size, and full step if delta_x <= alpha_for_y_tol
  // * dual-and-full: use the dual step size, and full step if delta_x <= alpha_for_y_tol
  // * acceptor: Call LSAcceptor to get step size for y
  // This option determines how the step size (alpha_y) will be calculated when updating the constraint multipliers.
  bOk = getStringInPList(pStructOptions, L"alpha_for_y", &pStrValue, "primal", iLog, 
		       CHECK_VALUES, 10, 
		       "primal",
		       "bound-mult",
		       "min",
		       "max",
		       "full",
		       "min-dual-infeas",
		       "safer-min-dual-infeas",
		       "primal-and-full",
		       "dual-and-full",
		       "acceptor");
  if (bOk) options->SetStringValue("alpha_for_y", pStrValue);
  FREE(pStrValue);
  // Tolerance for switching to full equality multiplier steps
  // This is only relevant if "alpha_for_y" is chosen "primal-and-full" or "dual-and-full". The step size for the equality constraint 
  // multipliers is taken to be one if the max-norm of the primal step is less than this tolerance
  bOk = getDoubleInPList(pStructOptions, L"alpha_for_y_tol", &dblValue, 10.0, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("alpha_for_y_tol", dblValue);
  // Tolerance for detecting numerically insignificant steps
  // If the search direction in the primal variables (x and s) is, in relative terms for each component, less than this value, the 
  // algorithm accepts the full step without line search. If this happens repeatedly, the algorithm will terminate with a corresponding exit 
  // message. The default value is 10 times machine precision.  
  bOk = getDoubleInPList(pStructOptions, L"tiny_step_tol", &dblValue, 10.0*std::numeric_limits<double>::epsilon(), iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("tiny_step_tol", dblValue);
  // Tolerance for quitting because of numerically insignificant steps.
  // If the search direction in the primal variables (x and s) is, in relative terms for each component, repeatedly less than tiny_step_tol,
  // and the step in the y variables is smaller than this threshold, the algorithm will terminate.
  bOk = getDoubleInPList(pStructOptions, L"tiny_step_y_tol", &dblValue, 1e-2, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("tiny_step_y_tol", dblValue);
  // Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates.
  // no: use the Newton step to update the multipliers
  // yes: use least-square multiplier estimates
  // This asks the algorithm to recompute the multipliers, whenever the current infeasibility is less than recalc_y_feas_tol. 
  // Choosing yes might be helpful in the quasi-Newton option.  However, each recalculation requires an extra factorization of the linear 
  // system.  If a limited memory quasi-Newton option is chosen, this is used by default.
  bOk = getStringInPList(pStructOptions, L"recalc_y", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("recalc_y", pStrValue);
  FREE(pStrValue);
  // Feasibility threshold for recomputation of multipliers.
  // If recalc_y is chosen and the current infeasibility is less than this value, then the multipliers are recomputed.
  bOk = getDoubleInPList(pStructOptions, L"recalc_y_feas_tol", &dblValue, 1e-6, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("recalc_y_feas_tol", dblValue);
  // Correction size for very small slacks.
  //   Due to numerical issues or the lack of an interior, the slack variables
  //   might become very small.  If a slack becomes very small compared to
  //   machine precision, the corresponding bound is moved slightly.  This
  //   parameter determines how large the move should be.  Its default value is
  //   mach_eps^{3/4}.  (See also end of Section 3.5 in implementation paper -
  //   but actual implementation might be somewhat different.)  
  bOk = getDoubleInPList(pStructOptions, L"slack_move", &dblValue, 1.81899e-12, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("slack_move", dblValue);

  /////////////////
  // Line Search //
  /////////////////
  // Maximal number of second order correction trial steps
  bOk = getIntInPList(pStructOptions, L"max_soc", &iValue, 4, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("max_soc", iValue);
  // Factor in the sufficient reduction rule for second order correction.
  //   This option determines how much a second order correction step must
  //   reduce the constraint violation so that further correction steps are
  //   attempted.  (See Step A-5.9 of Algorithm A in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"kappa_soc", &dblValue, 0.99, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("kappa_soc", dblValue);
  // Determines the upper bound on the acceptable increase of barrier objective function.
  //   Trial points are rejected if they lead to an increase in the barrier
  //   objective function by more than obj_max_inc orders of magnitude.
  bOk = getDoubleInPList(pStructOptions, L"obj_max_inc", &dblValue, 5, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("obj_max_inc", dblValue);
  // Maximal allowed number of filter resets
  //   A positive number enables a heuristic that resets the filter, whenever in
  //   more than "filter_reset_trigger" successive iterations the last rejected
  //   trial steps size was rejected because of the filter.  This option
  //   determine the maximal number of resets that are allowed to take place.
  bOk = getIntInPList(pStructOptions, L"max_filter_resets", &iValue, 5, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("max_filter_resets", iValue);
  // Number of iterations that trigger the filter reset.
  //   If the filter reset heuristic is active and the number of successive
  //   iterations in which the last rejected trial step size was rejected
  //   because of the filter, the filter is reset.
  bOk = getIntInPList(pStructOptions, L"filter_reset_trigger", &iValue, 5, iLog, CHECK_MIN, 1);
  if (bOk) options->SetIntegerValue("filter_reset_trigger", iValue);
  // Trigger counter for watchdog procedure
  // If the number of successive iterations in which the backtracking line search did not accept the first trial point exceeds this number, the 
  // watchdog procedure is activated.  Choosing "0" here disables the watchdog procedure.
  bOk = getIntInPList(pStructOptions, L"watchdog_shortened_iter_trigger", &iValue, 10, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("watchdog_shortened_iter_trigger", iValue);
  // Maximum number of watchdog iterations.
  // This option determines the number of trial iterations allowed before the watchdog procedure is aborted and the algorithm returns to the stored point.
  bOk = getIntInPList(pStructOptions, L"watchdog_trial_iter_max", &iValue, 3, iLog, CHECK_MIN, 1);
  if (bOk) options->SetIntegerValue("watchdog_trial_iter_max", iValue);
  // Determines upper bound for constraint violation in the filter.
  //   The algorithmic parameter theta_max is determined as theta_max_fact times
  //   the maximum of 1 and the constraint violation at initial point.  Any
  //   point with a constraint violation larger than theta_max is unacceptable
  //   to the filter (see Eqn. (21) in the implementation paper).
  bOk = getDoubleInPList(pStructOptions, L"theta_max_fact", &dblValue, 10000, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("theta_max_fact", dblValue);
  // Determines constraint violation threshold in the switching rule.
  //   The algorithmic parameter theta_min is determined as theta_min_fact times
  //   the maximum of 1 and the constraint violation at initial point.  The
  //   switching rules treats an iteration as an h-type iteration whenever the
  //   current constraint violation is larger than theta_min (see paragraph
  //   before Eqn. (19) in the implementation paper).
  bOk = getDoubleInPList(pStructOptions, L"theta_min_fact", &dblValue, 0.0001, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("theta_min_fact", dblValue);
  // Relaxation factor in the Armijo condition.
  //   (See Eqn. (20.) in the implementation paper)
  bOk = getDoubleInPList(pStructOptions, L"eta_phi", &dblValue, 1e-8, iLog, CHECK_BOTH, 0., 0.5);
  if (bOk) options->SetNumericValue("eta_phi", dblValue);
  // Multiplier for constraint violation in the switching rule.
  //   (See Eqn. (19) in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"delta", &dblValue, 1, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("delta", dblValue);
  // Exponent for linear barrier function model in the switching rule.
  //   (See Eqn. (19) in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"s_phi", &dblValue, 2.3, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("s_phi", dblValue);
  // Exponent for current constraint violation in the switching rule.
  //   (See Eqn. (19) in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"s_theta", &dblValue, 1.1, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("s_theta", dblValue);
  // Relaxation factor in the filter margin for the barrier function.
  //   (See Eqn. (18a) in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"gamma_phi", &dblValue, 1e-8, iLog, CHECK_BOTH, 0., 1.);
  if (bOk) options->SetNumericValue("gamma_phi", dblValue);
  // Relaxation factor in the filter margin for the constraint violation.
  //   (See Eqn. (18b) in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"gamma_theta", &dblValue, 1e-5, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("gamma_theta", dblValue);
  // Safety factor for the minimal step size (before switching to restoration phase).
  //   (This is gamma_alpha in Eqn. (20.) in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"alpha_min_frac", &dblValue, 0.05, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("alpha_min_frac", dblValue);
  // * none: no corrector
  // * affine: corrector step towards mu=0
  // * primal-dual: corrector step towards current mu
  bOk = getStringInPList(pStructOptions, L"corrector_type", &pStrValue, "none", iLog, 
		       CHECK_VALUES, 3,
		       "none",
		       "affine",
		       "primal-dual");
  if (bOk) options->SetStringValue("corrector_type", pStrValue);
  FREE(pStrValue);
  // Skip the corrector step in negative curvature iteration (unsupported!).
  //   The corrector step is not tried if negative curvature has been
  //   encountered during the computation of the search direction in the current
  //   iteration. This option is only used if "mu_strategy" is "adaptive".
  // Possible values:
  //  - no                      [don't skip]
  //  - yes                     [skip]
  bOk = getStringInPList(pStructOptions, L"skip_corr_if_neg_curv", &pStrValue, "yes", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("skip_corr_if_neg_curv", pStrValue);
  FREE(pStrValue);
  // Skip the corrector step during monotone barrier parameter mode (unsupported!).
  //   The corrector step is not tried if the algorithm is currently in the
  //   monotone mode (see also option "barrier_strategy").This option is only
  //   used if "mu_strategy" is "adaptive".
  // Possible values:
  //  - no                      [don't skip]
  //  - yes                     [skip]
  bOk = getStringInPList(pStructOptions, L"skip_corr_in_monotone_mode", &pStrValue, "yes", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("skip_corr_in_monotone_mode", pStrValue);
  FREE(pStrValue);
  // Complementarity tolerance factor for accepting corrector step (unsupported!).
  //   This option determines the factor by which complementarity is allowed to
  //   increase for a corrector step to be accepted.
  bOk = getDoubleInPList(pStructOptions, L"corrector_compl_avrg_red_fact", &dblValue, 1, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("corrector_compl_avrg_red_fact", dblValue);
  // Initial value of the penalty parameter.
  bOk = getDoubleInPList(pStructOptions, L"nu_init", &dblValue, 1e-6, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("nu_init", dblValue);
  // Increment of the penalty parameter.
  bOk = getDoubleInPList(pStructOptions, L"nu_inc", &dblValue, 0.0001, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("nu_inc", dblValue);
  // Value in penalty parameter update formula.
  bOk = getDoubleInPList(pStructOptions, L"rho", &dblValue, 0.1, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("rho", dblValue);

  ////////////////
  // Warm Start //
  ////////////////
  // Enables to specify bound multiplier values
  // * no: do not use the warm start initialization
  // * yes: use the warm start initialization
  bOk = getStringInPList(pStructOptions, L"warm_start_init_point", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("warm_start_init_point", pStrValue);
  FREE(pStrValue);
  // Indicates whether a problem with a structure identical to the previous one is to be solved.
  //   If "yes" is chosen, then the algorithm assumes that an NLP is now to be
  //   solved, whose structure is identical to one that already was considered
  //   (with the same NLP object).
  // Possible values:
  //  - no                      [Assume this is a new problem.]
  //  - yes                     [Assume this is problem has known structure]
  bOk = getStringInPList(pStructOptions, L"warm_start_same_structure", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("warm_start_same_structure", pStrValue);
  FREE(pStrValue);
  // Enables to specify how much should variables should be pushed inside the feasible region
  bOk = getDoubleInPList(pStructOptions, L"warm_start_bound_push", &dblValue, 1e-3, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("warm_start_bound_push", dblValue);
  bOk = getDoubleInPList(pStructOptions, L"warm_start_bound_frac", &dblValue, 1e-3, iLog, CHECK_BOTH, 0., 0.5);
  if (bOk) options->SetNumericValue("warm_start_bound_frac", dblValue);
  bOk = getDoubleInPList(pStructOptions, L"warm_start_slack_bound_frac", &dblValue, 1e-3, iLog, CHECK_BOTH, 0., 0.5);
  if (bOk) options->SetNumericValue("warm_start_slack_bound_frac", dblValue);
  bOk = getDoubleInPList(pStructOptions, L"warm_start_slack_bound_push", &dblValue, 1e-3, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("warm_start_slack_bound_push", dblValue);
  // Enables to specify how much should bound multipliers should be pushed inside the feasible region
  bOk = getDoubleInPList(pStructOptions, L"warm_start_mult_bound_push", &dblValue, 1e-3, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("warm_start_mult_bound_push", dblValue);
  bOk = getDoubleInPList(pStructOptions, L"warm_start_mult_init_max", &dblValue, 1e6, iLog, CHECK_NONE);
  if (bOk) options->SetNumericValue("warm_start_mult_init_max", dblValue);
  // Tells algorithm whether to use the GetWarmStartIterate method in the NLP.
  // Possible values:
  //  - no                      [call GetStartingPoint in the NLP]
  //  - yes                     [call GetWarmStartIterate in the NLP]
  bOk = getStringInPList(pStructOptions, L"warm_start_entire_iterate", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("warm_start_entire_iterate", pStrValue);
  FREE(pStrValue);

  bOk = getDoubleInPList(pStructOptions, L"warm_start_target_mu", &dblValue, 0, iLog, CHECK_NONE);
  if (bOk) options->SetNumericValue("warm_start_target_mu", dblValue);

  ///////////////////////
  // Restoration Phase //
  ///////////////////////
  // Enable heuristics to quickly detect an infeasible problem
  // * no: the problem probably be feasible
  // * yes: the problem has a good chance to be infeasible
  // This options is meant to activate heuristics that may speed up the infeasibility determination if you expect that there is a good chance for the problem to be
  // infeasible. In the filter line search procedure, the restoration phase is called more quickly than usually, and more reduction in 
  // the constraint violation is enforced before the restoration phase is left. If the problem is square, this option is enabled automatically.
  bOk = getStringInPList(pStructOptions, L"expect_infeasible_problem", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("expect_infeasible_problem", pStrValue);
  FREE(pStrValue);
  // Threshold for disabling "expect_infeasible_problem" option.
  // If the constraint violation becomes smaller than this threshold, the "expect_infeasible_problem" heuristics in the filter line 
  // search are disabled. If the problem is square, this options is set to 0.
  bOk = getDoubleInPList(pStructOptions, L"expect_infeasible_problem_ctol", &dblValue, 1e3, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("expect_infeasible_problem_ctol", dblValue);
  // Multiplier threshold for activating "expect_infeasible_problem" option.
  //   If the max norm of the constraint multipliers becomes larger than this
  //   value and "expect_infeasible_problem" is chosen, then the restoration phase is entered.
  bOk = getDoubleInPList(pStructOptions, L"expect_infeasible_problem_ytol", &dblValue, 1e8, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("expect_infeasible_problem_ytol", dblValue);
  // Tells algorithm to switch to restoration phase in first iteration.
  // * no: don't force start in restoration phase
  // * yes: force start in restoration phase
  // Setting this option to "yes" forces the algorithm to switch to the feasibility restoration phase in the first iteration. If the initial 
  // point is feasible, the algorithm will abort with a failure.
  bOk = getStringInPList(pStructOptions, L"start_with_resto", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("start_with_resto", pStrValue);
  FREE(pStrValue);
  // Required reduction in primal-dual error in the soft restoration phase.
  // The soft restoration phase attempts to reduce the primal-dual error with regular steps. If the damped 
  // primal-dual step (damped only to satisfy the fraction-to-the-boundary rule) is not decreasing the primal-dual error 
  // by at least this factor, then the regular restoration phase is called. Choosing "0" here disables the soft restoration phase.
  bOk = getDoubleInPList(pStructOptions, L"soft_resto_pderror_reduction_factor", &dblValue, 0.9999, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("soft_resto_pderror_reduction_factor", dblValue);
  // Maximum number of iterations performed successively in soft restoration phase.
  // If the soft restoration phase is performed for more than so many iterations in a row, the regular restoration phase is called.
  bOk = getIntInPList(pStructOptions, L"max_soft_resto_iters", &iValue, 10, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("max_soft_resto_iters", iValue);
  // Required infeasibility reduction in restoration phase
  // The restoration phase algorithm is performed, until a point is found that is acceptable to the filter and the infeasibility has been
  // reduced by at least the fraction given by this option.
  bOk = getDoubleInPList(pStructOptions, L"required_infeasibility_reduction", &dblValue, 0.9, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("required_infeasibility_reduction", dblValue);
  // Maximum number of successive iterations in restoration phase.
  // The algorithm terminates with an error message if the number of iterations successively taken in the restoration phase exceeds this number.
  bOk = getIntInPList(pStructOptions, L"max_resto_iter", &iValue, 3000000, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("max_resto_iter", iValue);
  // Threshold for resetting bound multipliers after the restoration phase.
  // After returning from the restoration phase, the bound multipliers are updated with a Newton step for complementarity. Here, the
  // change in the primal variables during the entire restoration phase is taken to be the corresponding primal Newton step. 
  // However, if after the update the largest bound multiplier exceeds the threshold specified by this option, the multipliers are all reset to 1.
  bOk = getDoubleInPList(pStructOptions, L"bound_mult_reset_threshold", &dblValue, 1000, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("bound_mult_reset_threshold", dblValue);
  // Threshold for resetting equality and inequality multipliers after restoration phase.
  // After returning from the restoration phase, the constraint multipliers are recomputed by a least square estimate.  This option triggers when
  // those least-square estimates should be ignored.
  bOk = getDoubleInPList(pStructOptions, L"constr_mult_reset_threshold", &dblValue, 0, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("constr_mult_reset_threshold", dblValue);
  // Determines if the original objective function should be evaluated at restoration phase trial points.
  // Setting this option to "yes" makes the restoration phase algorithm evaluate the objective function of the original problem at every trial 
  // point encountered during the restoration phase, even if this value is not required. In this way, it is guaranteed that the original 
  // objective function can be evaluated without error at all accepted iterates; otherwise the algorithm might fail at a point where the 
  // restoration phase accepts an iterate that is good for the restoration phase problem, but not the original problem. On the other hand, if 
  // the evaluation of the original objective is expensive, this might be costly.
  // * no: skip evaluation
  // * yes: evaluate at every trial point
  bOk = getStringInPList(pStructOptions, L"evaluate_orig_obj_at_resto_trial", &pStrValue, "yes", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("evaluate_orig_obj_at_resto_trial", pStrValue);
  FREE(pStrValue);
  // Penalty parameter in the restoration phase objective function.
  // This is the parameter rho in equation (31a) in the Ipopt implementation paper.
  bOk = getDoubleInPList(pStructOptions, L"resto_penalty_parameter", &dblValue, 1000.0, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("resto_penalty_parameter", dblValue);

  ///////////////////
  // Linear Solver //
  ///////////////////
  // Linear solver to be used for step calculation
  // * ma27: use the Harwell routine MA27
  // * ma57: use the Harwell routine MA57
  // * pardiso: use the Pardiso package
  // * wsmp: use WSMP package
  // * mumps: use MUMPS package
  // * custom: use custom linear solver
  // Determines which linear algebra package is to be used for the solution of the augmented linear system (for obtaining the search directions).
  // Note, the code must have been compiled with the linear solver you want to choose. Depending on your Ipopt installation, not all options are available.
  bOk = getStringInPList(pStructOptions, L"linear_solver", &pStrValue, "mumps", iLog, 
		       CHECK_VALUES, 6,
		       "ma27",
		       "ma57",
		       "pardiso",
		       "wsmp",
		       "mumps",
		       "custom");
  if (bOk) options->SetStringValue("linear_solver", pStrValue);
  FREE(pStrValue);
  // Method for scaling the linear systems
  // * none: no scaling will be performed
  // * mc19: use the Harwell routine MC19
  // "Determines the method used to compute symmetric scaling factors for the augmented system (see also the "linear_scaling_on_demand" option).
  // This scaling is independent of the NLP problem scaling.  By default, MC19 is only used if MA27 or MA57 are selected as linear solvers. 
  // This option is only available if Ipopt has been compiled with MC19.
  bOk = getStringInPList(pStructOptions, L"linear_system_scaling", &pStrValue, "none", iLog, 
		       CHECK_VALUES, 2, 
		       "none",
		       "mc19");
  if (bOk) options->SetStringValue("linear_system_scaling", pStrValue);
  FREE(pStrValue);
  // Enables heuristic for scaling only when seems required
  // * no: Always scale the linear system.
  // * yes: Start using linear system scaling if solutions seem not good.
  bOk = getStringInPList(pStructOptions, L"linear_scaling_on_demand", &pStrValue, "yes", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("linear_scaling_on_demand", pStrValue);
  FREE(pStrValue);
  // Maximal number of iterative refinement steps per linear system solve
  bOk = getIntInPList(pStructOptions, L"max_refinement_steps", &iValue, 10, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("max_refinement_steps", iValue);
  // Minimum number of iterative refinement steps per linear system solve
  bOk = getIntInPList(pStructOptions, L"min_refinement_steps", &iValue, 1, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("min_refinement_steps", iValue);
  // Iterative refinement tolerance
  //   Iterative refinement is performed until the residual test ratio is less
  //   than this tolerance (or until "max_refinement_steps" refinement steps are performed).
  bOk = getDoubleInPList(pStructOptions, L"residual_ratio_max", &dblValue, 1e-10, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("residual_ratio_max", dblValue);
  // Threshold for declaring linear system singular after failed iterative refinement.
  //   If the residual test ratio is larger than this value after failed
  //   iterative refinement, the algorithm pretends that the linear system is singular.
  bOk = getDoubleInPList(pStructOptions, L"residual_ratio_singular", &dblValue, 1e-5, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("residual_ratio_singular", dblValue);
  // Minimal required reduction of residual test ratio in iterative refinement.
  //   If the improvement of the residual test ratio made by one iterative
  //   refinement step is not better than this factor, iterative refinement is aborted.
  bOk = getDoubleInPList(pStructOptions, L"residual_improvement_factor", &dblValue, 1, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("residual_improvement_factor", dblValue);
  // Tolerance for heuristic to ignore wrong inertia.
  //   If positive, incorrect inertia in the augmented system is ignored, and we
  //   test if the direction is a direction of positive curvature.  This
  //   tolerance determines when the direction is considered to be sufficiently positive.
  bOk = getDoubleInPList(pStructOptions, L"neg_curv_test_tol", &dblValue, 0, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("neg_curv_test_tol", dblValue);

  //////////////////////////
  // Hessian Perturbation //
  //////////////////////////
  // Maximum value of regularization parameter for handling negative curvature.
  // In order to guarantee that the search directions are indeed proper descent directions, Ipopt requires that the inertia of the 
  // (augmented) linear system for the step computation has the correct number of negative and positive eigenvalues. The idea 
  // is that this guides the algorithm away from maximizers and makes Ipopt more likely converge to first order optimal points that 
  // are minimizers. If the inertia is not correct, a multiple of the identity matrix is added to the Hessian of the Lagrangian in the 
  // augmented system. This parameter gives the maximum value of the regularization parameter. If a regularization of that size is 
  // not enough, the algorithm skips this iteration and goes to the restoration phase. (This is delta_w^max in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"max_hessian_perturbation", &dblValue, 1e20, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("max_hessian_perturbation", dblValue);
  // Smallest perturbation of the Hessian block.
  // The size of the perturbation of the Hessian block is never selected smaller than this value, unless no perturbation is necessary. (This 
  // is delta_w^min in implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"min_hessian_perturbation", &dblValue, 1e-20, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("min_hessian_perturbation", dblValue);
  // Size of first x-s perturbation tried.
  // The first value tried for the x-s perturbation in the inertia correction scheme.
  // (This is delta_0 in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"first_hessian_perturbation", &dblValue, 1e-4, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("first_hessian_perturbation", dblValue);
  // Increase factor for x-s perturbation for very first perturbation.
  // The factor by which the perturbation is increased when a trial value was not sufficient - this value is used for the computation of the 
  // very first perturbation and allows a different value for for the first perturbation than that used for the remaining perturbations. 
  // (This is bar_kappa_w^+ in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"perturb_inc_fact_first", &dblValue, 100, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("perturb_inc_fact_first", dblValue);
  // Increase factor for x-s perturbation.
  // The factor by which the perturbation is increased when a trial value was not sufficient - this value is used for the computation of 
  // all perturbations except for the first. (This is kappa_w^+ in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"perturb_inc_fact", &dblValue, 8, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("perturb_inc_fact", dblValue);
  // Decrease factor for x-s perturbation.
  // The factor by which the perturbation is decreased when a trial value is deduced from the size of the most recent successful perturbation. 
  // (This is kappa_w^- in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"perturb_dec_fact", &dblValue, 0.333333, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("perturb_dec_fact", dblValue);
  // Size of the regularization for rank-deficient constraint Jacobians.
  // (This is bar delta_c in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"jacobian_regularization_value", &dblValue, 1e-8, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("jacobian_regularization_value", dblValue);
  // Exponent for mu in the regularization for rank-deficient constraint Jacobians.
  // (This is kappa_c in the implementation paper.)
  bOk = getDoubleInPList(pStructOptions, L"jacobian_regularization_exponent", &dblValue, 0.25, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("jacobian_regularization_exponent", dblValue);
  // Active permanent perturbation of constraint linearization.
  // no: perturbation only used when required
  // yes: always use perturbation
  // This options makes the delta_c and delta_d perturbation be used for the computation of every search direction.  Usually, it is only used 
  // when the iteration matrix is singular.
  bOk = getStringInPList(pStructOptions, L"perturb_always_cd", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("perturb_always_cd", pStrValue);
  FREE(pStrValue);

  //////////////////
  // Quasi-Newton //
  //////////////////
  // Can enable Quasi-Newton approximation of hessian
  // * exact: Use second derivatives provided by the NLP.
  // * limited-memory: Perform a limited-memory quasi-Newton approximation
  // This determines which kind of information for the Hessian of the Lagrangian function is used by the algorithm.
  bOk = getStringInPList(pStructOptions, L"hessian_approximation", &pStrValue, "limited-memory", iLog, 
		       CHECK_VALUES, 2,
		       "exact",
		       "limited-memory");
  if (bOk) options->SetStringValue("hessian_approximation", pStrValue);
  FREE(pStrValue);
  // Indicates in which subspace the Hessian information is to be approximated.
  // nonlinear-variables: only in space of nonlinear variables.
  // all-variables: in space of all variables (without slacks)
  bOk = getStringInPList(pStructOptions, L"hessian_approximation_space", &pStrValue, "nonlinear-variables", iLog, 
		       CHECK_VALUES, 2,
		       "nonlinear-variables",
		       "all-variables");
  if (bOk) options->SetStringValue("hessian_approximation_space", pStrValue);
  FREE(pStrValue);
  // Maximum size of the history for the limited quasi-Newton Hessian approximation.
  // This option determines the number of most recent iterations that are taken into account for the limited-memory quasi-Newton approximation.
  bOk = getIntInPList(pStructOptions, L"limited_memory_max_history", &iValue, 6, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("limited_memory_max_history", iValue);
  // Threshold for successive iterations where update is skipped.
  // If the update is skipped more than this number of successive iterations, we quasi-Newton approximation is reset.
  bOk = getIntInPList(pStructOptions, L"limited_memory_max_skipping", &iValue, 2, iLog, CHECK_MIN, 1);
  if (bOk) options->SetIntegerValue("limited_memory_max_skipping", iValue);
  // Quasi-Newton update formula for the limited memory approximation.
  // bfgs: BFGS update (with skipping)
  // sr1:  SR1 (not working well)
  // Determines which update formula is to be used for the limited-memory quasi-Newton approximation.
  bOk = getStringInPList(pStructOptions, L"limited_memory_update_type", &pStrValue, "bfgs", iLog, 
		       CHECK_VALUES, 2,
		       "bfgs",
		       "sr1");
  if (bOk) options->SetStringValue("limited_memory_update_type", pStrValue);
  FREE(pStrValue);
  // Initialization strategy for the limited memory quasi-Newton approximation.
  // scalar1: sigma = s^Ty/s^Ts
  // scalar2: sigma = y^Ty/s^Ty
  // constant: sigma = limited_memory_init_val
  // Determines how the diagonal Matrix B_0 as the first term in the limited memory approximation should be computed.
  bOk = getStringInPList(pStructOptions, L"limited_memory_initialization", &pStrValue, "scalar1", iLog, 
		       CHECK_VALUES, 3, 
		       "scalar1",
		       "scalar2",
		       "constant");
  if (bOk) options->SetStringValue("limited_memory_initialization", pStrValue);
  FREE(pStrValue);
  // Upper bound on value for B0 in low-rank update.
  // The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have 
  // been performed yet), and is constantly chosen as this value, if "limited_memory_initialization" is "constant".
  bOk = getDoubleInPList(pStructOptions, L"limited_memory_init_val_max", &dblValue, 1e8, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("limited_memory_init_val_max", dblValue);
  // Lower bound on value for B0 in low-rank update.
  // The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have 
  // been performed yet), and is constantly chosen as this value, if "limited_memory_initialization" is "constant".
  bOk = getDoubleInPList(pStructOptions, L"limited_memory_init_val_min", &dblValue, 1e-8, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("limited_memory_init_val_min", dblValue);
  // Value for B0 in low-rank update
  // The starting matrix in the low rank update, B0, is chosen to be this multiple of the identity in the first iteration (when no updates have
  // "been performed yet), and is constantly chosen as this value, if "limited_memory_initialization" is "constant".
  bOk = getDoubleInPList(pStructOptions, L"limited_memory_init_val", &dblValue, 1, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("limited_memory_init_val", dblValue);

  /////////////////////
  // Derivative Test //
  /////////////////////
  // Enable derivative checker
  // * none: do not perform derivative test
  // * first-order: perform test of first derivatives at starting point
  // * second-order: perform test of first and second derivatives at starting point
  // If this option is enabled, a (slow) derivative test will be performed before the optimization. The test is performed at the user provided 
  // starting point and marks derivative values that seem suspicious
  bOk = getStringInPList(pStructOptions, L"derivative_test", &pStrValue, "none", iLog, 
		       CHECK_VALUES, 4,
		       "none",
		       "first-order",
		       "second-order",
		       "only-second-order");
  if (bOk) options->SetStringValue("derivative_test", pStrValue);
  FREE(pStrValue);
  // Index of first quantity to be checked by derivative checker
  //   If this is set to -2, then all derivatives are checked.  Otherwise, for
  //   the first derivative test it specifies the first variable for which the
  //   test is done (counting starts at 0.).  For second derivatives, it
  //   specifies the first constraint for which the test is done; counting of
  //   constraint indices starts at 0, and -1 refers to the objective function
  //   Hessian.
  bOk = getIntInPList(pStructOptions, L"derivative_test_first_index", &iValue, -2, iLog, CHECK_MIN, -2);
  if (bOk) options->SetIntegerValue("derivative_test_first_index", iValue);
  // Size of the finite difference perturbation in derivative test.
  // This determines the relative perturbation of the variable entries.
  bOk = getDoubleInPList(pStructOptions, L"derivative_test_perturbation", &dblValue, 1e-8, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("derivative_test_perturbation", dblValue);
  // Threshold for indicating wrong derivative.
  // If the relative deviation of the estimated derivative from the given one is larger than this value, the corresponding derivative is marked as wrong.
  bOk = getDoubleInPList(pStructOptions, L"derivative_test_tol", &dblValue, 1e-4, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("derivative_test_tol", dblValue);
  // Indicates whether information for all estimated derivatives should be printed.
  // * no: Print only suspect derivatives
  // * yes: Print all derivatives
  bOk = getStringInPList(pStructOptions, L"derivative_test_print_all", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("derivative_test_print_all", pStrValue);
  FREE(pStrValue);
  // Maximal perturbation of an evaluation point.
  // If a random perturbation of a points is required, this number indicates the maximal perturbation. This is for example used when 
  // determining the center point at which the finite difference derivative test is executed.
  bOk = getDoubleInPList(pStructOptions, L"point_perturbation_radius", &dblValue, 10, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("point_perturbation_radius", dblValue);
  // Specifies technique to compute constraint Jacobian
  // exact: user-provided derivatives
  // finite-difference-values: user-provided structure, values by finite differences
  bOk = getStringInPList(pStructOptions, L"jacobian_approximation", &pStrValue, "exact", iLog, 
		       CHECK_VALUES, 2, 
		       "exact",
		       "finite-difference-values");
  if (bOk) options->SetStringValue("jacobian_approximation", pStrValue);
  FREE(pStrValue);
  // Size of the finite difference perturbation for derivative approximation.
  // This determines the relative perturbation of the variable entries.
  bOk = getDoubleInPList(pStructOptions, L"findiff_perturbation", &dblValue, 1e-7, iLog, CHECK_MIN, 0.);
  if (bOk) options->SetNumericValue("findiff_perturbation", dblValue);

#ifdef USE_MA27
  ////////////////////////
  // MA27 Linear Solver //
  ////////////////////////
  // Pivot tolerance for the linear solver MA27
  // A smaller number pivots for sparsity, a larger number pivots for stability.  This option is only available if Ipopt has been compiled with MA27.
  bOk = getDoubleInPList(pStructOptions, L"ma27_pivtol", &dblValue, 1e-8, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("ma27_pivtol", dblValue);
  // Maximal pivot tolerance for the linear solver MA27
  // Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system.  This option is only available if 
  // Ipopt has been compiled with MA27.
  bOk = getDoubleInPList(pStructOptions, L"ma27_pivtolmax", &dblValue, 1e-4, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("ma27_pivtolmax", dblValue);
  // Integer workspace memory for MA27.
  // The initial integer workspace memory = liw_init_factor * memory required by unfactored system. Ipopt will increase the workspace 
  // size by meminc_factor if required.  This option is only available if Ipopt has been compiled with MA27.
  bOk = getDoubleInPList(pStructOptions, L"ma27_liw_init_factor", &dblValue, 5, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("ma27_liw_init_factor", dblValue);
  // Real workspace memory for MA27.
  // The initial real workspace memory = la_init_factor * memory required by unfactored system. Ipopt will increase the workspace
  // size by meminc_factor if required.  This option is only available if Ipopt has been compiled with MA27.
  bOk = getDoubleInPList(pStructOptions, L"ma27_la_init_factor", &dblValue, 5, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("ma27_la_init_factor", dblValue);
  // Increment factor for workspace size for MA27.
  // If the integer or real workspace is not large enough, Ipopt will increase its size by this factor.  This option is only 
  // available if Ipopt has been compiled with MA27.
  bOk = getDoubleInPList(pStructOptions, L"ma27_meminc_factor", &dblValue, 5, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("ma27_meminc_factor", dblValue);
  // Always pretend inertia is correct.
  // no: check inertia
  // yes: skip inertia check
  // Setting this option to "yes" essentially disables inertia check. This option makes the algorithm non-robust and easily fail, but it 
  // might give some insight into the necessity of inertia control.
  bOk = getStringInPList(pStructOptions, L"ma27_skip_inertia_check", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("ma27_skip_inertia_check", pStrValue);
  FREE(pStrValue);
  // Enables MA27's ability to solve a linear system even if the matrix is singular.
  // no: Don't have MA27 solve singular systems
  // yes: Have MA27 solve singular systems
  // Setting this option to \"yes\" means that Ipopt will call MA27 to compute solutions for right hand sides, even if MA27 has detected that 
  // the matrix is singular (but is still able to solve the linear system). In some cases this might be better than using Ipopt's heuristic of 
  // small perturbation of the lower diagonal of the KKT matrix.
  bOk = getStringInPList(pStructOptions, L"ma27_ignore_singularity", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("ma27_ignore_singularity", pStrValue);
  FREE(pStrValue);
#endif

#ifdef USE_MA28
  ////////////////////////
  // MA28 Linear Solver //
  ////////////////////////
  // Pivot tolerance for linear solver MA28.
  //   This is used when MA28 tries to find the dependent constraints.
  bOk = getDoubleInPList(pStructOptions, L"ma28_pivtol", &dblValue, 0.01, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("ma27_pivtol", dblValue);
#endif

#ifdef USE_MA57
  ////////////////////////
  // MA57 Linear Solver //
  ////////////////////////
  // Pivot tolerance for the linear solver MA57
  // A smaller number pivots for sparsity, a larger number pivots for stability. This option is only available if Ipopt has been compiled with MA57.
  bOk = getDoubleInPList(pStructOptions, L"ma57_pivtol", &dblValue, 1e-8, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("ma57_pivtol", dblValue);
  // Maximal pivot tolerance for the linear solver MA57
  // Ipopt may increase pivtol as high as ma57_pivtolmax to get a more accurate solution to the linear system.  This option is only available 
  // if Ipopt has been compiled with MA57.
  bOk = getDoubleInPList(pStructOptions, L"ma57_pivtolmax", &dblValue, 1e-4, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("ma57_pivtolmax", dblValue);
  // Safety factor for work space memory allocation for the linear solver MA57.
  // If 1 is chosen, the suggested amount of work space is used.  However, choosing a larger number might avoid reallocation if the suggest values 
  // do not suffice.  This option is only available if Ipopt has been compiled with MA57.
  bOk = getIntInPList(pStructOptions, L"ma57_pre_alloc", &iValue, 3, iLog, CHECK_MIN, 1);
  if (bOk) options->SetIntegerValue("ma57_pre_alloc", iValue);
  // Controls pivot order in MA57
  //   This is INCTL(6) in MA57.
  bOk = getIntInPList(pStructOptions, L"ma57_pivot_order", &iValue, 5, iLog, CHECK_BOTH, 0, 5);
  if (bOk) options->SetIntegerValue("ma57_pivot_order", iValue);
#endif

#ifdef USE_MUMPS
  /////////////////////////
  // MUMPS Linear Solver //
  /////////////////////////
  // Pivot tolerance for the linear solver MUMPS.
  // A smaller number pivots for sparsity, a larger number pivots for stability.  This option is only available if Ipopt has been compiled with MUMPS.
  bOk = getDoubleInPList(pStructOptions, L"mumps_pivtol", &dblValue, 1e-6, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("mumps_pivtol", dblValue);
  // Maximum pivot tolerance for the linear solver MUMPS.
  // Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system.  This option is only available if 
  // Ipopt has been compiled with MUMPS.
  bOk = getDoubleInPList(pStructOptions, L"mumps_pivtolmax", &dblValue, 0.1, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("mumps_pivtolmax", dblValue);
  // Percentage increase in the estimated working space for MUMPS.
  // In MUMPS when significant extra fill-in is caused by numerical pivoting, larger values of mumps_mem_percent may help use the workspace more efficiently.
  bOk = getIntInPList(pStructOptions, L"mumps_mem_percent", &iValue, 1000, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("mumps_mem_percent", iValue);
  // Controls permuting and scaling in MUMPS
  // This is ICTL(6) in MUMPS.
  bOk = getIntInPList(pStructOptions, L"mumps_permuting_scaling", &iValue, 7, iLog, CHECK_BOTH, 0, 7);
  if (bOk) options->SetIntegerValue("mumps_permuting_scaling", iValue);
  // Controls pivot order in MUMPS
  // This is ICTL(8) in MUMPS.
  bOk = getIntInPList(pStructOptions, L"mumps_pivot_order", &iValue, 7, iLog, CHECK_BOTH, 0, 7);
  if (bOk) options->SetIntegerValue("mumps_pivot_order", iValue);
  // Controls scaling in MUMPS
  // This is ICTL(8) in MUMPS.
  bOk = getIntInPList(pStructOptions, L"mumps_scaling", &iValue, 77, iLog, CHECK_BOTH, -2, 77);
  if (bOk) options->SetIntegerValue("mumps_scaling", iValue);
  // Pivot threshold for detection of linearly dependent constraints in MUMPS.
  // When MUMPS is used to determine linearly dependent constraints, this is determines the threshold for a pivot to be considered zero. This is CNTL(3) in MUMPS.
  bOk = getDoubleInPList(pStructOptions, L"mumps_dep_tol", &dblValue, -1.0, iLog, CHECK_NONE);
  if (bOk) options->SetNumericValue("mumps_dep_tol", dblValue);
#endif

#ifdef USE_PARDISO
  ///////////////////////////
  // Pardiso Linear Solver //
  ///////////////////////////
  // Matching strategy for linear solver Pardiso
  // * complete: Match complete (IPAR(13)=1)
  // * complete+2x2: Match complete+2x2 (IPAR(13)=2)
  // * constraints: Match constraints (IPAR(13)=3)
  // This is IPAR(13) in Pardiso manual.  This option is only available if Ipopt has been compiled with Pardiso.
  bOk = getStringInPList(pStructOptions, L"pardiso_matching_strategy", &pStrValue, "complete+2x2", iLog, 
		       CHECK_VALUES, 3,
		       "complete",
		       "complete+2x2",
		       "constraints");
  if (bOk) options->SetStringValue("pardiso_matching_strategy", pStrValue);
  FREE(pStrValue);
  // Enables out-of-core version of linear solver Pardiso
  // Setting this option to a positive integer k makes Pardiso work in the out-of-core variant where the factor is split in 2^k subdomains. This 
  // is IPARM(50.) in the Pardiso manual.  This option is only available if Ipopt has been compiled with Pardiso.
  bOk = getIntInPList(pStructOptions, L"pardiso_out_of_core_power", &iValue, 0, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("pardiso_out_of_core_power", iValue);
  // Pardiso message level
  //   This determines the amount of analysis output from the Pardiso solver.
  //   This is MSGLVL in the Pardiso manual.
  bOk = getIntInPList(pStructOptions, L"pardiso_msglvl", &iValue, 0, iLog, CHECK_MIN, 0);
  if (bOk) options->SetIntegerValue("pardiso_msglvl", iValue);
  // Maximum number of Krylov-Subspace Iteration
  //   DPARM(1)
  bOk = getIntInPList(pStructOptions, L"pardiso_max_iter", &iValue, 500, iLog, CHECK_MIN, 1);
  if (bOk) options->SetIntegerValue("pardiso_max_iter", iValue);
  // Relative Residual Convergence
  //   DPARM(2)
  bOk = getDoubleInPList(pStructOptions, L"pardiso_iter_relative_tol", &dblValue, 1e-6, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("pardiso_iter_relative_tol", dblValue);
  // Maximum Size of Coarse Grid Matrix
  //   DPARM(3)
  bOk = getDoubleInPList(pStructOptions, L"pardiso_iter_coarse_size", &dblValue, 5000, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("pardiso_iter_coarse_size", dblValue);
  // Maximum Size of Grid Levels
  //   DPARM(4)
  bOk = getDoubleInPList(pStructOptions, L"pardiso_iter_max_levels", &dblValue, 10000, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("pardiso_iter_max_levels", dblValue);
  // dropping value for incomplete factor
  //   DPARM(5)
  bOk = getDoubleInPList(pStructOptions, L"pardiso_iter_dropping_factor", &dblValue, 0.5, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("pardiso_iter_dropping_factor", dblValue);
  // dropping value for sparsify schur complement factor
  //   DPARM(6)
  bOk = getDoubleInPList(pStructOptions, L"pardiso_iter_dropping_schur", &dblValue, 0.1, iLog, CHECK_BOTH, 0.,1.);
  if (bOk) options->SetNumericValue("pardiso_iter_dropping_schur", dblValue);
  // max fill for each row
  //   DPARM(7)
  bOk = getDoubleInPList(pStructOptions, L"pardiso_iter_max_row_fill", &dblValue, 10000000, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("pardiso_iter_max_row_fill", dblValue);
  //   DPARM(8)
  bOk = getDoubleInPList(pStructOptions, L"pardiso_iter_inverse_norm_factor", &dblValue, 5e6, iLog, CHECK_MIN, 1.);
  if (bOk) options->SetNumericValue("pardiso_iter_inverse_norm_factor", dblValue);
  // Switch on iterative solver in Pardiso library
  // Possible values:
  //  - no                      []
  //  - yes                     []
  bOk = getStringInPList(pStructOptions, L"pardiso_iterative", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("pardiso_iterative", pStrValue);
  FREE(pStrValue);
  // Maximal number of decreases of drop tolerance during one solve.
  //   This is relevant only for iterative Pardiso options.
  bOk = getIntInPList(pStructOptions, L"pardiso_max_droptol_corrections", &iValue, 4, iLog, CHECK_MIN, 1);
  if (bOk) options->SetIntegerValue("pardiso_max_droptol_corrections", iValue);

  // Toggle for handling case when elements were perturbed by Pardiso.
  // no: Always redo symbolic factorization when elements were perturbed
  // yes: Only redo symbolic factorization when elements were perturbed if also the inertia was wrong
  // This option is only available if Ipopt has been compiled with Pardiso.
  bOk = getStringInPList(pStructOptions, L"pardiso_redo_symbolic_fact_only_if_inertia_wrong", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("pardiso_redo_symbolic_fact_only_if_inertia_wrong", pStrValue);
  FREE(pStrValue);
  // Interpretation of perturbed elements
  // no: Don't assume that matrix is singular if elements were perturbed after recent symbolic factorization
  // yes: Assume that matrix is singular if elements were perturbed after recent symbolic factorization
  // This option is only available if Ipopt has been compiled with Pardiso.
  bOk = getStringInPList(pStructOptions, L"pardiso_repeated_perturbation_means_singular", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("pardiso_repeated_perturbation_means_singular", pStrValue);
  FREE(pStrValue);
  // Always pretent inertia is correct.
  // no: check inertia
  // yes: skip inertia check
  // Setting this option to "yes" essentially disables inertia check. This option makes the algorithm non-robust and easily fail, but it 
  // might give some insight into the necessity of inertia control.
  bOk = getStringInPList(pStructOptions, L"pardiso_skip_inertia_check", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");
  if (bOk) options->SetStringValue("pardiso_skip_inertia_check", pStrValue);
  FREE(pStrValue);
#endif

#ifdef USE_WSMP
  ////////////////////////
  // WSMP Linear Solver //
  ////////////////////////
  // Number of threads to be used in WSMP
  // This determines on how many processors WSMP is running on. This option is only available if Ipopt has been compiled with WSMP.
  bOk = getIntInPList(pStructOptions, L"wsmp_num_threads", &iValue, 1, iLog, CHECK_NONE);
  if (bOk) options->SetIntegerValue("wsmp_num_threads", iValue);
  // Determines how ordering is done in WSMP
  // This corresponds to the value of WSSMP's IPARM(16). This option is only available if Ipopt has been compiled with WSMP.
  bOk = getIntInPList(pStructOptions, L"wsmp_ordering_option", &iValue, 1, iLog, CHECK_NONE);
  if (bOk) options->SetIntegerValue("wsmp_ordering_option", iValue);
  // Pivot tolerance for the linear solver WSMP.
  // A smaller number pivots for sparsity, a larger number pivots for stability.  This option is only available if Ipopt has been compiled with WSMP.
  bOk = getDoubleInPList(pStructOptions, L"wsmp_pivtol", &iValue, 1e-4, iLog, CHECK_NONE);
  if (bOk) options->SetNumericValue("wsmp_pivtol", dblValue);
  // Maximum pivot tolerance for the linear solver WSMP.
  // Ipopt may increase pivtol as high as pivtolmax to get a more accurate solution to the linear system. This option is only available if Ipopt 
  // has been compiled with WSMP.
  bOk = getDoubleInPList(pStructOptions, L"wsmp_pivtolmax", &iValue, 0.1, iLog, CHECK_NONE);
  if (bOk) options->SetNumericValue("wsmp_pivtolmax", dblValue);
  // Determines how the matrix is scaled by WSMP.
  // This corresponds to the value of WSSMP's IPARM(10.). This option is only available if Ipopt has been compiled with WSMP.
  bOk = getIntInPList(pStructOptions, L"wsmp_scaling", &iValue, 0, iLog, CHECK_NONE);
  if (bOk) options->SetIntegerValue("wsmp_scaling", iValue);
#endif

  // wantsol: solution report without -AMPL: sum of
  // 1 ==> write .sol file
  // 2 ==> print primal variable values
  // 4 ==> print dual variable values
  // 8 ==> do not print solution message
  bOk = getIntInPList(pStructOptions, L"wantsol", &iValue, 8, iLog, CHECK_NONE);
  if (bOk) options->SetIntegerValue("wantsol", iValue);

  // print_user_options: If selected, the algorithm will print the list of all options set by
  //                     the user including their values and whether they have been used.  In
  //                     some cases this information might be incorrect, due to the internal program flow.
  // - no, don't print options (default value)
  // - yes, print options
  bOk = getStringInPList(pStructOptions, L"print_user_options", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");

  if (bOk) options->SetStringValue("print_user_options", pStrValue);
  FREE(pStrValue);

  // print_options_documentation: If selected, the algorithm will print the list of all available
  //                              algorithmic options with some documentation before solving the
  //                              optimization problem.
  // - no, don't print list (default value)
  // - yes, print list

  bOk = getStringInPList(pStructOptions, L"print_options_documentation", &pStrValue, "no", iLog, 
		       CHECK_VALUES, 2, "no", "yes");

  if (bOk) options->SetStringValue("print_options_documentation", pStrValue);
  FREE(pStrValue);

  return 0;
}
