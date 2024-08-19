import reaktoro as rkt

import numpy as np

from reaktoro_pse.core.reaktoro_state import reaktoroState
from reaktoro_pse.core.reaktoro_outputs import (
    reaktoroOutputSpec,
)
from reaktoro_pse.core.reaktoro_inputs import (
    reaktoroInputSpec,
)
from reaktoro_pse.core.reaktoro_jacobian import (
    reaktoroJacobianSpec,
)
import cyipopt
import idaes.logger as idaeslog
import time

__author__ = "Alexander Dudchenko, Ben Knueven, Ilayda Akkor"

_log = idaeslog.getLogger(__name__)
"""class to setup reaktor solver for reaktoro"""


class reaktoroSolver:
    def __init__(
        self,
        reaktoroBase,
        rktInputSpec,
        rktOutputSpec,
        rktJacobianSpec,
        block_name=None,
    ):
        self.blockName = block_name
        self.rktBase = reaktoroBase
        if isinstance(self.rktBase, reaktoroState) == False:
            raise TypeError("Reator jacobian require rektoroState class")
        self.rktInputSpec = rktInputSpec
        if isinstance(self.rktInputSpec, reaktoroInputSpec) == False:
            raise TypeError("Reator outputs require reaktoroOutputSpec class")
        self.rktOutputSpec = rktOutputSpec
        if isinstance(self.rktOutputSpec, reaktoroOutputSpec) == False:
            raise TypeError("Reator outputs require reaktoroOutputSpec class")
        self.rktJacobianSpec = rktJacobianSpec
        if isinstance(self.rktJacobianSpec, reaktoroJacobianSpec) == False:
            raise TypeError("Reator outputs require reaktoroOutputSpec class")
        self.rktSolverOptions = rkt.EquilibriumOptions()
        self.rktPresolveOptions = rkt.EquilibriumOptions()

        existing_constraints = self.rktInputSpec.rktEquilibriumSpecs.namesConstraints()

        existing_variables = (
            self.rktInputSpec.rktEquilibriumSpecs.namesControlVariables()
        )
        _log.info(f"rktSolver inputs: {existing_variables}")
        _log.info(f"rktSolver constraints: {existing_constraints}")
        self.rktSolver = rkt.EquilibriumSolver(self.rktInputSpec.rktEquilibriumSpecs)
        self.rktConditions = rkt.EquilibriumConditions(
            self.rktInputSpec.rktEquilibriumSpecs
        )

        self.rktSensitivity = rkt.EquilibriumSensitivity(
            self.rktInputSpec.rktEquilibriumSpecs
        )
        self.set_solver_options()
        self._sequential_fails = 0
        self._max_fails = 30

    def set_solver_options(
        self,
        tolerance=1e-32,
        presolve=False,
        presolve_tolerance=1e-12,
        max_iters=500,
        presolve_max_iters=500,
        hessian_type="J.tJ",
    ):
        """configuration for reaktro solver

        Keyword arguments:
        tolerance -- reaktoro solver tolerance (default: 1e-32)
        max_iters -- maxium iterations for reaktoro solver (default: 200)
        presolve -- presolve the model useing "presolve_tolerance" first (default: False)
        presolve_tolerance -- presolve tolerance if enabled (default: 1e-12)
        presolve_max_iters -- maximum iterations for presolve call if enabled (default: 200)
        """
        self.rktSolverOptions.epsilon = tolerance
        self.rktSolverOptions.optima.maxiters = max_iters
        self.rktPresolveOptions.epsilon = presolve_tolerance
        self.rktPresolveOptions.optima.maxiters = presolve_max_iters
        self.presolve = presolve
        self.rktSolver.setOptions(self.rktSolverOptions)
        self.hessian_type = hessian_type
        if self.rktInputSpec.assertChargeNeutrality:
            self.rktConditions.charge(0)

    def update_specs(self, params):
        for input_key in self.rktInputSpec.rktInputs.rktInputList:
            input_obj = self.rktInputSpec.rktInputs[input_key]

            if params is None:
                value = input_obj.get_value(update_temp=True)
            else:
                value = params.get(input_key)
                input_obj.currentValue = value
            unit = input_obj.mainUnit
            if input_key == "temperature":
                self.rktConditions.temperature(value, unit)
            elif input_key == "pressure":
                self.rktConditions.pressure(value, unit)
            else:
                # TODO figure out how deal with units...
                self.rktConditions.set(input_obj.get_rkt_input_name(), value)

    def get_outputs(self):
        output_arr = []
        for key, obj in self.rktOutputSpec.rktOutputs.items():
            output_arr.append(
                self.rktOutputSpec.evaluate_property(obj, update_values_in_object=True)
            )
        return output_arr

    def get_jacobian(self):
        self.tempJacobianMatrix = self.rktSensitivity.dudw()
        self.rktJacobianSpec.update_jacobian_absolute_values()
        jac_matrix = []
        for input_key in self.rktInputSpec.rktInputs.rktInputList:
            input_obj = self.rktInputSpec.rktInputs[input_key]
            jac_row = self.rktJacobianSpec.get_jacobian(
                self.tempJacobianMatrix, input_obj
            )
            jac_matrix.append(jac_row)

        return np.array(jac_matrix).T  # needs to transpose

    def solve_reaktoro_block(
        self,
        params=None,
        display=False,
        presolve=False,
    ):
        """here we solve reaktor model and return the jacobian matrix and solution, as
        well as update relevant reaktoroSpecs"""
        ts = time.time()
        self.update_specs(params)

        result = self.try_solve(presolve)
        self.outputs = self.get_outputs()
        self.jacobianMatrix = self.get_jacobian()
        if result.succeeded() == False or display:
            print(
                f"warning, solve was not successful for {self.blockName}, fail# {self._sequential_fails}"
            )
            self._sequential_fails += 1
            if self._sequential_fails > self._max_fails:
                assert False
            raise cyipopt.CyIpoptEvaluationError
        else:
            self._sequential_fails = 0
        return self.jacobianMatrix, self.outputs

    def try_solve(self, presolve=False):
        if self.presolve or presolve:
            """solve to loose tolerance first if selected"""
            self.rktSolver.setOptions(self.rktPresolveOptions)
            self.rktSolver.solve(
                self.rktBase.rktState,
                self.rktSensitivity,
                self.rktConditions,
            )
            self.rktSolver.setOptions(self.rktSolverOptions)

        result = self.rktSolver.solve(
            self.rktBase.rktState,
            self.rktSensitivity,
            self.rktConditions,
        )
        self.rktOutputSpec.update_supported_props()
        return result
