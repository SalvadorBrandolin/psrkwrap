from typing import Callable

import numpy as np

from scipy.optimize import minimize

from .psrk import PSRK


def mathias_copeman_optim(
    substance: str,
    sat_pressure_func: Callable,
    temperature_range: np.ndarray,
    temperature_intervals_num: int,
    c_guess: np.ndarray,
    method: str = "Nelder-Mead",
    tol=1e-4,
    max_func_eval=1000,
):
    """Optimize mathias copeman parameters with scipy.optimize.minimize.

    Parameters
    ----------
    substance : str
        Substance name in database.
    sat_pressure_func : Callable
        Function of temperature that returns saturation pressure in Pa at that
        temperature in K.
    temperature_range : np.ndarray
        Temperature range to perform the optimization.
    temperature_intervals_num : int
        Temperature interval number.
    c_guess : np.ndarray
        Mathias Copeman coefficients guess for optimization.
        np.array([c1, c2 ,c3])
    method : str, optional
        scipy minimize method, by default "Nelder-Mead"

    Returns
    -------
    OptimizeResult
        Output of scipy minimize.
    """

    temperatures = np.linspace(
        temperature_range[0], temperature_range[1], temperature_intervals_num
    )

    # =========================================================================
    # Function to optimize
    # =========================================================================
    def optim_func(variables):
        c1, c2, c3 = variables

        # =====================================================================
        # Build the "data base" por psrk with the adjustable variables
        # =====================================================================
        csv_data_mc = f"""Clapeyron Database File,,,,,
                        MathiasCopemanAlpha like Parameters,,,,,
                        species,c1,c2,c3
                        {substance},{c1},{c2},{c3}
        """

        # =====================================================================
        # Create PSRK object
        # =====================================================================
        model = PSRK(
            substances=[substance],
            pure_temperature=None,
            pure_pressure=None,
            alpha_userlocations=[csv_data_mc],
        )

        # =====================================================================
        # Eval the sum squared error
        # =====================================================================
        sse = 0
        for t in temperatures:
            experimental_sat_pressure = sat_pressure_func(t)
            psrk_sat_pressure = model.saturation_pressure(
                t, experimental_sat_pressure
            )

            abs_error = np.abs(
                experimental_sat_pressure - psrk_sat_pressure[0]
            )

            sse = sse + abs_error

        return sse

    # =========================================================================
    # Optimization
    # =========================================================================
    solution = minimize(
        fun=optim_func,
        x0=c_guess,
        method=method,
        tol=tol,
        options={"maxfev": max_func_eval},
    )

    return solution
