import os

from julia import Julia

import numpy as np

# Directorio de donde se instalo psrk.py
dirname = os.path.dirname(__file__)

# Ubicacion del compilado de julia sys.so
if os.path.exists(f"{dirname}/sys.so"):
    jl = Julia(sysimage=f"{dirname}/sys.so")

    from julia import Main

    # Uso a Clapeyron
    Main.using("Clapeyron")

    # Incluir Mathias Copeman
    Main.include(f"{dirname}/julia/mathiascopeman.jl")

    # Incluir PSRKFULL
    Main.include(f"{dirname}/julia/psrkfull.jl")


class PSRK:
    """PSRK model wrap of the PSRKFULL inspired in PSRK Clapeyron.

    Uses Julia with PyJulia and PyCall.

    Parameters
    ----------
    substances : list[str]
        Substance list.
    pure_temperature : float, optional
        Temperature to eval pure fugacity coefficients, by default 303.15 K.
    pure_pressure : int, optional
        Pressure to eval pure fugacity coefficients, by default 101325 Pa.
    phase : str, optional
        Phase to eval properties (liquid or vapour), by default "liquid".
    flash_nature : str, optional
        Option to define the kind of flash you want vle or lle. by default lle.
    userlocations : list[str], optional
        Location of data base: ["path/to/database/"]. by default library
        database.
    group_userlocations : list[str], optional
        Location of psrk groups data base: ["path/to/database/"]. by default
        library database.
    ideal_userlocations : list[str], optional
        Location of ReidIdeal data base: ["path/to/database"]. by default
        library database.
    alpha_userlocations : list[str], optional
        Location of MathiasCopeman data base: ["path/to/database"]. by default
        library database.
    mixing_userlocations : list[str], optional
        Location of data base: ["path/to/database"]. by default library
        database.
    translation_userlocations : list[str], optional
        Location of ConstantTranslation data base: ["path/to/database"].
        by default library database.

    Attributes
    ----------
    substances : list[str]
        Substance list.
    pure_temperature : float
        Temperature to eval pure fugacity coefficients, by default 303.15 K.
    pure_pressure : int
        Pressure to eval pure fugacity coefficients, by default 101325 Pa.
    phase : str
        Phase to eval properties (liquid or vapour), by default "liquid".
    flash_nature : str
        Option to define the kind of flash you want vle or lle. by default lle.
    userlocations : list[str]
        Location of data base: ["path/to/database/"]. by default library
        database.
    group_userlocations : list[str]
        Location of psrk groups data base: ["path/to/database/"]. by default
        library database.
    ideal_userlocations : list[str]
        Location of ReidIdeal data base: ["path/to/database"]. by default
        library database.
    alpha_userlocations : list[str]
        Location of MathiasCopeman data base: ["path/to/database"]. by default
        library database.
    mixing_userlocations : list[str]
        Location of data base: ["path/to/database"]. by default library
        database.
    translation_userlocations : list[str]
        Location of ConstantTranslation data base: ["path/to/database"].
        by default library database.
    critical_temperatures : np.ndarray
        Substances critical temperatures obtained from database.
    critical_pressures : np.ndarray
        Substances critical pressures obtained from database.
    c1 : np.ndarray
        Substances Mathias Copeman c1 coefficient.
    c2 : np.ndarray
        Substances Mathias Copeman c2 coefficient.
    c3 : np.ndarray
        Substances Mathias Copeman c3 coefficient.
    v_shift : np.ndarray
        Volume constant translation.
    functional_groups : np.ndarray
        UNIFAC functional groups.
    functional_groups_multiplicity : np.ndarray
        Multiplicity of each funcitonal group of each substance.
    r : np.ndarray
        Normalized group Van der Vals volume.
    q : np.ndarray
        Normalized group Surface Area.
    a_pair : np.ndarray
        UNIFAC binary A group Interaction Energy Parameter.
    b_pair : np.ndarray
        UNIFAC binary B group Interaction Energy Parameter.
    c_pair : np.ndarray
        UNIFAC binary C group Interaction Energy Parameter.
    """

    def __init__(
        self,
        substances: list[str],
        pure_temperature=303.15,
        pure_pressure=101325,
        phase="liquid",
        flash_nature="lle",
        userlocations=[f"{dirname}/julia/database/"],
        group_userlocations=[f"{dirname}/julia/database/"],
        ideal_userlocations=[f"{dirname}/julia/database/"],
        alpha_userlocations=[f"{dirname}/julia/database/"],
        mixing_userlocations=[f"{dirname}/julia/database/"],
        translation_userlocations=[f"{dirname}/julia/database/"],
    ) -> None:
        # =====================================================================
        # Constants
        # =====================================================================
        self.substances = substances
        self.pure_temperature = pure_temperature
        self.pure_pressure = pure_pressure
        self.userlocations = userlocations
        self.group_userlocations = group_userlocations
        self.ideal_userlocations = ideal_userlocations
        self.alpha_userlocations = alpha_userlocations
        self.mixing_userlocations = mixing_userlocations
        self.translation_userlocations = translation_userlocations

        if flash_nature == "lle":
            self.flash_nature = "lle"
        elif flash_nature == "vle":
            self.flash_nature = "vle"
        else:
            raise ValueError("Only lle and vle flash nature allowed")

        if phase == "liquid":
            self.phase = "liquid"
        elif phase == "vapour":
            self.phase = "vapour"
        else:
            raise ValueError("Wrong phase nature, use liquid or vapour")

        # =====================================================================
        # Making the PSRK clapeyron model (with personalized PSRK and alpha)
        # =====================================================================
        self.psrk = Main.PSRKFULL(
            self.substances,
            idealmodel=Main.ReidIdeal,
            userlocations=self.userlocations,
            group_userlocations=self.group_userlocations,
            ideal_userlocations=self.ideal_userlocations,
            alpha_userlocations=self.alpha_userlocations,
            mixing_userlocations=self.mixing_userlocations,
            translation_userlocations=self.translation_userlocations,
        )

        # =====================================================================
        # Model constants from the data base
        # =====================================================================
        self.molecular_weights = self.psrk.params.Mw.values

        # SRK parameters
        self.critical_temperatures = self.psrk.params.Tc.values
        self.critical_pressures = self.psrk.params.Pc.values
        self.c1 = self.psrk.alpha.params.c1.values
        self.c2 = self.psrk.alpha.params.c2.values
        self.c3 = self.psrk.alpha.params.c3.values
        self.v_shift = self.psrk.translation.params.v_shift.values

        # UNIFAC parameters
        self.functional_groups = self.psrk.mixing.activity.groups.groups
        self.functional_groups_multiplicity = (
            self.psrk.mixing.activity.groups.n_groups
        )
        self.r = self.psrk.mixing.activity.params.R.values
        self.q = self.psrk.mixing.activity.params.Q.values
        self.a_pair = self.psrk.mixing.activity.params.A.values
        self.b_pair = self.psrk.mixing.activity.params.B.values
        self.c_pair = self.psrk.mixing.activity.params.C.values

        # =====================================================================
        # Init of the pure fugacity coefficients
        # =====================================================================
        if (self.pure_pressure is None) and (self.pure_temperature is None):
            pass
        else:
            self.pure_fugacity_coefficients = (
                self.pure_fugacity_coefficients_calc(
                    self.pure_pressure, self.pure_temperature
                )
            )

    def mole_fractions(self, moles: np.ndarray) -> np.ndarray:
        """Return mole fractions of the mixture substances.

        Multiple mixtures compositions set can be added as a moles matrix. Each
        row represents each substance and each column represents each mixture
        composition.

        Parameters
        ----------
        moles : np.ndarray
            Moles of each substance. [mol]

        Returns
        -------
        np.ndarray
            Molar fraction of each substance.
        """
        sumatory = np.sum(moles, axis=0)
        return np.divide(moles, sumatory)

    def mass_fractions(self, mole_fractions: np.ndarray) -> np.ndarray:
        """Calculate the mass fractions.

        Use the mole fractions and the molecular weights to calculate the
        mass fractions.

        Parameters
        ----------
        mole_fractions : np.ndarray
            Mole fractions.

        Returns
        -------
        np.ndarray
            Mass fractions.
        """

        numerator = np.multiply(mole_fractions, self.molecular_weights)
        denominator = np.dot(mole_fractions, self.molecular_weights)

        return np.divide(numerator, denominator)

    def molar_density(
        self, pressure: float, temperature: float, mole_fractions: list[float]
    ) -> float:
        """Return the molar density in mol / l of the mixture at pressure,
        temperature and with mole fractions of each mix's substance at the
        chosen phase.

        Parameters
        ----------
        pressure : float
            Pressure. [Pa]
        temperature : float
            Temperature [K]
        moles : list[float]
            Moles of each substance. [mol]

        Returns
        -------
        float
            Molar density of mixture. [l]
        """
        molar_density = (
            Main.molar_density(
                self.psrk,
                pressure,
                temperature,
                mole_fractions,
                phase=self.phase,
            )
        ) / 1000
        return molar_density

    def molar_concentrations(
        self, pressure: float, temperature: float, mole_fractions: np.ndarray
    ) -> np.ndarray:
        """Return the molar concentration of each substance in mol / l of the
        mixture at pressure, temperature and with moles of each mix substance.

        Parameters
        ----------
        pressure : float
            Pressure. [Pa]
        temperature : float
            Temperature [K]
        mole_fractions : list[float]
            Moles of each substance. [mol]

        Returns
        -------
        np.array[float]
            Molar concentration of substances. [mol / l]
        """
        molar_density = self.molar_density(
            pressure, temperature, mole_fractions
        )

        concentrations = np.multiply(mole_fractions, molar_density)

        return concentrations

    def pure_fugacity_coefficients_calc(
        self, pressure: float, temperature: float
    ) -> np.ndarray:
        """Calculate the pure fugacity coefficients of the substances at
        pressure and temperature.

        Parameters
        ----------
        pressure : float
            Pressure. [Pa]
        temperature : float
            Temperature. [K]

        Returns
        -------
        ndarray
            Pure fugacity coefficients.
        """
        pure_fugacities_coeff = np.array([])

        for substance in self.substances:
            psrk_pure = Main.PSRKFULL(
                [substance],
                idealmodel=Main.ReidIdeal,
                userlocations=self.userlocations,
                group_userlocations=self.group_userlocations,
                ideal_userlocations=self.ideal_userlocations,
                alpha_userlocations=self.alpha_userlocations,
                mixing_userlocations=self.mixing_userlocations,
                translation_userlocations=self.translation_userlocations,
            )

            pure_fug = Main.fugacity_coefficient(
                psrk_pure, pressure, temperature
            )

            pure_fugacities_coeff = np.append(pure_fugacities_coeff, pure_fug)

        return pure_fugacities_coeff

    def fugacity_coefficients(
        self, pressure: float, temperature: float, mole_fractions: np.ndarray
    ) -> np.ndarray:
        fugacity_coefficients = Main.fugacity_coefficient(
            self.psrk,
            pressure,
            temperature,
            mole_fractions,
            phase=self.phase,
        )
        return fugacity_coefficients

    def activity(
        self, pressure: float, temperature: float, mole_fractions: list[float]
    ) -> np.ndarray:
        """Return the activity and activity coefficients of each substance of
        the mixture at pressure, temperature and mole fractions of each mix
        substance.

        Parameters
        ----------
        pressure : float
            Pressure. [Pa]
        temperature : float
            Temperature [K]
        moles : list[float]
            Moles of each substance. [mol]

        Returns
        -------
        np.ndarray, np.ndarray
            Activity of substances, activity coefficients.
        """
        # Calculate fugacity coefficients
        fug_coeff = self.fugacity_coefficients(
            pressure, temperature, mole_fractions
        )

        # Calculate activity coefficients
        activity_coeff = np.divide(fug_coeff, self.pure_fugacity_coefficients)

        # Calculate activity
        activity = np.multiply(activity_coeff, mole_fractions)

        return activity, activity_coeff

    def saturation_pressure(
        self, temperature: float, pressure_guess: float = 101325
    ) -> tuple:
        """Calculate the saturation pressure of a one substance PSRK.

        Only use whit 1 substance mixes.

        Parameters
        ----------
        temperature : float
            Temperature [K]
        pressure_guess : float
            Guess to find the saturation pressure [Pa]

        Returns
        -------
        tuple
            Array containing
            [saturation pressure [Pa], v_vapour [m3], v_liquid [m3]]
        """
        results = Main.saturation_pressure(
            self.psrk,
            temperature,
            Main.IsoFugacitySaturation(p0=pressure_guess),
        )

        return results

    def tp_flash(
        self,
        pressure: float,
        temperature: float,
        mole_fractions: np.ndarray,
        x0: np.ndarray,
        y0: np.ndarray,
    ) -> tuple:
        """flash at pressure and temperature.

        Parameters
        ----------
        pressure : float
            Pressure [Pa]
        temperature : float
            Temperature [K]
        mole_fractions : np.ndarray
            Mole fractions of substances.
        x0 : np.ndarray
            Guess of mole fractions of the phase x.
        y0 : np.ndarray
            Guess of mole fractions of the phase y.

        Returns
        -------
        tuple
            Returns a flash tuple with the format:
                (
                    np.array([[x1, x2, x3], [y1, y2, y3]]),
                    np.array([[nx1, nx2, nx3], [ny1, ny2, ny3]]),
                    G,
                )
            where:
            xi are the mole fractions of phase x.
            yi are the mole fractions of phase y.
            nxi are the mole numbers of phase x.
            nyi are the mole numbers of phase y.
            G is the Gibbs Free Energy of Equilibrium Mixture [J]
        """
        michelsen = Main.eval(
            f"""MichelsenTPFlash(
                equilibrium=:{self.flash_nature},
                x0={list(x0)},
                y0={list(y0)},
                ss_iters=50,
                nacc = 0,
                second_order = true
            )"""
        )

        flash = Main.tp_flash(
            self.psrk,
            pressure,
            temperature,
            mole_fractions,
            michelsen,
        )

        return flash
