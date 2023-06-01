import numpy as np

import psrk


# =============================================================================
# Test individual parameters loading
# =============================================================================
def test_database_ether():
    model = psrk.PSRK(["diisopropyl ether"])

    # SRK parameters
    assert (model.critical_temperatures == 5.00050e02).all()
    assert (model.critical_pressures == 2.88000e06).all()
    assert (model.c1 == 0.91199).all()
    assert (model.c2 == 0.82949).all()
    assert (model.c3 == -2.19890).all()
    assert (model.v_shift == 1.0456512672945569e-05).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["CH3", "CH", "CHO"]
    assert (model.functional_groups_multiplicity[0] == [4, 1, 1]).all()


def test_database_toluene():
    model = psrk.PSRK(["toluene"])

    # SRK parameters
    assert (model.critical_temperatures == 5.91750e02).all()
    assert (model.critical_pressures == 4.10800e06).all()
    assert (model.c1 == 0.94689).all()
    assert (model.c2 == -0.58961).all()
    assert (model.c3 == 1.21320).all()
    assert (model.v_shift == 1.4021522506795162e-5).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["ACH", "ACCH3"]
    assert (model.functional_groups_multiplicity[0] == [5, 1]).all()


def test_database_hexane():
    model = psrk.PSRK(["hexane"])

    # SRK parameters
    assert (model.critical_temperatures == 5.07600e02).all()
    assert (model.critical_pressures == 3.02500e06).all()
    assert (model.c1 == 0.93264).all()
    assert (model.c2 == 0).all()
    assert (model.c3 == 0).all()
    assert (model.v_shift == 1.569922941144466e-05).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["CH3", "CH2"]
    assert (model.functional_groups_multiplicity[0] == [2, 4]).all()


def test_database_heptane():
    model = psrk.PSRK(["heptane"])

    # SRK parameters
    assert (model.critical_temperatures == 5.40200e02).all()
    assert (model.critical_pressures == 2.74000e06).all()
    assert (model.c1 == 1.00313).all()
    assert (model.c2 == 0).all()
    assert (model.c3 == 0).all()
    assert (model.v_shift == 2.087190187010664e-05).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["CH3", "CH2"]
    assert (model.functional_groups_multiplicity[0] == [2, 5]).all()


def test_database_acetic_acid():
    model = psrk.PSRK(["acetic acid"])

    # SRK parameters
    assert (model.critical_temperatures == 5.91950e02).all()
    assert (model.critical_pressures == 5.78600e06).all()
    assert (model.c1 == 1.38764).all()
    assert (model.c2 == -1.58619).all()
    assert (model.c3 == 1.84943).all()
    assert (model.v_shift == 2.6814912786435057e-05).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["CH3", "COOH"]
    assert (model.functional_groups_multiplicity[0] == [1, 1]).all()


def test_database_octanoic_acid():
    model = psrk.PSRK(["octanoic acid"])

    # SRK parameters
    assert (model.critical_temperatures == 6.94260e02).all()
    assert (model.critical_pressures == 2.77900e06).all()
    assert (model.c1 == 1.60513).all()
    assert (model.c2 == -0.21195).all()
    assert (model.c3 == 1.15479).all()
    assert (model.v_shift == 3.4413410889015775e-05).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["CH3", "CH2", "COOH"]
    assert (model.functional_groups_multiplicity[0] == [1, 6, 1]).all()


def test_database_peracetic_acid():
    model = psrk.PSRK(["peracetic acid"])

    # SRK parameters
    assert (model.critical_temperatures == 5.52000e02).all()
    assert (model.critical_pressures == 6.40000e06).all()
    assert (model.c1 == 1.88407).all()
    assert (model.c2 == -2.36351).all()
    assert (model.c3 == 2.46947).all()
    assert (model.v_shift == 8.033984089651411e-06).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["CH3COO", "OH"]
    assert (model.functional_groups_multiplicity[0] == [1, 1]).all()


def test_database_peroctanoic_acid():
    model = psrk.PSRK(["peroctanoic acid"])

    # SRK parameters
    assert (model.critical_temperatures == 631.54436074476).all()
    assert (model.critical_pressures == 2639601.95207736).all()
    assert (model.c1 == 1.81245).all()
    assert (model.c2 == -0.06092).all()
    assert (model.c3 == 1.71043).all()
    assert (model.v_shift == 8.414016700e-06).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["CH3", "CH2", "CH2COO", "OH"]
    assert (model.functional_groups_multiplicity[0] == [1, 5, 1, 1]).all()


def test_database_water():
    model = psrk.PSRK(["water"])

    # SRK parameters
    assert (model.critical_temperatures == 6.47130e02).all()
    assert (model.critical_pressures == 2.20550e07).all()
    assert (model.c1 == 1.07830).all()
    assert (model.c2 == -0.58321).all()
    assert (model.c3 == 0.54619).all()
    assert (model.v_shift == 5.877439824702498e-06).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["H2O"]
    assert (model.functional_groups_multiplicity[0] == [1]).all()


def test_database_hydrogen_peroxide():
    model = psrk.PSRK(["hydrogen peroxide"])

    # SRK parameters
    assert (model.critical_temperatures == 7.30150e02).all()
    assert (model.critical_pressures == 2.17000e07).all()
    assert (model.c1 == 1.06615).all()
    assert (model.c2 == -0.39235).all()
    assert (model.c3 == 0.27125).all()
    assert (model.v_shift == 3.218736682603898e-06).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["OH"]
    assert (model.functional_groups_multiplicity[0] == [2]).all()


def test_database_limonene():
    model = psrk.PSRK(["limonene"])

    # SRK parameters
    assert (model.critical_temperatures == 6.53000e02).all()
    assert (model.critical_pressures == 2.82000e06).all()
    assert (model.c1 == 1.11438).all()
    assert (model.c2 == -0.62643).all()
    assert (model.c3 == 1.06066).all()
    assert (model.v_shift == 2.5292802413353628e-05).all()

    # UNIFAC parameters
    assert model.functional_groups[0] == ["CH3", "CH2", "CH", "CH2=C", "CH=C"]
    assert (model.functional_groups_multiplicity[0] == [2, 3, 1, 1, 1]).all()


# =============================================================================
# Test mix parameters loading
# =============================================================================
def test_pair_database_mix1():
    model = psrk.PSRK(["water", "hexane"])

    # SRK parameters
    assert (model.critical_temperatures == [6.47130e02, 5.07600e02]).all()
    assert (model.critical_pressures == [2.20550e07, 3.02500e06]).all()
    assert (model.c1 == [1.07830, 0.93264]).all()
    assert (model.c2 == [-0.58321, 0]).all()
    assert (model.c3 == [0.54619, 0]).all()
    assert (
        model.v_shift == [5.877439824702498e-06, 1.569922941144466e-05]
    ).all()

    # UNIFAC parameters
    assert model.functional_groups == [["H2O"], ["CH3", "CH2"]]
    assert (model.functional_groups_multiplicity[0] == [1]).all()
    assert (model.functional_groups_multiplicity[1] == [2, 4]).all()


def test_pair_database_mix2():
    model = psrk.PSRK(["acetic acid", "peracetic acid"])

    # SRK parameters
    assert (model.critical_temperatures == [5.91950e02, 5.52000e02]).all()
    assert (model.critical_pressures == [5.78600e06, 6.40000e06]).all()
    assert (model.c1 == [1.38764, 1.88407]).all()
    assert (model.c2 == [-1.58619, -2.36351]).all()
    assert (model.c3 == [1.84943, 2.46947]).all()
    assert (
        model.v_shift == [2.6814912786435057e-05, 8.033984089651411e-06]
    ).all()

    # UNIFAC parameters
    assert model.functional_groups == [["CH3", "COOH"], ["CH3COO", "OH"]]
    assert (model.functional_groups_multiplicity[0] == [1, 1]).all()
    assert (model.functional_groups_multiplicity[1] == [1, 1]).all()


def test_trio_database_mix():
    model = psrk.PSRK(["water", "hexane", "limonene"])

    # SRK parameters
    assert (
        model.critical_temperatures == [6.47130e02, 5.07600e02, 6.53000e02]
    ).all()
    assert (
        model.critical_pressures == [2.20550e07, 3.02500e06, 2.82000e06]
    ).all()
    assert (model.c1 == [1.07830, 0.93264, 1.11438]).all()
    assert (model.c2 == [-0.58321, 0, -0.62643]).all()
    assert (model.c3 == [0.54619, 0, 1.06066]).all()
    assert (
        model.v_shift
        == [
            5.877439824702498e-06,
            1.569922941144466e-05,
            2.5292802413353628e-05,
        ]
    ).all()

    # UNIFAC parameters
    assert model.functional_groups == [
        ["H2O"],
        ["CH3", "CH2"],
        ["CH3", "CH2", "CH", "CH2=C", "CH=C"],
    ]
    assert (model.functional_groups_multiplicity[0] == [1]).all()
    assert (model.functional_groups_multiplicity[1] == [2, 4]).all()
    assert (model.functional_groups_multiplicity[2] == [2, 3, 1, 1, 1]).all()


def test_quad_database_mix():
    model = psrk.PSRK(["water", "hexane", "limonene", "hydrogen peroxide"])

    # SRK parameters
    assert (
        model.critical_temperatures
        == [6.47130e02, 5.07600e02, 6.53000e02, 7.30150e02]
    ).all()
    assert (
        model.critical_pressures
        == [2.20550e07, 3.02500e06, 2.82000e06, 2.17000e07]
    ).all()
    assert (model.c1 == [1.07830, 0.93264, 1.11438, 1.06615]).all()
    assert (model.c2 == [-0.58321, 0, -0.62643, -0.39235]).all()
    assert (model.c3 == [0.54619, 0, 1.06066, 0.27125]).all()
    assert (
        model.v_shift
        == [
            5.877439824702498e-06,
            1.569922941144466e-05,
            2.5292802413353628e-05,
            3.218736682603898e-06,
        ]
    ).all()

    # UNIFAC parameters
    assert model.functional_groups == [
        ["H2O"],
        ["CH3", "CH2"],
        ["CH3", "CH2", "CH", "CH2=C", "CH=C"],
        ["OH"],
    ]
    assert (model.functional_groups_multiplicity[0] == [1]).all()
    assert (model.functional_groups_multiplicity[1] == [2, 4]).all()
    assert (model.functional_groups_multiplicity[2] == [2, 3, 1, 1, 1]).all()
    assert (model.functional_groups_multiplicity[3] == [2]).all()


def test_fifth_database_mix():
    model = psrk.PSRK(
        ["water", "hexane", "limonene", "hydrogen peroxide", "octanoic acid"]
    )

    # SRK parameters
    assert (
        model.critical_temperatures
        == [6.47130e02, 5.07600e02, 6.53000e02, 7.30150e02, 6.94260e02]
    ).all()
    assert (
        model.critical_pressures
        == [2.20550e07, 3.02500e06, 2.82000e06, 2.17000e07, 2.77900e06]
    ).all()
    assert (model.c1 == [1.07830, 0.93264, 1.11438, 1.06615, 1.60513]).all()
    assert (model.c2 == [-0.58321, 0, -0.62643, -0.39235, -0.21195]).all()
    assert (model.c3 == [0.54619, 0, 1.06066, 0.27125, 1.15479]).all()
    assert (
        model.v_shift
        == [
            5.877439824702498e-06,
            1.569922941144466e-05,
            2.5292802413353628e-05,
            3.218736682603898e-06,
            3.4413410889015775e-05,
        ]
    ).all()

    # UNIFAC parameters
    assert model.functional_groups == [
        ["H2O"],
        ["CH3", "CH2"],
        ["CH3", "CH2", "CH", "CH2=C", "CH=C"],
        ["OH"],
        ["CH3", "CH2", "COOH"],
    ]
    assert (model.functional_groups_multiplicity[0] == [1]).all()
    assert (model.functional_groups_multiplicity[1] == [2, 4]).all()
    assert (model.functional_groups_multiplicity[2] == [2, 3, 1, 1, 1]).all()
    assert (model.functional_groups_multiplicity[3] == [2]).all()
    assert (model.functional_groups_multiplicity[4] == [1, 6, 1]).all()


def test_sixth_database_mix():
    model = psrk.PSRK(
        [
            "water",
            "hexane",
            "limonene",
            "hydrogen peroxide",
            "octanoic acid",
            "peroctanoic acid",
        ]
    )

    # SRK parameters
    assert (
        model.critical_temperatures
        == [
            6.47130e02,
            5.07600e02,
            6.53000e02,
            7.30150e02,
            6.94260e02,
            631.54436074476,
        ]
    ).all()
    assert (
        model.critical_pressures
        == [
            2.20550e07,
            3.02500e06,
            2.82000e06,
            2.17000e07,
            2.77900e06,
            2639601.95207736,
        ]
    ).all()
    assert (
        model.c1 == [1.07830, 0.93264, 1.11438, 1.06615, 1.60513, 1.81245]
    ).all()
    assert (
        model.c2 == [-0.58321, 0, -0.62643, -0.39235, -0.21195, -0.06092]
    ).all()
    assert (model.c3 == [0.54619, 0, 1.06066, 0.27125, 1.15479, 1.71043]).all()
    assert (
        model.v_shift
        == [
            5.877439824702498e-06,
            1.569922941144466e-05,
            2.5292802413353628e-05,
            3.218736682603898e-06,
            3.4413410889015775e-05,
            8.414016700e-06,
        ]
    ).all()

    # UNIFAC parameters
    assert model.functional_groups == [
        ["H2O"],
        ["CH3", "CH2"],
        ["CH3", "CH2", "CH", "CH2=C", "CH=C"],
        ["OH"],
        ["CH3", "CH2", "COOH"],
        ["CH3", "CH2", "CH2COO", "OH"],
    ]
    assert (model.functional_groups_multiplicity[0] == [1]).all()
    assert (model.functional_groups_multiplicity[1] == [2, 4]).all()
    assert (model.functional_groups_multiplicity[2] == [2, 3, 1, 1, 1]).all()
    assert (model.functional_groups_multiplicity[3] == [2]).all()
    assert (model.functional_groups_multiplicity[4] == [1, 6, 1]).all()
    assert (model.functional_groups_multiplicity[5] == [1, 5, 1, 1]).all()


# =============================================================================
# Individual molar concentration
# =============================================================================
def test_density_ether():
    a, b, c, d = (6.9213e-01, 2.6974e-01, 5.0005e02, 2.8571e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["diisopropyl ether"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_toluene():
    a, b, c, d = (8.7920e-01, 2.7136e-01, 5.9175e02, 2.9241e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["toluene"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_hexane():
    a, b, c, d = (7.0824e-01, 2.6411e-01, 5.0760e02, 2.7537e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["hexane"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_heptane():
    a, b, c, d = (6.1259e-01, 2.6211e-01, 5.4020e02, 2.8141e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["heptane"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_acetic_acid():
    a, b, c, d = (1.4486e00, 2.5892e-01, 5.9195e02, 2.5290e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["acetic acid"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_octanoic_acid():
    a, b, c, d = (5.2887e-01, 2.6354e-01, 6.9426e02, 2.8099e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["octanoic acid"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_peracetic_acid():
    a, b, c, d = (1.2274e00, 2.4307e-01, 5.5200e02, 2.6690e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["peracetic acid"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_water():
    # Recordar que aca va la dippr116
    a, b, c, d, e = (1.7863e01, 5.8606e01, -9.5396e01, 2.1389e02, -1.4126e02)

    def dippr116(temperature):
        """this one is only for water"""
        tc = 6.47130e02
        tr = 1 - temperature / tc

        density = (
            a + b * tr**0.35 + c * tr ** (2 / 3) + d * tr + e * tr ** (4 / 3)
        )
        return density  # mol/L

    model = psrk.PSRK(["water"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr116(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_hydrogen_peroxide():
    a, b, c, d = (3.2151e00, 2.4982e-01, 7.3015e02, 2.8770e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["hydrogen peroxide"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


def test_density_limonene():
    a, b, c, d = (5.9659e-01, 2.8040e-01, 6.5300e02, 2.9257e-01)

    def dippr105(temperature):
        density = a / (b ** (1 + (1 - temperature / c) ** d))

        return density  # mol/L

    model = psrk.PSRK(["limonene"])

    density_model = model.molar_concentrations(101325, 303.15, [1])  # mol/L
    density_dippr = dippr105(303.15)  # mol/L

    assert np.allclose(density_model, density_dippr, atol=1e-10, rtol=1e-10)


# =============================================================================
# Individual saturation pressures
# =============================================================================
def test_psat_ether():
    """Datos del ajuste

    temperatura_min=187.65 temperatura_max=500.05 parametros=bibliografiapsrk
    """
    a, b, c, d, e = (4.1631e01, -4.6687e03, -2.8551e00, 6.3693e-04, 1.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["diisopropyl ether"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 0.27, atol=1e-3, rtol=1e-3)


def test_psat_toluene():
    """Datos del ajuste

    temperatura_min=178.18 temperatura_max=591.75 parametros=bibliografiapsrk
    """
    a, b, c, d, e = (7.6945e01, -6.7298e03, -8.1790e00, 5.3017e-06, 2.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["toluene"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 0.207845, atol=1e-6, rtol=1e-6)


def test_psat_hexane():
    """Datos del ajuste

    temperatura_min=177.83 temperatura_max=507.6 parametros=bibliografiapsrk
    """
    a, b, c, d, e = (1.0465e02, -6.9955e03, -1.2702e01, 1.2381e-05, 2.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["hexane"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 1.174614, atol=1e-6, rtol=1e-6)


def test_psat_heptane():
    """Datos del ajuste

    temperatura_min=182.57 temperatura_max=540.2 parametros=bibliografiapsrk
    """
    a, b, c, d, e = (8.7829e01, -6.9964e03, -9.8802e00, 7.2099e-06, 2.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["heptane"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 0.211479, atol=1e-6, rtol=1e-6)


def test_psat_acetic_acid():
    """Datos del ajuste

    temperatura_min=289.81 temperatura_max=591.95 parametros=ajustado
    """
    a, b, c, d, e = (5.3270e01, -6.3045e03, -4.2985e00, 8.8865e-18, 6.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["acetic acid"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 0.130854, atol=1e-6, rtol=1e-6)


def test_psat_octanoic_acid():
    """Datos del ajuste

    temperatura_min=289.65 temperatura_max=694.26 parametros=ajustado
    """
    a, b, c, d, e = (1.4016e02, -1.4813e04, -1.6004e01, 6.4239e-18, 6.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["octanoic acid"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 1.211635, atol=1e-6, rtol=1e-6)


def test_psat_peracetic_acid():
    """Datos del ajuste

    temperatura_min=272.95 temperatura_max=552 parametros=ajustado
    """
    a, b, c, d, e = (6.5788e01, -7.3120e03, -5.9242e00, 1.8660e-17, 6.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["peracetic acid"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 0.046691, atol=1e-6, rtol=1e-6)


def test_psat_water():
    """Datos del ajuste

    temperatura_min=273.16 temperatura_max=647.13 parametros=bibliografiapsrk
    """
    a, b, c, d, e = (7.3649e01, -7.2582e03, -7.3037e00, 4.1653e-06, 2.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["water"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 1.018524, atol=1e-6, rtol=1e-6)


def test_psat_hydrogen_peroxide():
    """Datos del ajuste

    temperatura_min=272.74 temperatura_max=730.15 parametros=ajustado
    """
    a, b, c, d, e = (1.1263e02, -9.4555e03, -1.3826e01, 1.1470e-02, 1.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["hydrogen peroxide"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 1.381962, atol=1e-6, rtol=1e-6)


def test_psat_limonene():
    """Datos del ajuste

    temperatura_min=198.8 temperatura_max=653 parametros=ajustado
    """
    a, b, c, d, e = (7.5574e01, -8.0797e03, -7.5596e00, 8.3872e-18, 6.0000e00)

    def dippr101(temperature):
        t = temperature
        p_sat = np.exp(a + b / t + c * np.log(t) + d * t**e)

        return p_sat

    model = psrk.PSRK(["limonene"])

    psat_model = model.saturation_pressure(303.15, 101325)[0]
    psat_dippr = dippr101(303.15)  # mol/L

    r_error = np.abs(psat_dippr - psat_model) / psat_dippr * 100

    assert np.allclose(r_error, 0.15714057766761058, atol=1e-6, rtol=1e-6)


# =============================================================================
# test consistency
# =============================================================================
def test_flash_consistency_1():
    model = psrk.PSRK(["water", "hexane"])

    x0 = [0.999, 0.001]
    y0 = [0.001, 0.999]

    f1, f2, f3, f4, f5 = (
        model.tp_flash(101325, 303.15, [0.5, 0.5], x0, y0),
        model.tp_flash(101325, 310.15, [0.5, 0.5], x0, y0),
        model.tp_flash(101325, 320.15, [0.5, 0.5], x0, y0),
        model.tp_flash(101325, 350.15, [0.5, 0.5], x0, y0),
        model.tp_flash(101325, 370.15, [0.5, 0.5], x0, y0),
    )

    flashes = [f1[0], f2[0], f3[0], f4[0], f5[0]]
    temperatures = [303.15, 310.15, 320.15, 350.15, 370.15]

    for flash, t in zip(flashes, temperatures):
        fug_coeff_phase_x = model.fugacity_coefficients(101325, t, flash[0])
        fug_coeff_phase_y = model.fugacity_coefficients(101325, t, flash[1])

        fug_phase_x = fug_coeff_phase_x * flash[0] * 101325
        fug_phase_y = fug_coeff_phase_y * flash[1] * 101325

        assert np.allclose(fug_phase_x, fug_phase_y, atol=0.0001)

        activity_x = model.activity(101325, t, flash[0])[0]
        activity_y = model.activity(101325, t, flash[1])[0]

        assert np.allclose(activity_x, activity_y, atol=0.00001)

        activity_coeff_x = model.activity(101325, t, flash[0])[1]
        activity_coeff_y = model.activity(101325, t, flash[1])[1]

        activity_hand_x = activity_coeff_x * flash[0]
        activity_hand_y = activity_coeff_y * flash[1]

        assert (activity_hand_x == activity_x).all()
        assert (activity_hand_y == activity_y).all()


def test_flash_consistency_2():
    model = psrk.PSRK(
        [
            "octanoic acid",
            "hydrogen peroxide",
            "water",
            "toluene",
            "peroctanoic acid",
        ]
    )

    x0 = [0.001, 0.3, 0.7, 0.000001, 0.001]
    y0 = [0.02, 0.003, 0.001, 0.8, 0.06]

    f1, f2, f3 = (
        model.tp_flash(101325, 303.15, [0.3, 0.1, 0.1, 0.4, 0.1], x0, y0),
        model.tp_flash(101325, 310.15, [0.3, 0.1, 0.1, 0.4, 0.1], x0, y0),
        model.tp_flash(101325, 320.15, [0.3, 0.1, 0.1, 0.4, 0.1], x0, y0),
    )

    flashes = [f1[0], f2[0], f3[0]]
    temperatures = [303.15, 310.15, 320.15, 350.15, 370.15]

    for flash, t in zip(flashes, temperatures):
        fug_coeff_phase_x = model.fugacity_coefficients(101325, t, flash[0])
        fug_coeff_phase_y = model.fugacity_coefficients(101325, t, flash[1])

        fug_phase_x = fug_coeff_phase_x * flash[0] * 101325
        fug_phase_y = fug_coeff_phase_y * flash[1] * 101325

        assert np.allclose(fug_phase_x, fug_phase_y, atol=0.0001)

        activity_x = model.activity(101325, t, flash[0])[0]
        activity_y = model.activity(101325, t, flash[1])[0]

        assert np.allclose(activity_x, activity_y, atol=0.00001)

        activity_coeff_x = model.activity(101325, t, flash[0])[1]
        activity_coeff_y = model.activity(101325, t, flash[1])[1]

        activity_hand_x = activity_coeff_x * flash[0]
        activity_hand_y = activity_coeff_y * flash[1]

        assert (activity_hand_x == activity_x).all()
        assert (activity_hand_y == activity_y).all()
