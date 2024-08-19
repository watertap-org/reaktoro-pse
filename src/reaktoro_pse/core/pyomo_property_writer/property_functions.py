from pyomo.environ import (
    log10,
    log,
)


def build_scaling_tendency_constraint(rktOutputObj):
    user_output_var = rktOutputObj.pyomoVar
    build_properties = rktOutputObj.pyomoBuildOptions.properties
    return (
        user_output_var
        == 10
        ** build_properties[("saturationIndex", rktOutputObj.propertyIndex)].pyomoVar
    )


def build_ph_constraint(rktOutputObj):
    user_output_var = rktOutputObj.pyomoVar
    build_properties = rktOutputObj.pyomoBuildOptions.properties
    return (
        -build_properties[("speciesActivityLn", "H+")].pyomoVar / log(10)
        == user_output_var
    )


def build_osmotic_constraint(rktOutputObj):
    user_output_var = rktOutputObj.pyomoVar
    build_properties = rktOutputObj.pyomoBuildOptions.properties
    return (
        user_output_var
        == -(8.31446261815324 * build_properties[("temperature", None)].pyomoVar)
        / build_properties[
            ("speciesStandardVolume", rktOutputObj.propertyIndex)
        ].pyomoVar
        * build_properties[("speciesActivityLn", rktOutputObj.propertyIndex)].pyomoVar
    )


def build_direct_scaling_tendency_constraint(rktOutputObj):
    # https://reaktoro.org/api/namespaceReaktoro.html#a55b9a29cdf35e98a6b07e67ed2edbc25
    user_output_var = rktOutputObj.pyomoVar
    build_properties = rktOutputObj.pyomoBuildOptions.properties
    build_options = rktOutputObj.pyomoBuildOptions.options
    # TODO: Needs to add temperature unit verificaiton and pressure unit verification
    temperature_var = build_properties[("temperature", None)].pyomoVar
    if build_options["logk_type"] == "Analytical":
        A_params = build_options["logk_paramters"]
        log_k = [A_params["A1"]]
        # temp dependence for phreeqc
        if A_params["A2"] != 0:
            log_k.append(A_params["A2"] * temperature_var)
        if A_params["A3"] != 0:
            log_k.append(A_params["A3"] * temperature_var**-1)
        if A_params["A4"] != 0:
            log_k.append(A_params["A4"] * log10(temperature_var))
        if A_params["A5"] != 0:
            log_k.append(A_params["A5"] * temperature_var**-2)
        if A_params["A6"] != 0:
            log_k.append(A_params["A6"] * temperature_var**2)

    if build_options["logk_type"] == "VantHoff":
        vfparams = build_options["logk_paramters"]
        log_k = [
            vfparams["lgKr"]
            - vfparams["dHr"]
            / 8.31446261815324
            * (1 / temperature_var - 1 / vfparams["Tr"])
        ]
    # pressure dependenance
    log_k.append(
        -(
            build_options["delta_V"]
            * (build_properties[("pressure", None)].pyomoVar - 101325)
            / (log(10) * 8.31446261815324 * temperature_var)
        )
    )

    activities = []
    for key, obj in build_properties.items():
        if "speciesActivityLn" in key:
            activities.append(obj.pyomoVar * obj.stoichiometricCoeff / log(10))
    return user_output_var == 10 ** (
        sum(activities)
        + sum(
            log_k
        )  # this is postive here and in log10 fom, so we add instead of subtract
    )
