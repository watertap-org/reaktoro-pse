from pyomo.environ import log10, log, exp


def build_scaling_tendency_constraint(rkt_output_object):
    user_output_var = rkt_output_object.pyomo_var
    build_properties = rkt_output_object.pyomo_build_options.properties
    return (
        user_output_var
        == 10
        ** build_properties[
            ("saturationIndex", rkt_output_object.property_index)
        ].pyomo_var
    )


def build_ph_constraint(rkt_output_object):
    user_output_var = rkt_output_object.pyomo_var
    build_properties = rkt_output_object.pyomo_build_options.properties
    return (
        -build_properties[("speciesActivityLn", "H+")].pyomo_var / log(10)
        == user_output_var
    )


def build_vapor_pressure_constraint(rkt_output_object):
    user_output_var = rkt_output_object.pyomo_var
    build_properties = rkt_output_object.pyomo_build_options.properties
    return (
        exp(
            build_properties[
                ("speciesActivityLn", rkt_output_object.property_index)
            ].pyomo_var
        )
        * 101325
        == user_output_var
    )


# work around provided in https://github.com/reaktoro/reaktoro/discussions/398
# can not be implemented as reference chemical potential is function of p and t
# def build_vapor_pressure_direct_constraint(rkt_output_object):
#     user_output_var = rkt_output_object.pyomo_var
#     build_properties = rkt_output_object.pyomo_build_options.properties
#     build_options = rkt_output_object.pyomo_build_options.options
#     return (
#         exp(
#             build_properties[
#                 ("speciesChemicalPotential", rkt_output_object.property_index)
#             ].pyomo_var
#             - build_options["reference_chemical_potential"]
#         )
#         / (
#             build_options["gas_constant"]
#             * build_properties[("temperature", None)].pyomo_var
#         )
#         * 1e5
#         == user_output_var
#     )


def build_osmotic_constraint(rkt_output_object):
    user_output_var = rkt_output_object.pyomo_var
    build_properties = rkt_output_object.pyomo_build_options.properties
    build_options = rkt_output_object.pyomo_build_options.options
    return (
        user_output_var
        == -(
            build_options["gas_constant"]
            * build_properties[("temperature", None)].pyomo_var
        )
        / build_properties[
            ("speciesStandardVolume", rkt_output_object.property_index)
        ].pyomo_var
        * build_properties[
            ("speciesActivityLn", rkt_output_object.property_index)
        ].pyomo_var
    )


def build_direct_scaling_tendency_constraint(rkt_output_object):
    # https://reaktoro.org/api/namespaceReaktoro.html#a55b9a29cdf35e98a6b07e67ed2edbc25
    user_output_var = rkt_output_object.pyomo_var
    build_properties = rkt_output_object.pyomo_build_options.properties
    build_options = rkt_output_object.pyomo_build_options.options
    # TODO: Needs to add temperature unit verificaiton and pressure unit verification
    temperature_var = build_properties[("temperature", None)].pyomo_var
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
            / build_options["gas_constant"]
            * (1 / temperature_var - 1 / vfparams["Tr"])
        ]
    # pressure dependenance
    log_k.append(
        -(
            build_options["delta_V"]
            * (build_properties[("pressure", None)].pyomo_var - 101325)
            / (log(10) * build_options["gas_constant"] * temperature_var)
        )
    )

    activities = []
    for key, obj in build_properties.items():
        if "speciesActivityLn" in key:
            activities.append(obj.pyomo_var * obj.stoichiometric_coeff / log(10))
    return user_output_var == 10 ** (
        sum(activities)
        + sum(
            log_k
        )  # this is postive here and in log10 fom, so we add instead of subtract
    )
