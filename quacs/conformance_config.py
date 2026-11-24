"""A QUACS configuration that checks whether CheMPAS runs QUACS parameterizations correctly.

Set `config_quacs_config_file` to the path of this file and `config_quacs_check_conformance`
to true in the `&quacs` namelist record. MPAS then writes the results of the checks to its
log at the end of a run.
"""

from quacs.conformance import VerticalOrderCheck, default_checks

checks = default_checks()

vertical_order = VerticalOrderCheck("base_state_pressure", "Pa", increases_upward=False)

parameterizations = [*checks, vertical_order]

hooks = {
    "initial_conditions": [check.name for check in checks],
    "before_dynamics": [check.name for check in parameterizations],
}
