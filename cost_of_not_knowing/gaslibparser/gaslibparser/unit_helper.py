import pint

try:
    import importlib.resources as pkg_resources
except ImportError:
    # python < 3.7
    import importlib_resources as pkg_resources

unit = pint.UnitRegistry(system="mks")
gaslib_units = pkg_resources.open_text(__package__, "gaslib_units.txt")
unit.load_definitions(gaslib_units)
