import math

pressure_unit = "bar"
flow_unit = "m_cube_per_hour"
mass_flow_unit = "kg_per_second"
length_unit = "meter"

speed_of_sound = 466  # todo taken from wikipedia for methan
norm_density = 0.87


def nikuradse(diameter, roughness):
    friction = math.pow(2 * math.log10(diameter / roughness) + 1.138, -2)
    assert friction > 0
    return friction


def norm_vol_per_hour_to_mass_flow(norm_vol_per_hour):
    return norm_density * norm_vol_per_hour / 3600.0


def mass_flow_to_norm_vol_per_hour(mass_flow):
    return 3600.0 * mass_flow / norm_density
