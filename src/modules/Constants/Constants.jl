# Constants.jl

export light_speed, electron_charge, reduced_planck_constant, electron_mass
export vacuum_permeability, electron_rest_energy, vacuum_permitticity
export electron_rest_energy_eV, electron_rest_energy_MeV, electron_rest_energy_GeV, electron_radius



"""
Speed of light in vacuum: ``c = $(light_speed)\\ [m/s]``
"""
const light_speed::Float64 = 299792458

"""
Elementary charge of an electron: ``e = $(electron_charge)\\ [C]``
"""
const electron_charge::Float64 = 1.602176634e-19

"""
Reduced Planck constant, also known as h-bar: ``ħ = $(reduced_planck_constant)\\ [J.s]``
"""
const reduced_planck_constant::Float64 = 1.054571817e-34

"""
Mass of an electron: ``mₑ = $(electron_mass)\\ [kg]``
"""
const electron_mass::Float64 = 9.1093837015e-31

"""
Vacuum permeability, also known as the magnetic constant: ``μ₀ = $(vacuum_permeability)\\ [T.m/A]``
"""
const vacuum_permeability::Float64 = 1.25663706212e-6

"""
Rest energy of an electron: ``E₀ = $(electron_rest_energy)\\ [kg.m^2/s^2]``
"""
const electron_rest_energy::Float64 = electron_mass * light_speed^2

"""
Vacuum permittivity, also known as the electric constant: ``ε₀ = $(vacuum_permitticity)\\ [V.s/(A.m)]``
"""
const vacuum_permitticity::Float64 = 1 / (vacuum_permeability * light_speed^2)

"""
Rest energy of an electron in electronvolts (eV): ``E₀ = $(electron_rest_energy_eV)\\ [eV]``
"""
const electron_rest_energy_eV::Float64 = electron_rest_energy / electron_charge

"""
Rest energy of an electron in megaelectronvolts (MeV): ``E₀ = $(electron_rest_energy_MeV)\\ [MeV]``
"""
const electron_rest_energy_MeV::Float64 = electron_rest_energy_eV / 1e6

"""
Rest energy of an electron in gigaelectronvolts (GeV): ``E₀ = $(electron_rest_energy_GeV)\\ [GeV]``
"""
const electron_rest_energy_GeV::Float64 = electron_rest_energy_eV / 1e9

"""
Classical electron radius, also known as the Thomson scattering radius: ``rₑ = $(electron_radius)\\ [m]``
"""
const electron_radius::Float64 = electron_charge^2 / (4 * pi * vacuum_permitticity * electron_rest_energy)
