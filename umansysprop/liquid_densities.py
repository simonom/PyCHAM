# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright (c) 2016 David Topping.
# All Rights Reserved.
# This file is part of umansysprop.
#
# umansysprop is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# umansysprop is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# umansysprop.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')


from . import groups
from . import molar_volumes
from . import boiling_points
import pybel



def girolami(compound):

    comp = groups.composition(compound)
    m = groups.girolami(compound)
    m['H'] = comp['H']
    volumes = {
        'H':              1.0,
        'short_period_1': 2.0,
        'short_period_2': 4.0,
        'long_period_1':  5.0,
        'long_period_2':  7.5,
        'long_period_3':  9.0,
        }
    result = groups.aggregate_matches(
        m, coefficients=volumes, groups=volumes.keys())
    # -----------SPO edit 24/04/2018-------------------------------------------
    # note that when compound is [H][H], the function here doesn't recognise it
    # leading to an error, so instead make up some values, which don't matter
    # because molecular hydrogen doesn't condense to the particle-phase anyway.
    # Testing this code on test_MCM.eqn inputs shows that only [H][H] suffers 
    # from this issue.  This issue also applies to [Na]
    if result==0 or comp['mass']==0:
    	
    	comp['mass'] = 1.0
    	result=5.0

    # -------------------------------------------------------------------------
    	
    result = comp['mass'] / (5.0 * result)
    # Make adjustments based on the presence of specific groups
    group_multipliers = {
        'hydroxyl_groups':                10.0,
        'hydroperoxide_groups':           10.0,
        'peroxyacid_groups':              10.0,
        'carboxylic_acid_groups':         10.0,
        'primary_secondary_amino_groups': 10.0,
        'amide_groups':                   10.0,
        'unfused_rings':                  10.0,
        'fused_rings':                    7.5,
        }
    percentage_increase = min(30.0, groups.aggregate_matches(
        m, coefficients=group_multipliers,groups=group_multipliers.keys()
        ))
    result += result * percentage_increase / 100.0
    return result


def _generic(
        compound, temperature, molar_volume, critical_properties):
    comp = groups.composition(compound)
    boiling_point = boiling_points.nannoolal(compound)
    phi = (
        (1 - (temperature / critical_properties.temperature)) ** (2/7) -
        (1 - (boiling_point / critical_properties.temperature)) ** (2/7)
        )
    result = molar_volume * critical_properties.compressibility ** phi
    return comp['mass'] / result


def schroeder(
        compound, temperature, critical_properties):
    return _generic(
        compound, temperature,
        molar_volumes.schroeder(compound),
        critical_properties)


def le_bas(
        compound, temperature, critical_properties):
    return _generic(
        compound, temperature,
        molar_volumes.le_bas(compound),
        critical_properties)


def tyn_and_calus(
        compound, temperature, critical_properties):
    return _generic(
        compound, temperature,
        molar_volumes.tyn_and_calus(critical_properties.volume),
        critical_properties)

