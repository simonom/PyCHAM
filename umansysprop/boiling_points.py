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


from . import data
from . import groups


def joback_and_reid(compound):
    return groups.aggregate_matches(
        groups.stein_and_brown(compound),
        data.JOBACK_BOILING_POINT
        ) + 198.2


def stein_and_brown(compound):
    result = groups.aggregate_matches(
        groups.stein_and_brown(compound),
        data.STEIN_AND_BROWN_BOILING_POINT
        ) + 198.2
    return (
        (0.4791 * result) + 282.7
        if result > 700 else
        (1.5577 * result) - (0.0007705 * result**2) - 94.84
        )


def nannoolal(compound):
    comp = groups.composition(compound)
    result1 = groups.aggregate_matches(
        groups.nannoolal_primary(compound),
        data.NANNOOLAL_BOILING_POINT_PRIMARY
        )
    result2 = groups.aggregate_matches(
        groups.nannoolal_secondary(compound),
        data.NANNOOLAL_BOILING_POINT_SECONDARY
        )
    result3 = groups.aggregate_interactions(
        compound,
        groups.nannoolal_interactions(compound),
        data.NANNOOLAL_BOILING_POINT_INTERACTIONS
        )
    num = result1 + result2 + result3
    denom = (comp['n-non-H'] ** 0.6583) + 1.6868
    return (num / denom) + 84.3395

