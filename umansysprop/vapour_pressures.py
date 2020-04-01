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


from math import log, sqrt

from . import data
from . import groups


def nannoolal(compound, temperature, boiling_point):
    result1 = groups.aggregate_matches(
        groups.nannoolal_primary(compound),
        data.NANNOOLAL_VAPOUR_PRESSURE_PRIMARY
        )
    result2 = groups.aggregate_matches(
        groups.nannoolal_secondary(compound),
        data.NANNOOLAL_VAPOUR_PRESSURE_SECONDARY
        )
    result3 = groups.aggregate_interactions(
        compound,
        groups.nannoolal_interactions(compound),
        data.NANNOOLAL_VAPOUR_PRESSURE_INTERACTIONS
        )
    # eqn. 6 of Nannoolal et al. (2008), with dB of that equation given by eq. eq. 7
    # log10(atm)
    return (4.1012 + result1 + result2 + result3 - 0.176055) * (
            (temperature / boiling_point - 1.0) /
            (temperature / boiling_point - 0.125)
            )


def myrdal_and_yalkowsky(compound, temperature, boiling_point):
    comp = groups.composition(compound)
    m = groups.stein_and_brown(compound)

    # NOTE that this method could be used to develop SIMPOL
    # CH2(nr)  group 2
    # CH(nr)   group 4
    # >C<(nr)  group 6
    C = m[2] + m[4] + m[6]

    # Count up the (non-ring) oxygen atoms in non-terminal positions (O in
    # ethers, hydroperoxides, esters, nitrate, peroxyacid, and PAN which
    # contributes 2)
    #
    # ether-O         group 23
    # hydroperoxides  group 25
    # ester-O         group 29
    # nitrate         group 101
    # peroxyacid      group 105
    # PAN             group 106
    O = m[23] + m[25] + m[29] + m[101] + m[105] + (2 * m[106])

    # Count up the (non-ring) nitrogens atoms in non-terminal positions (N in
    # secondary and tertiary amines and amides)
    #
    # amine (sec)     group 39
    # amine (tert)    group 41
    # amide (sec)     group 33
    # amide (tert)    group 35
    N = m[39] + m[41] + m[33] + m[35]

    # Count up the SP2 C and N atoms (=CH-, =C<, Aldehyde, Ketone, Ester,
    # Acid, Nitro group, Nitrate, Peroxyacid and PAN which contributes 2)
    #
    # =CH-             group 9
    # =C<              group 11
    # aldehyde         group 26
    # ketone           group 27
    # amide (prim)     group 32
    # amide (sec)      group 33
    # amide (tert)     group 35
    # ester-O          group 29
    # carboxylic acid  group 31
    # nitro            group 54
    # nitrate          group 101
    # peroxyacid       group 105
    # PAN              group 106
    sp2 = (
            m[9] + m[11] + m[26] + m[27] + m[32] + m[33] +
            m[35] + m[29] + m[31] + m[54] + m[101] + m[105] +
            (2 * m[106]))

    sp3 = (C + O + N)

    # Calculate tau using the number of rings provided by the composition;
    # ensure the value is at least 0.0
    tau = max(0.0, sp3 + (0.5 * sp2) + (0.5 * comp['rings']) - 1.0)

    # To calculate HBN we need to sum groups containing -OH groups, e.g.
    # alcohols (4 types), phenols, hydroperoxide, acid or peroxyacid. Also
    # requires molecular weight which can be obtained from the composition.
    # Summing hydroxyl and acid content (H-bonding groups)
    #
    # alcohols         groups 19-21
    # phenols          group 22
    # vinyl alcohols   group 102
    # hydroperoxide    group 25
    # carboxylic acid  group 31
    # peroxyacid       group 105
    hydroxyl_content = (
            m[19] + m[20] + m[21] + m[22] + m[102] +
            m[25] + m[31] + m[105])

    # Calculate the number of amines
    #
    # aliphatic amines    group 37
    # aromatic amines     group 38
    primary_amines = m[37] + m[38]

    hbn = ((hydroxyl_content ** 0.5) + (0.33 * primary_amines ** 0.5)) / comp['mass']

    return (
        -(86.0 + (0.4 * tau) + (1421 * hbn)) * (boiling_point - temperature) /
        (19.1 * temperature)
        ) + (
        ((-90.0 - 2.1 * tau) / 19.1) *
        (((boiling_point - temperature) / temperature) - log(boiling_point / temperature))
        )


def evaporation(compound, temperature):
    m = groups.evaporation(compound)

    # Calculate the sum of Carbonyl-Like groups, and Hydrogen Bonding groups
    CL_groups = m['5'] + m['6'] + m['7']
    HB_groups = m['8'] + m['9'] + m['10'] + m['11']
    HB_CL_groups = CL_groups + HB_groups

    if CL_groups > 1 or HB_groups > 1:
        # A correction is required if the compound is a diacid with at least
        # one extra HB or CL group
        if m['9'] > 1 and HB_CL_groups > 2:
            diacid_factor = 2.6 / HB_CL_groups
        else:
            diacid_factor = 1.0
        # For molecules with 2 or more CL/HB groups (but not 1xCL + 1xHB), a
        # more complex calculation for A is required. First, split the
        # contribution to the A parameter into lin, CL, and HB. Note that group
        # 12 can contribute to all three subsets.
        A_lin = groups.aggregate_matches(
            m, data.EVAPORATION_A,
            groups=('1', '2', '3', '4', '12lin')
            )
        # For the CL and HB contributions (if non-zero), the contribution is
        # reduced by dividing through by the square root of the number of
        # groups
        try:
            A_CL = groups.aggregate_matches(
                    m, data.EVAPORATION_A,
                    groups=('5', '6', '7', '13', '16', '17', '18', '12CL')
                    ) / sqrt(CL_groups)
        except ZeroDivisionError:
            A_CL = 0.0
        try:
            A_HB = groups.aggregate_matches(
                    m, data.EVAPORATION_A,
                    groups=('8', '9', '10', '11', '14', '15', '19', '20', '12HB')
                    ) / sqrt(HB_groups)
        except ZeroDivisionError:
            A_HB = 0.0

        A = A_lin + (A_CL + A_HB) * diacid_factor

        # The calculation of the B parameter also includes the diacid
        # correction for the terms associated with CL and HB groups. Note that
        # as this calculation stands it would give the wrong VP for a diacid
        # with a further CL or HB group and then a nitrate group on a ring.
        # This is because the contribution of group 12 to the B parameter is
        # assumed to be due to HB or CL groups.
        B_lin = groups.aggregate_matches(
                m, data.EVAPORATION_B,
                groups=('1', '2', '3', '4')
                )
        B_CL_HB = groups.aggregate_matches(
                m, data.EVAPORATION_B,
                groups=(
                    '5', '6', '7', '8', '9', '10', '11', '12',
                    '13', '14', '15', '16', '17', '18', '19', '20')
                )
        B = B_lin + B_CL_HB * diacid_factor

    else:
        # A simpler calculation can be used for monofunctionals and
        # hydrocarbons
        g = set(m.keys()) - {'12lin', '12CL', '12HB'}
        A = groups.aggregate_matches(
                m, data.EVAPORATION_A, groups=g)
        B = groups.aggregate_matches(
                m, data.EVAPORATION_B, groups=g)

    return A + (B / temperature ** 1.5)

