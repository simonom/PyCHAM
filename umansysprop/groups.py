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
try:
    range = xrange
except NameError:
    pass


from math import ceil
from itertools import chain, groupby
from operator import itemgetter

from . import data

def matches(patterns, compound):
    """
    Returns a mapping of group identifier to number of matches.

    The *patterns* parameter specifies a dictionary mapping a group number to
    OpenBabel Smarts objects.

    The *compound* parameter specifies an OpenBabel Molecule object (presumably
    derived from a SMILES string). The function calculates the number of
    matches of the Molecule with each Smarts string, returning the result as a
    mapping.

    :param patterns:
        The mapping of group numbers to Smarts objects

    :param compound:
        A Molecule objects to match against the SMARTS strings in *filename*

    :returns:
        A mapping of group number to match counts
    """
    return {
        group: len(smarts.findall(compound))
        for group, smarts in patterns.items()
        }


def all_matches(patterns, compounds, bias=True):
    """
    Returns a mapping of group identifier to a biased sum of matches.

    The *patterns* parameter specifies a dictionary mapping a group number to
    OpenBabel Smarts objects.

    The *compounds* parameter specifies a sequence of OpenBabel Molecule
    objects (presumably derived from SMILES strings) or, if *bias* is True
    (which it is by default), a mapping of OpenBabel Molecule objects to
    relative abundances.

    If *bias* is False, the function calculates the number of matches of each
    Molecule with each SMARTS string, returning a mapping of the sum of all
    matches with each SMARTS string. If *bias* is True, the function calculates
    the same sum but each match count will be multiplied by the value
    associated with the compound in the *compounds* mapping.

    :param patterns:
        The mapping of group numbers to Smarts objects

    :param compounds:
        A sequence of OpenBabel Molecule objects

    :param bias:
        A bool indicating whether or not *compounds* is a mapping (True) and
        that the value associated with each Molecule should be used as a bias
        in the resulting sum

    :returns:
        A mapping of group number to match counts
    """
    return {
        key: sum(count for (group, count) in items)
        for key, items in groupby(sorted(
            (group, count * (compounds[compound] if bias else 1))
            for compound in compounds
            for group, count in matches(patterns, compound).items()
            ), key=itemgetter(0))
        }


def aggregate_matches(matches, coefficients, groups=None):
    """
    Calculates the biased sum of a set of matches.

    The *matches* parameter specifies a mapping of groups to values to be
    summed.  The *coefficients* parameter specifies mapping of groups to
    coefficient values.

    The optional *groups* parameter specifies the subset of the keys of
    *matches* which are to be summed. If this is not specified, it defaults to
    all keys of *matches*.

    For each group in *groups*, the corresponding values from *matches* and
    *coefficients* will be multiplied. The result is the sum of all these
    multiplications.

    :param matches:
        A mapping of groups to values to be summed

    :param coefficients:
        A mapping of groups to coefficients

    :param groups:
        An optional subset of the keys of *matches* indicating the values to
        be summed

    :returns:
        A floating point value which is the sum of the multiplication of each
        value in *matches* against the corresponding value in *coefficients*
    """
    if groups is None:
        groups = matches.keys()
    return sum(
        matches[group] * coefficients[group]
        for group in groups
        )


def aggregate_interactions(compound, interactions, coefficients):
    """
    Calculates the biased sum of a matrix of interactions.

    The *compound* parameter specifies the compound that interactions are to be
    calculated for. The *interactions* parameter specifies a mapping of groups
    to values to be summed. The *coefficients* parameter specifies a mapping of
    groups to a mapping of groups to coefficient values, representing a
    2-dimensional matrix of coefficients.

    For example, if *interactions* is a mapping with the keys ``{A, B, C}``,
    then *coefficients* must have the structure::

        coefficients = {
            'A': {'A': 1.0, 'B': 1.0, 'C': 1.0},
            'B': {'A': 1.0, 'B': 1.0, 'C': 1.0},
            'C': {'A': 1.0, 'B': 1.0, 'C': 1.0},
            }

    The result is the sum of the multiplication of each group within
    *interactions* with every group within *interactions*, biased by
    *coefficients*.

    :param compound:
        The compound to calculate interactions for

    :param interactions:
        A mapping of groups to counts (e.g. from :func:`nannoolal_interactions`)

    :param coefficients:
        A mapping of groups to a mapping of groups to coefficients

    :returns:
        The biased sum of the multiplication of each group against each group
        in *interactions*
    """
    comp = composition(compound)
    # Calculate the cardinality of each interaction group
    n_interactions = sum(interactions.values())
    if n_interactions <= 1:
        return 0
    return sum(
        # Calculation for interaction of a group with itself; e.g. for glycerol
        # (3xOH): 3 * 2 * A:A-interaction-coeff / ((m - 1) * n)
        (
            interactions[g1] * (interactions[g2] - 1) *
            coefficients[g1][g2]
            ) / (comp['n-non-H'] * (n_interactions - 1))
        if g1 == g2 else
        # Calculation for interaction of grp1 with grp2:
        # 2 * #grp1 * #grp2 * grp1:grp2-interaction-coeff / ((m - 1) * n)
        # The multiplication by two is unnecessary as we will calculate this
        # for grp1:grp2, and grp2:grp1
        (
            interactions[g1] * interactions[g2] *
            coefficients[g1][g2]
            ) / (comp['n-non-H'] * (n_interactions - 1))
        for g1 in interactions
        for g2 in interactions
        )


def composition(compound):
    m = matches(data.COMPOSITION_SMARTS, compound)

    result = {}
    result['C'] = m[400]
    result['O'] = m[401]
    result['N'] = m[402]
    result['F'] = m[403]
    result['Cl'] = m[404]
    result['Br'] = m[405]
    result['I'] = m[406]
    result['H'] = sum((
        4 * m[419],            # CH4 groups
        3 * m[420],            # CH3 groups
        2 * (m[421] + m[422]), # CH2 groups
        1 * (m[423] + m[424] + m[425] + m[429]), # CH groups
        1 * m[426],            # OH groups
        2 * m[427],            # NH2 groups
        1 * m[428],            # NH groups
        ))
    result['n-non-H'] = sum((
        result['C'],
        result['O'],
        result['N'],
        result['F'],
        result['Cl'],
        result['Br'],
        result['I'],
        ))
    result['rings'] = sum((
        m[440] / 3, # 3-member rings
        m[441] / 4, # 4-member rings
        m[442] / 5, # 5-member rings
        m[443] / 6, # 6-member rings
        m[444] / 7, # 7-member rings
        m[445] / 8, # 8-member rings
        ))
    result['mass'] = sum((
         12.011 * result['C'],
          1.008 * result['H'],
         15.999 * result['O'],
         14.007 * result['N'],
         18.998 * result['F'],
         35.45  * result['Cl'],
         79.904 * result['Br'],
        126.904 * result['I'],
        ))
    return result


def stein_and_brown(compound):
    m = matches(data.STEIN_AND_BROWN_SMARTS, compound)
    # Copy groups 1-106
    result = {
        group: count
        for (group, count) in m.items()
        if 1 <= group <= 106
        }
    # Sum groups starting with CH< (non-ring)
    result[3] += m[704]
    # Sum ethers (non-ring)
    result[23] += m[700]
    # Sum ethers (ring)
    result[24] += m[306] + m[305]
    # Sum aldehydes
    result[26] += m[310] + m[702] + m[705] + m[706]
    # Sum ketones
    result[27] += m[304] + m[308] + m[309] + m[703]
    # Sum ketones (ring)
    result[28] += m[311]
    # Sum esters
    result[29] += m[303] + m[305] + (2 * m[307])
    result[29] += m[308] + m[309] + m[310] + (2 * m[705])
    # Sum esters (ring)
    result[30] += m[311] + m[306]
    # Sum carboxylic acids
    result[31] += m[701]
    # Sum nitrate count
    result[101] += m[700]
    # Sum vinyl chloride count
    result[103] += m[304] + (2 * m[703]) + m[706]
    return result


def nannoolal_primary(compound):
    comp = composition(compound)
    m = matches(data.NANNOOLAL_SMARTS_PRIMARY, compound)
    result = {}
    # (CH3) attached to carbon
    result[1] = m[1] + m[7000] + (2 * m[7001])
    # (CH3) attached to electronegative (non aromatic) atom
    result[2] = m[2] + m[7002]
    # (CH3) attached to an aromatic atom
    result[3] = m[3]
    # (CH2) attached to a (C) in a chain
    result[4] = m[4]
    # (CH) attached to a (C) in a chain
    result[5] = m[5]
    # (>C<) attached to a (C) in a chain
    result[6] = m[6]
    # (CHx) in a chain attached to one or more electronegative atoms
    result[7] = (
            m[7] + m[7007] + m[701] + m[702] + m[703] + m[704] +
            m[705] + m[706] + m[707] + m[708])
    # (CHx) attached to an aromatic atom
    result[8] = m[8] + m[801] + m[802]
    # (CH2) in a ring
    result[9] = m[9]
    # (>CH) in a ring
    result[10] = m[10]
    # (>C<) in a ring
    result[11] = m[11]
    # (CHx) in a ring attached to at least one electronegative atom that is
    # not part of a ring
    result[12] = m[12] + m[1201] + m[1202]
    # (CHx) in a ring attached to at least one N or O that is part of a ring
    result[13] = (
            m[13] + m[1301] + m[1302] + m[1303] + m[1304] +
            m[1305] + m[1306] + m[1307])
    # (CHx) in a ring attached to an aromatic atom
    result[14] = m[14] + m[1401] + m[1402]
    # Aromatic C-H
    result[15] = m[15]
    # C in an aromatic ring attached to aliphatic C
    result[16] = m[16]
    # C in an aromatic ring attached to electronegative atom
    result[17] = m[17]
    # 18: Fused aromatic ring compounds; not relevant to MCM
    # 19..23: Fluorine containing groups; not relevant to MCM
    # F attached to an aromatic ring
    result[24] = m[24]
    # A single Cl attached to non-aromatic C
    result[25] = m[25] + m[2501] + m[2502] + m[2503]
    # Pairs of Cl attached to C
    result[26] = 2 * (m[26] + m[2601] + m[2602])
    # Triplets of Cl attached to C
    result[27] = 3 * (m[27] + m[2701] + m[2702])
    # Cl attached to aromatic C
    result[28] = m[28]
    # Cl attached to double bonded C
    result[29] = m[29] + (2 * m[2901]) + m[2902] + m[7005]
    # Aliphatic bromides
    result[30] = m[30]
    # Aromatic bromides
    result[31] = m[31]
    # Iodides
    result[32] = m[32]
    # Tertiary alcohols
    result[33] = m[33]
    # Secondary alcohols
    result[34] = m[34] + m[3401] + (2 * m[7007])
    # Long chain (35) and short chain (36) alcohols. In Sep 2010, it was found
    # that 20 C5 primary alcohols were being wrongly assigned to short chain
    # (36) hence the composition check in the filter below...
    is_short_chain = (
            (m[36] == m[3601]) and
            (m[36] == m[3602]) and
            (m[36] == m[3603]) and
            (m[36] == m[3604])
            ) and (comp['C'] <= 4)
    result[35] = m[36] if not is_short_chain else 0
    result[36] = m[36] if is_short_chain else 0
    # Aromatic alcohols or phenols
    result[37] = m[37]
    # Ethers (group 38 hits all ethers, while 65 hits aromatic ethers of the
    # type found in Furans).
    result[38] = m[38] + m[7002] + m[7009] + m[7010] - m[65]
    # Epoxides
    result[39] = m[39] + m[3901] + m[3902] + m[3903] + m[3904] + m[3905]
    # Aliphatic primary amines
    result[40] = m[40]
    # Aromatic primary amines
    result[41] = m[41]
    # Secondary amines: aliphatic or aromatic
    result[42] = m[42]
    # Tertiary amines: aliphatic or aromatic
    result[43] = m[43]
    # Carboxylic acids
    result[44] = m[44] + m[7003] + m[7009]
    # Ester group in a chain
    result[45] = m[45]
    # Formate group
    result[46] = m[46]
    # Ester within a ring-lactone
    result[47] = m[47]
    # Tertiary amide C(=O)N<
    result[48] = m[48]
    # Secondary amide C(=O)NH-
    result[49] = m[49]
    # Primary amide C(=O)NH2
    result[50] = m[50]
    # Aliphatic ketone
    result[51] = m[51]
    # Aliphatic aldehyde
    result[52] = m[52] + m[304] + (2 * m[7008]) + m[7004]
    # 53..56: Sulphur groups; not relevant to the MCM
    # Nitrile group
    result[57] = m[57]
    # >C=C< in a chain with each C having at least one non-H neighbour
    result[58] = m[58] + m[5801] + m[5802]
    # ???
    result[59] = m[59] + m[5901] + m[5902] + m[5903]
    # >C=C< in a chain attached to at least one electronegative atom
    result[60] = (
            m[60] + m[6001] + m[6002] + m[6003] + m[6004] + m[6005] +
            m[6006] + m[6007] + m[6008] + m[6009])
    # CH2=C< in a chain
    result[61] = m[61] + m[6101] + m[6102]
    # >C=C< in a ring
    # NOTE: in the case of a hit on group 96 for a molecule such as maleic
    # anhydride, where group 96 covers the whole molecule the >C=C< in a
    # ring will be counted twice, hence the correction below
    result[62] = (
            m[62] + m[6201] + m[6202] + m[6203] + m[6204] +
            m[6205] + m[7006] - m[96])
    # -C#C- in a chain
    result[63] = m[63]
    # HC#C- at end of chain
    result[64] = m[64] + m[6401]
    # 65: Aromatic ethers; used in result[38] above
    # 66..67: Aromatic nitrogens in piperidine and pyridine; not relevant to MCM
    # Aliphatic nitro group
    result[68] = m[68]
    # Aromatic nitro group
    result[69] = m[69]
    # 70..71: Silicon containing groups; not relevant to MCM
    # Nitrate group
    result[72] = m[72] + m[7002]
    # 73: Phosphorus; not relevant to MCM
    # 74..75: Nitrogen (nitrites and oximes); not relevant to MCM
    # Acid anhydride
    # NOTE: Along with mixed anhydrides of formic acid plus a second acid
    # (7601), and the anhydrides of formic acid (7602). The groups will also
    # hit group 96 (cyclic anhydrides attached to double bonded Cs)
    result[76] = m[76] + m[7601] + m[7602] - m[96]
    # Acid chloride
    result[77] = m[77] + m[7701] + m[7005]
    # 78: Boron; not relevant to MCM
    # Carbonate in a chain
    result[79] = m[79]
    # 80: Isocyanate; not relevant to MCM
    # 81..82: Sulphur groups; not relevant to MCM
    # 83: Tin; not relevant to MCM
    # 84: Arsenic; not relevant to MCM
    # 85..86: Germanium; not relevant to MCM
    # 87: Cumulative double bonds (e.g. 1,2-Butadiene); not found in MCM
    # >C=C-C=C< in a ring
    result[88] = (
            m[88] + m[7006] + m[8801] + m[8802] + m[8803] + m[8804] +
            m[8805] + m[8806] + m[8807] + m[8808] + m[8809])
    # >C=C-C=C< in a chain
    result[89] = (
            m[89] + m[8901] + m[8902] + m[8903] + m[8904] + m[8905] +
            m[8906] + m[8907] + m[8908] + m[8909] + m[8910] + m[8911] +
            m[8912] + m[8913] + m[8914] + m[8915] + m[8916] + m[8917] +
            m[8918] + m[8919])
    # Aromatic aldehyde
    result[90] = m[90]
    # Nitrogen group; not relevant to MCM
    #result[91] = m[91]
    # Aromatic ketone
    result[92] = m[92]
    # 93: Silicon group; not relevant to MCM
    # Bridging peroxide
    result[94] = m[94]
    # 95: Conjugated triple bonds; not found in the MCM
    # 96: Cyclic anhydrides; used in result[76] above
    # Secondary amine in a ring
    result[97] = m[97]
    # 98..99: non-existent
    # 100..101: Nitrogen groups; not relveant to MCM
    # 102: Fluorine group; not relevant to MCM
    # Cyclic carbonate
    result[103] = m[103] + m[304] + m[7008]
    # 104..105: Sulphur groups; not relevant to MCM
    # 106: Nitrogen group; not relevant to MCM
    # 107: Sulphur group; not relevant to MCM
    # 108: Nitrogen/sulphur group; not relevant to MCM
    # 109: Nitrogen/oxygen group; not relevant to MCM
    # 110: non-existent
    # 111: Nitrogen group; not relevant to MCM
    # 112: non-existent
    # 113: Phosphorus group; not relevant to MCM
    # 114: non-existent
    # 115: Nitrogen group; not relevant to MCM
    # 116: non-existent
    # 117: Aluminum group; not relevant to MCM
    # 118: 1,2-diketone structure but no parameter is available for this group
    # 214: bridging bond between two aromatic rings; not found in the MCM
    # 215..216: Silicon compounds; not found in the MCM
    # Hydroperoxide group
    result[301] = m[301]
    # Peroxyacids
    result[302] = m[302]
    # PAN
    result[303] = m[303]
    return result


def nannoolal_secondary(compound):
    comp = composition(compound)
    m = matches(data.NANNOOLAL_SMARTS_SECONDARY, compound)
    result = {}
    # C=C-C=O. In Sep 2010, it was found that 12 quinone structures were double
    # counting the carbonyles. For this contribution a carbonyl alpha to one db
    # is the same as a carbonyl alpha to two db. Hence cases where the same
    # carbonyl was being double counted, C=C(C=O)C=C, needed to be subtracted
    result[118] = m[118] - m[1181]
    # C=O connected to C with two or more halogens
    result[119] = m[119]
    # C=O connected to two Cs with two or more halogens each
    result[120] = m[120]
    # Carbon with three halogens
    result[121] = m[121]
    # Secondary carbon with two halogens
    result[122] = m[122]
    # Correction for a Hydrogen free compound
    result[123] = int(comp['H'] == 0)
    # Correction for a compound with a single H
    result[124] = int(comp['H'] == 1)
    # ???
    result[125] = (m[125] // 3) + (m[1251] // 4)
    # ???
    result[126] = m[126] // 5
    # 127..129: identify ortho-, meta-, and para- structures respectively. Each
    # can only be counted once, and only if it is the only structure available.
    # However, the example of phloroglucinol highlighted that 1,3,5
    # tri-X-benzene obeys these rules and should be assigned one
    # meta-correction.
    result[127] = int(
            (m[127] == 1) and
            (m[128] == 0) and
            (m[129] == 0)
            )
    result[128] = int(
            (m[127] == 0) and
            ((m[128] == 1) or (m[128] == 3)) and
            (m[129] == 0)
            )
    result[129] = int(
            (m[127] == 0) and
            (m[128] == 0) and
            (m[129] == 2)
            )
    # ??? - need to add these, will consult with Mark Barley
    result[130] = m[130]
    # ???
    result[131] = m[131] + m[1311] + m[1312] + m[1313]
    # ???
    result[132] = m[132] + m[1321]
    # ???
    result[133] = m[133]
    return result


def nannoolal_interactions(compound):
    m = nannoolal_primary(compound)
    # Construct 12 groups for calculate interactions. Hydroperoxides are
    # assumed to be as effective as alcohols in group interactions and
    # similarly peroxyacids are treated the same as carboxylic acids
    result = {}
    # Alcohols (33-36) plus hydroperoxide (301)
    result['A'] = m[33] + m[34] + m[35] + m[36] + m[301]
    # Phenols
    result['B'] = m[37]
    # Carboxylic acids (44) plus peroxyacids (302)
    result['C'] = m[44] + m[302]
    # Ether
    result['D'] = m[38]
    # Epoxide
    result['E'] = m[39]
    # Esters (45-47) plus PAN (303)
    result['F'] = m[45] + m[46] + m[47] + m[303]
    # Ketones
    result['G'] = m[51] + m[92]
    # Aldehydes
    result['H'] = m[52] + m[90]
    # Primary amines
    result['M'] = m[40] + m[41]
    # Secondary amines
    result['N'] = m[42] + m[97]
    # Nitrile
    result['P'] = m[57]
    # Aromatic nitro
    result['Q'] = m[69]
    return result


def evaporation(compound):
    m = matches(data.EVAPORATION_SMARTS, compound)
    result = {}
    result['1'] = 1
    # Number of carbons (2) and "in-line" oxygen atoms in the compound. In-line
    # oxygens belong to groups such as ethers, esters, and peroxides. For use
    # with Gerard's compounds and other MCM derived datasets:- 5
    # (Methylperoxynitrate) is mis-assigned. This is probably best treated as
    # a Me group, an ether group, and a nitrate group.
    #
    # Ethers               group 201
    # Esters               group 202
    # Peroxides            group 203
    # Methylperoxynitrate  group 7002
    result['2'] = m[2] + m[201] + m[202] + (2 * m[203]) + m[7002]
    # Topological index "t" = number of branching C's
    #
    # CX4 with 4 C-C bonds                      group 301
    # CX4 with 3 C-C bonds                      group 302
    # CX4 with 3 C-C bonds + 1 C-O bond-not OH  group 303
    # CX4 with 2 C-C bonds + 1 C-O bond-not OH  group 304
    # CX4 with 2 C-C bonds + 2 C-O bond-not OH  group 305
    # CX4 with 1 C-C bonds + 2 C-O bond-not OH  group 306
    #
    # Also need to determine the number of rings. If this is an integer, then
    # unfused rings. If non-integer, then fused rings, the most common being
    # 2 fused rings
    #
    # 3-membered rings  group 307
    # 4-membered rings  group 308
    # 5-membered rings  group 309
    # 6-membered rings  group 310
    # 7-membered rings  group 311
    # 8-membered rings  group 312
    result['3'] = (
            2 * (m[301] + m[303] + m[305]) +
            (m[302] + m[304] + m[306]) -
            (
                m[307] / 3 + m[308] / 4 + m[309] / 5 +
                m[310] / 6 + m[311] / 7 + m[312] / 8)
            )
    # Nitrate groups. See notes above (in-line oxygen) regarding
    # Methylperoxynitrate.
    result['4'] = m[4] + m[7002]
    # Number of aldehyde and ketone groups.  For use with Gerard's compounds
    # and other MCM derived datasets:- 1288 (Formaldehyde) is mis-assigned
    #
    # Ketones in a chain     group 501
    # Ketones in a ring      group 502
    # Aldehydes              group 503
    # Formaldehyde           group 7001
    result['5'] = m[501] + m[502] + m[503] + m[7001]
    # Number of ester groups
    #
    # Esters in a chain     group 601
    # Esters in a ring      group 602
    result['6'] = m[601] + m[602]
    # Number of PAN groups
    result['7'] = m[7]
    # Alcohols
    #
    # Primary alcohols    group 801
    # Secondary alcohols  group 802
    # Tertiary alcohols   group 803
    # Vinyl alcohols      group 804
    result['8'] = m[801] + m[802] + m[803] + m[804]
    # Acid groups
    result['9'] = m[9]
    # Hydroperoxide groups
    result['10'] = m[10]
    # Peroxyacid groups
    result['11'] = m[11]
    # Correction for groups within (or attached to) a ring
    #
    # Ketone in a ring                              group 502
    # Esters in a ring                              group 602
    # Ethers in a ring                              group 1201
    # Peroxide in a ring                            group 1202
    # Alcohol attached to a carbon in a ring        group 1203
    # Nitrate attached to a carbon in a ring        group 1204
    # Hydroperoxide attached to a carbon in a ring  group 1205
    result['12'] = (
            m[502] + m[602] + m[1201] + m[1202] +
            m[1203] + m[1204] + m[1205])
    # Ring corrections also need to be split between Lin, Carbonyl-Like, and
    # Hydrogen Bonding types of functional group
    result['12lin'] = m[1201] + m[1202] + m[1204]
    result['12CL'] = m[1202] + m[502]
    result['12HB'] = m[1203] + m[1205]
    # Aldehyde or ketone conjugated to C=C
    result['13'] = m[13]
    # Correction for the type of alcohol. Primary and vinyl alcohols don't have
    # a correction. Just need to identify secondary (group 802) and tertiary
    # (group 803)
    result['14'] = m[802] + (2 * m[803])
    # Alcohol beta to a double bond (HO-CX4-C=c-). Limit this to 0 or 1
    result['15'] = min(1, m[15])
    # Ketone/aldehyde next to (i.e. O=C-C=O) a carbonyl containing group (ester
    # or PAN as well as ketone/aldehyde)
    #
    # Ketone/aldehyde next to an aldehyde  group 1601
    # Ketone/aldehyde next to a ketone     group 1602
    # Ketone/aldehyde next to an ester     group 1603
    # Ketone/aldehyde next to PAN          group 1604
    result['16'] = m[1601] + m[1602] + m[1603] + m[1604]
    # Ketone/aldehyde beta to (i.e. O=C-CX4-C=O) a carbonyl containing group
    # (ester or PAN as well as ketone/aldehyde)
    #
    # Ketone/aldehyde beta to an aldehyde  group 1701
    # Ketone/aldehyde beta to a ketone     group 1702
    # Ketone/aldehyde beta to an ester     group 1703
    # Ketone/aldehyde beta to a PAN        group 1704
    result['17'] = m[1701] + m[1702] + m[1703] + m[1704]
    # Identify features in which an aldehyde or ketone group is next to a C
    # atom bearing a non-CL group (alcohol, nitrate, hydroperoxide)
    #
    # Aldehyde or ketone next to a nitrate or hydroperoxide   group 1801
    # Aldehyde or ketone next to an alcohol                   group 1802
    result['18'] = m[1801] + m[1802]
    # Identify all features in which an alcohol is next to a functional
    # group (or a C atom bearing a functional group)
    #
    # Alcohol next to hydroperoxide or nitrate (e.g. OCCOO)  group 1901
    # Alcohol next to an alcohol (e.g. OCCO)                 group 1902
    # Alcohol next to an aldehyde (e.g. OCC=O)               group 1903
    # Alcohol next to a ketone (e.g. OCC(=O)CC)              group 1904
    # Alcohol next to an acid (e.g. OCC(=O)O)                group 1905
    # Alcohol next to an ester (e.g. OCC(=O)OCC)             group 1906
    # Alcohol next to a peroxyacid                           group 1907
    # Alcohol next to a PAN                                  group 1908
    result['19'] = (
            m[1901] + (2 * m[1902]) + m[1903] + m[1904] +
            m[1905] + m[1906] + m[1907] + m[1908])
    # Identify all features in which an acid group is next to a C atom bearing
    # a CL group (aldehyde, ketone, ester, and PAN)
    #
    # Acid next to an aldehyde     group 2001
    # Acid next to a ketone        group 2002
    # Acid next to an ester        group 2003
    # Acid next to PAN             group 2004
    result['20'] = m[2001] + m[2002] + m[2003] + m[2004]
    return result


def girolami(compound):
    m = matches(data.GIROLAMI_SMARTS, compound)
    result = {}
    result['short_period_1'] = m[2]
    result['short_period_2'] = m[3]
    result['long_period_1'] = m[4]
    result['long_period_2'] = m[5]
    result['long_period_3'] = m[6]
    result['hydroxyl_groups'] = m[7]
    result['hydroperoxide_groups'] = m[701]
    result['carboxylic_acid_groups'] = m[8]
    result['peroxyacid_groups'] = m[801]
    result['primary_secondary_amino_groups'] = m[9]
    result['amide_groups'] = m[10]
    all_rings = sum((
        m[1101] / 3,
        m[1102] / 4,
        m[1103] / 5,
        m[1104] / 6,
        m[1105] / 7,
        m[1106] / 8,
        ))
    all_rings = 0 if all_rings < 1 else ceil(all_rings)
    result['fused_rings'] = m[1201] + m[1202] + m[1203]
    result['unfused_rings'] = all_rings - result['fused_rings']
    return result


def schroeder(compound):
    comp = composition(compound)
    m = matches(data.SCHROEDER_SMARTS, compound)
    result = {}
    result[51] = comp['C']
    result[52] = comp['H']
    result[53] = comp['O']
    result[54] = comp['N']
    result[55] = comp['F']
    result[56] = comp['Cl']
    result[57] = comp['Br']
    result[58] = comp['I']
    result[59] = min(1.0, comp['rings'])
    result[60] = m[60] + (m[601] / 2)
    result[61] = m[61]
    return result


def le_bas(compound):
    comp = composition(compound)
    m = matches(data.LE_BAS_SMARTS, compound)
    result = {}
    result[71] = comp['C']
    result[72] = comp['H']
    result[73] = comp['O']
    result[75] = comp['F']
    result[76] = comp['Cl']
    result[77] = comp['Br']
    result[78] = comp['I']
    # Calculate the (integer) number of rings
    # 3-membered rings
    result[79] = ceil(m[79] / 3)
    # 4-membered rings
    result[80] = ceil(m[80] / 4)
    # 5-membered rings
    result[81] = ceil(m[81] / 5)
    # 6-membered rings
    result[82] = ceil(m[82] / 6)
    # Methyl esters
    result[731] = m[731] * 2
    # Ethyl esters
    result[732] = m[732] * 2
    # Propyl and higher esters are given by difference
    result[733] = (m[730] - m[731] - m[732]) * 2
    # Methyl ethers + dimethyl ether
    result[734] = m[734] + m[7341]
    # Ethyl ethers; this needs to deal with multiple hits on branched ethyl
    # ethers
    result[735] = (
        # Diethyl ether
        m[7351] +
        # Any ethyl ether with [CX4;H0,H1] on the other side of the oxygen
        # must be branched, or have heteroatom functionality (e.g. -OH or
        # -OOH)
        m[7352] +
        # The other side of the ether oxygen from the ethyl group can also
        # be an aromatic ring (phenyl ether) or a double bonded C (as in a
        # vinyl ether but can't be a carboxyl)
        m[7353] +
        # Finally, the first C may be a CH2 group which means that the
        # second carbon atom can be anything except a methyl group
        m[7354]
        )
    # Propyl and higher ethers are given by difference
    result[736] = m[733] - result[734] - result[735]
    # Acids (contribute two oxygens to the correction)
    result[737] = m[737] * 2
    # Need to count the number of Oxygens attached to Nitrogen; first with a
    # single bond - structures such as N-O-N are assumed to be too unstable to
    # exist in closed shell species
    result[738] = (
        m[7381] +
        # And then the more common double bond
        m[7382]
        )
    # Carbonates are esters which need to be corrected separately. Use the
    # intermediate value (ethyl ester) as a compromise. Carbonates detected
    # using the SMARTS for Extended Nannoolal
    result[739] = (
        # First cyclic carbonates
        m[739] +
        # Then the three types of non-cyclic carbonates
        m[7391] +
        m[7392] +
        m[7393]
        ) * 3
    # Nitrogen compounds containing the double bonded nitrogen need to be
    # identified along with nitrogen in amines. Identify double bonded nitrogen
    # (nitro, nitrate, PAN, etc. groups):- Nitrogen double bonded to another
    # atom. Also triple bonded N is assumed to have the same increment as
    # double bonded N.
    result[741] = (
        # Nitrogens double bonded to another atom
        m[741] +
        # Nitrogen triple bonded to carbon (nitriles)
        m[7411] -
        # Nitrogen double bonded to two other atoms (corrects double counting
        # of nitrogens m[741])
        m[7412]
        )
    # Nitrogen in primary amines/amides = NH2 groups
    result[742] = m[742]
    # Nitrogen in secondary amines/amides = NH groups
    result[743] = m[743]
    # Total number of nitrogen atoms (including nitrogen in tertiary
    # amines/amides = >N-groups; m[744])
    assert comp['N'] == result[741] + result[742] + result[743] + m[744]
    return result


def unifac(compound):
    comp = composition(compound)
    m = matches(data.AIOMFAC_UNIFAC_SMARTS, compound)
    result = {}
    # Ni for UNIFAC group 5 (CH2=CH-)
    result[5] = m[5]
    # Ni for UNIFAC group 6 (-CH=CH-)
    result[6] = m[6]
    # Ni for UNIFAC group 7 (CH2=C<)
    result[7] = m[7]
    # Ni for UNIFAC group 8 (CH=C<)
    result[8] = m[8]
    # Ni for UNIFAC group 70 (>C=C<)
    result[70] = m[70]
    # Ni for UNIFAC group 9 (C(Ar)-H)
    result[9]  = m[9]
    # Ni for UNIFAC group 10 (C(Ar)-)
    result[10] = m[10]
    # Ni for UNIFAC group 11 (C(Ar)-CH3)
    result[11] = m[11]
    # Ni for UNIFAC group 12 (C(Ar)-CH2-)
    result[12] = m[12]
    # Ni for UNIFAC group 13 (C(Ar)-CH<) and UNIFAC group (C(Ar)-C<). These
    # will be added into group 13 as T-butyl groups attached to aromatic C are
    # not defined in UNIFAC.
    result[13] = m[13] + m[1301]
    # Ni for UNIFAC group 14 (-OH) and ether plus alcohol (301)
    result[14] = m[14] + m[301]
    # Water
    result[16] = m[16]
    # Ni for UNIFAC group 17 (phenols)
    result[17] = m[17]
    # Identify compounds with a ketone structure (group 19), acid anyhydrides
    # split into an ester and a ketone (308); also identify ketone groups with
    # a CH3 attached (group 1901), with CH2 attached (group 1902), and two
    # carbonyl groups with a methylene group in between (group 1903). A CH3-C=O
    # structure is assigned to result[18].  A CH2-C=O is assigned to
    # result[19]. All other ketones are assigned to result[191], which also
    # includes a contribution from acid chlorides (group 5301).
    #
    # MCM_2171: [CX3;H0](Cl)(Cl)(=O) (Phosgene) needs to be assigned to 1 "bare
    # carbonyl" and 2 "bare chlorines" (703)
    result[18] = m[1901]
    result[19] = m[1902] + m[1903]
    result[191] = m[19] - m[1901] - m[1902] - m[1903] + m[5301] + m[308] + m[703]
    # Ni for UNIFAC group 20 (aldehydes), plus formyl acid anhydrides split
    # into an ester + aldehyde (309), plus diformyl acid anhydrides split into
    # a formyl-ester and an aldehyde (310)
    #
    # MCM_1288: Formaldehyde needs to be assigned to aldehyde (702)
    #
    # MCM_2737 (or DFC): C(OC(OC=O)=O)=O. This is split up as a formyl ester
    # with a normal ester and an aldehyde group (705).
    result[20] = m[20] + m[309] + m[310] + m[702] + m[705]
    # Identify compounds with an ester structure (group 22; formyl esters are
    # treated separately), ester plus ether plus nitro (303), carbonates split
    # into an ester and ether (306), formyl-carbonates split into a formyl
    # ester and a normal ester (307), acid anhydrides split into an ester and a
    # ketone (308), formyl acid anhydrides split into an ester + aldehyde
    # (309). Also identify ester groups with a CH3 attached (group 2201) and
    # with CH2 group attached (group 2202). A CH3-C(=O)O structure is assigned
    # to result[21]. A CH2-C(=O)O is assigned to result[22]. All other esters
    # are assigned to result[231].
    #
    # MCM_2737 (or DFC): C(OC(OC=O)=O)=O. This is split up as a formyl ester
    # with a normal ester and an aldehyde group (705).
    result[21] = m[2201]
    result[22] = m[2202]
    result[231] = (
            m[22] - m[2201] - m[2202] +
            m[303] + m[306] + m[307] + m[308] + m[309] +
            m[705])
    # Ni for UNIFAC group 23 (formyl ester), formyl-carbonates split into
    # a formyl ester and a normal ester (307), diformyl acid anhydrides
    # split into a formyl-ester and an aldehyde (310).
    #
    # MCM_2737 (or DFC): C(OC(OC=O)=O)=O. This is split up as a formyl ester
    # with a normal ester and an aldehyde group (705).
    result[23] = m[23] + m[307] + m[310] + m[705]
    # Identify compounds with a non-aromatic ether structure (group 25), ether
    # plus alcohol (301), ether plus carboxylic acid (302), ester plus ether
    # plus nitro (303), ether plus nitro (304), bridged peroxy groups split
    # into two ethers (305), carbonates split into an ester and ether (306);
    # also identify ether groups with a CH3 attached (group 2501), ether groups
    # with a CH2 attached to the oxygen but no CH3 groups (group 2502), and
    # with a CH group attached but no CH3 or CH2 groups (group 2503).
    #
    # MCM_0005: [CX4;H3][OX2][OX2][NX3](=O)=O. The CH3 group is correctly
    # assigned but the rest of the molecule needs to be assigned to 2xether +
    # 1xnitro (700).
    result[24] = m[2501]
    result[25] = m[2502]
    result[26] = m[2503]
    result[261] = (
            m[25] - m[2501] - m[2502] - m[2503] +
            m[301] + m[302] + m[303] + m[304] + (2 * m[305]) + m[306] +
            (2 * m[700]))
    # Compounds with an aromatic ether group (as found in THF)
    result[27] = m[27]
    # Amines: primary amines are easy to define. However, an NH2 group attached
    # to a tertiary C or a C=C is not defined within UNIFAC.
    # CH3-NH2- methyl amine
    result[28] = m[28]
    # CH2-NH2- methylene amine
    result[29] = m[29]
    # >CH-NH2- >CH amine
    result[30] = m[30]
    # For secondary and tertiary amines, the situation is more complex because
    # of the danger of multiple hits (on structures such as -CH2-N-CH2). First
    # need to identify how many >NH and >N- groups there are. This is done with
    # groups 3001 (secondary amines) and 3002 (tertiary amines).
    result[31] = min(m[31], m[3001])
    result[32] = min(m[32], m[3001] - result[31])
    result[33] = min(m[33], m[3001] - result[31] - result[32])
    result[34] = min(m[34], m[3002])
    result[35] = min(m[35], m[3002] - result[34])
    # Ni for UNIFAC group 36 (C(Ar)-NH2)
    result[36] = m[36]
    # Ni for UNIFAC group 40 (Acetonitrile)
    result[40] = m[40]
    # Ni for UNIFAC group 41 (Methylene nitrile- CH2-C#N)
    result[41] = m[41]
    # Ni for UNIFAC group 42 (Carboxylic acid), and ether plus carboxylic acid
    # (302)
    result[42] = m[42] + m[302]
    # Ni for UNIFAC group 43 (Formic acid)
    result[43] = m[43]
    # Ni for UNIFAC group 44 (-CH2Cl)
    result[44] = m[44]
    # Ni for UNIFAC group 45 (>CHCl)
    result[45] = m[45]
    # Ni for UNIFAC group 46 (>CCl)
    result[46] = m[46]
    # Ni for UNIFAC group 47 (Methylene chloride -CH2Cl2)
    result[47] = m[47]
    # Ni for UNIFAC group 48 (>CHCl2)
    result[48] = m[48]
    # Ni for UNIFAC group 49 (>CCl2)
    result[49] = m[49]
    # Ni for UNIFAC group 50 (Chloroform -CHCl3)
    result[50] = m[50]
    # Ni for UNIFAC group 51 (-CCl3)
    result[51] = m[51]
    # Ni for UNIFAC group 52 (Carbon Tetrachloride -CCl4)
    result[52] = m[52]
    # Ni for UNIFAC group 53 (C(Ar)-Cl)
    result[53] = m[53]
    # Acid chloride group (C=O)Cl. This group is split between a carbonyl and a
    # new group- a bare Cl atom (result[531]).  The contribution to the lone
    # carbonyl group (result[191] above) includes the contribution from acid
    # chlorides, and the bare Cl atom is assigned to 531.
    #
    # MCM_2171: [CX3;H0](Cl)(Cl)(=O) (Phosgene) needs to be assigned to 1 "bare
    # carbonyl" and 2 "bare chlorines" (703)
    #
    # MCM_2160: Methyl-chloride. The methyl group is assigned correctly but the
    # Cl atom is missed. Assign to a bare chlorine group (706).
    result[531] = m[5301] + (2 * m[703]) + m[706]
    # Vinyl chloride structure (C=C-Cl). Each hit with group 69 registers a
    # single Cl atom
    result[69] = m[69]
    # Identify compounds with a nitro structure (group 5401), ester plus ether
    # plus nitro (303), ether plus nitro (304); also identify nitro groups with
    # a CH3 attached (as in nitromethane:- group 54); CH2-N(=O)2 structure will
    # be assigned to result[55]. A CH-N(=O)2 is assigned to result[56]. All
    # other nitro groups are assigned to result[561].
    #
    # MCM_0005: [CX4;H3][OX2][OX2][NX3](=O)=O. The CH3 group is correctly
    # assigned but the rest of the molecule needs to be assigned to 2xether +
    # 1xnitro (700).
    result[54] = m[54]
    result[55] = m[55]
    result[56] = m[56]
    result[561] = m[5401] - m[54] - m[55] - m[56] + m[303] + m[304] + m[700]
    # Ni for UNIFAC group 57 (aromatic nitro)
    result[57] = m[57]
    # Iodo compounds are assigned to group 63
    result[63] = m[63]
    # Bromo compounds are assigned to group 64
    result[64] = m[64]
    # Ni for UNIFAC group 65 (CH#C-)
    result[65] = m[65]
    # Ni for UNIFAC group 66 (-C#C-)
    result[66] = m[66]
    # Ni for UNIFAC group 71 (aromatic fluoride- ACF)
    result[71] = m[71]
    # Ni for UNIFAC group 94 (primary acid amides C(=O)NH2)
    result[94] = m[94]
    # Group 95 (secondary acid amides CH3-NH-C(=O)-)
    result[95] = m[95]
    # Group 96 (secondary acid amides -CH2-NH-C(=O)-)
    result[96] = m[96]
    # Group 97 (tertiary acid amides (CH3)2=N-C(=O)-)
    result[97] = m[97]
    # Group 98 (tertiary acid amides (CH3)(CH2)=N-C(=O)-)
    result[98] = m[98]
    # Group 99 (tertiary acid amides (CH2)2=N-C(=O)-)
    result[99] = m[99]
    # Ni for UNIFAC group 1 (-CH3), compensating for alkyl (11), ketone (18),
    # ester (21), ether (24), amine (28, 31, 34), nitrile (40), nitro (54)
    # groups
    result[1] = (
        m[1] - result[11] - result[18] - result[21] - result[24] -
        result[28] - result[31] - result[34] - result[40] - result[54])
    # Ni for UNIFAC group 2 (>CH2), compensating for alkyl (12, 44, 47), ketone
    # (19), ester (22), ether (25), amine (29, 32, 35), nitrile (41), nitro(55)
    # groups
    result[2] = (
        m[2] - result[12] - result[19] - result[22] - result[25] -
        result[29] - result[32] - result[35] - result[41] -
        result[44] - result[47] - result[55])
    # Ni for UNIFAC group 3 (>CH-), compensating for alkyl (m13, 45, 48, 50),
    # ether (26), amine (30, 33), nitro (56) groups.
    #
    # MCM_1671: O=COCC([OX2;H1])[OX2;H1]. The CH group in between the two -OH
    # groups has not been assigned. The rest of the molecules have been
    # assigned correctly (704).
    result[3] = (
        m[3] - m[13] - result[26] - result[30] - result[33] - result[45] -
        result[48] - result[50] - result[56] + m[704])
    # Ni for UNIFAC group 4 (>C<), compensating for alkyl groups (m1301, 46,
    # 49, 51, 52)
    result[4] = (
        m[4] - m[1301] - result[46] - result[49] - result[51] - result[52])
    return result


def aiomfac(compound):
    comp = composition(compound)
    # The AIOMFAC result is largely based upon the UNIFAC groups, with the
    # exception that groups 1, 2, 3, and 4 (which provide the total number of
    # CH3, CH2, >CH- and >C< groups respectively) are split out to use in the
    # processing below.
    result = unifac(compound)
    ch_groups = {k: v for k, v in result.items() if k in (1, 2, 3, 4)}
    m = matches(data.AIOMFAC_HCTAIL_SMARTS, compound)

    # Assignment of aliphatic hydrocarbon groups to "hydrocarbon tail" and
    # -OH bearing groups as defined by Marcolli and Peters.
    #
    # Part 1: First, identify whether the compound contains -OH groups as the
    # only groups based upon heteroatoms, i.e. non-hydrocarbon groups. The
    # hydrocarbon groups are then assigned to 3 categories: those that hold an
    # -OH group, those that are part of a hydrocarbon tail, and "other"
    # including hydrocarbon rings and hydrocarbon groups between two -OH
    # bearing carbons.
    #
    # Identify whether the compound has one or more -OH groups and no other
    # non-hydrocarbon groups. This is accomplished by checking all UNIFAC
    # groups from 17 upwards, except 65 and 66 (alcohols), are zero.
    if all(result[g] == 0 for g in result if g >= 17 and g not in (65, 66, 70)):
        alcohols_alone = result[14]
    else:
        alcohols_alone = 0

    # Part 2: The above calculation has identified whether the compound only
    # contains alcohol groups. For such compounds we need to identify the
    # different types of CX4 (hydrocarbon groups) and assign them to one of
    # three classes: CX4 bearing -OH, CX4 in hydrocarbon chain, and other CX4.
    #
    # The CX4 bearing -OH is relatively easy to identify using SMARTS. To
    # identify the split between hydrocarbon chain and others it is easier to
    # identify the four types of CX4 that classify as "other" and assign the
    # rest to the hydrocarbon tail.
    #
    # "Other" CX4:-
    # 1. CX4 in a hydrocarbon ring
    # 2. Me groups attached to a CX4 group bearing an -OH group
    # 3. CX4 groups on a direct chain between two or more carbons bearing -OH
    #    groups
    #
    # Contribution to hydrocarbon tails = all CX4 minus CX4 (bearing -OH) minus
    # sum of "others" as defined above. This has to be done for each type of
    # CX4 (i.e. CH3, CH2, >CH- and >C<).

    # The CH2, >CH-, and >C< groups within a ring are assigned to the
    # appropriate groups in the result. NOTE: result groups 142-144 are
    # modified further on in this function.
    if (alcohols_alone > 0) and (m[1] > 0 or m[2] > 0 or m[3] > 0):
        result[142] = m[1] # ring CH2 groups
        result[143] = m[2] # ring CH< groups
        result[144] = m[3] # ring >C< groups
        ch_groups[2] -= result[142]
        ch_groups[3] -= result[143]
        ch_groups[4] -= result[144]
    else:
        result[142] = result[143] = result[144] = 0

    # The next step is to identify the type of CHn that bears the -OH. The
    # basic data for this is supplied by SMARTS 31-34. This correction is
    # applied to all alcohols, not just the "pure" alcohols. Note that
    # hydroperoxide groups that contribute an -OH group will not be picked up
    # by SMARTS 31-34 and hence will not assign any groups to result.
    #
    # The CHn bearing the -OH is assigned to the appropriate group in result
    # providing the group is available. It might not be available if the CH2,
    # CH, or C group bearing the -OH is already part of another UNIFAC group
    # (e.g. CH2(OH)Cl: m[44] + -OH).
    if result[14] > 0:
        result[149] = min(m[31], ch_groups[1]) # methanol
        result[150] = min(m[32], ch_groups[2]) # other primary alcohols (CH2-OH)
        result[151] = min(m[33], ch_groups[3]) # secondary alcohols (CH-OH)
        result[152] = min(m[34], ch_groups[4]) # tertiary alcohols (C-OH)
        ch_groups[1] -= result[149]
        ch_groups[2] -= result[150]
        ch_groups[3] -= result[151]
        ch_groups[4] -= result[152]
    else:
        result[149] = result[150] = result[151] = result[152] = 0

    # The next step is to identify the methyl groups that are attached to the
    # CHn that bear an -OH group. This is done using SMARTS 41-46.
    methyl_groups = sum((
        1 * m[41], # primary alcohol with 1 Me group (ethanol)
        2 * m[42], # secondary alcohol with 2 Me groups (iso-propyl alcohol)
        1 * m[43], # one Me group in secondary alcohol (one potential "tail")
        3 * m[44], # three Me groups in tertiary alcohol (tertiary butanol)
        2 * m[45], # two Me groups in tertiary alcohol (one potential "tail")
        1 * m[46], # one Me group in tertiary alcohol (two potential "tails")
        ))
    if (alcohols_alone > 0) and (methyl_groups > 0):
        result[141] = methyl_groups
        ch_groups[1] -= result[141]
    else:
        result[141] = 0

    # Part 3: The final part is the assignment of CHn to the "other" category
    # for the carbons on a chain between two or more -OH groups.
    #
    # Calculation for DIOLS. The analysis is similar to the above except that
    # now the carbons between the two -OH groups have to be treated as excluded
    # for a hydrocarbon tail. These are identified by SMARTS 100-110. For
    # diols only one of these SMARTS will be a hit.
    if alcohols_alone == 2:
        c_atoms_between_oh = sum(
            (key - 100) * value
            for key, value in m.items()
            if 100 <= key < 111
            )

    # Calculation for TRIOLS. The analysis is similar to above, except that now
    # the carbons between the three -OH groups have to be treated as excluded
    # for a hydrocarbon tail. These are identified by a series of SMARTS
    # 100-107.  For triols there will be three hits in this set of SMARTS.
    elif alcohols_alone == 3:
        n = sorted(
            key - 98
            for key, value in m.items()
            if 100 <= key < 111
            and value > 0
            )
        if len(n) > 0:
            if len(n) == 1:
                min_carbons = mid_carbons = max_carbons = n[0]
            elif len(n) == 2:
                min_carbons, max_carbons = n
                mid_carbons = max_carbons
            elif len(n) == 3:
                min_carbons, mid_carbons, max_carbons = n
            else:
                min_carbons, mid_carbons, max_carbons = n
            #WARNING: need appropriate capture here
            # assert False
            if (max_carbons + 1 - min_carbons - mid_carbons) == 0:
                # Number of C atoms between the two extreme alcohol groups minus
                # the three Cs carrying the -OH groups
                c_atoms_between_oh = max_carbons - 3
            else:
                # See Book 7 (XXX What's book 7?) p100 for working
                c_atoms_between_oh = ((sum(n) + 1) / 2) - 4

            # For all triols and diols it is assumed that the carbons between those
            # holding the -OH groups are CH2 unless we run out of CH2 groups, in which
            # case the >CH- groups will be used, or >C< if they run out.
            if alcohols_alone in (2, 3):
                n = min(c_atoms_between_oh, ch_groups[2]) # CH2 groups
                result[142] += n
                ch_groups[2] -= n
                c_atoms_between_oh -= n
                n = min(c_atoms_between_oh, ch_groups[3]) # >CH- groups
                result[143] += n
                ch_groups[3] -= n
                c_atoms_between_oh -= n
                n = min(c_atoms_between_oh, ch_groups[4]) # >C< groups
                result[144] += n
                ch_groups[4] -= n
                c_atoms_between_oh -= n

    # Part 4: The assignment of the groups in the hydrocarbon tail. This
    # consists of everything left over from the previous parts
    if alcohols_alone > 0:
        result[145] = ch_groups[1] # CH3 groups
        result[146] = ch_groups[2] # CH2 groups
        result[147] = ch_groups[3] # >CH- groups
        result[148] = ch_groups[4] # >C< groups
        ch_groups[1] = ch_groups[2] = ch_groups[3] = ch_groups[4] = 0
    else:
        result[145] = result[146] = result[147] = result[148] = 0
    result[1] = ch_groups[1]
    result[2] = ch_groups[2]
    result[3] = ch_groups[3]
    result[4] = ch_groups[4]

    # Rewrite several group numbers
    result[281] = result[191]
    result[282] = result[231]
    result[283] = result[261]
    result[284] = result[531]
    result[285] = result[561]
    del result[191]
    del result[231]
    del result[261]
    del result[531]
    del result[561]

    return result

def FP4(compound):
    #The SMARTS stored in FP4.smarts have been taken from the Pybel package and 
    #thus retain that copyright. Please read the file 'FP4.smarts.readme' for
    #more information
    comp = composition(compound)
    m = matches(data.FP4_SMARTS, compound)
    result = {}
    result[1]=m[1]
    result[2]=m[2]
    result[3]=m[3]
    result[4]=m[4]
    result[5]=m[5]
    result[6]=m[6]
    result[7]=m[7]
    result[8]=m[8]
    result[9]=m[9]
    result[10]=m[10]
    result[11]=m[11]
    result[12]=m[12]
    result[13]=m[13]
    result[14]=m[14]
    result[15]=m[15]
    result[16]=m[16]
    result[17]=m[17]
    result[18]=m[18]
    result[19]=m[19]
    result[20]=m[20]
    result[21]=m[21]
    result[22]=m[22]
    result[23]=m[23]
    result[24]=m[24]
    result[25]=m[25]
    result[26]=m[26]
    result[27]=m[27]
    result[28]=m[28]
    result[29]=m[29]
    result[30]=m[30]
    result[31]=m[31]
    result[32]=m[32]
    result[33]=m[33]
    result[34]=m[34]
    result[35]=m[35]
    result[36]=m[36]
    result[37]=m[37]
    result[38]=m[38]
    result[39]=m[39]
    result[40]=m[40]
    result[41]=m[41]
    result[42]=m[42]
    result[43]=m[43]
    result[44]=m[44]
    result[45]=m[45]
    result[46]=m[46]
    result[47]=m[47]
    result[48]=m[48]
    result[49]=m[49]
    result[50]=m[50]
    result[51]=m[51]
    result[52]=m[52]
    result[53]=m[53]
    result[54]=m[54]
    result[55]=m[55]
    result[56]=m[56]
    result[57]=m[57]
    result[58]=m[58]
    result[59]=m[59]
    result[60]=m[60]
    result[61]=m[61]
    result[62]=m[62]
    result[63]=m[63]
    result[64]=m[64]
    result[65]=m[65]
    result[66]=m[66]
    result[67]=m[67]
    result[68]=m[68]
    result[69]=m[69]
    result[70]=m[70]
    result[71]=m[71]
    result[72]=m[72]
    result[73]=m[73]
    result[74]=m[74]
    result[75]=m[75]
    result[76]=m[76]
    result[77]=m[77]
    result[78]=m[78]
    result[79]=m[79]
    result[80]=m[80]
    result[81]=m[81]
    result[82]=m[82]
    result[83]=m[83]
    result[84]=m[84]
    result[85]=m[85]
    result[86]=m[86]
    result[87]=m[87]
    result[88]=m[88]
    result[89]=m[89]
    result[90]=m[90]
    result[91]=m[91]
    result[92]=m[92]
    result[93]=m[93]
    result[94]=m[94]
    result[95]=m[95]
    result[96]=m[96]
    result[97]=m[97]
    result[98]=m[98]
    result[99]=m[99]
    result[100]=m[100]
    result[101]=m[101]
    result[102]=m[102]
    result[103]=m[103]
    result[104]=m[104]
    result[105]=m[105]
    result[106]=m[106]
    result[107]=m[107]
    result[108]=m[108]
    result[109]=m[109]
    result[110]=m[110]
    result[111]=m[111]
    result[112]=m[112]
    result[113]=m[113]
    result[114]=m[114]
    result[115]=m[115]
    result[116]=m[116]
    result[117]=m[117]
    result[118]=m[118]
    result[119]=m[119]
    result[120]=m[120]
    result[121]=m[121]
    result[122]=m[122]
    result[123]=m[123]
    result[124]=m[124]
    result[125]=m[125]
    result[126]=m[126]
    result[127]=m[127]
    result[128]=m[128]
    result[129]=m[129]
    result[130]=m[130]
    result[131]=m[131]
    result[132]=m[132]
    result[133]=m[133]
    result[134]=m[134]
    result[135]=m[135]
    result[136]=m[136]
    result[137]=m[137]
    result[138]=m[138]
    result[139]=m[139]
    result[140]=m[140]
    result[141]=m[141]
    result[142]=m[142]
    result[143]=m[143]
    result[144]=m[144]
    result[145]=m[145]
    result[146]=m[146]
    result[147]=m[147]
    result[148]=m[148]
    result[149]=m[149]
    result[150]=m[150]
    result[151]=m[151]
    result[152]=m[152]
    result[153]=m[153]
    result[154]=m[154]
    result[155]=m[155]
    result[156]=m[156]
    result[157]=m[157]
    result[158]=m[158]
    result[159]=m[159]
    result[160]=m[160]
    result[161]=m[161]
    result[162]=m[162]
    result[163]=m[163]
    result[164]=m[164]
    result[165]=m[165]
    result[166]=m[166]
    result[167]=m[167]
    result[168]=m[168]
    result[169]=m[169]
    result[170]=m[170]
    result[171]=m[171]
    result[172]=m[172]
    result[173]=m[173]
    result[174]=m[174]
    result[175]=m[175]
    result[176]=m[176]
    result[177]=m[177]
    result[178]=m[178]
    result[179]=m[179]
    result[180]=m[180]
    result[181]=m[181]
    result[182]=m[182]
    result[183]=m[183]
    result[184]=m[184]
    result[185]=m[185]
    result[186]=m[186]
    result[187]=m[187]
    result[188]=m[188]
    result[189]=m[189]
    result[190]=m[190]
    result[191]=m[191]
    result[192]=m[192]
    result[193]=m[193]
    result[194]=m[194]
    result[195]=m[195]
    result[196]=m[196]
    result[197]=m[197]
    result[198]=m[198]
    result[199]=m[199]
    result[200]=m[200]
    result[201]=m[201]
    result[202]=m[202]
    result[203]=m[203]
    result[204]=m[204]
    result[205]=m[205]
    result[206]=m[206]
    result[207]=m[207]
    result[208]=m[208]
    result[209]=m[209]
    result[210]=m[210]
    result[211]=m[211]
    result[212]=m[212]
    result[213]=m[213]
    result[214]=m[214]
    result[215]=m[215]
    result[216]=m[216]
    result[217]=m[217]
    result[218]=m[218]
    result[219]=m[219]
    result[220]=m[220]
    result[221]=m[221]
    result[222]=m[222]
    result[223]=m[223]
    result[224]=m[224]
    result[225]=m[225]
    result[226]=m[226]
    result[227]=m[227]
    result[228]=m[228]
    result[229]=m[229]
    result[230]=m[230]
    result[231]=m[231]
    result[232]=m[232]
    result[233]=m[233]
    result[234]=m[234]
    result[235]=m[235]
    result[236]=m[236]
    result[237]=m[237]
    result[238]=m[238]
    result[239]=m[239]
    result[240]=m[240]
    result[241]=m[241]
    result[242]=m[242]
    result[243]=m[243]
    result[244]=m[244]
    result[245]=m[245]
    result[246]=m[246]
    result[247]=m[247]
    result[248]=m[248]
    result[249]=m[249]
    result[250]=m[250]
    result[251]=m[251]
    result[252]=m[252]
    result[253]=m[253]
    result[254]=m[254]
    result[255]=m[255]
    result[256]=m[256]
    result[257]=m[257]
    result[258]=m[258]
    result[259]=m[259]
    result[260]=m[260]
    result[261]=m[261]
    result[262]=m[262]
    result[263]=m[263]
    result[264]=m[264]
    result[265]=m[265]
    result[266]=m[266]
    result[267]=m[267]
    result[268]=m[268]
    result[269]=m[269]
    result[270]=m[270]
    result[271]=m[271]
    result[272]=m[272]
    result[273]=m[273]
    result[274]=m[274]
    result[275]=m[275]
    result[276]=m[276]
    result[277]=m[277]
    result[278]=m[278]
    result[279]=m[279]
    result[280]=m[280]
    result[281]=m[281]
    result[282]=m[282]
    result[283]=m[283]
    result[284]=m[284]
    result[285]=m[285]
    result[286]=m[286]
    result[287]=m[287]
    result[288]=m[288]
    result[289]=m[289]
    result[290]=m[290]
    result[291]=m[291]
    result[292]=m[292]
    result[293]=m[293]
    result[294]=m[294]
    result[295]=m[295]
    result[296]=m[296]
    result[297]=m[297]
    result[298]=m[298]
    result[299]=m[299]
    result[300]=m[300]
    result[301]=m[301]
    result[302]=m[302]
    result[303]=m[303]
    result[304]=m[304]
    result[305]=m[305]
    result[306]=m[306]
    result[307]=m[307]
    result[308]=m[308]
    result[309]=m[309]    
    
    return result
    
    
def MACCS(compound):
    #The SMARTS stored in MACCS.smarts have been taken from the Pybel package and 
    #thus retain that copyright. Please read the file 'MACCS.smarts.readme' for
    #more information
    #Note also that in our implementation we dont force identification to 1
    comp = composition(compound)
    m = matches(data.MACCS_SMARTS, compound)
    result = {}
    result[1]=m[1]
    result[2]=m[2]
    result[3]=m[3]
    result[4]=m[4]
    result[5]=m[5]
    result[6]=m[6]
    result[7]=m[7]
    result[8]=m[8]
    result[9]=m[9]
    result[10]=m[10]
    result[11]=m[11]
    result[12]=m[12]
    result[13]=m[13]
    result[14]=m[14]
    result[15]=m[15]
    result[16]=m[16]
    result[17]=m[17]
    result[18]=m[18]
    result[19]=m[19]
    result[20]=m[20]
    result[21]=m[21]
    result[22]=m[22]
    result[23]=m[23]
    result[24]=m[24]
    result[25]=m[25]
    result[26]=m[26]
    result[27]=m[27]
    result[28]=m[28]
    result[29]=m[29]
    result[30]=m[30]
    result[31]=m[31]
    result[32]=m[32]
    result[33]=m[33]
    result[34]=m[34]
    result[35]=m[35]
    result[36]=m[36]
    result[37]=m[37]
    result[38]=m[38]
    result[39]=m[39]
    result[40]=m[40]
    result[41]=m[41]
    result[42]=m[42]
    result[43]=m[43]
    result[44]=m[44]
    result[45]=m[45]
    result[46]=m[46]
    result[47]=m[47]
    result[48]=m[48]
    result[49]=m[49]
    result[50]=m[50]
    result[51]=m[51]
    result[52]=m[52]
    result[53]=m[53]
    result[54]=m[54]
    result[55]=m[55]
    result[56]=m[56]
    result[57]=m[57]
    result[58]=m[58]
    result[59]=m[59]
    result[60]=m[60]
    result[61]=m[61]
    result[62]=m[62]
    result[63]=m[63]
    result[64]=m[64]
    result[65]=m[65]
    result[66]=m[66]
    result[67]=m[67]
    result[68]=m[68]
    result[69]=m[69]
    result[70]=m[70]
    result[71]=m[71]
    result[72]=m[72]
    result[73]=m[73]
    result[74]=m[74]
    result[75]=m[75]
    result[76]=m[76]
    result[77]=m[77]
    result[78]=m[78]
    result[79]=m[79]
    result[80]=m[80]
    result[81]=m[81]
    result[82]=m[82]
    result[83]=m[83]
    result[84]=m[84]
    result[85]=m[85]
    result[86]=m[86]
    result[87]=m[87]
    result[88]=m[88]
    result[89]=m[89]
    result[90]=m[90]
    result[91]=m[91]
    result[92]=m[92]
    result[93]=m[93]
    result[94]=m[94]
    result[95]=m[95]
    result[96]=m[96]
    result[97]=m[97]
    result[98]=m[98]
    result[99]=m[99]
    result[100]=m[100]
    result[101]=m[101]
    result[102]=m[102]
    result[103]=m[103]
    result[104]=m[104]
    result[105]=m[105]
    result[106]=m[106]
    result[107]=m[107]
    result[108]=m[108]
    result[109]=m[109]
    result[110]=m[110]
    result[111]=m[111]
    result[112]=m[112]
    result[113]=m[113]
    result[114]=m[114]
    result[115]=m[115]
    if m[116] > 1:
        result[116]=1
    result[117]=m[117]
    if m[118] > 1:
        result[118]=1
    result[119]=m[119]
    result[120]=m[120]
    result[121]=m[121]
    result[122]=m[122]
    if m[123] > 1:
        result[123]=1
    result[124]=m[124]
    if m[125] > 1:
        result[125]=m[125]
    result[126]=m[126]
    result[127]=m[127]
    if m[128] > 1:
        result[128]=1
    if m[129] > 1:
        result[129]=1
    result[130]=m[130]
    result[131]=m[131]
    result[132]=m[132]
    result[133]=m[133]
    if m[134] > 1:
        result[134]=1
    result[135]=m[135]
    result[136]=m[136]
    result[137]=m[137]
    if m[138] > 3:
        result[138]=1
    if m[139] > 2:
        result[139]=1
    if m[140] > 1:
        result[140]=1
    result[141]=m[141]
    result[142]=m[142]
    if m[143] > 1:
        result[143]=1
    if m[144] > 2:
        result[144]=1
    result[145]=m[145]
    result[146]=m[146]
    if m[147] > 1:
        result[147]=1
    result[148]=m[148]
    result[149]=m[149]
    result[150]=m[150]
    result[151]=m[151]
    result[152]=m[152]
    result[153]=m[153]
    result[154]=m[154]
    result[155]=m[155]
    result[156]=m[156]
    if m[157] > 1:
        result[157]=1
    result[158]=m[158]
    result[159]=m[159]
    result[160]=m[160]
    result[161]=m[161]
    result[162]=m[162]
    result[163]=m[163]

    return result
