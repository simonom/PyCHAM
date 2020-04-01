# vim: set et sw=4 sts=4 fileencoding=utf-8:
#
# Copyright (c) 2014 Manchester University.
# All Rights Reserved.

from __future__ import (
    unicode_literals,
    absolute_import,
    print_function,
    division,
    )
str = type('')
try:
    from itertools import izip
    zip = izip
except ImportError:
    pass


import re

import pybel
import pkg_resources as pkg


def smarts(s):
	# uncomment next two lines if using python2, but comment for python3
#     if not isinstance(s, bytes):
#         s = s.encode('ascii')
    try:
        return pybel.Smarts(s)
    except IOError as e:
        # Convert pybel's IOError (?!) into a ValueError
        raise ValueError(str(e))


_parse_re = re.compile(r'^\s*((?P<data>([^#]|\S#)*)(\s+#.*)?|#.*)$')
def _read_data(
        filename, key_conv=str, value_conv=float, key_col=0, value_col=1):
    result = {}
    cols = None
    for count, line in enumerate(
            pkg.resource_stream(__name__, filename), start=1):
        data = _parse_re.match(line.decode('utf-8')).group('data')
        if data:
            data = data.split()
            try:
                if cols is None:
                    cols = len(data)
                elif len(data) != cols:
                    raise ValueError(
                            'Unexpected number of values (expected %d)' % cols)
                key = key_conv(data[key_col])
                value = value_conv(data[value_col])
                if key in result:
                    raise ValueError(
                            'Duplicate definition for group %s' % key)
                result[key] = value
            except (IndexError, ValueError) as e:
                e.args += ('on line %d of %s' % (count, filename),)
                raise
    return result


def _read_smarts(filename):
    return _read_data(filename, key_conv=int, value_conv=smarts)


def _read_matrix(filename, key_conv=str, value_conv=float, symmetric=True):
    result = {}
    col_keys = []
    row_keys = []
    for count, line in enumerate(
            pkg.resource_stream(__name__, filename), start = 1):
        data = _parse_re.match(line.decode('utf-8')).group('data')
        if data:
            try:
                if not col_keys:
                    col_keys = [key_conv(key) for key in data.split()]
                else:
                    row_key, values = data.split(None, 1)
                    row_key = key_conv(row_key)
                    if row_key in result:
                        raise ValueError(
                                'Duplicate definition for row %s' % key)
                    values = [value_conv(value) for value in values.split()]
                    if not len(col_keys) == len(values):
                        raise ValueError(
                                'Expected %d values but found %d' % (
                                    len(col_keys), len(values)))
                    row_keys.append(row_key)
                    result[row_key] = {
                        col_key: value
                        for col_key, value in zip(col_keys, values)
                        }
            except ValueError as e:
                e.args += ('on line %d of %s' % (count, filename))
                raise
    if symmetric:
        if sorted(row_keys) != sorted(col_keys):
            raise ValueError('Column and row keys are not identical')
        for row_key in row_keys:
            for col_key in col_keys:
                if result[row_key][col_key] != result[col_key][row_key]:
                    raise ValueError(
                            'Value %f in row %s, column %s does not match '
                            '%f in row %s, column %s' % (
                                result[row_key][col_key],
                                row_key, col_key,
                                result[col_key][row_key],
                                col_key, row_key,
                                ))
    return result


STEIN_AND_BROWN_BOILING_POINT          = _read_data('joback.data', key_conv=int, value_col=1)
JOBACK_BOILING_POINT                   = _read_data('joback.data', key_conv=int, value_col=2)
JOBACK_TEMPERATURE                     = _read_data('joback.data', key_conv=int, value_col=3)
JOBACK_PRESSURE                        = _read_data('joback.data', key_conv=int, value_col=4)
JOBACK_VOLUME                          = _read_data('joback.data', key_conv=int, value_col=5)
SCHROEDER_DENSITY                      = _read_data('schroeder.data', key_conv=int)
LE_BAS_DENSITY                         = _read_data('le_bas.data', key_conv=int)
NANNOOLAL_BOILING_POINT_PRIMARY        = _read_data('nannoolal_primary.data', key_conv=int, value_col=1)
NANNOOLAL_VAPOUR_PRESSURE_PRIMARY      = _read_data('nannoolal_primary.data', key_conv=int, value_col=2)
NANNOOLAL_TEMPERATURE_PRIMARY          = _read_data('nannoolal_primary.data', key_conv=int, value_col=3)
NANNOOLAL_PRESSURE_PRIMARY             = _read_data('nannoolal_primary.data', key_conv=int, value_col=4)
NANNOOLAL_VOLUME_PRIMARY               = _read_data('nannoolal_primary.data', key_conv=int, value_col=5)
NANNOOLAL_BOILING_POINT_SECONDARY      = _read_data('nannoolal_secondary.data', key_conv=int, value_col=1)
NANNOOLAL_VAPOUR_PRESSURE_SECONDARY    = _read_data('nannoolal_secondary.data', key_conv=int, value_col=2)
NANNOOLAL_TEMPERATURE_SECONDARY        = _read_data('nannoolal_secondary.data', key_conv=int, value_col=3)
NANNOOLAL_PRESSURE_SECONDARY           = _read_data('nannoolal_secondary.data', key_conv=int, value_col=4)
NANNOOLAL_VOLUME_SECONDARY             = _read_data('nannoolal_secondary.data', key_conv=int, value_col=5)
EVAPORATION_A                          = _read_data('evaporation.data', value_col=1)
EVAPORATION_B                          = _read_data('evaporation.data', value_col=2)
AIOMFAC_SALT_CATION_GROUP              = _read_data('aiomfac_salts.data', key_conv=int, value_col=1, value_conv=int)
AIOMFAC_SALT_CATION_STOICH             = _read_data('aiomfac_salts.data', key_conv=int, value_col=2)
AIOMFAC_SALT_ANION_GROUP               = _read_data('aiomfac_salts.data', key_conv=int, value_col=3, value_conv=int)
AIOMFAC_SALT_ANION_STOICH              = _read_data('aiomfac_salts.data', key_conv=int, value_col=4)
AIOMFAC_SALT_DENSITY                   = _read_data('aiomfac_salts.data', key_conv=int, value_col=5)
AIOMFAC_SALT_MASS                      = _read_data('aiomfac_salts.data', key_conv=int, value_col=6)
AIOMFAC_SALT_NAME                      = _read_data('aiomfac_salts.data', key_conv=int, value_col=7, value_conv=str)
AIOMFAC_ION_CHARGE                     = _read_data('aiomfac_ions.data', key_conv=int)
AIOMFAC_MAIN_GROUP                     = _read_data('aiomfac_main.data', key_conv=int, value_col=1, value_conv=int)
AIOMFAC_MASS                           = _read_data('aiomfac_main.data', key_conv=int, value_col=2)
AIOMFAC_RI                             = _read_data('aiomfac_main.data', key_conv=int, value_col=3)
AIOMFAC_QI                             = _read_data('aiomfac_main.data', key_conv=int, value_col=4)

AIOMFAC_ION_CHARGE_ABS = {
    group: abs(charge)
    for group, charge in AIOMFAC_ION_CHARGE.items()
    }

AIOMFAC_ION_SALT = {
    (cation, AIOMFAC_SALT_ANION_GROUP[group]): group
    for group, cation in AIOMFAC_SALT_CATION_GROUP.items()
    }

# NOTE: 2012-12-12 - The following extensions are ones suggsted by Mark Barley,
# November 2012 These are to make sure the MCM compounds are parsed correctly.
# For ANY official AIOMFAC extensions then these can be removed if required.
AIOMFAC_MAIN_GROUP[281] = AIOMFAC_MAIN_GROUP[19]
AIOMFAC_MAIN_GROUP[282] = AIOMFAC_MAIN_GROUP[23]
AIOMFAC_MAIN_GROUP[283] = AIOMFAC_MAIN_GROUP[26]
AIOMFAC_MAIN_GROUP[284] = AIOMFAC_MAIN_GROUP[53]
AIOMFAC_MAIN_GROUP[285] = AIOMFAC_MAIN_GROUP[56]
AIOMFAC_MASS[281] = AIOMFAC_MASS[19]
AIOMFAC_MASS[282] = AIOMFAC_MASS[23]
AIOMFAC_MASS[283] = AIOMFAC_MASS[26]
AIOMFAC_MASS[284] = AIOMFAC_MASS[53]
AIOMFAC_MASS[285] = AIOMFAC_MASS[56]
AIOMFAC_RI[281] = AIOMFAC_RI[19]
AIOMFAC_RI[282] = AIOMFAC_RI[23]
AIOMFAC_RI[283] = AIOMFAC_RI[26]
AIOMFAC_RI[284] = AIOMFAC_RI[53]
AIOMFAC_RI[285] = AIOMFAC_RI[56]
AIOMFAC_QI[281] = AIOMFAC_QI[19]
AIOMFAC_QI[282] = AIOMFAC_QI[23]
AIOMFAC_QI[283] = AIOMFAC_QI[26]
AIOMFAC_QI[284] = AIOMFAC_QI[53]
AIOMFAC_QI[285] = AIOMFAC_QI[56]

COMPOSITION_SMARTS                     = _read_smarts('composition.smarts')
EVAPORATION_SMARTS                     = _read_smarts('evaporation.smarts')
GIROLAMI_SMARTS                        = _read_smarts('girolami.smarts')
LE_BAS_SMARTS                          = _read_smarts('le_bas.smarts')
MOLLER_SMARTS                          = _read_smarts('moller.smarts')
NANNOOLAL_SMARTS_PRIMARY               = _read_smarts('nannoolal_primary.smarts')
NANNOOLAL_SMARTS_SECONDARY             = _read_smarts('nannoolal_secondary.smarts')
SCHROEDER_SMARTS                       = _read_smarts('schroeder.smarts')
STEIN_AND_BROWN_SMARTS                 = _read_smarts('stein_and_brown.smarts')
AIOMFAC_UNIFAC_SMARTS                  = _read_smarts('aiomfac_unifac.smarts')
AIOMFAC_HCTAIL_SMARTS                  = _read_smarts('aiomfac_hctail.smarts')
AIOMFAC_ION_SMARTS                     = _read_smarts('aiomfac_ions.smarts')
FP4_SMARTS                             = _read_smarts('FP4.smarts')
MACCS_SMARTS                           = _read_smarts('MACCS.smarts')

NANNOOLAL_BOILING_POINT_INTERACTIONS   = _read_matrix('nannoolal_bp_interactions.matrix')
NANNOOLAL_VAPOUR_PRESSURE_INTERACTIONS = _read_matrix('nannoolal_vp_interactions.matrix')
NANNOOLAL_TEMPERATURE_INTERACTIONS     = _read_matrix('nannoolal_tc_interactions.matrix')
NANNOOLAL_PRESSURE_INTERACTIONS        = _read_matrix('nannoolal_pc_interactions.matrix')
NANNOOLAL_VOLUME_INTERACTIONS          = _read_matrix('nannoolal_vc_interactions.matrix')
AIOMFAC_SR_INTERACTIONS                = _read_matrix('aiomfac_sr_interactions.matrix', key_conv=int, symmetric=False)
# Now adding interaction matrices from the org-ion and ion-ion
AIOMFAC_MR_ORG_ION_INTERACTIONS_b1     = _read_matrix('aiomfac_mr_org_ion_interactions_b1.matrix', key_conv=int, symmetric=False)
AIOMFAC_MR_ORG_ION_INTERACTIONS_b2     = _read_matrix('aiomfac_mr_org_ion_interactions_b2.matrix', key_conv=int, symmetric=False)
#AIOMFAC_MR_ORG_ION_INTERACTIONS_b3     = _read_matrix('aiomfac_mr_org_ion_interactions_b3.matrix', key_conv=int, symmetric=False)
AIOMFAC_MR_ION_ION_INTERACTIONS_b1     = _read_matrix('aiomfac_mr_ion_ion_interactions_b1.matrix', key_conv=int, symmetric=False)
AIOMFAC_MR_ION_ION_INTERACTIONS_b2     = _read_matrix('aiomfac_mr_ion_ion_interactions_b2.matrix', key_conv=int, symmetric=False)
AIOMFAC_MR_ION_ION_INTERACTIONS_b3     = _read_matrix('aiomfac_mr_ion_ion_interactions_b3.matrix', key_conv=int, symmetric=False)
#AIOMFAC_MR_ION_ION_INTERACTIONS_c1     = _read_matrix('aiomfac_mr_ion_ion_interactions_c1.matrix', key_conv=int, symmetric=False)
AIOMFAC_MR_ION_ION_INTERACTIONS_c2     = _read_matrix('aiomfac_mr_ion_ion_interactions_c2.matrix', key_conv=int, symmetric=False)
AIOMFAC_MR_ION_ION_INTERACTIONS_c3     = _read_matrix('aiomfac_mr_ion_ion_interactions_c3.matrix', key_conv=int, symmetric=False)
#AIOMFAC_MR_ION_ION_INTERACTIONS_Rc     = _read_matrix('aiomfac_mr_ion_ion_interactions_Rc.matrix', key_conv=int, symmetric=False)
#AIOMFAC_MR_ION_ION_INTERACTIONS_Qc     = _read_matrix('aiomfac_mr_ion_ion_interactions_Qc.matrix', key_conv=int, symmetric=False)


