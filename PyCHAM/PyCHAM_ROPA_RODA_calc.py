# module to estimate rates of production and rates of destruction,
# can accept one or more components to return results for

def main(output_by_sim, comp_names_to_plot):
    import netCDF4 as nc
    import numpy as np

    ds = nc.Dataset(output_by_sim) # open file
            
    timehr = ds['time'][:]/3600. # get time (seconds to hours)

    # get total number of components
    num_comp = np.array((ds['number_of_components'][:]))[0]

    # get chemical scheme names of components
    comp_names = ds['component_chemical_scheme_name'][:].tolist()

    # prepare to hold concentrations of components (molecules/cm^3)
    yrec = np.zeros((len(timehr), len(comp_names)))

    # loop through these components to get their gas-phase concentrations (molecules/cm^3)
    for ci in range(len(comp_names)):
                
        # set variable name
        var_name = str('number_concentration_of_' +
            'gas_phase_' + comp_names[ci] + '_molecules_in_air')

        # get values in molecules/cm^3
        try:
            yrec[:, ci] = ds[str('/concentrations_g/' + var_name)][:]
        except: # in case this component not saved, then continue to next component
            continue

    # rates of reactions (molecules/cm^3/s), time in rows
    # and reactions in columns
    reac_rate = ds['rates_of_reaction_of_all_reactions'][:, :]

    # indices of components acting as reactants (columns) per reaction (rows)
    rindx_g = ds['reactant_indices'][:]

    # reactant stoichiometries per reaction (row) per reactant (column)
    rstoi_g = ds['reactant_stoichiometries'][:]

    # number of reactants per reaction
    nreac_g = ds['number_of_reactants_per_reaction'][:]

    # indices of components acting as products (columns) per reaction (rows)
    pindx_g = ds['product_indices'][:]

    # product stoichiometries per reaction (row) per reactant (column)
    pstoi_g = ds['product_stoichiometries'][:]

    # number of products per reaction
    nprod_g = ds['number_of_products_per_reaction'][:]

    # strings of gas-phase chemical reaction equation per reaction
    eq_str = ds['equations_per_reaction'][:]


    # -----------------------------------------------
    # now calculate reaction rates
    # get the gas-phase concentrations(molecules/cm^3)
    # in preparation for calculation of rates
    y = np.zeros((yrec.shape[0], yrec.shape[1]))
    y[:, 0:num_comp] = yrec[:, 0:num_comp]

    # prepare dictionaries to hold production and loss rates
    # and equation strings for all components
    lr_all = {}
    leq_all = {}
    pr_all = {}
    peq_all = {}

    for comp_name in (comp_names_to_plot): # loop through components to plot

        ci = comp_names.index(comp_name) # get index of this component

        # prepare to hold loss rates (molecules/cm^3/s)
        lr_all[comp_name] = np.zeros((y.shape[0], 0))
        # prepare to hold loss equation strings
        leq_all[comp_name] = []

        # prepare to hold production rates (molecules/cm^3/s),
        # times in rows and reactions in columns, considered 
        # components in third dimension
        pr_all[comp_name] = np.zeros((y.shape[0], 0))
        # prepare to hold production equations strings
        peq_all[comp_name] = []

        # loop through reactions to
        # identify the reactions where this component occurs as a reactant 
        # (loss process)
        # and to identify where it occurs as a product 
        # (production process), note ci is index of component of interest
        for reaci in range(rindx_g.shape[0]):

            # if a reactant
            if ci in rindx_g[reaci, 0:nreac_g[reaci]]:
                # get loss rate at all times
                lr = rstoi_g[reaci, 0:nreac_g[reaci]][ci == rindx_g[reaci, 0:nreac_g[reaci]]]*reac_rate[:, reaci]
                
                # include with all gas-phase chemistry loss rates
                lr_all[comp_name] = np.concatenate((lr_all[comp_name], lr.reshape(-1, 1)), axis=1)

                # include with all gas-phase chemistry loss reactions
                leq_all[comp_name].append(eq_str[reaci])
            
            # if a product
            if ci in pindx_g[reaci, 0:nprod_g[reaci]]:
                
                # get production rate of component of interest
                # (molecules/cm^3/s)
                pr = pstoi_g[reaci, 0:nprod_g[reaci]][ci == pindx_g[reaci, 0:nprod_g[reaci]]]*reac_rate[:, reaci]

                # include with all gas-phase chemistry production rates
                pr_all[comp_name] = np.concatenate((pr_all[comp_name], pr.reshape(-1, 1)), axis=1)

                # include with all gas-phase chemistry production reactions
                peq_all[comp_name].append(eq_str[reaci])

    return(lr_all, leq_all, pr_all, peq_all, timehr)