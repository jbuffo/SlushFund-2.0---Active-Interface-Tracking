# SlushFund-2.0---Active-Interface-Tracking
Multiphase Reactive Transport codes

Running the one-dimensional active interface tracking model is accomplished by calling the Active_Track_Zolotov_mod_int function with the following input arguments, for example if you type:

[Temperature, Salinity, Liquid_Fraction, Ice_Properties]=Active_Track_Zolotov_mod_int(TBOTTOM, SBOTTOM, INTERFACE_DEPTH, T_END, RHO_BR, DT)

Where,  
TBOTTOM is the temperature of the underlying fluid reservoir, in Kelvin.  
SBOTTOM is the salinity (ppt) of the underlying fluid reservoir.  
INTERFACE_DEPTH is the depth below the surface that the run is started at (in meters). This is used to calculate the heat flux at the top of the domain when using a Neumann boundary condition. i.e. for a surface temperature of 100 K and an ocean/reservoir temperature of 273 K, if you are interested in calculating the solidification of a water body 100 m below the surface, a thermal gradient of 2.73 K/m is applied to the top of the domain. This assumes a conductive equilibrium profile exists above the fluid, which is typical in ice systems that are not undergoing solid state convection.  
T_END is the end time of the run (in seconds). All runs start at zero seconds and will iterate until T_END is reached.  
RHO_BR is the density of the ambient ocean/reservoir (in kg/m^3). This is used to calculate heat fluxes as well the one dimensional gravity drainage parameterization.  
DT is the temporal discretization step (in seconds). This will determine the number of iterations and the temporal accuracy of the simulation.

The output of the model is 4 arrays:
Temperature is the current temperatures within the domain (in Kelvin). 1D vertical array with top cell of domain as first value in array.  
Salinity is the current salinities within the domain (in ppt). 1D vertical array with top cell of domain as first value in array.  
Liquid_Fraction is the current liquid fraction in the domain.1D vertical array with top cell of domain as first value in array.  
Ice_Properties contains information on the ice that has been frozen out during the simulation and cataloged, it would be the ice that exists above the active layer of the domain which has reached an impermeability limit and is assumed to freeze in the composition of the ice it had when reaching this limit. There are 4 columns in this array [Temperature Salinity Liquid_Fraction Darcy_Velocity]  
Temperature is the temperature at the time the ice solidified.  
Salinity is the salinity of the liquid portion of the cell when it solidified.  
Liquid Fraction is the liquid fraction of the cell when it solidified.  
Darcy Velocity has been set to zero in this model as gravity drainage is taken to be the main mode of desalination.  

Each element in these lists corresponds to a vertical cell of 1 cm unless modified within the code.

Bulk Ice salinity can be calculated by multiplying Salinity and Liquid Fraction.

NOTE!!! - COMMENT OUT OR ALTER LINES 303-305 IF YOU WISH TO NOT SAVE OR MODIFY THE SAVE FILE NAME

————————————————————
Optional Manual Modification
The code can easily me edited and modified to include alternate ice/ocean/brine properties and additional physics. The code was written functionally for easy application, but manually editing lines within the code will drastically extend the models capabilities.

Select terms that can easily be modified (property - variable in code):  
Surface Temperature - Ttop  
Melting temperature of pure ice - Tm  
Thermal diffusion of ice - k_i  
Thermal diffusion of brine - k_br  
Molecular diffusion coefficient (salt/impurity diffusion) - k_s  
Grid size - dz  
Domain size - H  
Specific heat of ocean/brine - c_br  
Specific heat of ice - c_i  
Latent heat of fusion - L  
Density of ice - rho_i  
Error tolerances (salt, liquid fraction, temperature) - STol, PhiTol, TTol  
Gravity - g  
ocean/brine viscosity - mu  
Critical porosity cutoff - phi_c  

Boundary condition selection for upper boundary - select between dirilecht/neumann at line 59

To modify liquidus curves - in the MATLAB file one_D_adv_ARD_FREZ_track.m find line 37, either select a preexisting form of the ‘delta_T’ variable by leaving your selection uncommented or introduce and alternative ‘delta_T’ variable. A similar form is recommended, where delta_T is calculated by subtracting the salinity/impurity dependent freezing point (here typically a quadratic function) from the freezing point of pure ice (Tm).

