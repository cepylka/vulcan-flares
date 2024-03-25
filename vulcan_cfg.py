import os

# ---

# sadly, it isn't trivial to implement proper CLI arguments in this project,
# so we'll have to do that with environment variables

vulcanRunsAreChainedVarName = "VULCAN_RUNS_ARE_CHAINED"
vulcanRunOrdinalNumberVarName = "VULCAN_RUN_ORDINAL_NUMBER"

# if this environment variable is set (to any value), then VULCAN runs are meant to be "chained"
vulcanRunsAreChained = False
if vulcanRunsAreChainedVarName in os.environ:
    vulcanRunsAreChained = True

vulcanRunOrdinalNumber = None
if vulcanRunsAreChained:
    # ordinal number of the current VULCAN run
    if vulcanRunOrdinalNumberVarName in os.environ:
        try:
            vulcanRunOrdinalNumber = int(os.environ[vulcanRunOrdinalNumberVarName])
        except ValueError:
            raise SystemExit(
                f"[ERROR] The value of {vulcanRunOrdinalNumberVarName} is not an integer"
            )
        # one more check, just in case
        if not isinstance(vulcanRunOrdinalNumber, int):
            raise SystemExit(
                f"[ERROR] The value of {vulcanRunOrdinalNumberVarName} is not an integer"
            )
    else:
        raise SystemExit(
            f"[ERROR] Missing {vulcanRunOrdinalNumberVarName} environment variable"
        )

# ---

# =============================================================================
# Configuration file of VULCAN:
# =============================================================================
#This is Mamonova et al. 2025 parameters for rumnning VULCAN
# Stellar flux is variable and the light curve with flares is for 10 days
# TP profiles prodused by HELIOS with FASTCHEM for abundansies inferred
# for const_mix = {'H2':6.44000882e-09,'H2O':9.99999794e-01,'CH4': 5.55219897e-09, 'He':1.93839964e-07,'NH3':1.63705970e-09} #H2O%100 for 100% H2O atmo
# therefore the run is for 'EQ' with corrected solar abundancies file: /fastchem_vulcan/input/solar_element_abundances.dat,
# which needed to be updated for the next atmosphere in the grid
#As we change the solar abundancies file, we choose use_solar = True
# The initial frequency to update the actinic flux and optical depth is 50 steps
# Atmospheric escape limited to diffusion and calculates for H,H2, H2O, He
#This is NCNO full photonetwork

# ====== Setting up the elements included in the network ======
atom_list = ['H', 'O', 'C', 'N', 'He']

# water percentage in atmosphere
#
# - 00001 is 0.01%
# - 00010 is 0.10%
# - 00015 is 0.15%
# - 00100 is 1.00%
# - 00150 is 1.50%
# - 01000 is 10.00%
# - 05000 is 50.00%
# - 10000 is 100.00%
#
# and so on
waterPercentage = "10000H2O"
#
simulationID = f"ME1MS03_{waterPercentage}"

# ====== Setting up paths and filenames for the input and output files  ======
# input:
network = 'thermo/NCHO_full_photo_network.txt'
use_lowT_limit_rates = True
# all the nasa9 files must be placed in the folder: thermo/NASA9/
gibbs_text = 'thermo/gibbs_text.txt'
cross_folder = 'thermo/photo_cross/'
com_file = 'thermo/all_compose.txt'
# TP and Kzz (optional) file
atm_file = f"./atm/atm_{simulationID}_Kzz_{waterPercentage}.txt"

# there might be a need to take stellar flux files in a custom order
sfluxFileCustomOrderMap = {
    1: 1,
    2: 2,
    3: 3,
    4: 9,
    5: 5,
    6: 6,
    7: 7,
    8: 8,
    9: 4
}
sfluxFileIsCustomOrder = True
sfluxFileNumber = (
    sfluxFileCustomOrderMap.get(vulcanRunOrdinalNumber)
    if sfluxFileIsCustomOrder
    else vulcanRunOrdinalNumber
)
if sfluxFileNumber is None:
    raise SystemExit(
        f"[ERROR] There is no mapped number for the VULCAN run ordinal number #{vulcanRunOrdinalNumber}"
    )
# the flux density at the stellar surface
sflux_file = (
    f"./atm/stellar_flux/ME1M03_sflux_timesteps_60secH-{sfluxFileNumber}.pkl"
    if vulcanRunsAreChained
    else "./atm/stellar_flux/ME1M03_sflux_timesteps_60secH-1.pkl"
)
# whether sflux_file contains fluxes just for one moment of time,
# or is it a set of fluxes data with a time component,
# like in ./atm/stellar_flux/gj876_sflux_timesteps_60sec.pkl
fluxWithTime = True
# is it a plain-text file or a binary pickle
sflux_file_is_plaintext = False

# the file for the top boundary conditions
top_BC_flux_file = 'atm/BC_top_GJ.txt'
# the file for the lower boundary conditions
bot_BC_flux_file = 'atm/BC_bot_mars.txt'
# the file to initialize the abundances for ini_mix = 'vulcan_ini'
vul_ini = (
    f"./output/{simulationID}-{vulcanRunOrdinalNumber-1}.vul"  # takes results of the previous run
    if vulcanRunsAreChained
    else f"./output/{simulationID}.vul"
)
# output:
output_dir = 'output/'
plot_dir = 'plot/'
movie_dir = 'plot/movie/'
# output file name
out_name = (
    f"{simulationID}-{vulcanRunOrdinalNumber}.vul"
    if vulcanRunsAreChained
    else f"{simulationID}.vul"
)

# ====== Setting up the elemental abundance ======
use_solar = True  # True: using the solar abundance from Table 10. K.Lodders 2009; False: using the customized elemental abundance.
solar_ele = f"./fastchem_vulcan/input/element_abundances_{waterPercentage}.dat"
# customized elemental abundance (only read when use_solar = False)
O_H = 0.13#0.21#0.21# 0.5  #
C_H = 2.7761E-04
N_H = 8.1853E-05
# S_H = 1.3183E-05
He_H = 0.09692
# O_H = 0.5#100% H2O  =(1/3)/(2/3)
# C_H = 6.2e-3#9.39E-7#6.2e-3 * 1e-6  *0.5 #(6.224628e-3 for solar atm) 9.39392980501438e-09 1.969694973739396e-09 4.92423501923749e-07
# N_H = 1.3e-3#1.97E-7#1.3e-3# * 1e-6 * 0.5#(1.299e-3 for solar atm)
# S_H = 1.3183E-7
# He_H = 0.09#4.92E-5#0.325 * 1e-6  * 0.5#(0.325 for solar atm)

# options: 'EQ', 'const_mix', 'vulcan_ini', 'table' (for 'vulcan_ini, the T-P grids have to be exactly the same)
ini_mix = "EQ"
if vulcanRunsAreChained:
    ini_mix = "EQ" if vulcanRunOrdinalNumber == 1 else "vulcan_ini"

fastchem_met_scale = 0.1 # scaling factor for other elements in fastchem (e.g., if fastchem_met_scale = 0.1, other elements such as Si and Mg will take 0.1 solar values)

# Initialsing uniform (constant with pressure) mixing ratios (only reads when ini_mix = const_mix)

const_mix = {'H2':6.44000882e-09,'H2O':9.99999794e-01,'CH4': 5.55219897e-09, 'He':1.93839964e-07,'NH3':1.63705970e-09} #H2O%100
# const_mix = {'H2':3.68707469e-01,'H2O':6.31098173e-01,'CH4': 5.55112841e-07, 'He':1.93802588e-04,'NH3':1.63674404e-07} #H2O%60
# const_mix = {'H2':4.37900975e-01,'H2O':3.98112018e-01,'CH4': 6.65013658e-04, 'He':1.63321994e-01,'NH3':1.37111745e-04} #H2O%40
# const_mix = {'H2':5.84832141e-01,'H2O':2.51180851e-01,'CH4': 6.65013658e-04, 'He':1.63321994e-01,'NH3':1.37111745e-04} #H2O%25
# const_mix = {'H2':6.77515279e-01,'H2O':1.58497714e-01,'CH4': 6.65013658e-04, 'He':1.63321994e-01,'NH3':1.37111745e-04} #H2O%15
# const_mix = {'H2':7.36009679e-01,'H2O':1.00003313e-01,'CH4': 6.65013658e-04, 'He':1.63321994e-01,'NH3':1.37111745e-04} #H2O%10
# const_mix = {'H2':7.72912074e-01,'H2O':6.31009181e-02,'CH4': 6.65013658e-04, 'He':1.63321994e-01,'NH3':1.37111745e-04} #H2O%6
# const_mix = {'H2':7.72912074e-01,'H2O':6.31009181e-02,'CH4': 6.65013658e-04, 'He':1.63321994e-01,'NH3':1.37111745e-04} #H2O%6

# const_mix = {'CH4':2.7761E-4*2, 'O2':4.807e-4, 'He':0.09691, 'N2':8.1853E-5, 'H2':1. -2.7761E-4*2*4/2}
# const_mix = {'CH4':7.2E-4, 'O2':4.807e-4, 'He':0.16, 'N2':8.1853E-5, 'H2O':0.84 , 'H2':1.7E-3,'NH3':2.2E-4,'CO':1.0E-9}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-7,'H2O':1.00000 - 1.E-6, 'H':1.E-6 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*1.E-6, 'H2': 8.5e-07, 'CH4':4.2E-10, 'N2':6.2E-11, 'CO':8.4E-16, 'CO2':8.40E-17}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':3.24E-07,'H2O':1.0 - 1.E-06, 'H':1E-06 - 3.25E-07 - 1.5E-07 - 6.2E-08 - 1.3E-08, 'H2': 1.5E-07, 'CH4':6.2E-08, 'N2':1.3E-08}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.36904,'H2O':1.00000 - 0.36904, 'H':0.36904 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.36904, 'H2': 8.5e-01*0.36904, 'CH4':4.2E-4*0.36904, 'N2':6.2E-7*0.36904, 'CO':8.4E-10*0.36904, 'CO2':8.40E-11*0.36904}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.60189,'H2O':1.00000 - 0.60189, 'H':0.60189 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.60189, 'H2': 8.5e-01*0.60189, 'CH4':4.2E-4*0.60189, 'N2':6.2E-7*0.60189, 'CO':8.4E-10*0.60189, 'CO2':8.40E-11*0.60189}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.74881,'H2O':1.00000 - 0.74881, 'H':0.74881 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.74881, 'H2': 8.5e-01*0.74881, 'CH4':4.2E-4*0.74881, 'N2':6.2E-7*0.74881, 'CO':8.4E-10*0.74881, 'CO2':8.40E-11*0.74881}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.84151,'H2O':1.00000 - 0.84151, 'H':0.84151 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.84151, 'H2': 8.5e-01*0.84151, 'CH4':4.2E-4*0.84151, 'N2':6.2E-7*0.84151, 'CO':8.4E-10*0.84151, 'CO2':8.40E-11*0.84151}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.9,'H2O':1.00000 - 0.9, 'H':0.9 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.9, 'H2': 8.5e-01*0.9, 'CH4':4.2E-4*0.9, 'N2':6.2E-7*0.9, 'CO':8.4E-10*0.9, 'CO2':8.40E-11*0.9}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.9369,'H2O':1.00000 - 0.9369, 'H':0.9369 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.9369, 'H2': 8.5e-01*0.9369, 'CH4':4.2E-4*0.9369, 'N2':6.2E-7*0.9369, 'CO':8.4E-10*0.9369, 'CO2':8.40E-11*0.9369}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.96019,'H2O':1.00000 - 0.96019, 'H':0.96019 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.96019, 'H2': 8.5e-01*0.96019, 'CH4':4.2E-4*0.96019, 'N2':6.2E-7*0.96019, 'CO':8.4E-10*0.96019, 'CO2':8.40E-11*0.96019}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.97488,'H2O':1.00000 - 0.97488, 'H':0.97488 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.97488, 'H2': 8.5e-01*0.97488, 'CH4':4.2E-4*0.97488, 'N2':6.2E-7*0.97488, 'CO':8.4E-10*0.97488, 'CO2':8.40E-11*0.97488}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.98415,'H2O':1.00000 - 0.98415, 'H':0.98415 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.98415, 'H2': 8.5e-01*0.98415, 'CH4':4.2E-4*0.98415, 'N2':6.2E-7*0.98415, 'CO':8.4E-10*0.98415, 'CO2':8.40E-11*0.98415}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.99,'H2O':1.00000 - 0.99, 'H':0.99 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.99, 'H2': 8.5e-01*0.99, 'CH4':4.2E-4*0.99, 'N2':6.2E-7*0.99, 'CO':8.4E-10*0.99, 'CO2':8.40E-11*0.99}#'PH3':7.5E-7, 'H2S':3.7E-5,

# const_mix = {'He':1.4e-1*0.99369,'H2O':1.00000 - 0.99369, 'H':0.99369 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.99369, 'H2': 8.5e-01*0.99369, 'CH4':4.2E-4*0.99369, 'N2':6.2E-7*0.99369, 'CO':8.4E-10*0.99369, 'CO2':8.40E-11*0.99369}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.99602,'H2O':1.00000 - 0.99602, 'H':0.99602 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.99602, 'H2': 8.5e-01*0.99602, 'CH4':4.2E-4*0.99602, 'N2':6.2E-7*0.99602, 'CO':8.4E-10*0.99602, 'CO2':8.40E-11*0.99602}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.999,'H2O':1.00000 - 0.999, 'H':0.999 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.999, 'H2': 8.5e-01*0.999, 'CH4':4.2E-4*0.999, 'N2':6.2E-7*0.999, 'CO':8.4E-10*0.999, 'CO2':8.40E-11*0.999}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.99749,'H2O':1.00000 - 0.99749, 'H':0.99749 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.99749, 'H2': 8.5e-01*0.99749, 'CH4':4.2E-4*0.99749, 'N2':6.2E-7*0.99749, 'CO':8.4E-10*0.99749, 'CO2':8.40E-11*0.99749}#'PH3':7.5E-7, 'H2S':3.7E-5,
# const_mix = {'He':1.4e-1*0.99842,'H2O':1.00000 - 0.99842, 'H':0.99842 - (8.5e-01+1.4e-01+4.2E-4+6.2E-05+8.4E-10+8.4E-11)*0.99842, 'H2': 8.5e-01*0.99842, 'CH4':4.2E-4*0.99842, 'N2':6.2E-7*0.99842, 'CO':8.4E-10*0.99842, 'CO2':8.40E-11*0.99842}#'PH3':7.5E-7, 'H2S':3.7E-5,

# ====== Setting up photochemistry ======
use_photo = True
# astronomy input
r_star = 0.303 # stellar radius in solar radius
Rp = 0.0892141781*7.1492E9 # Planetary radius (cm) (for computing gravity) (7.1492E9 is R_jup)
orbit_radius =0.0304 # planet-star distance in A.U.
sl_angle = 58 /180.*3.14159 # the zenith angle of the star in degree (usually 58 deg for the dayside average)
f_diurnal = 0.5 # to account for the diurnal average of solar flux (i.e. 0.5 for Earth; 1 for tidally-locked planets)
scat_sp = ['H2', 'He'] # the bulk gases that contribute to Rayleigh scattering
T_cross_sp = [] # warning: slower start! available atm: 'CO2','H2O','NH3', 'SH','H2S','SO2', 'S2', 'COS', 'CS2'

edd = 0.5 # the Eddington coefficient 
dbin1 = 0.1  # the uniform bin width < dbin_12trans (nm)
dbin2 = 2.   # the uniform bin width > dbin_12trans (nm)
dbin_12trans = 240. # the wavelength switching from dbin1 to dbin2 (nm)

# the frequency to update the actinic flux and optical depth
ini_update_photo_frq = 100
final_update_photo_frq = 500

# ====== Setting up ionchemistry ======
use_ion = False
if use_photo == False and use_ion == True:
    print ('Warning: use_ion = True but use_photo = False')
# photoionization needs to run together with photochemistry


# ====== Setting up parameters for the atmosphere ======
atm_base = 'H2O' #Options: 'H2','H2O', 'N2', 'O2', 'CO2 -- the bulk gas of the atmosphere: changes the molecular diffsion, thermal diffusion factor, and settling velocity
rocky = True # for the surface gravity
nz = 106   # number of vertical layers
P_b = 1e7  # pressure at the bottom (dyne/cm^2)
P_t = 1e1 # pressure at the top (dyne/cm^2)
use_Kzz = True
use_moldiff = True
use_vz = False
atm_type = 'file'  # Options: 'isothermal', 'analytical', 'file', or 'vulcan_ini' 'table'
Kzz_prof = 'file' # Options: 'const','file' or 'Pfunc' (Kzz increased with P^-0.4)
K_max = 1e5        # for Kzz_prof = 'Pfunc'
K_p_lev = 0.1      # for Kzz_prof = 'Pfunc'
vz_prof = 'const'  # Options: 'const' or 'file'
gs = 1320.         # surface gravity (cm/s^2)  (HD189:2140  HD209:936)
Tiso = 1000 # only read when atm_type = 'isothermal'
# setting the parameters for the analytical T-P from (126)in Heng et al. 2014. Only reads when atm_type = 'analytical' 
# T_int, T_irr, ka_L, ka_S, beta_S, beta_L
para_warm = [120., 1500., 0.1, 0.02, 1., 1.]
para_anaTP = para_warm
const_Kzz = 1.E10 # (cm^2/s) Only reads when use_Kzz = True and Kzz_prof = 'const'
const_vz = 0 # (cm/s) Only reads when use_vz = True and vz_prof = 'const'

# frequency for updating dz and dzi due to change of mu
update_frq = 100

# ====== Setting up the boundary conditions ======
# Boundary Conditions:
use_topflux = True
use_botflux = False
use_fix_sp_bot = {} # fixed mixing ratios at the lower boundary
diff_esc = ["H","H2","H2O","He"] # species for diffusion-limit escape at TOA
max_flux = 1e16  # upper limit for the diffusion-limit fluxes

# ====== Reactions to be switched off  ======
remove_list = [] # in pairs e.g. [1,2]

# ====== Species escape list from the upper level  ======
dflux_sp = ["H","H2","H2O","He"]

# # == Condensation ======
# use_condense = False
# use_settling = False
# use_relax = ['H2O']
# humidity = 0.25 # only for water
# r_p = {'H2O_l_s': 0.01}  # particle radius in cm (1e-4 = 1 micron)
# rho_p = {'H2O_l_s': 0.9} # particle density in g cm^-3
# start_conden_time = 1e10
# condense_sp = ["H2O"]
# non_gas_sp = ["H2O_l_s"]
# fix_species = ['H2O','H2O_l_s']    # fixed the condensable species after condensation-evapoation EQ has reached
# fix_species_time = 0  # after this time to fix the condensable species
# use_ini_cold_trap = False

# == Condensation ======
use_condense = False
use_settling = False
start_conden_time = 1e10
condense_sp = []
non_gas_sp = []
fix_species = []      # fixed the condensable species after condensation-evapoation EQ has reached
fix_species_time = 0  # after this time to fix the condensable species


# ====== steady state check ======
st_factor = 0.5
conv_step = 500

# ====== Setting up numerical parameters for the ODE solver ====== 
ode_solver = 'Ros2' # case sensitive
use_print_prog = True
use_print_delta = False
print_prog_num = 500  # print the progress every x steps 
dttry = 1.E-10
trun_min = 1e2
runtime = 865000.0#1.E22
dt_min = 1.E-14
dt_max = 40#runtime*1e-5
dt_var_max = 2.
dt_var_min = 0.5
count_min = 120
count_max = int(5E6)
atol = 1.E-1 # Try decreasing this if the solutions are not stable
mtol = 1.E-22
mtol_conv = 1.E-20
pos_cut = 0
nega_cut = -1.
loss_eps = 1e-1
yconv_cri = 0.01 # for checking steady-state
slope_cri = 1.e-4
yconv_min = 0.1
flux_cri = 0.1
flux_atol = 1. # the tol for actinc flux (# photons cm-2 s-1 nm-1)
conver_ignore = [] # added 2023. to get rid off non-convergent species, e.g. HC3N without sinks 

# ====== Setting up numerical parameters for Ros2 ODE solver ====== 
rtol = 0.2             # relative tolerence for adjusting the stepsize 
post_conden_rtol = 0.1 # switched to this value after fix_species_time

# ====== Setting up for ouwtput and plotting ======
# plotting:
plot_TP = False
use_live_plot = False
use_live_flux = False
use_plot_end = False
use_plot_evo = False
use_save_movie = False
use_flux_movie = False
plot_height = False
use_PIL = False
live_plot_frq = 1
save_movie_rate = 5000 #live_plot_frq
y_time_freq = 1  #  storing data for every 'y_time_freq' step
plot_spec = ['H2O', 'H2', 'CH4', 'CO', 'CO2', 'CN', 'N2', 'OH', 'H', 'O', 'O2', 'O3','NH3']
# output:
output_humanread = False
use_shark = False
save_evolution = True   # save the evolution of chemistry (y_time and t_time) for every save_evo_frq step
save_evo_frq = 1
