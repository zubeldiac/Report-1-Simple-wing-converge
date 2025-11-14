

import FLOWUnsteady as uns
import FLOWVLM as vlm

run_name        = "wing-example-convergence"            # Name of this simulation

save_path       = ""#run_name                  # Where to save this simulation



# ----------------- SIMULATION PARAMETERS --------------------------------------
AOA             = 4.2                       # (deg) angle of attack
magVinf         = 49.7                      # (m/s) freestream velocity
rho             = 0.93                      # (kg/m^3) air density
qinf            = 0.5*rho*magVinf^2         # (Pa) static pressure

Vinf(X, t)      = magVinf*[cosd(AOA), 0.0, sind(AOA)]  # Freestream function


# ----------------- GEOMETRY PARAMETERS ----------------------------------------
b               = 2.489                     # (m) span length
ar              = 5.0                       # Aspect ratio b/c_tip
tr              = 1.0                       # Taper ratio c_tip/c_root
twist_root      = 0.0                       # (deg) twist at root
twist_tip       = 0.0                       # (deg) twist at tip
lambda          = 45.0                      # (deg) sweep
gamma           = 0.0                       # (deg) dihedral

# Discretization
n               = 40                        #was 50 start with 10 of spanwise elements per                          #CHANGE THIS first
r               = 10.0                      # Geometric expansion of elements
central         = false                     # Whether expansion is central


# ----------------- SOLVER PARAMETERS ------------------------------------------
# Time parameters
wakelength      = 2.75*b                    # (m) length of wake to be resolved
ttot            = wakelength/magVinf        # (s) total simulation time
nsteps          = 100 #was50 start with 10 ,  was 200 originally                     # Number of time steps            #CHANGE THIS second

# VLM and VPM parameters
p_per_step      = 6 #was 1 Start with 1  when changing p_per_step manually            # Number of particle sheds per time step  #Change this third
lambda_vpm      = 2.0                       # VPM core overlap
#p_per_step = ceil(Int, (lambda_vpm*magVinf*(ttot/nsteps))/(b/n)) #when changing nsteps and n. This edits it automatically to help with overlapping
sigma_vpm_overwrite = lambda_vpm * magVinf * (ttot/nsteps)/p_per_step # Smoothing core size
sigma_vlm_solver= -1                        # VLM-on-VLM smoothing radius (deactivated with <0)
sigma_vlm_surf  = 0.05*b                    # VLM-on-VPM smoothing radius

shed_starting   = true                      # Whether to shed starting vortex
vlm_rlx         = 0.7                       # VLM relaxation


# ----------------- 1) VEHICLE DEFINITION --------------------------------------
println("Generating geometry...")

# Generate wing
wing = vlm.simpleWing(b, ar, tr, twist_root, lambda, gamma;
                        twist_tip=twist_tip, n=n, r=r, central=central);


println("Generating vehicle...")

# Generate vehicle
system = vlm.WingSystem()                   # System of all FLOWVLM objects
vlm.addwing(system, "Wing", wing)

vlm_system = system;                        # System solved through VLM solver
wake_system = system;                       # System that will shed a VPM wake

vehicle = uns.VLMVehicle(   system;
                            vlm_system=vlm_system,
                            wake_system=wake_system
                         );


# ------------- 2) MANEUVER DEFINITION -----------------------------------------

Vvehicle(t) = zeros(3)                      # Translational velocity of vehicle over time
anglevehicle(t) = zeros(3)                  # Angle of the vehicle over time

angle = ()                                  # Angle of each tilting system (none)
RPM = ()                                    # RPM of each rotor system (none)

maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)


# ------------- 3) SIMULATION DEFINITION ---------------------------------------

Vref = 0.0                                  # Reference velocity to scale maneuver by
RPMref = 0.0                                # Reference RPM to scale maneuver by
Vinit = Vref*Vvehicle(0)                    # Initial vehicle velocity
Winit = pi/180*(anglevehicle(1e-6) - anglevehicle(0))/(1e-6*ttot)  # Initial angular velocity

                                            # Maximum number of particles
max_particles = (nsteps+1)*(vlm.get_m(vehicle.vlm_system)*(p_per_step+1) + p_per_step)

simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot;
                                                    Vinit=Vinit, Winit=Winit);




# ------------- 4) MONITORS DEFINITIONS ----------------------------------------


#Added by me/////////////////////////////
figs, figaxs = [], []

time_vec = Float64[]
Cl_vec = Float64[]
Cd_vec = Float64[]
#/////////////////////////////////////////


calc_aerodynamicforce_fun = uns.generate_calc_aerodynamicforce(;
                                    add_parasiticdrag=true,
                                    add_skinfriction=false,
                                    airfoilpolar="xf-rae101-il-1000000.csv"
                                    )

D_dir = [cosd(AOA), 0.0, sind(AOA)]         # Direction of drag
L_dir = uns.cross(D_dir, [0,1,0])           # Direction of lift

figs, figaxs = [], []                       # Figures generated by monitor

# Generate wing monitor
monitor_wing = uns.generate_monitor_wing(wing, Vinf, b, ar,
                                            rho, qinf, nsteps;
                                            calc_aerodynamicforce_fun=calc_aerodynamicforce_fun,
                                            L_dir=L_dir,
                                            D_dir=D_dir,
                                            out_CLwing=Cl_vec, #Added by me
                                            out_CDwing=Cd_vec, #ADDED BY ME
                                            #out_figs=figs,
                                            #out_figaxs=figaxs,
                                            save_path="",#was save_path
                                            run_name="", #was run_name
                                            figname="", #was wing monitor
                                            disp_plot = false, #Was true, this change hopefully gets rid of pngs
                                            );




# ------------- 5) RUN SIMULATION ----------------------------------------------
println("Running simulation...")

uns.run_simulation(simulation, nsteps;
                    # ----- SIMULATION OPTIONS -------------
                    Vinf=Vinf,
                    rho=rho,
                    # ----- SOLVERS OPTIONS ----------------
                    p_per_step=p_per_step,
                    max_particles=max_particles,
                    sigma_vlm_solver=sigma_vlm_solver,
                    sigma_vlm_surf=sigma_vlm_surf,
                    sigma_rotor_surf=sigma_vlm_surf,
                    sigma_vpm_overwrite=sigma_vpm_overwrite,
                    shed_starting=shed_starting,
                    vlm_rlx=vlm_rlx,
                    extra_runtime_function=monitor_wing,
                    # ----- OUTPUT OPTIONS ------------------
                    #save_path=save_path,
                    #run_name=run_name
                    );




import Plots

# Construct the time vector: total time / number of steps, times the index
time_vec = LinRange(0, ttot, length(Cl_vec)) #Make sure it keeps into account the changing number of N steps when changing nsteps

# 1. Create the Lift Coefficient (Cl) Plot (using Plots.plot)
cl_plot = Plots.plot(time_vec, Cl_vec,
                     label="Lift Coefficient (\$C_l\$)",
                     xlabel="Time (s)",
                     ylabel="\$C_l\$",
                     title="Lift Coefficient over Time",
                     #legend=:topright
                    )

# 2. Create the Drag Coefficient (Cd) Plot (using Plots.plot again)
cd_plot = Plots.plot(time_vec, Cd_vec,
                     label="Drag Coefficient (\$C_d\$)",
                     xlabel="Time (s)",
                     ylabel="\$C_d\$",
                     title="Drag Coefficient over Time",
                     legend=:topright
                    )

# Display both plots (this step is optional but often needed in non-interactive environments)
Plots.display(cl_plot)
Plots.display(cd_plot)