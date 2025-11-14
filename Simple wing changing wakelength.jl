
import FLOWUnsteady as uns
import FLOWVLM as vlm
#import Plots

#run_name        = "wing-example_practice_graph"            # Name of this simulation

#save_path       = ""#run_name                  # Where to save this simulation
#paraview        = false                      # Whether to visualize with Paraview



function changing_wakelength(w_k)
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
    n               = 50                        # Number of spanwise elements per side
    r               = 10.0                      # Geometric expansion of elements
    central         = false                     # Whether expansion is central

    # ----------------- SOLVER PARAMETERS ------------------------------------------
    # Time parameters
    wakelength      = w_k*b    #was 2.75                # (m) length of wake to be resolved
    ttot            = wakelength/magVinf        # (s) total simulation time
    nsteps          = 40 #was 50,  200 originally                     # Number of time steps

    # VLM and VPM parameters
    p_per_step      = 1                         # Number of particle sheds per time step

    lambda_vpm      = 2.0                       # VPM core overlap
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
    return Cl_vec, Cd_vec, ttot
end 



Cl_vec1, Cd_vec1, ttot = changing_wakelength(1.0) 
time_vec1 = LinRange(0, ttot, length(Cl_vec1))
Cl_vec2, Cd_vec2, ttot = changing_wakelength(1.5) 
time_vec2 = LinRange(0, ttot, length(Cl_vec2))
Cl_vec3, Cd_vec3, ttot = changing_wakelength(2.0) 
time_vec3 = LinRange(0, ttot, length(Cl_vec3))
Cl_vec4, Cd_vec4, ttot = changing_wakelength(2.5) 
time_vec4 = LinRange(0, ttot, length(Cl_vec4))
Cl_vec5, Cd_vec5, ttot = changing_wakelength(2.75) 
time_vec5 = LinRange(0, ttot, length(Cl_vec5))


using PyPlot
PyPlot.matplotlib[:use]("qt5agg")

figure()
plot(time_vec1, Cl_vec1, label = "wakelength = 1")
plot(time_vec2, Cl_vec2, label = "wakelength = 1.5")
plot(time_vec3, Cl_vec3, label = "wakelength = 2.0")
plot(time_vec4, Cl_vec4, label = "wakelength = 2.5")
plot(time_vec5, Cl_vec5, label = "wakelength = 2.75")
legend()
ylabel("Cl")
xlabel("time (s)")
show()

figure()
plot(time_vec1, Cd_vec1, label = "wakelength = 1")
plot(time_vec2, Cd_vec2, label = "wakelength = 1.5")
plot(time_vec3, Cd_vec3, label = "wakelength = 2.0")
plot(time_vec4, Cd_vec4, label = "wakelength = 2.5")
plot(time_vec5, Cd_vec5, label = "wakelength = 2.75")
legend()
ylabel("Cd")
xlabel("time (s)")
show()




#=
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


=#


#=
wake_lengths = [1.0, 1.5, 2.0, 2.5, 2.75]

# Corrected Lift Coefficient (Cl) Plot
cl_plot = Plots.plot(
    time_vec1, Cl_vec1, #label="Wake length $(wake_lengths[1])",
    time_vec2, Cl_vec2, #label="Wake length $(wake_lengths[2])",
    time_vec3, Cl_vec3, #label="Wake length $(wake_lengths[3])",
    time_vec4, Cl_vec4, #label="Wake length $(wake_lengths[4])",
    time_vec5, Cl_vec5, #label="Wake length $(wake_lengths[5])",
    xlabel="Time (s)",
    ylabel="\$C_l\$",
    title="Lift Coefficient over Time",
    legend=:topright
)

# Corrected Drag Coefficient (Cd) Plot
cd_plot = Plots.plot(
    time_vec1, Cd_vec1, #label="Wake length $(wake_lengths[1])",
    time_vec2, Cd_vec2, #label="Wake length $(wake_lengths[2])",
    time_vec3, Cd_vec3, #label="Wake length $(wake_lengths[3])",
    time_vec4, Cd_vec4, #label="Wake length $(wake_lengths[4])",
    time_vec5, Cd_vec5, #label="Wake length $(wake_lengths[5])",
    xlabel="Time (s)",
    ylabel="\$C_d\$", # You had $C_l$ here, corrected to $C_d$
    title="Drag Coefficient over Time",
    legend=:topright
)

Plots.display(cl_plot)
Plots.display(cd_plot)

=#


#=
wake_lengths = [1.0, 1.5, 2.0, 2.5, 2.75]

# Create an array of (time_vec, Cl_vec, label) tuples for convenience
cl_data = [
    (time_vec1, Cl_vec1, "Wake length $(wake_lengths[1])"),
    (time_vec2, Cl_vec2, "Wake length $(wake_lengths[2])"),
    (time_vec3, Cl_vec3, "Wake length $(wake_lengths[3])"),
    (time_vec4, Cl_vec4, "Wake length $(wake_lengths[4])"),
    (time_vec5, Cl_vec5, "Wake length $(wake_lengths[5])"),
]

# Corrected Lift Coefficient (Cl) Plot
cl_plot = Plots.plot(
    # The splatting operator (...) is often needed when passing an array of args
    cl_data...; 
    xlabel="Time (s)",
    ylabel="\$C_l\$",
    title="Lift Coefficient over Time",
    legend=:topright,
    label=["Wake length $(wl)" for wl in wake_lengths] # Pass labels as a vector
)


# Create an array of (time_vec, Cd_vec, label) tuples for convenience
cd_data = [
    (time_vec1, Cd_vec1, "Wake length $(wake_lengths[1])"),
    (time_vec2, Cd_vec2, "Wake length $(wake_lengths[2])"),
    (time_vec3, Cd_vec3, "Wake length $(wake_lengths[3])"),
    (time_vec4, Cd_vec4, "Wake length $(wake_lengths[4])"),
    (time_vec5, Cd_vec5, "Wake length $(wake_lengths[5])"),
]

# Corrected Drag Coefficient (Cd) Plot
cd_plot = Plots.plot(
    cd_data...;
    xlabel="Time (s)",
    ylabel="\$C_d\$",
    title="Drag Coefficient over Time",
    legend=:topright,
    label=["Wake length $(wl)" for wl in wake_lengths] # Pass labels as a vector
)

Plots.display(cl_plot)
Plots.display(cd_plot)
=#