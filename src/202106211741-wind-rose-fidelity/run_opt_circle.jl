using SNOW
using Snopt
using PyPlot
using OrderedCollections
using Distributed
using ClusterManagers

plotresults = true
verbose = true

##########################################################
# SET UP FLOWFARM MODELS
##########################################################

### parse input arguments ###

ntasks = parse(Int, ENV["SLURM_NTASKS"])
if ntasks>1
    addprocs(SlurmManager(ntasks - 1))
end

if length(ARGS)<1
    # defaults
    append!(ARGS, ["001", "020", "01", "Gaussian", "9", "10"])
end

layout_number = lpad(ARGS[1], 3, "0")
ndirs = lpad(ARGS[2], 3, "0")
nspeeds = lpad(ARGS[3], 2, "0")
wake_model = ARGS[4]
nturbines = parse(Int64,ARGS[5])
diameter_spacing = parse(Int,ARGS[6])

### set input file path/names ###

# farm name
@everywhere farm_name = "circle-$($nturbines)turb-$($diameter_spacing)diam"
# farm definition
@everywhere path_to_farm_definition_directory = "../inputfiles/farms/random-layouts/$($farm_name)/"
@everywhere farm_definition_filename = "initial-design-$($layout_number).yaml"
# wind resource
@everywhere path_to_wind_resource_directory = "../inputfiles/wind/wind-rose-fidelity/horns-rev/"
if parse(Int64,nspeeds)==1
    @everywhere wind_resource_filename = "hornsrev-windresource-$($ndirs)dirs-averagespeeds.yaml"
else
    @everywhere wind_resource_filename = "hornsrev-windresource-$($ndirs)dirs-$($nspeeds)speeds.yaml"
end
# flow models
@everywhere path_to_flow_models_directory = "../inputfiles/model-sets/flow-models/"
@everywhere flow_models_filename = "$($wake_model)-NoLocalTI.yaml"

### import problem setup file ###

@everywhere include("problem_setup_circle.jl")


##########################################################
# RUN OPTIMIZATION
##########################################################

### set up optimization ###

# get SNOW arguments
@everywhere func!, x0, ng, lx, ux, lg, ug = ff.get_SNOW_optimization_inputs(wind_farm_opt_problem)

# test objective and constraint functions
nx = length(x0)
g = zeros(ng)
df = zeros(nx)
dg = zeros(ng,nx)
f = func!(g,df,dg,x0)

# separate initial design variable values into x and y
turbine_x = copy(x0[1:nturbines])
turbine_y = copy(x0[nturbines+1:end])
rotor_diameter = wind_farm_opt_problem.wind_farm_problem.farm_description.turbine_definitions[1].rotor_diameter

if plotresults
    plot(0,0)
    # add initial turbine location to plot
    for i = 1:nturbines
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C0"))
    end
end

# set WEC options
wec_values = [3.0, 2.6, 2.2, 1.8, 1.4, 1.0, 1.0]
with_TI = [false, false, false, false, false, false, true]
opt_tolerance = [1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-6, 1e-6]

# set SNOPT options
snopt_options = Dict{String, Any}()
snopt_options["Derivative option"] = 1
snopt_options["Verify level"] = -1
snopt_options["Major iteration limit"] = 50000

### run optimizations, once at each wec value ###

startx = []
starty = []
optx = []
opty = []
project_directory = pwd()

for i = 1:length(wec_values)

    # set up SNOW/SNOPT options
    log_directory = "results/log/$farm_name/$(ndirs)dirs/$(nspeeds)speeds/$(wake_model)/"
    if !isdir(log_directory); mkpath(log_directory); end
    cd(log_directory)
    snopt_options["Summary file"] = "snopt_summary_wec$(round(wec_values[i],digits=1)).out"
    snopt_options["Print file"] = "snopt_print_wec$(round(wec_values[i],digits=1)).out"
    snopt_options["Major optimality tolerance"] = opt_tolerance[i]
    solver = SNOPT(options=snopt_options)
    options = Options(;solver, derivatives=UserDeriv())

    # set wec value
    wind_farm_opt_problem.flow_models.wake_deficit_model.wec_factor[1] = wec_values[i]
    
    # initialize and save starting point
    x0_i = [turbine_x; turbine_y]
    push!(startx,turbine_x)
    push!(starty,turbine_y)
    if verbose
        println("\nnow running with WEC = ", wind_farm_opt_problem.flow_models.wake_deficit_model.wec_factor[1])
        println("initial turbine coordinates: ", x0_i[1:5])
    end

    # run optimization
    xopt_i, fopt_i, info_i, out_i = minimize(func!, x0_i, ng, lx, ux, lg, ug, options)
    
    # extract and save optimized point
    turbine_x[:] = xopt_i[1:nturbines]
    turbine_y[:] = xopt_i[nturbines+1:end]
    push!(optx,turbine_x)
    push!(opty,turbine_y)

    if verbose
        println("final turbine coordinates: ", xopt_i[1:5])
    end

    # navigate back to project directory
    cd(project_directory)

end



##########################################################
# CALCULATE AND SAVE FINAL RESULTS
##########################################################

### define support functions ###

function save_final_aep!(annual_energy_production_Dict, dict_entry_name, wind_farm_problem, flow_models)
    
    # calculate AEP
    aep, state_aep = ff.calculate_aep!(wind_farm_problem, flow_models, return_state_aeps=true)

    # create binned AEP vector
    wd = wind_farm_problem.site_description.wind_resource.wind_directions
    wd_bins = unique(wd)
    ws = wind_farm_problem.site_description.wind_resource.wind_speeds
    ws_bins_or_avgs = unique(ws)
    if sum([ws[i] in [ws[1:i-1]; ws[i+1:end]] for i=1:length(ws)-1])==0
        ws_dist_type = :averaged
    else
        ws_dist_type = :binned
    end
    aep_binned = [state_aep[wd.==wd_bins[i]] for i=1:length(wd_bins)]

    # save AEP to Dict
    annual_energy_production_Dict[dict_entry_name] = OrderedDict{String, Any}()
    annual_energy_production_Dict[dict_entry_name]["units"] = :W
    annual_energy_production_Dict[dict_entry_name]["default"] = aep
    annual_energy_production_Dict[dict_entry_name]["wind_direction"] = OrderedDict{String, Any}()
    annual_energy_production_Dict[dict_entry_name]["wind_direction"]["distribution_type"] = :binned
    annual_energy_production_Dict[dict_entry_name]["wind_direction"]["units"] = :radians
    annual_energy_production_Dict[dict_entry_name]["wind_direction"]["bins"] = wd_bins
    annual_energy_production_Dict[dict_entry_name]["wind_speed"] = OrderedDict{String, Any}()
    annual_energy_production_Dict[dict_entry_name]["wind_speed"]["distribution_type"] = ws_dist_type
    annual_energy_production_Dict[dict_entry_name]["wind_speed"]["units"] = :(m/s)
    if ws_dist_type==:averaged
        annual_energy_production_Dict[dict_entry_name]["wind_speed"]["averages_by_direction"] = ws_bins_or_avgs
    else
        annual_energy_production_Dict[dict_entry_name]["wind_speed"]["bins"] = ws_bins_or_avgs
    end
    annual_energy_production_Dict[dict_entry_name]["binned"] = aep_binned

end

### split design variables vector into x and y ###

turbine_x = copy(optx[end])
turbine_y = copy(opty[end])
turbine_xy = [[turbine_x[i],turbine_y[i]] for i=1:nturbines]

### define additional FLOWFarm structs ###

# wind farm problem (final)
wind_farm_problem_final = deepcopy(wind_farm_problem)
wind_farm_problem_final.farm_description.turbine_x[:] = turbine_x
wind_farm_problem_final.farm_description.turbine_y[:] = turbine_y

# wind farm problem with higher fidelity (initial and final)
site_definition_360dirs_20speeds = ff.set_site_definition_YAML(path_to_site_definition_directory*site_definition_filename, 
    path_to_wind_resource_directory="../inputfiles/wind/wind-rose-fidelity/horns-rev/", 
    wind_resource_filename="hornsrev-windresource-360dirs-20speeds.yaml")
wfp_360dirs_20speeds_final = ff.WindFarmProblem(wind_farm_problem_final.farm_description, site_definition_360dirs_20speeds)

# flow models (Gaussian and Jensen)
fm_gaussianwakemodel = ff.set_flow_models_YAML("../inputfiles/model-sets/flow-models/Gaussian-MaxTI.yaml")
fm_jensenwakemodel = ff.set_flow_models_YAML("../inputfiles/model-sets/flow-models/JensenCosine-MaxTI.yaml")

### calculate AEPs and store in structs ###
turbine_z = zeros(nturbines)
turbine_definition_ids = ones(Int64,nturbines)
turbine_definition_references = "VestasV80_2MW.yaml"
farm_design_Dict_final = ff._create_farm_definition_YAML_data_structure(turbine_xy, turbine_z, turbine_definition_ids, turbine_definition_references;
    title="Circular Wind Farm - Final Design")
annual_energy_production_Dict = farm_design_Dict_final["definitions"]["plant_energy"]["properties"]["annual_energy_production"] = OrderedDict{String, Any}()

## reduced fidelity 
    # Gaussian wake model
    wfp_gaussianwakemodel_reduceddirs_reducedspeeds = deepcopy(wind_farm_problem)
    save_final_aep!(annual_energy_production_Dict, "gaussianwakemodel_reduceddirs_reducedspeeds", wfp_gaussianwakemodel_reduceddirs_reducedspeeds, fm_gaussianwakemodel)
    # Jensen wake model
    wfp_jensenwakemodel_reduceddirs_reducedspeeds = deepcopy(wind_farm_problem)
    save_final_aep!(annual_energy_production_Dict, "jensenwakemodel_reduceddirs_reducedspeeds", wfp_jensenwakemodel_reduceddirs_reducedspeeds, fm_jensenwakemodel)

## higher fidelity
    # Gaussian wake model
    wfp_gaussianwakemodel_360dirs_20speeds = deepcopy(wfp_360dirs_20speeds_final)
    save_final_aep!(annual_energy_production_Dict, "gaussianwakemodel_360dirs_20speeds", wfp_gaussianwakemodel_360dirs_20speeds, fm_gaussianwakemodel)
    # Jensen wake model
    wfp_jensenwakemodel_360dirs_20speeds = deepcopy(wfp_360dirs_20speeds_final)
    save_final_aep!(annual_energy_production_Dict, "jensenwakemodel_360dirs_20speeds", wfp_jensenwakemodel_360dirs_20speeds, fm_jensenwakemodel)

### save results ###

path_to_final_designs_directory = "results/final-designs/$farm_name/$(ndirs)dirs/$(nspeeds)speeds/$(wake_model)/"
if !isdir(path_to_final_designs_directory); mkpath(path_to_final_designs_directory); end
final_design_filename = "final-design-$layout_number.yaml"
YAML.write_file(path_to_final_designs_directory*final_design_filename, farm_design_Dict_final)


##########################################################
# PLOT FINAL RESULTS
##########################################################

if plotresults
    # add final turbine locations to plot
    for i = 1:nturbines
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1", linestyle="--")) 
    end
    
    # add wind farm boundary to plot
    boundary_radius = wind_farm_opt_problem.wind_farm_problem.site_description.boundaries[1].radius
    boundary_center = wind_farm_opt_problem.wind_farm_problem.site_description.boundaries[1].center
    plt.gcf().gca().add_artist(plt.Circle((boundary_center[1],boundary_center[2]), boundary_radius, fill=false,color="C2"))

    # set up and show plot
    axis("square")
    xlim(-boundary_radius-200,boundary_radius+200)
    ylim(-boundary_radius-200,boundary_radius+200)
    plt.show()
end
