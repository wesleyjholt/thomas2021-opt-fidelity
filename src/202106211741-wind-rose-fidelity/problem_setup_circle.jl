import FLOWFarm; const ff=FLOWFarm
import YAML

function split_directory_and_filename(path_to_file)
    path_to_file_segmented = split(path_to_file,"/")
    path_to_directory = ""
    if length(path_to_file_segmented) > 1
        for segment = path_to_file_segmented[1:end-1]; path_to_directory *= segment*"/"; end
    end
    filename = path_to_file[end]

    return path_to_directory, filename
end

#### SET DEFAULTS (IF FILE PATHS ARE NOT ALREADY DEFINED) ####

# farm definition
if !(@isdefined farm_definition_filename)
    path_to_farm_definition_directory = "../inputfiles/farms/random-layouts/circle-9turb-10diam/"
    farm_definition_filename = "initial-design-001.yaml"
end

# turbine definition
farm_definition_Dict = YAML.load(open(path_to_farm_definition_directory*farm_definition_filename))
try path_to_turbine_definition = farm_definition_Dict["definition"]["wind_resource"]["items"]["ref"]; catch; end
if (@isdefined path_to_turbine_definition) && path_to_wind_resource !== nothing
    path_to_turbine_definition_directory, turbine_definition_filename = split_directory_and_filename(path_to_turbine_definition)
elseif !(@isdefined turbine_definition_filename) 
    path_to_turbine_definition_directory = "../inputfiles/turbines/vestas-v80/"
    turbine_definition_filename = "VestasV80_2MW.yaml"
end

# site definition
if !(@isdefined site_definition_filename)
    path_to_site_definition_directory = "../inputfiles/sites/"
    site_definition_filename = "circle_site_withoutwindresource.yaml"
end

# wind resource
site_definition_Dict = YAML.load(open(path_to_site_definition_directory*site_definition_filename))
try path_to_wind_resource = site_definition_Dict["definition"]["wind_resource"]["items"]["ref"]; catch; end
if (@isdefined path_to_wind_resource) && path_to_wind_resource !== nothing
    path_to_wind_resource_directory, wind_resource_filename = split_directory_and_filename(path_to_wind_resource)
elseif !(@isdefined wind_resource_filename)
    path_to_wind_resource_directory = "../inputfiles/wind/wind-rose-fidelity/horns-rev/"
    wind_resource_filename = "hornsrev-windresource-020dirs-averagespeeds.yaml"
end

# flow models
if !(@isdefined flow_models_filename)
    path_to_flow_models_directory = "../inputfiles/model-sets/flow-models/"
    flow_models_filename = "Gaussian-NoLocalTI.yaml"
end

#### SET FLOWFARM PARAMETERS ####

# farm definition
farm_definition = ff.set_farm_definition_YAML(path_to_farm_definition_directory*farm_definition_filename, 
    path_to_turbine_files_directory=path_to_turbine_definition_directory, 
    turbine_definition_filenames=turbine_definition_filename)

# site definition
site_definition = ff.set_site_definition_YAML(path_to_site_definition_directory*site_definition_filename, 
    path_to_wind_resource_directory=path_to_wind_resource_directory, 
    wind_resource_filename=wind_resource_filename, 
    average_wind_speeds=true)

# combine farm definition and wind resource into a single wind farm problem
wind_farm_problem = ff.WindFarmProblem(farm_definition, site_definition)

# flow analysis models
flow_models = ff.set_flow_models_YAML(path_to_flow_models_directory*flow_models_filename)

# define optimization problem
objective = ff.MaximizeAEP(1e-6)
design_variables = [ff.TurbineXYPositions()]
constraints = [ff.TurbineSpacingConstraint(), 
               ff.BoundaryConstraint()]

# combine everything into a single wind farm optimization problem data structure
wind_farm_opt_problem = ff.WindFarmOptimizationProblem(objective, design_variables, constraints, wind_farm_problem, flow_models)
