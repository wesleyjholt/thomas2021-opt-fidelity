import FLOWFarm; const ff=FLOWFarm
using DelimitedFiles
using FLOWMath
using Distributions

function generate_all_wind_resource_files()
    # create Horns Rev wind resource files
    ndirs_horns_rev = [8,10,11,12,13,15,16,18,20,24,30,36,45,60,90,180,360]
    nspeeds_horns_rev = [1,20]
    for ndirs = ndirs_horns_rev
        for nspeeds = nspeeds_horns_rev
            wind_resource = HornsRevWindResource(ndirs, nspeeds)
            resample_wind_resource(wind_resource)
        end
    end
end

struct NantucketWindResource
    ndirs
    nspeeds
    measurement_height
    air_density
    ambient_ti
    shear_exponent
end
NantucketWindResource(ndirs, nspeeds) = NantucketWindResource(ndirs, nspeeds, 70.0, 1.225, 0.108, 0.31)

struct HornsRevWindResource
    ndirs
    nspeeds
    measurement_height
    air_density
    ambient_ti
    shear_exponent
end
HornsRevWindResource(ndirs, nspeeds) = HornsRevWindResource(ndirs, nspeeds, 70.0, 1.225, 0.108, 0.31)

"""
    resample_wind_resource(wind_resource::NantucketWindResource)

Resamples the Nantucket wind resource for the specified number of directions. (average speeds)

# Arguments
-`wind_resource::NantucketWindResource`: container holding parameters for the Nantucket wind rose
"""
function resample_wind_resource(wind_resource::NantucketWindResource)

    # load original wind resource
    input_wind_resource_file_path = "../data/wind-resource-files/nantucket/nantucket_windresource_036dirs_averagespeeds.txt"
    wind_resource_data = readdlm(input_wind_resource_file_path, skipstart=1)
    directions_orig = wind_resource_data[:,1]
    speeds_orig = wind_resource_data[:,2]
    direction_probabilities_orig = wind_resource_data[:,3]

    # resample discretized wind resource
    directions, direction_probabilities, average_speeds = resample_discretized_wind_resource_average_speeds(directions_orig, direction_probabilities_orig, speeds_orig, wind_resource.ndirs)

    # write to a yaml file
    output_wind_resource_filepath = "wind_resource_files/wind-resource-files/nantucket/nantucket_windresource_$(lpad(wind_resource.ndirs, 3, "0"))dirs_averagespeeds.yaml"
    title = Symbol("Nantucket Wind Rose, $(wind_resource.ndirs) wind directions, average wind speeds")
    description = Symbol("Wind resource conditions")
    ff.write_wind_resource_YAML(output_wind_resource_filepath, directions, direction_probabilities, average_speeds, wind_resource.measurement_height, wind_resource.air_density, wind_resource.ambient_ti, wind_resource.shear_exponent,
        title=title, description=description)
end

"""
    resample_wind_resource(wind_resource::HornsRevWindResource)

Resamples the Horns Rev 1 wind resource for the specified number of directions and speeds.
If only one speed is asked for, then a wind rose with average speeds for direction is created.

# Arguments
- `wind_resource::HornsRevWindResource`: container holding parameters for the Horns Rev 1 wind rose
"""
function resample_wind_resource(wind_resource::HornsRevWindResource)
    # set original values for direction and speed joint distribution (from paper by Ju Feng and Wen Zhong Shen, 2015, "Modelling Wind for Wind Farm Layout Optimization Using Joint Distribution of Wind Speed and Wind Direction")
    directions_orig = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]
    direction_probabilities_orig = [0.0482, 0.0406, 0.0359, 0.0527, 0.0912, 0.0697, 0.0917, 0.1184, 0.1241, 0.1134, 0.117, 0.0969]
    speed_weibull_params_orig = [2.09 8.89; 2.13 9.27; 2.29 8.23; 2.30 9.78; 2.67 11.64; 2.45 11.03; 2.51 11.50; 2.40 11.92; 2.35 11.49; 2.27 11.08; 2.24 11.34; 2.19 10.76]

    # resample discretized wind rose
    if wind_resource.nspeeds != 1
        # wind speed discretized distributions
        average_speeds = false
        directions, direction_probabilities, speeds, speed_probabilities = resample_discretized_wind_resource_multiple_speeds_weibull(directions_orig, direction_probabilities_orig, speed_weibull_params_orig, wind_resource.ndirs, wind_resource.nspeeds)
        output_wind_resource_filepath = "wind_resource_files/horns-rev/hornsrev-windresource-$(lpad(wind_resource.ndirs,3,"0"))dirs-$(lpad(wind_resource.nspeeds,2,"0"))speeds.yaml"
    else
        # average wind speeds
        average_speeds = true
        directions, direction_probabilities, speeds = resample_discretized_wind_resource_average_speeds_weibull(directions_orig, direction_probabilities_orig, speed_weibull_params_orig, wind_resource.ndirs)
        speed_probabilities = []
        output_wind_resource_filepath = "wind_resource_files/horns-rev/hornsrev-windresource-$(lpad(wind_resource.ndirs,3,"0"))dirs-averagespeeds.yaml"
    end

    # write to a yaml file
    title = Symbol("Horns Rev 1 Wind Rose, $(wind_resource.ndirs) wind directions, $(wind_resource.nspeeds) wind speeds")
    description = "Wind resource conditions, using direction and speed distribution from the paper by Ju Feng and Wen Zhong Shen: \nModelling Wind for Wind Farm Layout Optimization Using Joint Distribution of Wind Speed and Wind Direction \nhttps://backend.orbit.dtu.dk/ws/portalfiles/portal/110639016/Modelling_Wind_for_Wind_Farm_Layout.pdf"
    ff.write_wind_resource_YAML(output_wind_resource_filepath, directions, direction_probabilities, speeds, speed_probabilities, wind_resource.measurement_height, wind_resource.air_density, wind_resource.ambient_ti, wind_resource.shear_exponent,
        title=title, description=description)
end

function resample_discretized_wind_resource_average_speeds(directions_orig, direction_probabilities_orig, speeds_orig, ndirs)

    # create cyclical data
    directions_orig_cyclical = [directions_orig[end-1:end] .- 360.0; directions_orig; directions_orig[1:2] .+ 360.0]
    direction_probabilities_orig_cyclical = [direction_probabilities_orig[end-1:end]; direction_probabilities_orig; direction_probabilities_orig[1:2]]
    speeds_orig_cyclical = [speeds_orig[end-1:end]; speeds_orig; speeds_orig[1:2]]

    # get direction vector
    directions = range(0, stop=360*(ndirs-1)/ndirs, length=ndirs)

    # use akima spline to interpolate the pmf to the specified number of wind directions
    direction_probabilities = akima(directions_orig_cyclical, direction_probabilities_orig_cyclical, directions)
    direction_probabilities /= sum(direction_probabilities)

    # use akima spline to interpolate speeds
    speeds = akima(directions_orig_cyclical, speeds_orig_cyclical, directions)

    return directions, direction_probabilities, speeds
end

function resample_discretized_wind_resource_average_speeds_weibull(directions_orig, direction_probabilities_orig, speed_weibull_params_orig, ndirs)

    # create cyclical data
    directions_orig_cyclical = [directions_orig[end-1:end] .- 360.0; directions_orig; directions_orig[1:2] .+ 360.0]
    direction_probabilities_orig_cyclical = [direction_probabilities_orig[end-1:end]; direction_probabilities_orig; direction_probabilities_orig[1:2]]
    speed_weibull_params_orig_cyclical = [speed_weibull_params_orig[end-1:end,:]; speed_weibull_params_orig; speed_weibull_params_orig[1:2,:]]

    # get direction vector
    directions = range(0, stop=360*(ndirs-1)/ndirs, length=ndirs)

    # use akima spline to interpolate the direction probabilities and normalize the interpolated probabilities to make them an actual probability mass function 
    direction_probabilities = akima(directions_orig_cyclical, direction_probabilities_orig_cyclical, directions)
    direction_probabilities /= sum(direction_probabilities)

    # use akima spline to interpolate the speed weibull parameters to the specified number of wind directions
    speed_weibull_params = zeros(ndirs, 2)
    speed_weibull_params[:,1] = akima(directions_orig_cyclical, speed_weibull_params_orig_cyclical[:,1], directions)
    speed_weibull_params[:,2] = akima(directions_orig_cyclical, speed_weibull_params_orig_cyclical[:,2], directions)

    # get average speeds
    average_speeds = zeros(ndirs)
    n_integration_bins = 200
    delta_speed_integration = 25.0/n_integration_bins
    speed_integration_bins = range(0.0 + delta_speed_integration/2, stop = 25.0 - delta_speed_integration/2, step = delta_speed_integration)
    for i = 1:ndirs
        dist = Weibull(speed_weibull_params[i,:]...)
        average_speeds[i] = sum(delta_speed_integration .* speed_integration_bins .* pdf.(dist, speed_integration_bins))
    end

    return directions, direction_probabilities, average_speeds
end

function resample_discretized_wind_resource_multiple_speeds_weibull(directions_orig, direction_probabilities_orig, speed_weibull_params_orig, ndirs, nspeeds)

    # create cyclical data
    directions_orig_cyclical = [directions_orig[end-1:end] .- 360.0; directions_orig; directions_orig[1:2] .+ 360.0]
    direction_probabilities_orig_cyclical = [direction_probabilities_orig[end-1:end]; direction_probabilities_orig; direction_probabilities_orig[1:2]]
    speed_weibull_params_orig_cyclical = [speed_weibull_params_orig[end-1:end,:]; speed_weibull_params_orig; speed_weibull_params_orig[1:2,:]]

    # get direction and speed vectors
    directions = range(0, stop=360*(ndirs-1)/ndirs, length=ndirs)
    delta_speed = 25.0/nspeeds
    speeds = range(0.0 + delta_speed/2, stop = 25.0 - delta_speed/2, step = delta_speed)
    speed_bin_edges = range(0.0, stop=25.0, length=nspeeds+1)

    # use akima spline to interpolate the probabilities and weibull parameters to the specified number of wind directions
    direction_probabilities = akima(directions_orig_cyclical, direction_probabilities_orig_cyclical, directions)
    speed_weibull_params = zeros(ndirs, 2)
    speed_weibull_params[:,1] = akima(directions_orig_cyclical, speed_weibull_params_orig_cyclical[:,1], directions)
    speed_weibull_params[:,2] = akima(directions_orig_cyclical, speed_weibull_params_orig_cyclical[:,2], directions)

    # normalize the interpolated probabilities to make them an actual probability mass function 
    direction_probabilities /= sum(direction_probabilities)

    # get the speed probabilities for each direction
    speed_probabilities = fill(zeros(length(speeds)), ndirs)
    for i = 1:ndirs
        dist = Weibull(speed_weibull_params[i,:]...)
        speed_probabilities[i] = cdf.(dist, speed_bin_edges[2:end]) - cdf.(dist, speed_bin_edges[1:end-1])
        # normalize
        speed_probabilities[i] ./= sum(speed_probabilities[i])
    end

    return directions, direction_probabilities, speeds, speed_probabilities
end
