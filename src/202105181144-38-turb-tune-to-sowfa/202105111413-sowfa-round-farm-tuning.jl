using FLOWFarm; const ff = FLOWFarm
using DelimitedFiles
using Statistics
using DataFrames
using CSV
import PyPlot; const plt = PyPlot

function calculate_state_turbine_powers(turbine_x, turbine_y, turbine_z, rotor_diameter,
    hub_height, turbine_yaw, ct_model, generator_efficiency, cut_in_speed,
    cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set::ff.AbstractModelSet;
    rotor_sample_points_y=[0.0], rotor_sample_points_z=[0.0], hours_per_year=365.25*24.0, weighted=true)

    nturbines = length(turbine_x)

    wind_probabilities = wind_resource.wind_probabilities

    nstates = length(wind_probabilities)

    arr_type = promote_type(typeof(turbine_x[1]),typeof(turbine_y[1]),typeof(turbine_z[1]),typeof(rotor_diameter[1]),typeof(hub_height[1]),typeof(turbine_yaw[1]),
            typeof(generator_efficiency[1]),typeof(cut_in_speed[1]),typeof(cut_out_speed[1]),typeof(rated_speed[1]),typeof(rated_power[1]))
    
    turbine_powers_by_direction = zeros(arr_type,(nstates, nturbines))

    for i = 1:nstates

        rot_x, rot_y = ff.rotate_to_wind_direction(turbine_x, turbine_y, wind_resource.wind_directions[i])

        sorted_turbine_index = sortperm(rot_x)

        turbine_velocities = ff.turbine_velocities_one_direction(rot_x, rot_y, turbine_z, rotor_diameter, hub_height, turbine_yaw,
                            sorted_turbine_index, ct_model, rotor_sample_points_y, rotor_sample_points_z, wind_resource,
                            model_set, wind_farm_state_id=i, velocity_only=true)

        wt_power = ff.turbine_powers_one_direction(generator_efficiency, cut_in_speed, cut_out_speed, rated_speed,
                            rated_power, rotor_diameter, turbine_velocities, turbine_yaw, wind_resource.air_density, power_models)
        turbine_powers_by_direction[i,:] = wt_power
    end
    
    return turbine_powers_by_direction
end

# function to compare directions 
function sowfa_base_comparison(nsamplepoints=1)

    # load Niayifar LES data 
    lesfile = "../inputfiles/results-niayifar-2016/thomas2019-FinalDirectionalGeneratorPowerOutputBaseline.txt"
    sowfa_les_data = readdlm(lesfile, skipstart=0) 
    sowfa_les_data = sowfa_les_data[:,1:5]
    println(size(sowfa_les_data))

    # load FLOWFarm modelset
    include("../inputfiles/model-sets/round-farm-38-turbs-12-dirs.jl")

    # rotor swept area sample points (normalized by rotor radius)
    rotor_sample_points_y, rotor_sample_points_z = ff.rotor_sample_points(nsamplepoints)

    # run FLOWFarm with local ti
    turbine_powers_by_direction_ff = calculate_state_turbine_powers(turbine_x, turbine_y, turbine_z, rotor_diameter,
        hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
        cut_out_speed, rated_speed, rated_power, wind_resource, power_models, model_set,
        rotor_sample_points_y=rotor_sample_points_y, rotor_sample_points_z=rotor_sample_points_z)

    state_powers_ff = sum(turbine_powers_by_direction_ff, 2)

    # format sowfa data for plotting
    turbine_powers_by_direction_sowfa = zeros((nstates, nturbines))
    # state_powers_sowfa = sum()
    for i in 1:nstates
        for j in 1:nturbines
            turbine_powers_by_direction_sowfa[i, j] = sowfa_les_data[(i-1)*nturbines + j, 5]
        end
        turbine_powers_by_direction_sowfa[i, :] /= maximum(turbine_powers_by_direction_sowfa[i, :])
        turbine_powers_by_direction_ff[i, :] /= maximum(turbine_powers_by_direction_ff[i,:])
    end
    
    # println(turbine_powers_by_direction_ff)
    # println("BREAK")
    # println(turbine_powers_by_direction_sowfa)
    difference_turbines = turbine_powers_by_direction_sowfa - turbine_powers_by_direction_ff

    heatmap(difference, xticks=1:2:nturbines, yticks=1:nstates, c=:cividis)

end

function find_upstream_turbines(turbinex, turbiney, winddirection, diameter)

    # find wake count for all turbines in given wind direction 
    wake_count = ff.number_of_wakes_iec(turbinex, turbiney, winddirection, diameter)

    # return unwaked turbines 
    return collect(1:length(turbinex))[wake_count .== 0]

end

function tune_flowfarm_to_sowfa(;case="low-ti")

    # include functions from other files 
    include("../202105111413-SOWFA-comparison/202105111413-sowfa-round-farm-comparison.jl")

    # load sowfa data
    turbine_powers_by_direction_sowfa, _ = get_data(journal=true, case=case)

    # set path to model set 
    model_set_file = "../inputfiles/model-sets/round-farm-38-turbs-12-dirs-$(case).jl"

    # load flowfarm set up
    include(model_set_file)

    # get the number of wind directions and wind turbines 
    ndirections = length(winddirections)
    nturbines = length(turbine_x)

    # initialize arrays for tuned wind speed and TI values
    opt_speeds = zeros(ndirections)
    opt_tis = zeros(ndirections)

    # tune inflow wind speed for each direction using front turbines
    for i = 1:ndirections
        # find upstream turbines 
        upstream_turbines = find_upstream_turbines(turbine_x, turbine_y, winddirections[i], rotor_diameter[1])
        # println("dir = $(winddirections[i]), upstream turbines = $(upstream_turbines)")

        # get sowfa data for turbines of interest 
        upstream_turbines_power = turbine_powers_by_direction_sowfa[i, upstream_turbines]

        # generate appropriate objective function 
        obj_func_windspeed_local(x1, x2) = obj_func_windspeed(x1, x2, ti=0.05589339140297106, case="low-ti", wd=i)

        # run least squared fit to the sowfa data 
        fit = curve_fit(obj_func_windspeed_local, upstream_turbines, upstream_turbines_power, [8.0])

        # store result
        opt_speeds[i] = fit.param[1]
        println("dir: $(round(winddirections[i]*180.0/pi,digits=0)), opt speed: $(opt_speeds[i])")

    end

    # tune local TI for each direction using all turbines 
    for i = 1:ndirections 
        # select all turbines 
        downstream_turbines = collect(1:nturbines)

        # get sowfa data for turbines of interest 
        downstream_turbines_power = turbine_powers_by_direction_sowfa[i, downstream_turbines]

        # generate appropriate objective function
        obj_func_ti_local(x1, x2) = obj_func_ti(x1, x2, windspeed=opt_speeds[i], case="low-ti", wd=i)

        # run least squared fit to the sowfa data 
        fit = curve_fit(obj_func_ti_local, downstream_turbines, downstream_turbines_power, [0.04])

        # store result
        opt_tis[i] = fit.param[1]
        println("dir: $(round(winddirections[i]*180.0/pi, digits=0)), opt ti: $(opt_tis[i])")

    end

    # save results 
    df = DataFrame(dir=(round.(winddirections*180.0./pi, digits=0)), speed=opt_speeds, ti=opt_tis)
    CSV.write("tuned-parameters-$(case).csv", df)

end