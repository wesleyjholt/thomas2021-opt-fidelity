import YAML
import CSV
using PyPlot
using DataFrames
using Statistics

function create_windrosefidelity_results_dataframe(farm_directory,ndirs_vec,nspeeds_vec,wake_model_vec,simulation_number_vec)

    results = DataFrame(
                ndirs = Int[], 
                nspeeds = Int[], 
                wake_model = String[], 
                simulation_number = Int[], 
                aep_Gaussian = Float64[], 
                aep_JensenCosine = Float64[], 
                aep_Gaussian_binned = Vector{Vector{Float64}}[], 
                aep_JensenCosine_binned = Vector{Vector{Float64}}[]
              )

    # wind directions
    for ndirs in ndirs_vec

        # wind speeds
        for nspeeds in nspeeds_vec
            
            # wake model
            for wake_model in wake_model_vec

                # simulation ID
                for simulation_number in simulation_number_vec
                    
                    path_to_results_file = "$farm_directory/$(lpad(ndirs,3,"0"))dirs/$(lpad(nspeeds,2,"0"))speeds/$wake_model/final-design-$(lpad(simulation_number,3,"0")).yaml"
                    if isfile(path_to_results_file)
                        dat = YAML.load_file(path_to_results_file)
                        aep_data = dat["definitions"]["plant_energy"]["properties"]["annual_energy_production"]
                        push!(results, 
                           (ndirs = ndirs, 
                            nspeeds = nspeeds, 
                            wake_model = wake_model,
                            simulation_number = simulation_number,
                            aep_Gaussian = aep_data["gaussianwakemodel_360dirs_20speeds"]["default"],
                            aep_JensenCosine = aep_data["jensenwakemodel_360dirs_20speeds"]["default"],
                            aep_Gaussian_binned = aep_data["gaussianwakemodel_360dirs_20speeds"]["binned"],
                            aep_JensenCosine_binned = aep_data["jensenwakemodel_360dirs_20speeds"]["binned"]
                           )
                        )
                    else
                        push!(results, 
                           (ndirs = ndirs, 
                            nspeeds = nspeeds, 
                            wake_model = wake_model,
                            simulation_number = simulation_number,
                            aep_Gaussian = NaN,
                            aep_JensenCosine = NaN,
                            aep_Gaussian_binned = [fill(NaN,20) for i=1:360],
                            aep_JensenCosine_binned = [fill(NaN,20) for i=1:360]
                           )
                        )
                    end

                end
            end
        end
    end

    return results
end

function is_outlier(points;thresh=3.5)
    med = median(points)
    diff = abs.(points .- med)
    med_abs_deviation = median(diff)
    deviation = diff / med_abs_deviation

    return deviation .> thresh
end

function plot_histogram_comparing_wake_models(results; ndirs, nspeeds, wake_models=["Gaussian","JensenCosine"])
    # create overlayed histograms comparing the AEP when optimized using different wake models

    plt.figure()
    for wake_model in wake_models
        subset_indices = (results[:,:ndirs].==ndirs) .& (results[:,:nspeeds].==nspeeds) .& (results[:,:wake_model].==wake_model)
        histogram_data = eval(Meta.parse("results[$subset_indices,:aep_Gaussian]"))
        histogram_data = histogram_data[.!isnan.(histogram_data)]
        if !isempty(histogram_data)
            plt.hist(histogram_data[.!is_outlier(histogram_data)], alpha=0.5)
        end
    end
end

# """
#     plot_aeps(data_file_names, labels, plot_type::ConfidenceIntervalScatterPlot; show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_boxplots.png", title="AEP Values for Optimized Layouts")
# """
# function plot_aeps(data_file_names::Array{String,1}, x_values, plot_type::ConfidenceIntervalScatterPlot; max_aep = 1.0, fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_confidence_interval_plot_with_scatter.png", _title="AEP Values for Optimized Layouts", _xlabel="")

#     # get figure and axes handles
#     fig, ax = get_fig_ax_handles(fig_handle, ax_handle)
#     # get the AEP values from the files
#     aeps_raw = get_aep_values_from_file_names(data_file_names)
#     # filter out NaN values
#     aeps = filter_aeps(aeps_raw)/max_aep
#     # create scatter plots
#     create_aep_scatter_plots(x_values, aeps)
#     # create mean line with confidence interval
#     create_confidence_interval_plot(x_values, aeps)
#     # add plot labels
#     title(_title)
#     xlabel(_xlabel)
#     if max_aep == 1.0
#         ylabel("AEP (MWh)")
#     else
#         ylabel("Normalized AEP")
#         ylim([0.979, 1.0])
#         yticks([0.98, 0.985, 0.99, 0.995, 1.0])
#     end
#     # save figure
#     if save_fig
#         mkpath(path_to_fig_directory)
#         savefig(path_to_fig_directory * fig_file_name , dpi=600)
#     end
#     # show figure
#     if show_fig
#         plt.show()
#     end
# end


function create_confidence_interval_plot(results; nspeeds, wake_model, confidence_interval_width=1.0, upper_10percent=false)
    x_values = unique(results[:,:ndirs])
    aeps = [results[(results[:,:ndirs].==x_value) .& (results[:,:nspeeds].==nspeeds) .& (results[:,:wake_model].==wake_model),"aep_Gaussian"] for x_value in x_values]
    if !upper_10percent
        aeps_means = mean.(aeps)
        aeps_std = std.(aeps)
    else
        aeps_means = calc_mean_upper10percent.(aeps)
        aeps_std = calc_std_upper10percent.(aeps)
    end
    aeps_lower_CI = aeps_means - aeps_std*confidence_interval_width
    aeps_upper_CI = aeps_means + aeps_std*confidence_interval_width
    plt.plot(360 ./ x_values, aeps_means, alpha=0.5, linewidth=1)
    aeps_lower_CI = reshape(aeps_lower_CI, length(aeps_lower_CI))
    aeps_upper_CI = reshape(aeps_upper_CI, length(aeps_upper_CI))
    plt.fill_between(360 ./ x_values, aeps_lower_CI, aeps_upper_CI, alpha=0.1, linewidth=0.0)
end

function create_line_plot(x_values, aeps)
    for i = 1:length(aeps[1])
        plt.plot(x_values, aeps[i], color="C0", alpha=0.5)
    end
end

function calc_mean_upper10percent(x)
    n_x = length(x)
    index_90percentile = Int(round(0.8*n_x))
    sorted_indices = sortperm(x)
    mean_upper10percent_aeps = mean(x[sorted_indices[index_90percentile:end]])
    return mean_upper10percent_aeps
end

function calc_std_upper10percent(x)
    n_x = length(x)
    index_90percentile = Int(round(0.8*n_x))
    sorted_indices = sortperm(x)
    std_upper10percent_aeps = std(x[sorted_indices[index_90percentile:end]])
    return std_upper10percent_aeps
end

function filter_nans(df; equal_number_of_simulations=true)
    hasnan = isnan.(df[:,:aep_Gaussian])
    if equal_number_of_simulations
        bad_simulation_numbers = unique(df[hasnan,:simulation_number])
        good_indices = findall((x) -> !(x in bad_simulation_numbers), df[:,:simulation_number])
        df_filtered = df[good_indices,:]
    else
        df_filtered = df[.!hasnan,:]
    end
    return df_filtered
end

function aep_polar_bar_plot(results; ndirs, nspeeds, wake_model, simulation_number)

    ax = plt.subplot(projection="polar")
    wind_direction = range(0.0, stop=2*pi, length=360) .- pi/360

    if typeof(wake_model)==String
        subset_indices = 
            (results[:,:ndirs].==ndirs) .& 
            (results[:,:nspeeds].==nspeeds) .& 
            (results[:,:wake_model].==wake_model) .& 
            (results[:,:simulation_number].==simulation_number)
        
        ax.bar(wind_direction, aep, bottom=0.0, alpha=0.4)

    elseif (typeof(wake_model)==Vector{String}) && (length(wake_model)==2)
        aep = [Vector{Float64}(),Vector{Float64}()]
        for i = 1:2
            subset_indices[i] = 
                (results[:,:ndirs].==ndirs) .& 
                (results[:,:nspeeds].==nspeeds) .& 
                (results[:,:wake_model].==wake_model[i]) .& 
                (results[:,:simulation_number].==simulation_number)
        end
        
    else
        error("Invalid wake model entry")
    end

    subset_indices_Gaussian = (results[:,:ndirs].==ndirs) .& (results[:,:nspeeds].==nspeeds) .& (results[:,:wake_model].=="Gaussian") .& (results[:,:simulation_number].==simulation_number)
    subset_indices_JensenCosine = (results[:,:ndirs].==ndirs) .& (results[:,:nspeeds].==nspeeds) .& (results[:,:wake_model].=="JensenCosine") .& (results[:,:simulation_number].==simulation_number)
    aep_Gaussian = eval(Meta.parse((results[subset_indices_Gaussian,:aep_Gaussian_binned][1])))
    aep_Gaussian = reshape(sum.(aep_Gaussian), length(aep_Gaussian))
    aep_JensenCosine = eval(Meta.parse((results[subset_indices_JensenCosine,:aep_JensenCosine_binned][1])))
    aep_JensenCosine = reshape(sum.(aep_JensenCosine), length(aep_JensenCosine))
    width = pi/360
    # colors = plt.cm.viridis(aep ./ maximum(maximum.(aep)))

    
    ax.bar(wind_direction, aep_JensenCosine, width=width, bottom=0.0, alpha=0.4)
    ax.bar(wind_direction, aep_Gaussian, width=width, bottom=0.0, alpha=0.4)

    return aep_Gaussian-aep_JensenCosine
end

function extract_binned_aep(results; ndirs, nspeeds, wake_model, simulation_number)
    # get indices of rows of interest
    subset_indices = 
        (results[:,:ndirs].==ndirs) .& 
        (results[:,:nspeeds].==nspeeds) .& 
        (results[:,:wake_model].==wake_model) .& 
        (results[:,:simulation_number].==simulation_number)
    # get the binned AEPs
    aep = eval(Meta.parse((results[subset_indices,:aep_Gaussian_binned][1])))
    return reshape(sum.(aep), length(aep))
end


# set farm directory
farm_directory = "results/final-designs/circle-38turb-7diam"
ndirs_vec = [8,10,11,12,13,15,16,20,24,30,36,45,60,90,180,360]
nspeeds_vec = [1,20]
wake_model_vec = ["Gaussian","JensenCosine"]
simulation_number_vec = 1:100

if !isfile("$farm_directory/combined_results.csv")
    results = create_windrosefidelity_results_dataframe(farm_directory,ndirs_vec,nspeeds_vec,wake_model_vec,simulation_number_vec)
    CSV.write("$farm_directory/combined_results.csv",results)
else
    results=DataFrame(CSV.File("$farm_directory/combined_results.csv"))
end
results = DataFrame(CSV.File("$farm_directory/combined_results_full.csv"))
results = results[results[:,:ndirs].!=360,:]
results_filtered = filter_nans(results,equal_number_of_simulations=false)

# plot_histogram_comparing_wake_models(results_filtered,ndirs=16,nspeeds=20,wake_models=["Gaussian","JensenCosine"])
figure()
create_confidence_interval_plot(results_filtered,nspeeds=1,wake_model="Gaussian")
create_confidence_interval_plot(results_filtered,nspeeds=1,wake_model="JensenCosine")

figure()
diff = aep_polar_bar_plot(results; ndirs=12, nspeeds=1, wake_model="Gaussian", simulation_number=1)
aep_polar_bar_plot(results; ndirs=12, nspeeds=1, wake_model="JensenCosine", simulation_number=1)

