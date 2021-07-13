import FLOWFarm; const ff=FLOWFarm
import YAML

function generate_all_initial_farm_design_files()
    for i = 1:100
        generate_initial_farm_design_circle(i,9,5,80.0)
        generate_initial_farm_design_circle(i,9,7,80.0)
        generate_initial_farm_design_circle(i,9,10,80.0)
        generate_initial_farm_design_hornsrev(i)
    end
end


###################################################################################
# SUPPORT FUNCTIONS
###################################################################################

function generate_initial_farm_design_circle(layout_number, nturbines, diameter_spacing, rotor_diameter)

    #### CIRCLE FARM ####

    # file name, path, title, etc.
    farmname = "circle-$(nturbines)turb-$(diameter_spacing)diam"
    filepath = "../inputfiles/farms/random-layouts/$farmname/"
    filename = "initial-design-$(lpad(layout_number,3,"0")).yaml"
    title = "Circular Wind Farm, Randomly Positioned Turbines"

    # generate a random farm layout with uniform turbine type and height
    boundary_radius = sqrt(((sqrt(nturbines)-1)*diameter_spacing*rotor_diameter)^2/pi)
    boundary = ff.CircleBoundary([0.0,0.0], boundary_radius)
    min_spacing = 4.0*rotor_diameter
    param = get_farm_definition_parameters(nturbines, boundary, min_spacing)
    
    # save to yaml file
    description = "$nturbines turbines, $diameter_spacing-diameter spacing, $(round(boundary_radius))-meter boundary radius"
    mkpath(filepath)
    ff.write_farm_definition_YAML(filepath*filename, param..., title=title, description=description)
end

function generate_initial_farm_design_hornsrev(layout_number)

    #### HORNS REV FARM ####

    title = "Horns Rev Wind Farm, Randomly Positioned Turbines"
    description = "80 turbines, parallelogram shaped boundary"
    farmname = "horns-rev"
    horns_rev_site = YAML.load(open("../inputfiles/sites/hornsrev_site_withoutwindresource.yaml"))
    boundary_vertices = horns_rev_site["definitions"]["boundaries"]["items"]["I"]["positions"]
    boundary_vertices = [boundary_vertices[i][j] for i=1:length(boundary_vertices),j=1:2]

    generate_initial_farm_design_polygon(layout_number, 80, boundary_vertices, 80.0; title=title, description=description, farmname=farmname)
end

function generate_initial_farm_design_polygon(layout_number, nturbines, boundary_vertices, rotor_diameter; title="", description="", farmname="polygon")

    #### POLYGON FARM ####

    # file name, path, title, etc.
    farmname = farmname
    filepath = "../inputfiles/farms/random-layouts/$farmname/"
    filename = "initial-design-$(lpad(layout_number,3,"0")).yaml"
    title = title
    description = description

    # generate random farm layout with uniform type and height
    boundary = ff.PolygonBoundary(boundary_vertices)
    min_spacing = 4.0*rotor_diameter
    param = get_farm_definition_parameters(nturbines, boundary, min_spacing)

    # save to yaml file
    mkpath(filepath)
    ff.write_farm_definition_YAML(filepath*filename, param..., title=title, description=description)
end

function get_farm_definition_parameters(nturbines, boundary, min_spacing)
    
    # set and return the fields for FLOWFarm's farm_definition struct
    turbine_x, turbine_y = ff.generate_random_layout(nturbines, boundary, min_spacing)
    turbine_xy = [[turbine_x[i], turbine_y[i]] for i=1:nturbines]
    turbine_z = zeros(nturbines)
    turbine_definition_ids = ones(Int64,nturbines)
    turbine_definition_references = ["VestasV80_2MW.yaml"]

    return turbine_xy, turbine_z, turbine_definition_ids, turbine_definition_references
end
