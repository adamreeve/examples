# read in solution, which may be split into multiple files
@exnodes=<./uniaxial_stress.part*.exnode>;
@exelems=<./uniaxial_stress.part*.exelem>;
foreach $filename (@exnodes) {
    print "Reading $filename\n";
    gfx read node "$filename";
}
foreach $filename (@exelems) {
    print "Reading $filename\n";
    gfx read elem "$filename";
}

# define deformed geometry and pressure
gfx define field "deformed_geom" component Dependent.1 Dependent.2 Dependent.3
gfx define field "hydrostatic_pressure" component Dependent.4

gfx create window 1

gfx create spectrum pressure
gfx modify spectrum pressure clear overwrite_colour
gfx modify spectrum pressure linear reverse range 0.0 1.0 extend_above extend_below rainbow colour_range 0 1

# display deformed geometry
gfx define faces egroup "example_region"
gfx modify g_element "example_region" lines coordinate deformed_geom select_on material default selected_material default_selected
gfx modify g_element "example_region" surfaces coordinate deformed_geom select_on material tissue selected_material default_selected data hydrostatic_pressure spectrum pressure render_shaded
gfx modify g_element "example_region" node_points coordinate deformed_geom glyph sphere General size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected

# display undeformed lines
gfx modify g_element "example_region" lines select_on material green selected_material default_selected

gfx modify g_element "/" point  glyph axes general size "1*1*1" centre 0,0,0 font default select_on material default selected_material default_selected;

gfx modify spectrum pressure autorange

gfx modify window 1 set antialias 2
