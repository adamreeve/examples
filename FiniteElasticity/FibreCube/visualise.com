gfx read node FibreCube.part0.exnode;
gfx read elem FibreCube.part0.exelem;

gfx create window 1;

gfx define faces egroup "Region"

# View undeformed geometry
gfx modify g_element "Region" lines coordinate Geometry select_on material green;

# View deformed geometry
gfx modify g_element "Region" lines coordinate DeformedGeometry select_on material white;
gfx modify g_element "Region" node_points coordinate DeformedGeometry glyph sphere General size "0.05*0.05*0.05" centre 0,0,0 font default select_on material white;

# Add axis
gfx modify g_element "/" point  glyph axes_solid general size "0.3*0.3*0.3" centre 0,0,0 font default select_on material default;

# View fibre orientations
gfx modify g_element "Region" element_points coordinate DeformedGeometry discretization "3*3*3" use_elements glyph cylinder_solid size "0.20*0.05*0.05" scale_factors "0*0*0" centre "0.5,0.0,0.0" orientation Fibre material gold
gfx modify g_element "Region" element_points coordinate DeformedGeometry discretization "3*3*3" use_elements glyph sheet size "0.20*0.12*0.12" scale_factors "0*0*0" centre "0.0,0.0,0.0" orientation Fibre material blue

gfx modify window 1 set antialias 4
