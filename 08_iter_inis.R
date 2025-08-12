### generate .ini files needed to run Circuitscape in Julia
### Hester Bronnvik
### 2025-05-28
### hbronnvik@ab.mpg.de

## define the common content
template <- "[Options for advanced mode]
ground_file_is_resistances = True
remove_src_or_gnd = keepall
ground_file = {{param1}}
use_unit_currents = False
source_file = {{param2}}
use_direct_grounds = False

[Mask file]
mask_file = None
use_mask = False

[Calculation options]
low_memory_mode = False
parallelize = False
solver =  cg+amg
print_timings = False
preemptive_memory_release = False
print_rusages = False
max_parallel = 0

[Short circuit regions (aka polygons)]
polygon_file = (Browse for a short-circuit region file)
use_polygons = False

[Options for one-to-all and all-to-one modes]
use_variable_source_strengths = False
variable_source_file = None

[Output options]
set_null_currents_to_nodata = False
set_focal_node_currents_to_zero = False
set_null_voltages_to_nodata = False
compress_grids = False
write_cur_maps = 1
write_volt_maps = False
output_file = Current/{{param3}}
write_cum_cur_map_only = False
log_transform_maps = False
write_max_cur_maps = True

[Version]
version = 4.0.7

[Options for reclassification of habitat data]
reclass_file = (Browse for file with reclassification data)
use_reclass_table = False

[Logging Options]
log_level = INFO
log_file = curconn_out.log
profiler_log_file = None
screenprint_log = False

[Options for pairwise and one-to-all and all-to-one modes]
included_pairs_file = (Browse for a file with pairs to include or exclude)
use_included_pairs = False
point_file = None

[Connection scheme for raster habitat data]
connect_using_avg_resistances = False
connect_four_neighbors_only = False

[Habitat raster or graph]
habitat_map_is_resistances = True
habitat_file = resist/{{param4}}

[Circuitscape mode]
data_type = raster
scenario = advanced"

# the three settings that need to change are the sources/grounds
# ground_file = dem_grounds_east.txt
# source_file = dem_sources_east.txt
# the habitat file
# habitat_file = tile_out/both_prob_rast_dem.tif
# and the output file
# output_file = Current/curconn_demE.out

# Define parameters that change for each file
# grounds change west to east
param1 <- c(rep("dem_grounds_east.txt", 7),
                   rep("dem_grounds.txt", 7))
# sources change west to east
param2 <- c(rep("dem_sources_east.txt", 7),
                   rep("dem_sources.txt", 7))
# each of the resistance layers needs to be used once for each flyway
param4 <- list.files("/home/hbronnvik/Documents/chapter3/Raven/resist", pattern = ".tif")
param4 <- param4[!grepl("aux", param4)]
param4 <- rep(param4, 2)
# and write out the correct file
param3 <- paste0("curconn_", 
                 sub(".tif", "", sub("both_prob_rast_", "", param4)),
                 c(rep("E", 7), rep("W", 7)),
                 ".out")

# Loop to generate the files
for (i in 1:length(param3)) {
  content <- gsub("\\{\\{param1\\}\\}", param1[i], template)
  content <- gsub("\\{\\{param2\\}\\}", param2[i], content)
  content <- gsub("\\{\\{param3\\}\\}", param3[i], content)
  content <- gsub("\\{\\{param4\\}\\}", param4[i], content)
  
  filename <- sprintf("config_%02d.ini", i)
  writeLines(content, paste0("/home/hbronnvik/Documents/chapter3/Raven/iter_inis/", filename))
}
