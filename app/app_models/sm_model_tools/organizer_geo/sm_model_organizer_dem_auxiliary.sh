DEM="/home/fabio/Desktop/tmp/regione_marche/model.dem.txt"
OUT="/home/fabio/Desktop/tmp/regione_marche/"

mkdir -p "$OUT"

# If DEM is lon/lat degrees and elevation is meters, use -s 111120
SCALE="-s 111120"

gdaldem aspect    "$DEM" "$OUT/model.aspect.txt"    -of AAIGrid $SCALE -compute_edges
gdaldem slope     "$DEM" "$OUT/model.slope.txt"     -of AAIGrid $SCALE -compute_edges
gdaldem hillshade "$DEM" "$OUT/model.hillshade.txt" -of AAIGrid $SCALE -compute_edges
gdaldem roughness "$DEM" "$OUT/model.roughness.txt" -of AAIGrid -compute_edges
