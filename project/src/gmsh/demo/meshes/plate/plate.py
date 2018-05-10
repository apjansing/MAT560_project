import pygmsh

geom = pygmsh.opencascade.Geometry(
  characteristic_length_min=0.1,
  characteristic_length_max=0.1,
  )
outerDisk = geom.add_disk([ 0.0,  0.0, 0.0], 0.5)
print(geom.get_code())
points, cells, point_data, cell_data, field_data = pygmsh.generate_mesh(geom, geo_filename = "plate.geo")