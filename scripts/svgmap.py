#! /usr/bin/env python3 -B

import click, json
from math import log, tan, pi, radians

def project(coordinate, center, zoom):
  # https://en.wikipedia.org/wiki/Web_Mercator_projection
  coordinate = (radians(coordinate[0]), radians(coordinate[1]))
  center = (radians(center[0]), radians(center[1]))
  lon, lat = coordinate
  x = (lon + pi) / (2 * pi)
  y = (pi - log(tan(pi/4 +  lat/2))) / (2 * pi)
  clon, clat = center
  cx = (clon + pi) / (2 * pi)
  cy = (pi - log(tan(pi/4 +  clat/2))) / (2 * pi)
  # print('c', x, y)
  x = (x - cx)
  y = (y - cy)
  # print('o', x, y)
  x = x * pow(2.0, zoom)
  y = y * pow(2.0, zoom)
  # print('s', x, y)
  x = (x + 1) * 512
  y = (y + 1) * 512
  return (x, y)

def geojson_to_svg(js, center=[10.465, 44.80], zoom = 15):
  svg = '<?xml version="1.0"?>\n'
  svg += '<svg xmlns="http://www.w3.org/2000/svg" width="1024" height="1024">\n'
  center_min = [+1000, +1000]
  center_max = [-1000, -1000]
  for feature in js['features']:
    geometry = feature['geometry']
    if geometry['type'] == 'Polygon':
      for outline in geometry['coordinates']:
        for coordinate in outline:
          center_min[0] = min(coordinate[0], center_min[0])
          center_min[1] = min(coordinate[1], center_min[1])
          center_max[0] = max(coordinate[0], center_min[0])
          center_max[1] = max(coordinate[1], center_min[1])
    elif geometry['type'] == 'LineString':
      pass
    else:
      pass
  print(center)
  print(center_min, center_max)
  center[0] = (center_max[0] + center_min[0]) / 2
  center[1] = (center_max[1] + center_min[1]) / 2
  print(center)
  for feature in js['features']:
    geometry = feature['geometry']
    if geometry['type'] == 'Polygon':
      points = ''
      for outline in geometry['coordinates']:
        for coordinate in outline:
          x, y = project(coordinate, center, zoom)
          points += f'{x},{y} '
      svg += f'<polygon points="{points}" style="fill:lime;stroke:black;stroke-width:1" />\n'
    elif geometry['type'] == 'LineString':
      points = ''
      for coordinate in geometry['coordinates']:
        x, y = project(coordinate, center, zoom)
        points += f'{x},{y} '
      svg += f'<polyline points="{points}" style="fill:none;stroke:black;stroke-width:2" />\n'
    else:
      print('unsupported ' + geometry['type'])
  svg += '</svg>\n'
  return svg

@click.command()
@click.argument('mapname')
@click.argument('svgname')
def run(mapname='map.json',svgname='map.svg'):
  with open(mapname) as f:
    js = json.load(f)
  if mapname.endswith('.geojson'):
    svg = geojson_to_svg(js)
  else:
    pass
  with open(svgname, 'w') as f:
    f.write(svg)

run()
