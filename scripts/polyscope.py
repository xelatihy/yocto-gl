#! /usr/bin/env python3 -B

import polyscope as ps
import potpourri3d as pp
import click

@click.command()
@click.argument('meshes', nargs=-1)
def run(meshes=[]):
  # Initialize polyscope
  ps.init()

  # load meshes and register them
  for mesh in meshes:
    verts, faces = pp.read_mesh(mesh)
    ps.register_surface_mesh(mesh, verts, faces, smooth_shade=True)

  # View the point cloud and mesh we just registered in the 3D UI
  ps.show()

if __name__ == '__main__':
  run()
