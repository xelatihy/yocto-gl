# Yocto/Sym

Simple rigib body simulator with collision support for convex and concave
triangle meshes.

The rigib body code performs collision detection and response in gravity.
For collision detection, we use only mesh vertices, so increase object
tesselation to make the simulation more accurate. This allows us to support
convex and concave objects and make the simulation very stable compared to
convex collision detection such as GJK or MPR.

The solver is based on the sequential impulse techniques, more correctly
known as Projected Guass-Sidel. Friction is grossly approximated now,
waiting for a refactoring before getting better.

This library depends in yocto_math.h Optionally depend on yocto_bvh.h/.cpp
for internal acceleration. Disable this by setting YSYM_NO_BVH.


## Usage for Scene Creation

1. create a scene with `make_scene()`
2. add materials with `add_rigid_material()`
3. add shapes with `add_rigid_shape()`
4. add instances with `add_instance()`

## Usage for Simulation

1. either build the point-overlap acceleration structure with
  `init_overlap()` or supply your own with `set_overlap_callbacks()`
2. prepare simuation internal data `init_simulation()`
3. define simulation params with the `simulation_params` structure
4. advance the simiulation with `advance_simulation()`
5. set or get rigid body frames and velocities with `set_XXX()` and
`get_XXX()` methods

## History

- v 0.16: simpler logging
- v 0.15: removal of group overlap
- v 0.14: use yocto_math in the interface and remove inline compilation
- v 0.13: move to add api
- v 1.12: internally use yocto_bvh if desired
- v 0.11: added instances
- v 0.10: switch to .h/.cpp pair
- v 0.9: doxygen comments
- v 0.8: opaque API (allows for changing internals without altering API)
- v 0.7: internally use pointers for performance transparency
- v 0.6: new formulation for moment computation (and bug fixes)
- v 0.5: faster collision detection
- v 0.4: [major API change] move to modern C++ interface
- v 0.3: removal of C interface
- v 0.2: use of STL containers
- v 0.1: C++ implementation
- v 0.0: initial release in C99

## Namespace ysym

Simple rigid body simulator

### Struct scene

~~~ .cpp
struct scene;
~~~

Simulation scene.

### Function make_scene()

~~~ .cpp
scene* make_scene();
~~~

Initialze a scene.

### Function free_scene()

~~~ .cpp
void free_scene(scene*& scn);
~~~

Free a scene.

- Parameters:
    - scn: rigid body scene

### Function add_rigid_shape()

~~~ .cpp
int add_rigid_shape(scene* scn, int ntriangles, const ym::vec3i* triangles,
    int nverts, const ym::vec3f* pos);
~~~

Sets a rigid shape.

- Parameters:
    - scn: scene
    - ntriangles: number of triangles
    - triangles: triangles indices
    - nverts: number of vertices
    - pos: vertex positions

### Function add_rigid_material()

~~~ .cpp
int add_rigid_material(scene* scn, float density);
~~~

Sets a material.

- Parameters:
    - scn: scene
    - density: body density (zero for not-simulated)

### Function add_rigid_body()

~~~ .cpp
int add_rigid_body(scene* scn, const ym::frame3f& frame, int sid, int mid,
    const ym::vec3f& lin_vel =;
~~~

Set a rigid body.

- Parameters:
    - scn: scene
    - frame: body frame
    - sid: shape id
    - mid: material id
    - lin_vel: linear velocity
    - ang_vel: angular velocity

### Function get_rigid_body_frame()

~~~ .cpp
ym::frame3f get_rigid_body_frame(const scene* scn, int bid);
~~~

Get a rigid body frame.

- Parameters:
    - scn: scene
    - bid: body index
- Returns:
    - body frame

### Function get_rigid_body_velocity()

~~~ .cpp
std::pair<ym::vec3f, ym::vec3f> get_rigid_body_velocity(
    const scene* scn, int bid);
~~~

Get a rigid body linear and angular velocity.

- Parameters:
    - scn: scene
    - bid: body index
- Returns:
    - linear and angular velocity

### Function set_rigid_body_frame()

~~~ .cpp
void set_rigid_body_frame(scene* scn, int bid, const ym::frame3f& frame);
~~~

Sets a rigid body frame.

- Parameters:
    - scn: scene
    - bid: body index
    - frame: rigid body frame

### Function set_rigid_body_velocity()

~~~ .cpp
void set_rigid_body_velocity(
    scene* scn, int bid, const ym::vec3f& lin_vel, const ym::vec3f& ang_vel);
~~~

Sets a rigid body linear and angular velocity.

- Parameters:
    - scn: scene
    - bid: body index
    - lin_vel: linear velocity
    - ang_vel: angular velocity

### Struct overlap_point

~~~ .cpp
struct overlap_point {
    float dist = 0;
    int iid = -1;
    int sid = -1;
    int eid = -1;
    ym::vec4f euv = {0, 0, 0, 0};
    operator bool() const; 
}
~~~

Point-scene overlap.

- Members:
    - dist:      overlap distance
    - iid:      instance index
    - sid:      shape index
    - eid:      elements index
    - euv:      element baricentric coordinates
    - operator bool():      check whether it was a hit


### Typedef overlap_shapes_cb

~~~ .cpp
using overlap_shapes_cb = std::function<void(std::vector<ym::vec2i>*)>;
~~~

Shape-shape intersection (conservative)

- Out Parameters:
    - overlaps: overlaps array

### Function Alias function <overlap_point(int sid, const ym::vec3f& pt, float max_dist) \>()

~~~ .cpp
using overlap_shape_cb =
    std::function<overlap_point(int sid, const ym::vec3f& pt, float max_dist)>;
~~~

Closest element intersection callback

- Parameters:
    - sid: shape to check
    - pt: point
    - max_dist: maximum distance
- Return:
    - overlap point

### Typedef overlap_refit_cb

~~~ .cpp
using overlap_refit_cb = std::function<void(const scene* scn, int nshapes)>;
~~~

Refit data structure after transform updates

- Parameters:
    - scn: scene

### Function set_overlap_callbacks()

~~~ .cpp
void set_overlap_callbacks(scene* scn, overlap_shapes_cb overlap_shapes,
    overlap_shape_cb overlap_shape, overlap_refit_cb overlap_refit);
~~~

Set overlap functions

### Function init_overlap()

~~~ .cpp
void init_overlap(scene* scn);
~~~

Initialize overlap functions using internal structures.

### Function compute_moments()

~~~ .cpp
void compute_moments(int ntriangles, const ym::vec3i* triangles, int nverts,
    const ym::vec3f* pos, float* volume, ym::vec3f* center, ym::mat3f* inertia);
~~~

Computes the moments of a shape.

- Parameters:
    - triangles: triangle indices
    - pos: vertex positions
- Output Parameters:
    - volume: volume
    - center: center of mass
    - inertia: inertia tensore (wrt center of mass)

### Function compute_moments()

~~~ .cpp
void compute_moments(int ntetra, const ym::vec4i* tetra, int nverts,
    const ym::vec3f* pos, float* volume, ym::vec3f* center, ym::mat3f* inertia);
~~~

Computes the moments of a shape.

- Parameters:
    - triangles: triangle indices
    - pos: vertex positions
- Output Parameters:
    - volume: volume
    - center: center of mass
    - inertia: inertia tensore (wrt center of mass)

### Function init_simulation()

~~~ .cpp
void init_simulation(scene* scn);
~~~

Initialize the simulation

- Paramaters:
    - scene: rigib body scene

### Struct simulation_params

~~~ .cpp
struct simulation_params {
    float dt = 1 / 60.0f;
    ym::vec3f gravity = {0.f, -9.82f, 0.f};
    int solver_iterations = 20;
    float lin_drag = 0.01;
    float ang_drag = 0.01;
}
~~~

Simulation Parameters.

- Members:
    - dt:      delta time
    - gravity:      gravity
    - solver_iterations:      solver iterations
    - lin_drag:      global linear velocity drag
    - ang_drag:      global angular velocity drag


### Function advance_simulation()

~~~ .cpp
void advance_simulation(scene* scn, const simulation_params& params);
~~~

Advance the simulation one step at a time.

- Paramaters:
    - scene: rigib body scene
    - dt: time step

