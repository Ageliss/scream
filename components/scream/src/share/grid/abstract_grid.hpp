#ifndef SCREAM_ABSTRACT_GRID_HPP
#define SCREAM_ABSTRACT_GRID_HPP

#include "share/grid/grid_utils.hpp"
#include "ekat/scream_types.hpp"

namespace scream
{

class AbstractGrid
{
public:
  using gid_type         = long;          // TODO: template class on gid type
  using device_type      = DefaultDevice; // TODO: template class on device type
  using kokkos_types     = KokkosTypes<device_type>;
  using dofs_list_type   = kokkos_types::view_1d<gid_type>;
  using dofs_coords_type = kokkos_types::view_2d<Real>;

  virtual ~AbstractGrid () = default;

  // Basic info: name and grid type
  virtual GridType type () const = 0;
  virtual const std::string& name () const = 0;

  // Dof counters
  virtual int get_num_local_dofs () const = 0;
  virtual int get_num_global_dofs () const = 0;

  // Dofs info: gids and coords
  virtual dofs_list_type    get_dofs_gids () const = 0;
  virtual dofs_coords_type  get_dofs_coords (const LatLonType) const = 0;
};

} // namespace scream

#endif // SCREAM_ABSTRACT_GRID_HPP
