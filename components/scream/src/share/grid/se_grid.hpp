#ifndef SCREAM_SE_GRID_HPP
#define SCREAM_SE_GRID_HPP

#include "share/grid/abstract_grid.hpp"

#include "ekat/mpi/scream_comm.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/scream_assert.hpp"

namespace scream
{

class SEGrid : public AbstractGrid
{
public:
  using base_type     = AbstractGrid;
  using dofs_map_type = kokkos_types::view<int*[3]>; // elem, igp, jgp

  SEGrid (const std::string& grid_name,
          const GridType type,
          const Comm& comm);

  SEGrid (const std::string& grid_name,
          const GridType type,
          const int num_local_dofs,
          const Comm& comm);

  SEGrid (const dofs_list_type& dofs_gids,
          const std::string& grid_name,
          const GridType type,
          const Comm& comm);

  SEGrid (const dofs_map_type&  dof_to_elgp,
          const dofs_list_type& dofs_gids,
          const std::string& grid_name,
          const GridType type,
          const Comm& comm);

  virtual ~SEGrid () = default;

  // Basic info: name and grid type
  GridType type () const override { return m_type; }
  const std::string& name () const override { return m_grid_name; }

  // Dof counters
  int get_num_local_dofs  () const override { return m_num_local_dofs; }
  int get_num_global_dofs () const override { return m_num_global_dofs; }

  // Dofs info: gids and coords
  dofs_list_type    get_dofs_gids () const override { return m_dofs_gids; }
  dofs_coords_type  get_dofs_coords (const LatLonType) const override ;

  // Method specific to SEGrid
  dofs_map_type get_dofs_map () const { return m_dof_to_elgp; }
protected:

  const std::string   m_grid_name;
  const GridType      m_type;

  int                 m_num_local_dofs;
  int                 m_num_global_dofs;
  dofs_list_type      m_dofs_gids;
  dofs_map_type       m_dof_to_elgp;

  Comm                m_comm;
};

} // namespace scream

#endif // SCREAM_SE_GRID_HPP
