#include "share/grid/se_grid.hpp"

#include "ekat/mpi/scream_comm.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/scream_assert.hpp"

namespace scream
{

SEGrid::
SEGrid (const std::string& grid_name,
        const GridType type,
        const Comm& comm)
 : SEGrid(grid_name,type,0,comm)
{
  // Nothing to do here
}

SEGrid::
SEGrid (const std::string& grid_name,
         const GridType type,
         const int num_local_dofs,
         const Comm& comm)
 : SEGrid(dofs_list_type("",num_local_dofs),grid_name,type,comm)
{
  // Put some invalid number in the list of dofs
  Kokkos::deep_copy(m_dofs_gids,-1);
}

SEGrid::
SEGrid (const dofs_list_type& dofs_gids,
        const std::string& grid_name,
        const GridType type,
        const Comm& comm)
 : SEGrid (dofs_map_type("dof_to_elgp",dofs_gids.extent_int(0)),
           dofs_gids,grid_name,type,comm)
{
  // Put some invalid number in the dof_to_elgp view
  Kokkos::deep_copy(m_dof_to_elgp,-1);
}

SEGrid::
SEGrid (const dofs_map_type&  dof_to_elgp,
        const dofs_list_type& dofs_gids,
        const std::string& grid_name,
        const GridType type,
        const Comm& comm)
 : m_grid_name      (grid_name)
 , m_type           (type)
 , m_num_local_dofs (dofs_gids.extent_int(0))
 , m_dofs_gids      (dofs_gids)
 , m_dof_to_elgp    (dof_to_elgp)
 , m_comm           (comm)
{
  scream_require_msg(type==GridType::SE_CellBased || type==GridType::SE_NodeBased,
                     "Error! Grid type not (yet) supported by SEGrid.\n");
  scream_require_msg(dofs_gids.extent_int(0)==dof_to_elgp.extent_int(0),
                     "Error! Dofs gids and dofs map views have mismatching dimensions.\n");

  // Set the number of global dofs to the sum across all ranks of the number of local dofs
  MPI_Allreduce(&m_num_local_dofs,&m_num_global_dofs,1,MPI_INT,MPI_SUM,m_comm.mpi_comm());
}

SEGrid::dofs_coords_type
SEGrid::get_dofs_coords(const LatLonType lat_lon_type) const {
  
}

} // namespace scream
