//
// 文件名:     PatchStrategy.C
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:45:28 2014 $
// 描述　:     网格片策略类派生类实现, RRAM 多物理场 (电-热-漂移扩散) 2D 版本.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "BaseElement.h"
#include "CSRMatrixData.h"
#include "CSRMatrixVariable.h"
#include "CellData.h"
#include "CellVariable.h"
#include "EdgeData.h"
#include "EdgeVariable.h"
#include "ElementManager.h"
#include "FaceData.h"
#include "FaceVariable.h"
#include "IntegratorManager.h"
#include "LinearTetrahedron.h"
#include "LinearTriangle.h"
#include "Matrix.h"
#include "NodeData.h"
#include "NodeVariable.h"
#include "Patch.h"
#include "PatchGeometry.h"
#include "PatchStrategy.h"
#include "PatchTopology.h"
#include "ShapeFunctionManager.h"
#include "TetrahedronIntegrator.h"
#include "TetrahedronShapeFunction.h"
#include "TriangleIntegrator.h"
#include "TriangleShapeFunction.h"
#include "Vector.h"
#include "VectorData.h"
#include "VectorVariable.h"
#include "fem/GridInfo.h"
#include "fem/Nedelec2D.h"
#include "fem/TriQuad.h"
#include "math.h"

#include "JAUMIN_Macros.h"

///////////////写入文件//////////////
#include <fstream>
#include <sstream>

// 材料参数查询表: [材料ID][温度系数]
double K_T[5][5] = {{166.246, 5.12e-1, -1.22e-3, 1.1574e-6, -4.1667e-10},
                    {420.33208, -6.809e-2, 0, 0, 0},
                    {0.54335, 1.05e-3, 0, 0, 0},
                    {332.14097, -1.07848e0, 1.58e-3, -1.08505e-6, 2.81425e-10},
                    {20, 0, 0, 0, 0}};
double Sigma_T[5][5] = {{1.31430000e8, -5.59730e5, 1071, -9.582e-1, 3.2504e-4},
                        {2.91115000e8, -1.564300e6, 3700, -3.9347e0, 1.5644e-3},
                        {0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0},
                        {0, 0, 0, 0, 0}};

/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
PatchStrategy::PatchStrategy(const string &object_name, bool is_from_restart,
                             tbox::Pointer<tbox::Database> input_db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!object_name.empty());
#endif

  /// (形函数, 积分器, 单元)管理器
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<ShapeFunctionManager<NDIM> > func_manager =
      ShapeFunctionManager<NDIM>::getManager();
  tbox::Pointer<IntegratorManager<NDIM> > integrator_manager =
      IntegratorManager<NDIM>::getManager();

  /// 创建单元.
#if (NDIM == 2) // 三角形
  tbox::Pointer<LinearTriangle> linear_ele = new LinearTriangle("LinearTriangle");
  tbox::Pointer<TriangleIntegrator> linear_int = new TriangleIntegrator(1, "Triangle");
  tbox::Pointer<TriangleShapeFunction> linear_sha = new TriangleShapeFunction("Triangle");
#else // 四面体
  tbox::Pointer<LinearTetrahedron> linear_ele = new LinearTetrahedron("LinearTetrahedron");
  tbox::Pointer<TetrahedronIntegrator> linear_int = new TetrahedronIntegrator(1, "Tetrahedron");
  tbox::Pointer<TetrahedronShapeFunction> linear_sha = new TetrahedronShapeFunction("Tetrahedron");
#endif
  /// 将单元添加到单元管理器.
  ele_manager->addElement(linear_ele);
  integrator_manager->addIntegrator(linear_int);
  func_manager->addShapeFunction(linear_sha);
  d_object_name = object_name;

  // 读取从输入文件或重启动文件读入数据.
  if (is_from_restart) {
    getFromRestart(input_db);
  } else {
    getFromInput(input_db);
  }

  registerModelVariable();
}

/*************************************************************************
 *
 * 析构函数.
 *
 ************************************************************************/
PatchStrategy::~PatchStrategy() {
  tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
}

/*************************************************************************
 *
 * 注册变量和数据片.
 *
 ************************************************************************/
void PatchStrategy::registerModelVariable() {
  /*************************************************************************
   * 自由度
   ************************************************************************/
  d_dof_info = new solv::DOFInfo<NDIM>(true, false, false, false);
  d_dof_info_DD = new solv::DOFInfo<NDIM>(true, false, false, false);

  /*************************************************************************
   * 用于静态电场求解变量和数据片（电流连续性方程）.
   ************************************************************************/
  DECLARE_MATVEC_VARIABLE(solution, Vector, double, d_dof_info);
  DECLARE_MATVEC_VARIABLE(rhs, Vector, double, d_dof_info);
  DECLARE_MATVEC_VARIABLE(matrix, CSRMatrix, double, d_dof_info);
  DECLARE_VARIABLE(Ex, Cell, double, NDIM, 1);
  DECLARE_VARIABLE(plotdata, Node, double, 1, 1);
  DECLARE_MATVEC_VARIABLE(error, Vector, double, d_dof_info);

  REGISTER_VARIABLE(d_E_solution_id, solution, CURRENT, 1);
  REGISTER_VARIABLE(d_E_rhs_id, rhs, CURRENT, 1);
  REGISTER_VARIABLE(d_E_matrix_id, matrix, CURRENT, 1);
  REGISTER_VARIABLE(d_potential_plot_id, plotdata, CURRENT, 1);
  REGISTER_VARIABLE(d_E_error_id, error, CURRENT, 1);
  REGISTER_VARIABLE(d_Ex_id, Ex, CURRENT, 1);

  /*************************************************************************
   * 矢量有限元信息数据片（2D: 三角单元边元信息）
   ************************************************************************/
  DECLARE_VARIABLE(Cell_jacobian, Cell, double, (NDIM + 1) * (NDIM + 1), 1);
  DECLARE_VARIABLE(Cell_volume, Cell, double, 1, 1);
  DECLARE_VARIABLE(Edge_order, Edge, bool, 1, 1);
  DECLARE_VARIABLE(Edge_flag, Edge, int, 1, 1);

  REGISTER_VARIABLE(d_Edge_order_id, Edge_order, CURRENT, 1);
  REGISTER_VARIABLE(d_Edge_flag_id, Edge_flag, CURRENT, 1);
  REGISTER_VARIABLE(d_Cell_jacobian_id, Cell_jacobian, CURRENT, 1);
  REGISTER_VARIABLE(d_Cell_volume_id, Cell_volume, CURRENT, 1);

  /*************************************************************************
   * 漂移扩散方程变量数据片
   ************************************************************************/
  DECLARE_MATVEC_VARIABLE(DD_solution, Vector, double, d_dof_info_DD);
  DECLARE_MATVEC_VARIABLE(DD_rhs, Vector, double, d_dof_info_DD);
  DECLARE_MATVEC_VARIABLE(DD_matrix, CSRMatrix, double, d_dof_info_DD);
  DECLARE_VARIABLE(DD_nD, Node, double, 1, 1);
  DECLARE_VARIABLE(DD_nD_temp, Node, double, 1, 1);
  DECLARE_MATVEC_VARIABLE(DD_error, Vector, double, d_dof_info_DD);
  DECLARE_VARIABLE(DD_J, Edge, double, 2, 1);

  REGISTER_VARIABLE(d_DD_solution_id, DD_solution, CURRENT, 1);
  REGISTER_VARIABLE(d_DD_rhs_id, DD_rhs, CURRENT, 1);
  REGISTER_VARIABLE(d_DD_matrix_id, DD_matrix, CURRENT, 1);
  REGISTER_VARIABLE(d_nD_plot_id, DD_nD, CURRENT, 1);
  REGISTER_VARIABLE(d_DD_nD_temp_id, DD_nD_temp, CURRENT, 1);
  REGISTER_VARIABLE(d_DD_J_id, DD_J, CURRENT, 1);
  REGISTER_VARIABLE(d_DD_error_id, DD_error, CURRENT, 1);
  d_DD_mobility_id = -1;
  d_DD_diff_id = -1;

  /*************************************************************************
   * 用于热场变量和数据片.
   ************************************************************************/
  DECLARE_MATVEC_VARIABLE(thermal_matrix, CSRMatrix, double, d_dof_info);
  DECLARE_MATVEC_VARIABLE(thermal_rhs, Vector, double, d_dof_info);
  DECLARE_MATVEC_VARIABLE(thermal_solution, Vector, double, d_dof_info);
  DECLARE_VARIABLE(thermal_plotdata, Node, double, 1, 1);
  DECLARE_VARIABLE(thermal_old, Node, double, 1, 1);
  DECLARE_MATVEC_VARIABLE(thermal_error, Vector, double, d_dof_info);

  REGISTER_VARIABLE(d_thermal_solution_id, thermal_solution, CURRENT, 1);
  REGISTER_VARIABLE(d_thermal_rhs_id, thermal_rhs, CURRENT, 1);
  REGISTER_VARIABLE(d_thermal_matrix_id, thermal_matrix, CURRENT, 1);
  REGISTER_VARIABLE(d_temperature_plot_id, thermal_plotdata, CURRENT, 1);
  REGISTER_VARIABLE(d_thermal_old_id, thermal_old, CURRENT, 1);
  REGISTER_VARIABLE(d_thermal_error_id, thermal_error, CURRENT, 1);

  /*************************************************************************
   * 用于材料变量和数据片.
   ************************************************************************/
  DECLARE_VARIABLE(epsilonr, Cell, double, 1, 1);
  DECLARE_VARIABLE(Sigma, Cell, double, 1, 1);
  DECLARE_VARIABLE(Density, Cell, double, 1, 1);
  DECLARE_VARIABLE(Cp, Cell, double, 1, 1);
  DECLARE_VARIABLE(K, Cell, double, 1, 1);

  REGISTER_VARIABLE(d_epsilonr_id, epsilonr, CURRENT, 1);
  REGISTER_VARIABLE(d_Sigma_id, Sigma, CURRENT, 1);
  REGISTER_VARIABLE(d_K_id, K, CURRENT, 1);
  REGISTER_VARIABLE(d_Cp_id, Cp, CURRENT, 1);
  REGISTER_VARIABLE(d_Density_id, Density, CURRENT, 1);

  /*************************************************************************
   * 其他变量和数据片.
   ************************************************************************/
  DECLARE_VARIABLE(Cell_flag, Cell, int, 1, 1);
  REGISTER_VARIABLE(d_Cell_flag_id, Cell_flag, CURRENT, 1);

  /// Yin-Da Wang
  /// SG电流的侧风效应误差估计（2D版本使用棱边变量）
  DECLARE_VARIABLE(Recovery_basis, Edge, double, NDIM, 1);
  REGISTER_VARIABLE(d_Recovery_basis_id, Recovery_basis, CURRENT, 1);
  /// 体插值的电流误差
  DECLARE_VARIABLE(J_error, Cell, double, 1, 1);
  REGISTER_VARIABLE(d_J_error_id, J_error, CURRENT, 1);
}

/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 ************************************************************************/
void PatchStrategy::initializeComponent(algs::IntegratorComponent<NDIM> *component) const {
  const string &component_name = component->getName();

  if (component_name == "INIT") {
    component->registerInitPatchData(d_Edge_order_id);
    component->registerInitPatchData(d_Edge_flag_id);
    component->registerInitPatchData(d_Cell_jacobian_id);
    component->registerInitPatchData(d_Cell_volume_id);
    component->registerInitPatchData(d_nD_plot_id);
    component->registerInitPatchData(d_DD_nD_temp_id);
    component->registerInitPatchData(d_DD_J_id);
    component->registerInitPatchData(d_epsilonr_id);
    component->registerInitPatchData(d_Sigma_id);
    component->registerInitPatchData(d_K_id);
    component->registerInitPatchData(d_Density_id);
    component->registerInitPatchData(d_Cp_id);
    component->registerInitPatchData(d_temperature_plot_id);
    component->registerInitPatchData(d_thermal_old_id);
    component->registerInitPatchData(d_potential_plot_id);
    component->registerInitPatchData(d_Cell_flag_id);
    component->registerInitPatchData(d_Recovery_basis_id);
    component->registerInitPatchData(d_J_error_id);
    d_dof_info->registerToInitComponent(component);
    d_dof_info_DD->registerToInitComponent(component);

  } else if (component_name == "ALLOC") {
    component->registerPatchData(d_E_matrix_id);
    component->registerPatchData(d_E_solution_id);
    component->registerPatchData(d_E_rhs_id);
    component->registerPatchData(d_Ex_id);
    component->registerPatchData(d_DD_matrix_id);
    component->registerPatchData(d_DD_solution_id);
    component->registerPatchData(d_DD_rhs_id);
    component->registerPatchData(d_thermal_matrix_id);
    component->registerPatchData(d_thermal_solution_id);
    component->registerPatchData(d_thermal_rhs_id);
    component->registerPatchData(d_DD_error_id);
    component->registerPatchData(d_E_error_id);
    component->registerPatchData(d_thermal_error_id);

  } else if (component_name == "ALLOC_MULTIPHYSICS") {
    component->registerPatchData(d_E_matrix_id);
    component->registerPatchData(d_E_solution_id);
    component->registerPatchData(d_E_rhs_id);
    component->registerPatchData(d_Ex_id);
    component->registerPatchData(d_DD_matrix_id);
    component->registerPatchData(d_DD_solution_id);
    component->registerPatchData(d_DD_rhs_id);
    component->registerPatchData(d_thermal_matrix_id);
    component->registerPatchData(d_thermal_solution_id);
    component->registerPatchData(d_thermal_rhs_id);
    component->registerPatchData(d_DD_error_id);
    component->registerPatchData(d_E_error_id);
    component->registerPatchData(d_thermal_error_id);

  } else if (component_name == "EDGE_ORDER") {
    component->registerPatchData(d_Edge_order_id);
  } else if (component_name == "COMM_EDGE_FLAG") {
  } else if (component_name == "EDGE_FLAG") {
  } else if (component_name == "INIT_GEOM") {
  } else if (component_name == "RHS") {
  } else if (component_name == "MAT") {
    component->registerCommunicationPatchData(d_DD_nD_temp_id, d_DD_nD_temp_id);
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
  } else if (component_name == "CONS") {
  } else if (component_name == "E_RHS") {
    component->registerCommunicationPatchData(d_Cell_flag_id, d_Cell_flag_id);
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
    component->registerCommunicationPatchData(d_DD_nD_temp_id, d_DD_nD_temp_id);
    component->registerCommunicationPatchData(d_potential_plot_id, d_potential_plot_id);
  } else if (component_name == "E_MAT") {
    component->registerCommunicationPatchData(d_DD_nD_temp_id, d_DD_nD_temp_id);
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
  } else if (component_name == "E_CONS") {
  } else if (component_name == "calculate_Ex") {
    component->registerCommunicationPatchData(d_potential_plot_id, d_potential_plot_id);
  } else if (component_name == "DD_RHS") {
    component->registerCommunicationPatchData(d_Cell_flag_id, d_Cell_flag_id);
  } else if (component_name == "DD_MAT") {
    component->registerCommunicationPatchData(d_Cell_flag_id, d_Cell_flag_id);
    component->registerCommunicationPatchData(d_potential_plot_id, d_potential_plot_id);
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
    component->registerCommunicationPatchData(d_nD_plot_id, d_nD_plot_id);
  } else if (component_name == "DD_CONS") {
  } else if (component_name == "DD_J_SG") {
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
    component->registerCommunicationPatchData(d_potential_plot_id, d_potential_plot_id);
  } else if (component_name == "thermal_RHS") {
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
    component->registerCommunicationPatchData(d_Cell_flag_id, d_Cell_flag_id);
    component->registerCommunicationPatchData(d_DD_nD_temp_id, d_DD_nD_temp_id);
    component->registerCommunicationPatchData(d_potential_plot_id, d_potential_plot_id);
    component->registerCommunicationPatchData(d_Ex_id, d_Ex_id);
  } else if (component_name == "thermal_MAT") {
    component->registerCommunicationPatchData(d_DD_nD_temp_id, d_DD_nD_temp_id);
  } else if (component_name == "thermal_CONS") {
  } else if (component_name == "T_RHS") {
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
    component->registerCommunicationPatchData(d_Cell_flag_id, d_Cell_flag_id);
    component->registerCommunicationPatchData(d_DD_nD_temp_id, d_DD_nD_temp_id);
    component->registerCommunicationPatchData(d_potential_plot_id, d_potential_plot_id);
    component->registerCommunicationPatchData(d_Ex_id, d_Ex_id);
  } else if (component_name == "T_MAT") {
    component->registerCommunicationPatchData(d_DD_nD_temp_id, d_DD_nD_temp_id);
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
  } else if (component_name == "T_CONS") {
  } else if (component_name == "POST") {
  } else if (component_name == "E_POST") {
    component->registerCommunicationPatchData(d_potential_plot_id, d_potential_plot_id);
  } else if (component_name == "DD_POST") {
  } else if (component_name == "ERROR_EST") {
    component->registerCommunicationPatchData(d_Ex_id, d_Ex_id);
    component->registerCommunicationPatchData(d_nD_plot_id, d_nD_plot_id);
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
  } else if (component_name == "DD_ITER_POST") {
  } else if (component_name == "thermal_POST") {
    component->registerCommunicationPatchData(d_thermal_solution_id, d_thermal_solution_id);
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
    component->registerCommunicationPatchData(d_thermal_old_id, d_thermal_old_id);
  } else if (component_name == "RED") {
  } else if (component_name == "calculate_Sigma") {
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
  } else if (component_name == "calculate_K") {
    component->registerCommunicationPatchData(d_temperature_plot_id, d_temperature_plot_id);
  } else if (component_name == "CELL_FLAG") {
  } else if (component_name == "Maxdelta") {
  } else if (component_name == "Dt") {
  } else {
    TBOX_ERROR("\n::initializeComponent() : component " << component_name << " is not matched. "
                                                        << endl);
  }
}

/*************************************************************************
 *  初始化数据片（支持初值构件）.
 *
 *  包含：电、漂移扩散浓度、热、材料属性的完整初始化.
 ************************************************************************/
void PatchStrategy::initializePatchData(hier::Patch<NDIM> &patch, const double time,
                                        const bool initial_time, const string &component_name) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(component_name == "INIT");
#endif
  NULL_USE(time);

  if (initial_time) {
    DECLARE_ADJACENCY(patch, cell, node, Cell, Node);

    GET_PATCH_DATA(patch, plot_data, d_potential_plot_id, Node, double);
    GET_PATCH_DATA(patch, thermal_plot_data, d_temperature_plot_id, Node, double);
    GET_PATCH_DATA(patch, thermal_old_data, d_thermal_old_id, Node, double);

    int num_nodes = patch.getNumberOfNodes(1);
    int local_num_nodes = patch.getNumberOfNodes();

    int *dis_ptr = d_dof_info->getDOFDistribution(patch, hier::EntityUtilities::NODE);
    for (int i = 0; i < num_nodes; ++i) {
      dis_ptr[i] = 1;
    }
    d_dof_info->buildPatchDOFMapping(patch);

    int *dis_ptr_DD = d_dof_info_DD->getDOFDistribution(patch, hier::EntityUtilities::NODE);
    for (int i = 0; i < num_nodes; ++i) {
      dis_ptr_DD[i] = 0;
    }
    for (int imaterial = 0; imaterial < 2; imaterial++) {
      for (int i = 0; i < material_entity[imaterial].size(); i++) {
        if (!HAS_ENTITY_SET(patch, material_entity[imaterial][i], CELL, 1))
          continue;
        DECLARE_ENTITY_SET(patch, cells, material_entity[imaterial][i], CELL, 1);

        for (int jcount = 0; jcount < cells.size(); jcount++) {
          int c = cells[jcount];
          int n_vertex_local = (NDIM == 2) ? 3 : 4;
          for (int inode = 0; inode < n_vertex_local; inode++) {
            int nodeid = cell_node_idx[cell_node_ext[c] + inode];
            dis_ptr_DD[nodeid] = 1;
          }
        }
      }
    }
    d_dof_info_DD->buildPatchDOFMapping(patch);

    for (int i = 0; i < local_num_nodes; ++i) {
      (*plot_data)(0, i) = 0;
      (*thermal_plot_data)(0, i) = BC.T_initial;
      (*thermal_old_data)(0, i) = BC.T_initial;
    }

    /*************************************************************************
     * 材料数据片初始化.
     ************************************************************************/
    GET_PATCH_DATA(patch, epsilonr_data, d_epsilonr_id, Cell, double);
    GET_PATCH_DATA(patch, Sigma_data, d_Sigma_id, Cell, double);
    GET_PATCH_DATA(patch, K_data, d_K_id, Cell, double);
    GET_PATCH_DATA(patch, Density_data, d_Density_id, Cell, double);
    GET_PATCH_DATA(patch, Cp_data, d_Cp_id, Cell, double);

    for (int imaterial = 0; imaterial < NumberOfMaterial; imaterial++) {
      unsigned int id = imaterial + 1;
      for (int i = 0; i < material_entity[imaterial].size(); i++) {
        if (!HAS_ENTITY_SET(patch, material_entity[imaterial][i], CELL, 1))
          continue;
        DECLARE_ENTITY_SET(patch, cells, material_entity[imaterial][i], CELL, 1);

        for (int jcount = 0; jcount < cells.size(); jcount++) {
          int c = cells[jcount];
          (*epsilonr_data)(0, c) = d_region_table[id].d_epsilonr;
          (*Sigma_data)(0, c) = d_region_table[id].d_sigma;
          (*K_data)(0, c) = d_region_table[id].d_K;
          (*Density_data)(0, c) = d_region_table[id].d_density;
          (*Cp_data)(0, c) = d_region_table[id].d_Cp;
        }
      }
    }

    /*************************************************************************
     * 载流子浓度数据片初始化.
     ************************************************************************/
    GET_PATCH_DATA(patch, nD_data, d_nD_plot_id, Node, double);
    GET_PATCH_DATA(patch, nD_temp_data, d_DD_nD_temp_id, Node, double);
    for (int i = 0; i < num_nodes; ++i) {
      (*nD_data)(0, i) = 0;
      (*nD_temp_data)(0, i) = 0;
    }

    for (int imaterial = 0; imaterial < 2; imaterial++) {
      unsigned int id = imaterial;
      for (int i = 0; i < material_entity[imaterial].size(); i++) {
        if (!HAS_ENTITY_SET(patch, material_entity[imaterial][i], CELL, 1))
          continue;
        DECLARE_ENTITY_SET(patch, cells, material_entity[imaterial][i], CELL, 1);

        for (int jcount = 0; jcount < cells.size(); jcount++) {
          int c = cells[jcount];
          int n_vertex_local = (NDIM == 2) ? 3 : 4;
          for (int inode = 0; inode < n_vertex_local; inode++) {
            int nodeid = cell_node_idx[cell_node_ext[c] + inode];
            (*nD_data)(0, nodeid) = nd_int[id];
            (*nD_temp_data)(0, nodeid) = nd_int[id];
          }
        }
      }
    }

    for (int con_id = 0; con_id < BC.num_concentration; con_id++) {
      for (int con_face = 0; con_face < BC.concentration_bc[con_id].con_face.size(); con_face++) {
        if (HAS_ENTITY_SET(patch, BC.concentration_bc[con_id].con_face[con_face], NODE, 1)) {
          DECLARE_ENTITY_SET(patch, entity_idx, BC.concentration_bc[con_id].con_face[con_face],
                             NODE, 1);
          for (int i = 0; i < entity_idx.size(); ++i) {
            int index = entity_idx[i];
            (*nD_data)(0, index) = BC.concentration_bc[con_id].concentration;
          }
        }
      }
    }
  }
}

/*************************************************************************
 *  注册可视化数据.
 ************************************************************************/
void PatchStrategy::registerPlotData(
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer) {
  javis_data_writer->registerPlotQuantity("plot", "SCALAR", d_potential_plot_id);
  javis_data_writer->registerPlotQuantity("thermal_plot", "SCALAR", d_temperature_plot_id);
  javis_data_writer->registerPlotQuantity("nD", "SCALAR", d_nD_plot_id);
  javis_data_writer->registerPlotQuantity("DD_error", "SCALAR", d_J_error_id);
  javis_data_writer->registerPlotQuantity("cell_identity", "SCALAR", d_Cell_flag_id);
}

/*************************************************************************
 *  输出数据成员到重启动数据库.
 ************************************************************************/
void PatchStrategy::putToDatabase(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif
}

/*************************************************************************
 *  取得时间步长.
 ************************************************************************/
double PatchStrategy::getPatchDt(hier::Patch<NDIM> &patch, const double time,
                                 const bool initial_time, const int flag_last_dt,
                                 const double last_dt, const string &component_name) {
  return Dt;
}

/*************************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 *
 * 调度所有物理模块的矩阵组装、右端项、约束加载、后处理.
 ************************************************************************/
void PatchStrategy::computeOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt,
                                   const bool initial_time, const string &component_name) {
  if (component_name == "INIT_GEOM") {
    initGeometryOnPatch(patch);
  } else if (component_name == "EDGE_FLAG") {
  } else if (component_name == "COMM_EDGE_FLAG") {
  } else if (component_name == "MAT") {
    buildMatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "RHS") {
    buildRHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "CONS") {
    applyConstraint(patch, time, dt, component_name);
  } else if (component_name == "E_MAT") {
    buildEMatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "E_RHS") {
    buildERHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "E_CONS") {
    applyEConstraint(patch, time, dt, component_name);
  } else if (component_name == "calculate_Ex") {
    calculateEx(patch, time, dt, component_name);
  } else if (component_name == "DD_MAT") {
    buildDDMatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "DD_RHS") {
    buildDDRHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "DD_CONS") {
    applyDDConstraint(patch, time, dt, component_name);
  } else if (component_name == "DD_J_SG") {
    calculateDD_J(patch, time, dt, component_name);
  } else if (component_name == "thermal_MAT") {
    buildthermalMatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "thermal_RHS") {
    buildthermalRHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "thermal_CONS") {
    applythermalConstraint(patch, time, dt, component_name);
  } else if (component_name == "T_MAT") {
    buildTMatrixOnPatch(patch, time, dt, component_name);
  } else if (component_name == "T_RHS") {
    buildTRHSOnPatch(patch, time, dt, component_name);
  } else if (component_name == "T_CONS") {
    applyTConstraint(patch, time, dt, component_name);
  } else if (component_name == "calculate_Sigma") {
    calculateSigma(patch, time, dt, component_name);
  } else if (component_name == "calculate_K") {
    calculateK(patch, time, dt, component_name);
  } else if (component_name == "CELL_FLAG") {
    computeCellFlagOnPatch(patch);
  } else if (component_name == "POST") {
    postProcess(patch, time, dt, component_name);
  } else if (component_name == "E_POST") {
    postEProcess(patch, time, dt, component_name);
  } else if (component_name == "DD_POST") {
    postDDProcess(patch, time, dt, component_name);
  } else if (component_name == "ERROR_EST") {
    calculateErrorOnPatch(patch, time, dt, component_name);
  } else if (component_name == "DD_ITER_POST") {
    postDDIterError(patch, time, dt, component_name);
  } else if (component_name == "thermal_POST") {
    postthermalProcess(patch, time, dt, component_name);
  } else {
    TBOX_ERROR(" PatchStrategy :: component name \"" << component_name << "\" is not matched! ");
  }
}

/*************************************************************************
 *  标记单元（材料区域标识）.
 ************************************************************************/
void PatchStrategy::computeCellFlagOnPatch(hier::Patch<NDIM> &patch) {
  GET_PATCH_DATA(patch, Cell_flag, d_Cell_flag_id, Cell, int);
  Cell_flag->fillAll(0);

  for (int imaterial = 0; imaterial < NumberOfMaterial; imaterial++) {
    unsigned int id = imaterial + 1;
    for (int i = 0; i < material_entity[imaterial].size(); i++) {
      if (!HAS_ENTITY_SET(patch, material_entity[imaterial][i], CELL, 1))
        continue;
      DECLARE_ENTITY_SET(patch, cells, material_entity[imaterial][i], CELL, 1);

      for (int jcount = 0; jcount < cells.size(); jcount++) {
        int c = cells[jcount];
        (*Cell_flag)(0, c) = id;
      }
    }
  }
}

/*************************************************************************
 *  计算单元体积、Jacobian矩阵（2D: 使用 TriGeom 替代 TetGeom）.
 ************************************************************************/
void PatchStrategy::initGeometryOnPatch(hier::Patch<NDIM> &patch) {
  int num_cell_ghost = patch.getNumberOfCells(1);
  GET_PATCH_DATA(patch, jacobian, d_Cell_jacobian_id, Cell, double);
  GET_PATCH_DATA(patch, volume, d_Cell_volume_id, Cell, double);

#if (NDIM == 2)
  GET_COORD_DATA(patch, node_coord, Node);
  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);

  for (int c = 0; c < num_cell_ghost; c++) {
    int n0 = cell_node_idx[cell_node_ext[c] + 0];
    int n1 = cell_node_idx[cell_node_ext[c] + 1];
    int n2 = cell_node_idx[cell_node_ext[c] + 2];
    double x0 = (*node_coord)(0, n0), y0 = (*node_coord)(1, n0);
    double x1 = (*node_coord)(0, n1), y1 = (*node_coord)(1, n1);
    double x2 = (*node_coord)(0, n2), y2 = (*node_coord)(1, n2);

    double detJ = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    double area = 0.5 * fabs(detJ);
    (*volume)(0, c) = area;

    double inv_detJ = 1.0 / detJ;
    double *jac = &((*jacobian)(0, c));
    // nabla[0] = (∇N_0.x, ∇N_0.y, 0)
    jac[0] = (y1 - y2) * inv_detJ;
    jac[1] = (x2 - x1) * inv_detJ;
    jac[2] = 0.0;
    // nabla[1] = (∇N_1.x, ∇N_1.y, 0)
    jac[3] = (y2 - y0) * inv_detJ;
    jac[4] = (x0 - x2) * inv_detJ;
    jac[5] = 0.0;
    // nabla[2] = (∇N_2.x, ∇N_2.y, 0)
    jac[6] = (y0 - y1) * inv_detJ;
    jac[7] = (x1 - x0) * inv_detJ;
    jac[8] = 0.0;
  }
#else
  JAUMIN::appu::TetGeom tetrahedron(patch);
  for (int c = 0; c < num_cell_ghost; c++) {
    (*volume)(0, c) = tetrahedron.volume(c);
    tetrahedron.jacobian(c, &((*jacobian)(0, c)));
  }
#endif
}

// ===== 电学模块 =====

void PatchStrategy::applyConstraint(hier::Patch<NDIM> &patch, const double time, const double dt,
                                    const string &component_name) {
  GET_PATCH_DATA(patch, mat_data, d_E_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, vec_data, d_E_rhs_id, Vector, double);

  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  int *row_start = mat_data->getRowStartPointer();
  int *col_idx = mat_data->getColumnIndicesPointer();
  double *mat_val = mat_data->getValuePointer();
  double *vec_val = vec_data->getPointer();

  // 固定电压边界
  for (int vol_id = 0; vol_id < BC.num_voltage; vol_id++) {
    for (int vol_face = 0; vol_face < BC.voltage_bc[vol_id].vol_face.size(); vol_face++) {
      if (HAS_ENTITY_SET(patch, BC.voltage_bc[vol_id].vol_face[vol_face], NODE, 1)) {
        DECLARE_ENTITY_SET(patch, entity_idx, BC.voltage_bc[vol_id].vol_face[vol_face], NODE, 1);
        for (int i = 0; i < entity_idx.size(); ++i) {
          int index = dof_map[entity_idx[i]];
          vec_val[index] = BC.voltage_bc[vol_id].voltage;
          for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
            if (col_idx[j] == index) {
              mat_val[j] = 1.0;
            } else {
              mat_val[j] = 0.0;
              vec_val[col_idx[j]] -= (*mat_data)(col_idx[j], index) * vec_val[index];
              (*mat_data)(col_idx[j], index) = 0.0;
            }
          }
        }
      }
    }
  }

  // 时变电压边界
  for (int vol_id = 0; vol_id < BC.num_timevoltage; vol_id++) {
    for (int vol_face = 0; vol_face < BC.timevoltage_bc[vol_id].vol_face.size(); vol_face++) {
      if (HAS_ENTITY_SET(patch, BC.timevoltage_bc[vol_id].vol_face[vol_face], NODE, 1)) {
        DECLARE_ENTITY_SET(patch, entity_idx, BC.timevoltage_bc[vol_id].vol_face[vol_face], NODE,
                           1);
        for (int i = 0; i < entity_idx.size(); ++i) {
          int index = dof_map[entity_idx[i]];
          vec_val[index] = voltagePulse(time);
          for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
            if (col_idx[j] == index) {
              mat_val[j] = 1.0;
            } else {
              mat_val[j] = 0.0;
              vec_val[col_idx[j]] -= (*mat_data)(col_idx[j], index) * vec_val[index];
              (*mat_data)(col_idx[j], index) = 0.0;
            }
          }
        }
      }
    }
  }
}

void PatchStrategy::applyEConstraint(hier::Patch<NDIM> &patch, const double time, const double dt,
                                     const string &component_name) {
  /// 电学约束加载 (E_CONS 构件).
  /// 与 applyConstraint 同构: 固定电压 + 时变电压 Dirichlet 处理.
  GET_PATCH_DATA(patch, mat_data, d_E_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, vec_data, d_E_rhs_id, Vector, double);

  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  int *row_start = mat_data->getRowStartPointer();
  int *col_idx = mat_data->getColumnIndicesPointer();
  double *mat_val = mat_data->getValuePointer();
  double *vec_val = vec_data->getPointer();

  // 固定电压边界
  for (int vol_id = 0; vol_id < BC.num_voltage; vol_id++) {
    for (int vol_face = 0; vol_face < BC.voltage_bc[vol_id].vol_face.size(); vol_face++) {
      if (HAS_ENTITY_SET(patch, BC.voltage_bc[vol_id].vol_face[vol_face], NODE, 1)) {
        DECLARE_ENTITY_SET(patch, entity_idx, BC.voltage_bc[vol_id].vol_face[vol_face], NODE, 1);
        for (int i = 0; i < entity_idx.size(); ++i) {
          int index = dof_map[entity_idx[i]];
          vec_val[index] = BC.voltage_bc[vol_id].voltage;
          for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
            if (col_idx[j] == index) {
              mat_val[j] = 1.0;
            } else {
              mat_val[j] = 0.0;
              vec_val[col_idx[j]] -= (*mat_data)(col_idx[j], index) * vec_val[index];
              (*mat_data)(col_idx[j], index) = 0.0;
            }
          }
        }
      }
    }
  }

  // 时变电压边界
  for (int vol_id = 0; vol_id < BC.num_timevoltage; vol_id++) {
    for (int vol_face = 0; vol_face < BC.timevoltage_bc[vol_id].vol_face.size(); vol_face++) {
      if (HAS_ENTITY_SET(patch, BC.timevoltage_bc[vol_id].vol_face[vol_face], NODE, 1)) {
        DECLARE_ENTITY_SET(patch, entity_idx, BC.timevoltage_bc[vol_id].vol_face[vol_face], NODE,
                           1);
        for (int i = 0; i < entity_idx.size(); ++i) {
          int index = dof_map[entity_idx[i]];
          vec_val[index] = voltagePulse(time);
          for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
            if (col_idx[j] == index) {
              mat_val[j] = 1.0;
            } else {
              mat_val[j] = 0.0;
              vec_val[col_idx[j]] -= (*mat_data)(col_idx[j], index) * vec_val[index];
              (*mat_data)(col_idx[j], index) = 0.0;
            }
          }
        }
      }
    }
  }
}

void PatchStrategy::buildRHSOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt,
                                    const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  GET_PATCH_DATA(patch, vec_data, d_E_rhs_id, Vector, double);
  GET_PATCH_DATA(patch, cell_flag, d_Cell_flag_id, Cell, int);

  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
      }
    }
    tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
    ele_vec->resize(n_vertex);
    for (int ii = 0; ii < n_vertex; ++ii)
      (*ele_vec)[ii] = 0.0;
    ele->buildEElementRHS(vertex, dt, time, (*cell_flag)(0, i), ele_vec);
    for (int i2 = 0; i2 < n_vertex; ++i2)
      vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
  }
}

void PatchStrategy::buildERHSOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt,
                                     const string &component_name) {
  /// 电学右端项组装 (E_RHS 构件).
  /// 与 3D OneOrderPatchStrategy::buildRHSOnPatch 同构:
  /// 遍历单元, 由 buildEElementRHS 计算单元右端项并组装到全局向量.
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  GET_PATCH_DATA(patch, vec_data, d_E_rhs_id, Vector, double);
  GET_PATCH_DATA(patch, cell_flag, d_Cell_flag_id, Cell, int);

  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
      }
    }
    tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
    ele_vec->resize(n_vertex);
    for (int ii = 0; ii < n_vertex; ++ii)
      (*ele_vec)[ii] = 0.0;
    ele->buildEElementRHS(vertex, dt, time, (*cell_flag)(0, i), ele_vec);
    for (int i2 = 0; i2 < n_vertex; ++i2)
      vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
  }
}

void PatchStrategy::buildMatrixOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt,
                                       const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, mat_data, d_E_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, Sigma_data, d_Sigma_id, Cell, double);
  GET_PATCH_DATA(patch, nD_data, d_DD_nD_temp_id, Node, double);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, cell_flag, d_Cell_flag_id, Cell, int);

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<double> T(n_vertex);
    tbox::Array<double> nD(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
      }
      T[i1] = (*T_data)(0, cell_node_idx[j]);
      nD[i1] = nD_data->getPointer()[cell_node_idx[j]];
    }

    tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
    ele_mat->resize(n_vertex, n_vertex);
    for (int ii = 0; ii < n_vertex; ++ii)
      for (int jj = 0; jj < n_vertex; ++jj)
        (*ele_mat)(ii, jj) = 0.0;

    ele->buildEStiffElementMatrix(vertex, dt, time, T, nD, (*cell_flag)(0, i), ele_mat);

    for (int ii = 0; ii < n_vertex; ++ii) {
      int row = mapping[ii];
      for (int jj = 0; jj < n_vertex; ++jj)
        mat_data->addMatrixValue(row, mapping[jj], (*ele_mat)(ii, jj));
    }
  }
  mat_data->assemble();
}

void PatchStrategy::calculateEx(hier::Patch<NDIM> &patch, const double time, const double dt,
                                const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, Ex_data, d_Ex_id, Cell, double);
  GET_PATCH_DATA(patch, Phi_data, d_potential_plot_id, Node, double);

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<double> phi(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k)
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
      phi[i1] = (*Phi_data)(0, cell_node_idx[j]);
    }
    tbox::Pointer<tbox::Vector<double> > ele_Ex = new tbox::Vector<double>();
    ele_Ex->resize(NDIM);
    for (int k = 0; k < NDIM; ++k)
      (*ele_Ex)[k] = 0.0;
    ele->calculateElementEx(vertex, dt, time, phi, ele_Ex);
    for (int k = 0; k < NDIM; ++k)
      (*Ex_data)(k, i) = (*ele_Ex)[k];
  }
}

// ===== 漂移扩散模块 =====

void PatchStrategy::applyDDConstraint(hier::Patch<NDIM> &patch, const double time, const double dt,
                                      const string &component_name) {
  GET_PATCH_DATA(patch, DDmat_data, d_DD_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, DDvec_data, d_DD_rhs_id, Vector, double);

  // 注意: DD 矩阵/向量按 d_dof_info_DD 编号(只对半导体材料分配 DOF),
  // 必须用 d_dof_info_DD 的映射, 用 d_dof_info(E 的全节点恒等映射)会越界写.
  // 与 3D OneOrderPatchStrategy::applyDDConstraint 一致.
  int *dof_map = d_dof_info_DD->getDOFMapping(patch, hier::EntityUtilities::NODE);
  int *row_start = DDmat_data->getRowStartPointer();
  int *col_idx = DDmat_data->getColumnIndicesPointer();
  double *mat_val = DDmat_data->getValuePointer();
  double *vec_val = DDvec_data->getPointer();

  for (int con_id = 0; con_id < BC.num_concentration; con_id++) {
    for (int con_face = 0; con_face < BC.concentration_bc[con_id].con_face.size(); con_face++) {
      if (HAS_ENTITY_SET(patch, BC.concentration_bc[con_id].con_face[con_face], NODE, 1)) {
        DECLARE_ENTITY_SET(patch, entity_idx, BC.concentration_bc[con_id].con_face[con_face], NODE,
                           1);
        for (int i = 0; i < entity_idx.size(); ++i) {
          int index = dof_map[entity_idx[i]];
          vec_val[index] = BC.concentration_bc[con_id].concentration;
          for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
            if (col_idx[j] == index) {
              mat_val[j] = 1.0;
            } else {
              mat_val[j] = 0.0;
              vec_val[col_idx[j]] -= (*DDmat_data)(col_idx[j], index) * vec_val[index];
              (*DDmat_data)(col_idx[j], index) = 0.0;
            }
          }
        }
      }
    }
  }
}

void PatchStrategy::buildDDRHSOnPatch(hier::Patch<NDIM> &patch, const double time, const double dt,
                                      const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  DECLARE_ADJACENCY(patch, cell, edge, Cell, Edge);
  GET_PATCH_DATA(patch, DDvec_data, d_DD_rhs_id, Vector, double);
  GET_PATCH_DATA(patch, nD_data, d_nD_plot_id, Node, double);
  GET_PATCH_DATA(patch, Cell_flag, d_Cell_flag_id, Cell, int);

  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info_DD->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    bool is_DD_solve = (*Cell_flag)(0, i) == 1 || (*Cell_flag)(0, i) == 2;
    if (is_DD_solve) {
      int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
      tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
      tbox::Array<int> mapping(n_vertex);
      tbox::Array<double> nD(n_vertex);
      for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
        mapping[i1] = dof_map[cell_node_idx[j]];
        for (int k = 0; k < NDIM; ++k)
          vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
        nD[i1] = (*nD_data)(0, cell_node_idx[j]);
      }
      tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
      ele_vec->resize(n_vertex);
      for (int ii = 0; ii < n_vertex; ++ii)
        (*ele_vec)[ii] = 0.0;
      ele->buildDDElementRHS(vertex, dt, time, nD, ele_vec);
      for (int i2 = 0; i2 < n_vertex; ++i2)
        DDvec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
    }
  }
}

void PatchStrategy::buildDDMatrixOnPatch(hier::Patch<NDIM> &patch, const double time,
                                         const double dt, const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, DDmat_data, d_DD_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, nD_data, d_nD_plot_id, Node, double);
  GET_PATCH_DATA(patch, Phi_data, d_potential_plot_id, Node, double);
  GET_PATCH_DATA(patch, J_data, d_DD_J_id, Edge, double);
  GET_PATCH_DATA(patch, E_data, d_Ex_id, Cell, double);
  GET_PATCH_DATA(patch, edgeorder_data, d_Edge_order_id, Edge, bool);
  GET_PATCH_DATA(patch, Cell_flag, d_Cell_flag_id, Cell, int);

#if (NDIM == 2)
  appu::Nedelec2D shapefunc(patch, patch.getPatchData(d_Edge_order_id),
                            patch.getPatchData(d_Cell_jacobian_id));
  appu::TriQuad quad(patch, patch.getPatchData(d_Cell_volume_id),
                     patch.getPatchData(d_Cell_jacobian_id));
#else
  appu::Nedelec shapefunc(patch, patch.getPatchData(d_Edge_order_id),
                          patch.getPatchData(d_Cell_jacobian_id));
  appu::TetQuad quad(patch, patch.getPatchData(d_Cell_volume_id),
                     patch.getPatchData(d_Cell_jacobian_id));
#endif

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  DECLARE_ADJACENCY(patch, cell, edge, Cell, Edge);
  DECLARE_ADJACENCY(patch, edge, node, Edge, Node);

  int num_cells = patch.getNumberOfCells(1);
  int *dof_map_DD = d_dof_info_DD->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    bool is_DD_solve = (*Cell_flag)(0, i) == 1 || (*Cell_flag)(0, i) == 2;
    // 仅对半导体区域进行漂移扩散求解
    if (is_DD_solve) {
      // JAUMIN cell-node 邻接表对三角形也存 4 个槽(第 4 个是填充), 按 3 截断
      int n_vertex_local = (NDIM == 2) ? 3
                         : cell_node_ext[i + 1] - cell_node_ext[i];
      int n_edge_local   = (NDIM == 2) ? 3
                         : cell_edge_ext[i + 1] - cell_edge_ext[i];
      /// 单元矩阵
      tbox::Pointer<tbox::Matrix<double> > ele_K = new tbox::Matrix<double>();
      ele_K->resize(n_vertex_local, n_vertex_local);
      for (int ii = 0; ii < n_vertex_local; ++ii)
        for (int jj = 0; jj < n_vertex_local; ++jj)
          (*ele_K)(ii, jj) = 0.0;

      for (int ii = 0; ii < n_vertex_local; ++ii) {
        for (int jj = 0; jj < n_vertex_local; ++jj) {
          int nodejj = cell_node_idx[cell_node_ext[i] + jj];
          for (int le = 0; le < n_edge_local; le++) {
            int edge_id = cell_edge_idx[cell_edge_ext[i] + le];
            int node_id1 = edge_node_idx[edge_node_ext[edge_id] + 0];
            int node_id2 = edge_node_idx[edge_node_ext[edge_id] + 1];
            if (nodejj == node_id1 || nodejj == node_id2) {
              /// 结点上的电流密度值
              double edgeJ = (nodejj == node_id1) ? (*J_data)(0, edge_id) : (*J_data)(1, edge_id);
              double Kij[1] = {0};
              quad.quadcalculateKij(i, &shapefunc, 0, le, ii, edgeJ, &Kij[0]);
              (*ele_K)(ii, jj) += Kij[0];
            }
          }
        }
      }

      int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
      tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
      tbox::Array<double> Tn(n_vertex);
      tbox::Array<double> Phi(n_vertex);
      tbox::Array<int> mapping(n_vertex);
      for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
        mapping[i1] = dof_map_DD[cell_node_idx[j]];
        for (int k = 0; k < NDIM; ++k)
          vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
        Tn[i1] = (*T_data)(0, cell_node_idx[j]);
        Phi[i1] = (*Phi_data)(0, cell_node_idx[j]);
      }

      tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
      ele_mat->resize(n_vertex, n_vertex);
      for (int ii = 0; ii < n_vertex; ++ii)
        for (int jj = 0; jj < n_vertex; ++jj)
          (*ele_mat)(ii, jj) = 0.0;

      ele->buildDDStiffElementMatrix(vertex, dt, time, ele_K, Phi, Tn, ele_mat);

      for (int ii = 0; ii < n_vertex; ++ii) {
        int row = mapping[ii];
        for (int jj = 0; jj < n_vertex; ++jj)
          DDmat_data->addMatrixValue(row, mapping[jj], (*ele_mat)(ii, jj));
      }
    }
  }
  DDmat_data->assemble();
}

void PatchStrategy::calculateDD_J(hier::Patch<NDIM> &patch, const double time, const double dt,
                                  const string &component_name) {
  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, Phi_data, d_potential_plot_id, Node, double);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, nD_data, d_nD_plot_id, Node, double);
  GET_PATCH_DATA(patch, J_data, d_DD_J_id, Edge, double);
  GET_PATCH_DATA(patch, edgeorder_data, d_Edge_order_id, Edge, bool);

  DECLARE_ADJACENCY(patch, edge, node, Edge, Node);
  int num_edges = patch.getNumberOfEdges(1);

  for (int edge = 0; edge < num_edges; edge++) {
    int node_id1 = edge_node_idx[edge_node_ext[edge]];
    int node_id2 = edge_node_idx[edge_node_ext[edge] + 1];
    double x1 = (*node_coord)(0, node_id1);
    double x2 = (*node_coord)(0, node_id2);
    double y1 = (*node_coord)(1, node_id1);
    double y2 = (*node_coord)(1, node_id2);
    /// 二维情况下不存在z轴
    double z1 = (NDIM > 2) ? (*node_coord)(2, node_id1) : 0.0;
    double z2 = (NDIM > 2) ? (*node_coord)(2, node_id2) : 0.0;
    /// 全局边方向
    int deltaij = -1;
    bool in_order = (*edgeorder_data)(0, edge);
    if (!in_order)
      deltaij = 1;
    /// 边长度
    double l_edge = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
    /// 取中点温度
    double T = 0.5 * ((*T_data)(0, node_id2) + (*T_data)(0, node_id1));
    /// 扩散系数
    double Diff_coef = Diff0 * exp(-EA / (K_b * T));
    /// 氧空穴迁移率
    double Mobility_coef = Echarge * Diff_coef / (K_b * T);
    /// 电势
    double phi_i = (*Phi_data)(0, node_id1);
    double phi_j = (*Phi_data)(0, node_id2);
    /// 电场强度
    double E_ij = -(phi_i - phi_j) / l_edge;
    double Velocity_ij = Mobility_coef * E_ij * deltaij;
    /// Peclet数
    double beta_ij = Velocity_ij * l_edge / (2.0 * Diff_coef);

    if (abs(beta_ij) < 1e-7) {
      (*J_data)(0, edge) = deltaij * Diff_coef;
      (*J_data)(1, edge) = -deltaij * Diff_coef;
    } else {
      (*J_data)(0, edge) = deltaij * (l_edge * Velocity_ij * 0.5) * (1.0 / tanh(beta_ij) - deltaij);
      (*J_data)(1, edge) =
          -deltaij * (l_edge * Velocity_ij * 0.5) * (1.0 / tanh(beta_ij) + deltaij);
    }
  }
}

// ===== 热传导模块 =====

void PatchStrategy::applythermalConstraint(hier::Patch<NDIM> &patch, const double time,
                                           const double dt, const string &component_name) {
  GET_PATCH_DATA(patch, mat_data, d_thermal_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, vec_data, d_thermal_rhs_id, Vector, double);

  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  int *row_start = mat_data->getRowStartPointer();
  int *col_idx = mat_data->getColumnIndicesPointer();
  double *mat_val = mat_data->getValuePointer();
  double *vec_val = vec_data->getPointer();

  for (int temp_id = 0; temp_id < BC.num_temperature; temp_id++) {
    for (int temp_face = 0; temp_face < BC.FixTemperature_bc[temp_id].temp_face.size();
         temp_face++) {
      if (HAS_ENTITY_SET(patch, BC.FixTemperature_bc[temp_id].temp_face[temp_face], NODE, 1)) {
        DECLARE_ENTITY_SET(patch, entity_idx, BC.FixTemperature_bc[temp_id].temp_face[temp_face],
                           NODE, 1);
        for (int i = 0; i < entity_idx.size(); ++i) {
          int index = dof_map[entity_idx[i]];
          vec_val[index] = BC.FixTemperature_bc[temp_id].temperature;
          for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
            if (col_idx[j] == index) {
              mat_val[j] = 1.0;
            } else {
              mat_val[j] = 0.0;
              vec_val[col_idx[j]] -= (*mat_data)(col_idx[j], index) * vec_val[index];
              (*mat_data)(col_idx[j], index) = 0.0;
            }
          }
        }
      }
    }
  }
}

void PatchStrategy::buildthermalRHSOnPatch(hier::Patch<NDIM> &patch, const double time,
                                           const double dt, const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  GET_PATCH_DATA(patch, vec_data, d_thermal_rhs_id, Vector, double);
  GET_PATCH_DATA(patch, Ex_data, d_Ex_id, Cell, double);
  GET_PATCH_DATA(patch, Sigma_data, d_Sigma_id, Cell, double);
  GET_PATCH_DATA(patch, nD_data, d_DD_nD_temp_id, Node, double);
  GET_PATCH_DATA(patch, phi_data, d_potential_plot_id, Node, double);
  GET_PATCH_DATA(patch, cell_flag, d_Cell_flag_id, Cell, int);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);

  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    tbox::Array<double> T_val(n_vertex);
    tbox::Array<double> nD(n_vertex);
    tbox::Array<double> phi(n_vertex);
    tbox::Array<double> E(NDIM);
    for (int k = 0; k < NDIM; k++)
      E[k] = (*Ex_data)(k, i);

    double u_val = 0;
    for (int k = 0; k < NDIM; k++)
      u_val += (*Ex_data)(k, i) * (*Ex_data)(k, i);
    u_val *= (*Sigma_data)(0, i);

    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      T_val[i1] = T_data->getPointer()[cell_node_idx[j]];
      nD[i1] = nD_data->getPointer()[cell_node_idx[j]];
      phi[i1] = phi_data->getPointer()[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k)
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
    }
    tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
    ele_vec->resize(n_vertex);
    for (int ii = 0; ii < n_vertex; ++ii)
      (*ele_vec)[ii] = 0.0;

    ele->thermal_buildElementRHS(vertex, dt, time, u_val, phi, E, T_val, nD, (*cell_flag)(0, i),
                                 ele_vec);
    for (int i2 = 0; i2 < n_vertex; ++i2)
      vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
  }
}

void PatchStrategy::buildthermalMatrixOnPatch(hier::Patch<NDIM> &patch, const double time,
                                              const double dt, const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, mat_data, d_thermal_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, Density_data, d_Density_id, Cell, double);
  GET_PATCH_DATA(patch, K_data, d_K_id, Cell, double);
  GET_PATCH_DATA(patch, Cp_data, d_Cp_id, Cell, double);
  GET_PATCH_DATA(patch, nD_data, d_DD_nD_temp_id, Node, double);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, cell_flag, d_Cell_flag_id, Cell, int);

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<double> T(n_vertex);
    tbox::Array<double> nD(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k)
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
      T[i1] = (*T_data)(0, cell_node_idx[j]);
      nD[i1] = nD_data->getPointer()[cell_node_idx[j]];
    }
    tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
    ele_mat->resize(n_vertex, n_vertex);
    for (int ii = 0; ii < n_vertex; ++ii)
      for (int jj = 0; jj < n_vertex; ++jj)
        (*ele_mat)(ii, jj) = 0.0;

    ele->thermal_buildStiffElementMatrix(vertex, dt, time, T, nD, (*cell_flag)(0, i), ele_mat);

    for (int ii = 0; ii < n_vertex; ++ii) {
      int row = mapping[ii];
      for (int jj = 0; jj < n_vertex; ++jj)
        mat_data->addMatrixValue(row, mapping[jj], (*ele_mat)(ii, jj));
    }
  }
  mat_data->assemble();
}

void PatchStrategy::buildTMatrixOnPatch(hier::Patch<NDIM> &patch, const double time,
                                        const double dt, const string &component_name) {
  /// 温度场矩阵组装 (T_MAT 构件): 热传导方程 ∫K(T,nD) ∇N_i·∇N_j dΩ.
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, mat_data, d_thermal_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, nD_data, d_DD_nD_temp_id, Node, double);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, cell_flag, d_Cell_flag_id, Cell, int);

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<double> T(n_vertex);
    tbox::Array<double> nD(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k)
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
      T[i1] = (*T_data)(0, cell_node_idx[j]);
      nD[i1] = nD_data->getPointer()[cell_node_idx[j]];
    }
    tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
    ele_mat->resize(n_vertex, n_vertex);
    for (int ii = 0; ii < n_vertex; ++ii)
      for (int jj = 0; jj < n_vertex; ++jj)
        (*ele_mat)(ii, jj) = 0.0;

    ele->buildTStiffElementMatrix(vertex, dt, time, T, nD, (*cell_flag)(0, i), ele_mat);

    for (int ii = 0; ii < n_vertex; ++ii) {
      int row = mapping[ii];
      for (int jj = 0; jj < n_vertex; ++jj)
        mat_data->addMatrixValue(row, mapping[jj], (*ele_mat)(ii, jj));
    }
  }
  mat_data->assemble();
}

void PatchStrategy::buildTRHSOnPatch(hier::Patch<NDIM> &patch, const double time,
                                     const double dt, const string &component_name) {
  /// 温度场右端项组装 (T_RHS 构件): 焦耳热源 ∫(|E|²σ)·N_i dΩ.
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  GET_PATCH_DATA(patch, vec_data, d_thermal_rhs_id, Vector, double);
  GET_PATCH_DATA(patch, Ex_data, d_Ex_id, Cell, double);
  GET_PATCH_DATA(patch, Sigma_data, d_Sigma_id, Cell, double);
  GET_PATCH_DATA(patch, nD_data, d_DD_nD_temp_id, Node, double);
  GET_PATCH_DATA(patch, phi_data, d_potential_plot_id, Node, double);
  GET_PATCH_DATA(patch, cell_flag, d_Cell_flag_id, Cell, int);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);

  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    tbox::Array<double> T_val(n_vertex);
    tbox::Array<double> nD(n_vertex);
    tbox::Array<double> phi(n_vertex);
    tbox::Array<double> E(NDIM);
    for (int k = 0; k < NDIM; k++)
      E[k] = (*Ex_data)(k, i);

    double u_val = 0;
    for (int k = 0; k < NDIM; k++)
      u_val += (*Ex_data)(k, i) * (*Ex_data)(k, i);
    u_val *= (*Sigma_data)(0, i);

    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      T_val[i1] = T_data->getPointer()[cell_node_idx[j]];
      nD[i1] = nD_data->getPointer()[cell_node_idx[j]];
      phi[i1] = phi_data->getPointer()[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k)
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
    }
    tbox::Pointer<tbox::Vector<double> > ele_vec = new tbox::Vector<double>();
    ele_vec->resize(n_vertex);
    for (int ii = 0; ii < n_vertex; ++ii)
      (*ele_vec)[ii] = 0.0;

    ele->buildTElementRHS(vertex, dt, time, u_val, phi, E, T_val, nD, (*cell_flag)(0, i),
                          ele_vec);
    for (int i2 = 0; i2 < n_vertex; ++i2)
      vec_data->addVectorValue(mapping[i2], (*ele_vec)[i2]);
  }
}

void PatchStrategy::applyTConstraint(hier::Patch<NDIM> &patch, const double time,
                                     const double dt, const string &component_name) {
  /// 温度场约束加载 (T_CONS 构件): 固定温度 Dirichlet 处理.
  (void)dt;
  GET_PATCH_DATA(patch, mat_data, d_thermal_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, vec_data, d_thermal_rhs_id, Vector, double);

  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  int *row_start = mat_data->getRowStartPointer();
  int *col_idx = mat_data->getColumnIndicesPointer();
  double *mat_val = mat_data->getValuePointer();
  double *vec_val = vec_data->getPointer();

  for (int temp_id = 0; temp_id < BC.num_temperature; temp_id++) {
    for (int temp_face = 0; temp_face < BC.FixTemperature_bc[temp_id].temp_face.size();
         temp_face++) {
      if (HAS_ENTITY_SET(patch, BC.FixTemperature_bc[temp_id].temp_face[temp_face], NODE, 1)) {
        DECLARE_ENTITY_SET(patch, entity_idx, BC.FixTemperature_bc[temp_id].temp_face[temp_face],
                           NODE, 1);
        for (int i = 0; i < entity_idx.size(); ++i) {
          int index = dof_map[entity_idx[i]];
          vec_val[index] = BC.FixTemperature_bc[temp_id].temperature;
          for (int j = row_start[index]; j < row_start[index + 1]; ++j) {
            if (col_idx[j] == index) {
              mat_val[j] = 1.0;
            } else {
              mat_val[j] = 0.0;
              vec_val[col_idx[j]] -= (*mat_data)(col_idx[j], index) * vec_val[index];
              (*mat_data)(col_idx[j], index) = 0.0;
            }
          }
        }
      }
    }
  }
}

// ===== 材料参数更新 =====

void PatchStrategy::calculateSigma(hier::Patch<NDIM> &patch, const double time, const double dt,
                                   const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, Sigma_data, d_Sigma_id, Cell, double);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, nD_data, d_nD_plot_id, Node, double);

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int imaterial = 0; imaterial < NumberOfMaterial; imaterial++) {
    for (int i = 0; i < material_entity[imaterial].size(); i++) {
      if (!HAS_ENTITY_SET(patch, material_entity[imaterial][i], CELL, 1))
        continue;
      DECLARE_ENTITY_SET(patch, cells, material_entity[imaterial][i], CELL, 1);

      for (int jcount = 0; jcount < cells.size(); jcount++) {
        int c = cells[jcount];
        int n_vertex_raw = cell_node_ext[c + 1] - cell_node_ext[c];
        int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<double> T(n_vertex);
        tbox::Array<double> nD(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        for (int i1 = 0, j = cell_node_ext[c]; i1 < n_vertex; ++i1, ++j) {
          mapping[i1] = dof_map[cell_node_idx[j]];
          for (int k = 0; k < NDIM; ++k)
            vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
          T[i1] = (*T_data)(0, cell_node_idx[j]);
          nD[i1] = (*nD_data)(0, cell_node_idx[j]);
        }
        tbox::Pointer<double> ele_Sigma = new double();
        (*ele_Sigma) = (*Sigma_data)(0, c);
        ele->calculateElementSigma(vertex, time, dt, T, nD, ele_Sigma);
        (*Sigma_data)(0, c) = (*ele_Sigma);
      }
    }
  }
}

void PatchStrategy::calculateK(hier::Patch<NDIM> &patch, const double time, const double dt,
                               const string &component_name) {
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, K_data, d_K_id, Cell, double);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, nD_data, d_nD_plot_id, Node, double);

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int imaterial = 0; imaterial < NumberOfMaterial; imaterial++) {
    for (int i = 0; i < material_entity[imaterial].size(); i++) {
      if (!HAS_ENTITY_SET(patch, material_entity[imaterial][i], CELL, 1))
        continue;
      DECLARE_ENTITY_SET(patch, cells, material_entity[imaterial][i], CELL, 1);

      for (int jcount = 0; jcount < cells.size(); jcount++) {
        int c = cells[jcount];
        int n_vertex_raw = cell_node_ext[c + 1] - cell_node_ext[c];
        int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
        tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
        tbox::Array<double> T(n_vertex);
        tbox::Array<double> nD(n_vertex);
        tbox::Array<int> mapping(n_vertex);
        for (int i1 = 0, j = cell_node_ext[c]; i1 < n_vertex; ++i1, ++j) {
          mapping[i1] = dof_map[cell_node_idx[j]];
          for (int k = 0; k < NDIM; ++k)
            vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
          nD[i1] = (*nD_data)(0, cell_node_idx[j]);
          T[i1] = (*T_data)(0, cell_node_idx[j]);
        }
        tbox::Pointer<double> ele_K = new double();
        (*ele_K) = (*K_data)(0, c);
        ele->calculateElementK(vertex, time, dt, T, nD, ele_K);
        (*K_data)(0, c) = (*ele_K);
      }
    }
  }
}

// ===== 后处理 =====

void PatchStrategy::postProcess(hier::Patch<NDIM> &patch, const double time, const double dt,
                                const string &component_name) {
  GET_PATCH_DATA(patch, plot_nd, d_potential_plot_id, Node, double);
  GET_PATCH_DATA(patch, vec, d_E_solution_id, Vector, double);
  GET_PATCH_DATA(patch, err, d_E_error_id, Vector, double);

  int num_nodes = patch.getNumberOfNodes(0);
  // 2D: DOF 编号稀疏(影像区/悬空节点无自由度), 必须经 dof_map 查表取值
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  double *vec_ptr = vec->getPointer();
  double *plot_ptr = plot_nd->getPointer();
  double *error_ptr = err->getPointer();
  for (int i = 0; i < num_nodes; ++i) {
    int mapping = dof_map[i];
    if (mapping == -1) continue;  // 无自由度节点不写解向量
    error_ptr[i] = abs(vec_ptr[mapping] - plot_ptr[i]);
    plot_ptr[i] = vec_ptr[mapping];
  }
}

void PatchStrategy::postEProcess(hier::Patch<NDIM> &patch, const double time, const double dt,
                                 const string &component_name) {
  /// 电场后处理 (E_POST 构件).
  /// 与 3D OneOrderPatchStrategy::postProcess 同构:
  /// 将解向量写入节点可视化数据片, 并计算新旧解之差作为误差.
  GET_PATCH_DATA(patch, plot_nd, d_potential_plot_id, Node, double);
  GET_PATCH_DATA(patch, vec, d_E_solution_id, Vector, double);
  GET_PATCH_DATA(patch, err, d_E_error_id, Vector, double);

  int num_nodes = patch.getNumberOfNodes(0);
  // 2D: DOF 编号稀疏(影像区/悬空节点无自由度), 必须经 dof_map 查表取值
  // 3D 中 DOF 编号恒等于节点索引, 直接 vec_ptr[i] 即可; 2D 用 vec_ptr[mapping].
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  double *vec_ptr = vec->getPointer();
  double *plot_ptr = plot_nd->getPointer();
  double *error_ptr = err->getPointer();
  for (int i = 0; i < num_nodes; ++i) {
    int mapping = dof_map[i];
    if (mapping == -1) continue;  // 无自由度节点不写解向量
    error_ptr[i] = abs(vec_ptr[mapping] - plot_ptr[i]);
    plot_ptr[i] = vec_ptr[mapping];
  }
}

void PatchStrategy::postthermalProcess(hier::Patch<NDIM> &patch, const double time, const double dt,
                                       const string &component_name) {
  GET_PATCH_DATA(patch, thermal_plot_nd, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, thermal_old_nd, d_thermal_old_id, Node, double);
  GET_PATCH_DATA(patch, thermal_vec, d_thermal_solution_id, Vector, double);
  GET_PATCH_DATA(patch, thermal_err, d_thermal_error_id, Vector, double);

  int num_nodes = patch.getNumberOfNodes(1);
  // 2D: DOF 编号稀疏, 解向量必须经 dof_map 查表; 无自由度节点不写解向量
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);
  double *thermal_vec_ptr = thermal_vec->getPointer();
  double *thermal_plot_ptr = thermal_plot_nd->getPointer();
  double *thermal_old_ptr = thermal_old_nd->getPointer();
  double *thermal_error_ptr = thermal_err->getPointer();
  for (int i = 0; i < num_nodes; ++i) {
    thermal_old_ptr[i] = thermal_plot_ptr[i];
    int mapping = dof_map[i];
    if (mapping == -1) continue;  // 无自由度节点不写解向量
    thermal_plot_ptr[i] = thermal_vec_ptr[mapping];
    thermal_error_ptr[i] = abs(thermal_plot_ptr[i] - thermal_old_ptr[i]);
  }
}

void PatchStrategy::postDDProcess(hier::Patch<NDIM> &patch, const double time, const double dt,
                                  const string &component_name) {
  GET_PATCH_DATA(patch, DD_plot_nd, d_nD_plot_id, Node, double);
  GET_PATCH_DATA(patch, DD_temp_nd, d_DD_nD_temp_id, Node, double);

  int num_nodes = patch.getNumberOfNodes(1);
  double *DD_temp_ptr = DD_temp_nd->getPointer();
  double *DD_plot_ptr = DD_plot_nd->getPointer();
  int *dis_ptr_DD = d_dof_info_DD->getDOFDistribution(patch, hier::EntityUtilities::NODE);
  for (int i = 0; i < num_nodes; ++i) {
    DD_plot_ptr[i] = DD_temp_ptr[i];
    if (dis_ptr_DD[i] < 1)
      DD_plot_ptr[i] = 0.;
  }
}

void PatchStrategy::postDDIterError(hier::Patch<NDIM> &patch, const double time, const double dt,
                                    const string &component_name) {
  /// 氧空穴数据片
  GET_PATCH_DATA(patch, DD_temp_nd, d_DD_nD_temp_id, Node, double);
  /// 氧空穴解向量数据片
  GET_PATCH_DATA(patch, DD_vec, d_DD_solution_id, Vector, double);
  /// 误差向量数据片
  GET_PATCH_DATA(patch, DD_err, d_DD_error_id, Vector, double);

  int *dof_map = d_dof_info_DD->getDOFMapping(patch, hier::EntityUtilities::NODE);
  int num_nodes = patch.getNumberOfNodes(1);
  double *DD_vec_ptr = DD_vec->getPointer();
  double *DD_temp_ptr = DD_temp_nd->getPointer();
  double *DD_error_ptr = DD_err->getPointer();
  for (int i = 0; i < num_nodes; ++i) {
    int mapping = dof_map[i];
    if (mapping == -1) {
      // 影像区/悬空节点无自由度: 解向量槽位不存在, 只更新节点数据片
      DD_temp_ptr[i] = 0;
      continue;
    }
    // 误差向量按 DOF 槽写入(mapping), 不能按节点索引 i:
    // DD 的 DOF 数(2880) < 节点数(4034), 按 i 写会越界损坏堆.
    DD_error_ptr[mapping] = abs(DD_vec_ptr[mapping] - DD_temp_ptr[i]);
    // 更新上一时刻氧空穴数据片
    DD_temp_ptr[i] = DD_vec_ptr[mapping];
  }
}

void PatchStrategy::buildEMatrixOnPatch(hier::Patch<NDIM> &patch, const double time,
                                        const double dt, const string &component_name) {
  /// 电场矩阵组装 (电流连续性方程 ∇·(σ∇V) = 0).
  /// 与 buildMatrixOnPatch 同构, 供 E_MAT 构件调用.
  tbox::Pointer<ElementManager<NDIM> > ele_manager = ElementManager<NDIM>::getManager();
  tbox::Pointer<BaseElement<NDIM> > ele = ele_manager->getElement(d_element_type);

  GET_COORD_DATA(patch, node_coord, Node);
  GET_PATCH_DATA(patch, mat_data, d_E_matrix_id, CSRMatrix, double);
  GET_PATCH_DATA(patch, Sigma_data, d_Sigma_id, Cell, double);
  GET_PATCH_DATA(patch, nD_data, d_DD_nD_temp_id, Node, double);
  GET_PATCH_DATA(patch, T_data, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, cell_flag, d_Cell_flag_id, Cell, int);

  DECLARE_ADJACENCY(patch, cell, node, Cell, Node);
  int num_cells = patch.getNumberOfCells(1);
  int *dof_map = d_dof_info->getDOFMapping(patch, hier::EntityUtilities::NODE);

  for (int i = 0; i < num_cells; ++i) {
    int n_vertex_raw = cell_node_ext[i + 1] - cell_node_ext[i];
      int n_vertex = (NDIM == 2) ? 3 : n_vertex_raw;
    tbox::Array<hier::DoubleVector<NDIM> > vertex(n_vertex);
    tbox::Array<double> T(n_vertex);
    tbox::Array<double> nD(n_vertex);
    tbox::Array<int> mapping(n_vertex);
    for (int i1 = 0, j = cell_node_ext[i]; i1 < n_vertex; ++i1, ++j) {
      mapping[i1] = dof_map[cell_node_idx[j]];
      for (int k = 0; k < NDIM; ++k) {
        vertex[i1][k] = (*node_coord)(k, cell_node_idx[j]);
      }
      T[i1] = (*T_data)(0, cell_node_idx[j]);
      nD[i1] = nD_data->getPointer()[cell_node_idx[j]];
    }

    tbox::Pointer<tbox::Matrix<double> > ele_mat = new tbox::Matrix<double>();
    ele_mat->resize(n_vertex, n_vertex);
    for (int ii = 0; ii < n_vertex; ++ii)
      for (int jj = 0; jj < n_vertex; ++jj)
        (*ele_mat)(ii, jj) = 0.0;

    ele->buildEStiffElementMatrix(vertex, dt, time, T, nD, (*cell_flag)(0, i), ele_mat);

    for (int ii = 0; ii < n_vertex; ++ii) {
      int row = mapping[ii];
      for (int jj = 0; jj < n_vertex; ++jj)
        mat_data->addMatrixValue(row, mapping[jj], (*ele_mat)(ii, jj));
    }
  }
  mat_data->assemble();
}

void PatchStrategy::calculateErrorOnPatch(hier::Patch<NDIM> &patch, const double time,
                                          const double dt, const string &component_name) {
  TBOX_WARNING("calculateErrorOnPatch: 2D version not yet fully implemented.");
}

// ===== 归约 =====

void PatchStrategy::reduceOnPatch(double *vector, int len, hier::Patch<NDIM> &patch,
                                  const double time, const double dt,
                                  const string &component_name) {
  GET_PATCH_DATA(patch, T_new_nd, d_temperature_plot_id, Node, double);
  GET_PATCH_DATA(patch, T_old_nd, d_thermal_old_id, Node, double);

  double maxdeltaT = 0;
  double *T_new_ptr = T_new_nd->getPointer();
  double *T_old_ptr = T_old_nd->getPointer();
  int num_nodes = patch.getNumberOfNodes(0);
  for (int i = 0; i < num_nodes; ++i) {
    double deltaT = abs(T_new_ptr[i] - T_old_ptr[i]);
    if (deltaT > maxdeltaT)
      maxdeltaT = deltaT;
  }
  vector[0] = maxdeltaT;
}

// ===== 输入/重启动 =====

void PatchStrategy::setParameter(tbox::Pointer<tbox::Database> input_db) {}

void PatchStrategy::getFromInput(tbox::Pointer<tbox::Database> db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!db.isNull());
#endif

  if (db->keyExists("element_type")) {
    d_element_type = db->getString("element_type");
  } else {
    TBOX_ERROR(d_object_name << ": No key `element_type' found in data.");
  }
  if (db->keyExists("integrator_type")) {
    d_integrator_type = db->getString("integrator_type");
  } else {
    TBOX_ERROR(d_object_name << ": No key `integrator_type' found in data.");
  }
  if (db->keyExists("shape_func_type")) {
    d_shape_func_type = db->getString("shape_func_type");
  } else {
    TBOX_ERROR(d_object_name << ": No key `shape_func_type' found in data.");
  }
  if (db->keyExists("constraint_types")) {
    d_constraint_types = db->getStringArray("constraint_types");
  } else {
    TBOX_ERROR(d_object_name << ": No key `constraint_types' found in data.");
  }
  if (db->keyExists("constraint_marks")) {
    d_constraint_marks = db->getIntegerArray("constraint_marks");
  } else {
    TBOX_ERROR(d_object_name << ": No key `constraint_marks' found in data.");
  }
  if (db->keyExists("Dt")) {
    Dt = db->getDouble("Dt");
  } else {
    TBOX_ERROR(d_object_name << ": No key `Dt' found in data.");
  }

  /// 读入材料参数数据
  if (db->isDatabase("Material")) {
    tbox::Pointer<tbox::Database> material = db->getDatabase("Material");
    NumberOfMaterial = material->getInteger("numofmaterial");
    tbox::Array<tbox::Array<int> > temp(NumberOfMaterial);
    for (int i = 0; i < NumberOfMaterial; i++) {
      stringstream materialname;
      materialname << "material" << i + 1;
      temp[i] = material->getIntegerArray(materialname.str());
    }
    material_entity = temp;

    for (int region = 1;; region++) {
      stringstream region_db_name;
      region_db_name << "region_" << region;
      if (!material->isDatabase(region_db_name.str()))
        break;
      Region r;
      tbox::Pointer<tbox::Database> region_db = material->getDatabase(region_db_name.str());
      r.d_id = region_db->getInteger("set_id");
      r.d_epsilonr = region_db->getDouble("epsilonr");
      r.d_mur = region_db->getDouble("mur");
      r.d_sigma = region_db->getDouble("sigma");
      r.d_K = region_db->getDouble("K");
      r.d_density = region_db->getDouble("density");
      r.d_Cp = region_db->getDouble("Cp");
      d_region_table[region] = r;
    }
  } else {
    tbox::pout << "No matching material database..." << endl;
    TBOX_ASSERT(-1);
  }

  if (db->isDatabase("BoundaryCondition")) {
    tbox::Pointer<tbox::Database> BC_db = db->getDatabase("BoundaryCondition");
    BC.T_initial = BC_db->getDouble("temperature_initial");
    BC.num_temperature = BC_db->getInteger("num_temperature");
    tbox::Array<FixTemperature> temp1(BC.num_temperature);
    for (int i = 0; i < BC.num_temperature; i++) {
      stringstream temperature_name;
      temperature_name << "temperature" << i + 1;
      tbox::Pointer<tbox::Database> temperature_name_db =
          BC_db->getDatabase(temperature_name.str());
      temp1[i].temp_face = temperature_name_db->getIntegerArray("temp_face");
      temp1[i].temperature = temperature_name_db->getDouble("temperature");
    }
    BC.FixTemperature_bc = temp1;

    BC.num_voltage = BC_db->getInteger("num_voltage");
    tbox::Array<Voltage> temp2(BC.num_voltage);
    for (int i = 0; i < BC.num_voltage; i++) {
      stringstream voltage_name;
      voltage_name << "voltage" << i + 1;
      tbox::Pointer<tbox::Database> voltage_name_db = BC_db->getDatabase(voltage_name.str());
      temp2[i].vol_face = voltage_name_db->getIntegerArray("vol_face");
      temp2[i].voltage = voltage_name_db->getDouble("voltage");
    }
    BC.voltage_bc = temp2;

    BC.num_timevoltage = BC_db->getInteger("num_timevoltage");
    tbox::Array<timeVoltage> temp3(BC.num_timevoltage);
    for (int i = 0; i < BC.num_timevoltage; i++) {
      stringstream timevoltage_name;
      timevoltage_name << "timevoltage" << i + 1;
      tbox::Pointer<tbox::Database> timevoltage_name_db =
          BC_db->getDatabase(timevoltage_name.str());
      temp3[i].vol_face = timevoltage_name_db->getIntegerArray("vol_face");
    }
    BC.timevoltage_bc = temp3;

    BC.num_concentration = BC_db->getInteger("num_concentration");
    tbox::Array<FixConcentration> temp4(BC.num_concentration);
    for (int i = 0; i < BC.num_concentration; i++) {
      stringstream concentration_name;
      concentration_name << "concentration" << i + 1;
      tbox::Pointer<tbox::Database> concentration_name_db =
          BC_db->getDatabase(concentration_name.str());
      temp4[i].con_face = concentration_name_db->getIntegerArray("con_face");
      temp4[i].concentration = concentration_name_db->getDouble("concentration");
    }
    BC.concentration_bc = temp4;
  } else {
    tbox::pout << "No matching BoundaryCondition database..." << endl;
    TBOX_ASSERT(-1);
  }
}

void PatchStrategy::getFromRestart(tbox::Pointer<tbox::Database> db) {
  getFromInput(db);
  tbox::Pointer<tbox::Database> root_db = tbox::RestartManager::getManager()->getRootDatabase();
  tbox::Pointer<tbox::Database> sub_db = root_db->getDatabase(d_object_name);
}
