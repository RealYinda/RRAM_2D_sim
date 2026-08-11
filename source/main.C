//
// 文件名: main.C
// 软件包: 2D FEM application
// 版权  : (c) 2004-2010 北京应用物理与计算数学研究所
// 版本号: $Revision: 136 $
// 修改  : $Date: 2011-07-22 08:41:30 +0800 (五, 2011-07-22) $
// 描述  : 主控程序(非结构网格上, 求解Heat方程).
//

#include "GridGeometry.h"
#include "PatchHierarchy.h"
#include "HierarchyTimeIntegrator.h"
#include "JaVisDataWriter.h"
#include "JAUMINManager.h"
#include "InputManager.h"
#include "RestartManager.h"
#include "TimerManager.h"
#include "HeatLevelStrategy.h"
#include "PatchStrategy.h"

using namespace JAUMIN;

/*!
*************************************************************************
*
* @brief 2D RRAM 电热多物理场有限元求解器.
*
* 该程序分以下几个步骤:
* -# 预处理: 初始化计算环境, 解析输入文件, 读取主程序控制参数;
* -# 创建网格层时间积分算法类对象, 主要包括:
*    网格层(单块) hier::PatchLevel<NDIM>
*    -# 有限元算法 Heat
*    -# 网格层积分算法 FEMLevelIntegrator
*    -# 网格层时间积分算法 algs::HierarchyTimeIntegrator<NDIM>
* -# 初始化网格片层次结构和物理量数据片;
* -# 组装矩阵和右端项并求解线性系统;
* -# 后处理: 释放应用类对象, 释放计算资源.
*
************************************************************************
*/

static void prefixInputDirName(const std::string& input_filename,
                               tbox::Pointer<tbox::Database> input_db) {
  std::string path_name = "";
  std::string::size_type slash_pos = input_filename.find_last_of('/');
  if (slash_pos != std::string::npos)
    path_name = input_filename.substr(0, slash_pos + 1);

  std::string mesh_file = input_db->getDatabase("GridGeometry")
                         ->getDatabase("MeshImportationParameter")
                         ->getString("file_name");

  slash_pos = mesh_file.find_first_of('/');
  if (slash_pos != 0) {
    input_db->getDatabase("GridGeometry")
        ->getDatabase("MeshImportationParameter")
        ->putString("file_name", path_name + mesh_file);
  }
}

int main(int argc, char* argv[]) {
  // 初始化计算环境.
  tbox::JAUMINManager::startup();
  {
    /*******************************************************************************
     *                               预  处  理 *
     *******************************************************************************/

    // 解析命令行参数:
    bool is_from_restart = false;  // 此次运行是否为重启动.
    std::string input_filename;         // 输入文件名.
    std::string restart_read_dirname;   // 重启动文件所在路径
    int restore_num = 0;           // 重启动时间序号
    if ((argc != 2) && (argc != 4)) {
      tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
                 << "<restart dir> <restore number> [options]\n"
                 << "  options:\n"
                 << "  none at this time" << std::endl;
      tbox::MPI::abort();
      return (-1);
    } else {
      input_filename = argv[1];

      if (argc == 4) {  // 重启动情形
        is_from_restart = true;
        restart_read_dirname = argv[2];  // 重启动文件所在路径
        restore_num = atoi(argv[3]);     // 重启动时间序号
      }
    }

    /// 把信息输出到log文件
    tbox::plog << "input_filename = " << input_filename << std::endl;
    if (is_from_restart) {
      tbox::plog << "restart_read_dirname = " << restart_read_dirname << std::endl;
      tbox::plog << "restore_num = " << restore_num << std::endl;
    }

    // 解析输入文件的计算参数到输入数据库, 称之为根数据库.
    tbox::Pointer<tbox::Database> input_db =
        new tbox::InputDatabase("input_db");
    tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);
    prefixInputDirName(input_filename, input_db);

    // 从根数据库中获得名称为"Main"的子数据库.
    tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

    // 从"Main"子数据库中获取日志文件控制参数.
    std::string log_file_name = main_db->getString("log_file_name");
    bool log_all_nodes = main_db->getBoolWithDefault("log_all_nodes", false);
    if (log_all_nodes) {
      tbox::PIO::logAllNodes(log_file_name);
    } else {
      tbox::PIO::logOnlyNodeZero(log_file_name);
    }

    // 从"Main"子数据库中获取可视化输出控制参数.
    int javis_dump_interval =
        main_db->getIntegerWithDefault("javis_dump_interval", 0);
    std::string javis_dump_dirname;
    int javis_number_procs_per_file = 0;
    if (javis_dump_interval > 0) {
      javis_dump_dirname = main_db->getString("javis_dump_dirname");
      javis_number_procs_per_file =
          main_db->getIntegerWithDefault("javis_number_procs_per_file", 1);
    }

    // 从"Main"子数据库中获取重启动输出控制参数.
    int restart_dump_interval =
        main_db->getIntegerWithDefault("restart_dump_interval", 0);
    std::string restart_dump_dirname;
    if (restart_dump_interval > 0) {
      restart_dump_dirname = main_db->getString("restart_dump_dirname");
    }

    // 打开重启动输入文件
    tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
    if (is_from_restart) {
      restart_manager->openRestartFile(restart_read_dirname, restore_num,
                                       tbox::MPI::getNodes());
    }

    /// 创建网格几何对象.
    tbox::Pointer<hier::GridGeometry<NDIM> > grid_geometry =
        new hier::GridGeometry<NDIM>("GridGeometry",
                                     input_db->getDatabase("GridGeometry"));

    /// 创建网格拓扑对象.
    tbox::Pointer<hier::GridTopology<NDIM> > grid_topology =
        new hier::GridTopology<NDIM>("GridTopology",
                                     input_db->getDatabase("GridTopology"));

    /// 创建网格片层次结构.
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
        new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry,
                                       grid_topology, true);

    /*******************************************************************************
     *                     创建网格层和积分算法类对象 *
     *******************************************************************************/

    /// 有限元计算输入数据库.
    tbox::Pointer<tbox::Database> fem_db = input_db->getDatabase("FEM");

    /// 创建网格片算法策略类.
    tbox::Pointer<PatchStrategy> poisson_patch_strategy =
        new PatchStrategy("HeatPatchIntegrator", is_from_restart, fem_db);

    /// 创建网格层时间积分算法类（应用级: 提供有限元求解流程定制支撑）.
    tbox::Pointer<HeatLevelStrategy> poisson_level_integrator =
        new HeatLevelStrategy("HeatLevelIntegrator", poisson_patch_strategy,
                              fem_db->getDatabase("Solver"));

    /// 创建网格层时间积分算法类.
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            input_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            poisson_level_integrator, true);

    /*******************************************************************************
     *                初 始 化 网 格 层 和 物 理 量 *
     *******************************************************************************/
    time_integrator->initializeHierarchy();

    /*******************************************************************************
     *  JaVis可视化输出器类. *
     *******************************************************************************/
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer;
    if (javis_dump_interval > 0) {
      javis_data_writer = new appu::JaVisDataWriter<NDIM>(
          "Heat_JaVis_Writer", javis_dump_dirname, javis_number_procs_per_file);

      poisson_patch_strategy->registerPlotData(
          javis_data_writer); /**< 注册待输出的绘图量.  */
    }

    // 输出初始化信息到日志文件.
    tbox::plog << "\nCheck input data and variables before simulation:" << std::endl;
    tbox::plog << "Input database..." << std::endl;
    input_db->printClassData(tbox::plog);
    tbox::plog << "\nVariable database..." << std::endl;
    hier::VariableDatabase<NDIM>::getDatabase()->printClassData(tbox::plog);
    if (javis_dump_interval > 0) {
      javis_data_writer->writePlotData(time_integrator->getPatchHierarchy(), 0,
                                       1.0);
    }

    // 关闭重启动输入文件
    if (is_from_restart) {
      restart_manager->closeRestartFile();
    }

    double loop_time = time_integrator->getIntegratorTime();
    double loop_time_end = time_integrator->getEndTime();
    int iteration_num = time_integrator->getIntegratorStep();

    while ((loop_time < loop_time_end) && time_integrator->stepsRemaining()) {
      iteration_num = time_integrator->getIntegratorStep() + 1;

      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      tbox::pout << "At begining of solving  step :   " << iteration_num
                 << std::endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;
      double dt_actual;

      // 将数值解推进一个时间步, 返回其时间步长.
      dt_actual = time_integrator->advanceHierarchy();

      loop_time += dt_actual;

      tbox::pout << "Dt = " << dt_actual << ", Simulation time is " << loop_time
                 << std::endl;
      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
      tbox::pout << "At end of solving # " << std::endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << std::endl;

      if ((restart_dump_interval > 0) &&
          ((iteration_num % restart_dump_interval) == 0)) {
        restart_manager->writeRestartFile(restart_dump_dirname, iteration_num);
      }

      // 输出可视化数据.

      if ((javis_dump_interval > 0) &&
          ((iteration_num % javis_dump_interval) == 0 ||
           !((loop_time < loop_time_end) &&
             time_integrator->stepsRemaining()))) {
        javis_data_writer->writePlotData(time_integrator->getPatchHierarchy(),
                                         iteration_num, loop_time);
      }
    }

    /************************************************************************************
     *                               模  拟  结  束 *
     ************************************************************************************/
    // 输出计时器统计的时间数据.
    tbox::TimerManager::getManager()->print(tbox::plog);
  }
  // 释放计算资源.
  tbox::JAUMINManager::shutdown();

  return (0);
}
