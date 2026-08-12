//
// 文件名:     HeatLevelStrategy.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元网格层策略类派生类.
// 类别　:     %Internal File% ( Don't delete this line )
//
//

#ifndef included_HeatLevelStrategy
#define included_HeatLevelStrategy

#include "BaseLinearSolver.h"
#include "DtIntegratorComponent.h"
#include "EdgeOrderIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "LinearSolverManager.h"
#include "MemoryIntegratorComponent.h"
#include "NumericalIntegratorComponent.h"
#include "PatchLevel.h"
#include "ReductionIntegratorComponent.h"
#include "StandardComponentPatchStrategy.h"
#include "TimeIntegratorLevelStrategy.h"
using namespace JAUMIN;

/**
 * @brief 该类从网格层时间积分算法策略类 algs::TimeIntegratorLevelStrategy 派生,
 *        实现有限元方法的流程控制.
 */
class HeatLevelStrategy : public algs::TimeIntegratorLevelStrategy<NDIM> {
public:
  /**
   * @brief 构造函数.
   * @param object_name   输入参数, 字符串, 表示对象名称.
   * @param strategy      输入参数, 指针, 网格片策略类派生类对象.
   * @param input_db      输入参数, 数据库指针, 指向解法器的数据库.
   *
   */
  HeatLevelStrategy(const std::string &object_name,
                    tbox::Pointer<algs::StandardComponentPatchStrategy<NDIM> > strategy,
                    tbox::Pointer<tbox::Database> input_db);

  /**
   * @brief 析构函数.
   */
  virtual ~HeatLevelStrategy();

  ///@name 重载基类algs::TimeIntegratorLevelStrategy<NDIM>的函数
  //@{
  /**
   * @brief 初始化该积分算法: 创建所有计算需要的积分构件.
   * @param manager 输入参数, 指针, 指向积分构件管理器.
   */
  void initializeLevelIntegrator(tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager);

  /**
   * @brief 初始化指定网格层的数据片.
   *
   * 具体地，该函数完成以下操作：
   * - 若输入参数 initial_time 为真，则根据初始条件，
   *   为当前网格层上的所有数据片 <uval,current> 赋初值;
   * - 若输入参数 initial_time 为假（此时刚完成负载调整），
   *   则将数据片 <uval,current> 的值从旧网格层复制到新网格层.
   *
   * @param level 输入参数, 指针, 指向待初始化网格层.
   * @param init_data_time 输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   *
   * @note
   * 该函数调用了1个初值构件对象，该对象又进一步自动调用函数
   * Euler::initializePatchData(), 逐个网格片地完成初始时刻的数据初始化.
   */
  void initializeLevelData(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                           const double init_data_time, const bool initial_time = true);

  /**
   * @brief 返回指定网格层的时间步长.
   *
   * @param level 输入参数, 指针, 指向网格层.
   * @param dt_time 输入参数, 双精度浮点型, 表示计算时间步长的当前时刻.
   * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param flag_last_dt 输入参数, 整型, 表示上个时间步积分返回的状态.
   * @param last_dt 输入参数, 双精度浮点型, 表示上个时间步长.
   * @return 双精度浮点型, 表示网格层的时间步长.
   *
   * @note
   * 该函数调用了1个时间步长构件对象，该对象又进一步自动调用函数
   * FEM::getPatchDt(), 逐个网格片地计算稳定时间步长.
   */
  double getLevelDt(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level, const double dt_time,
                    const bool initial_time, const int flag_last_dt, const double last_dt);

  /**
   * @brief 网格层向前积分一个时间步.
   *
   * @param level 输入参数, 指针, 指向待积分的网格层.
   * @param current_time 输入参数, 双精度浮点型, 表示时间步的起始时刻.
   * @param predict_dt 输入参数, 双精度浮点型, 表示为该时间步预测的时间步长.
   * @param max_dt 输入参数, 双精度浮点型, 表示时间步允许的最大时间步长.
   * @param min_dt 输入参数, 双精度浮点型, 表示时间步允许的最小时间步长.
   * @param first_step 输入参数, 逻辑型, 真值当前为重构后或时间步序列的第1步.
   * @param step_number, 输入参数, 整型, 表示积分步数.
   * @param actual_dt 输出参数, 双精度浮点型, 表示时间步实际采用的时间步长.
   *
   * @return 整型, 表示该时间步积分的状态.
   *
   */
  int advanceLevel(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                   const double current_time, const double predict_dt, const double max_dt,
                   const double min_dt, const bool first_step, const int step_number,
                   double &actual_dt);

  ///@name 全局后验误差 ‖e‖ = (Σ_e η²_e)^{1/2} (每个时间步收敛后经 MPI_SUM 归约得到)
  //@{
  double getGlobalDDError() const { return d_global_DD_err; }
  double getGlobalTError() const { return d_global_T_err; }
  double getGlobalEError() const { return d_global_E_err; }
  //@}
private:
  /*!@brief 对象名称. */
  std::string d_object_name;

  /*!@brief 全局后验误差: 漂移扩散/温度/电场 (每个时间步更新) */
  double d_global_DD_err, d_global_T_err, d_global_E_err;

  /// 解法器数据库.
  tbox::Pointer<tbox::Database> d_solver_db;

  /*!@brief 有限元网格片策略类. */
  tbox::Pointer<algs::StandardComponentPatchStrategy<NDIM> > d_patch_strategy;

  /*!@brief 解法器: 求解矩阵系统 */
  tbox::Pointer<solv::LinearSolverManager<NDIM> > d_solver_manager;
  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_solver;

  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_E_solver;
  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_DD_solver;
  tbox::Pointer<solv::BaseLinearSolver<NDIM> > d_T_solver;

  /*!@brief 初值构件: 管理物理量当前值数据片的内存以及初始化 */
  tbox::Pointer<algs::InitializeIntegratorComponent<NDIM> > d_init_set_value;
  /*!@brief 矩阵构件: 半导体单元标识数据片 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Cell_flag_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_com_edge_flag_intc;
  /*!@brief 矩阵构件: 组装矩阵 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_mat_intc;

  /*!@brief 右端项构件: 组装右端项 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_rhs_intc;

  /*!@brief 内存构件: 管理多物理场矩阵/向量数据片 */
  tbox::Pointer<algs::MemoryIntegratorComponent<NDIM> > d_alloc_multiphysics_data;

  /*!@brief 数值构件 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_num_intc;
  /*!@brief 数值构件 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_constraint_intc;

  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_DD_J_SG_intc;

  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_DD_mat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_DD_rhs_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_DD_cons_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_DD_iter_post_intc;
  /*!@brief 数值构件: 循环结束后统一更新浓度场 */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_DD_post_intc;

  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_E_mat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_E_rhs_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_E_cons_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_E_post_intc;
  /*!@brief 数值构件: 计算电场模值 Ex */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_Ex_intc;

  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_T_mat_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_T_rhs_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_T_cons_intc;
  /*!@brief 数值构件: 温度场后处理 (误差向量 + plot 更新) */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_T_post_intc;
  /*!@brief 数值构件: 多物理场后验误差估计 (ZZ 梯度恢复) */
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_error_est_intc;

  /*!@brief 归约构件: 全局后验误差 (MPI_SUM, 与 3D PossionLevelStrategy 一致) */
  tbox::Pointer<algs::ReductionIntegratorComponent<NDIM> > d_reduction_intc;



  tbox::Pointer<algs::EdgeOrderIntegratorComponent<NDIM> > d_edge_order_intc;
  tbox::Pointer<algs::NumericalIntegratorComponent<NDIM> > d_int_geom_intc;

  /*!@brief 步长构件 */
  tbox::Pointer<algs::DtIntegratorComponent<NDIM> > d_dt_intc;

  /*!@brief JPSOL向量 */
  tbox::Pointer<JPSOL::JVector<NDIM, double> > d_DD_sol;
  tbox::Pointer<JPSOL::JVector<NDIM, double> > d_DD_error;

  tbox::Pointer<JPSOL::JVector<NDIM, double> > d_E_sol;
  tbox::Pointer<JPSOL::JVector<NDIM, double> > d_E_error;

  tbox::Pointer<JPSOL::JVector<NDIM, double> > d_T_sol;
  tbox::Pointer<JPSOL::JVector<NDIM, double> > d_T_error;
};
#endif
