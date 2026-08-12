//
// 文件名:     HeatLevelStrategy.C
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元网格层策略类.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "HeatLevelStrategy.h"
#include "PatchStrategy.h"
#ifdef DEBUG_CHECK_ASSERTIONS
#include <assert.h>
#endif
#include "JAUMINVector.h"
#include <cmath>
#include <iomanip>
#include <mpi.h>

/*************************************************************************
 * 构造函数.
 *************************************************************************/
HeatLevelStrategy::HeatLevelStrategy(
    const std::string &object_name,
    tbox::Pointer<algs::StandardComponentPatchStrategy<NDIM> > strategy,
    tbox::Pointer<tbox::Database> input_db) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!object_name.empty());
  assert(strategy.getPointer() != NULL);
#endif
  /// 创建网格片积分算法.
  d_patch_strategy = strategy;

  /// 从数据库中获取解法器类型.
  d_solver_db = input_db;
  d_solver_manager = solv::LinearSolverManager<NDIM>::getManager();
  d_E_solver = d_solver_manager->lookupLinearSolver(
      d_solver_db->getDatabase("SolverE")->getString("solver_name"));
  d_DD_solver = d_solver_manager->lookupLinearSolver(
      d_solver_db->getDatabase("SolverDD")->getString("solver_name"));
  d_T_solver = d_solver_manager->lookupLinearSolver(
      d_solver_db->getDatabase("SolverT")->getString("solver_name"));
  d_object_name = object_name;

  /// 全局后验误差初值
  d_global_DD_err = 0.;
  d_global_T_err = 0.;
  d_global_E_err = 0.;
}

/*************************************************************************
 * 析构函数.
 ************************************************************************/
HeatLevelStrategy::~HeatLevelStrategy() {}

/*************************************************************************
 * 初始化该积分算法: 创建所有计算需要的积分构件.
 *
 * 该函数创建了1个内存构件, 1个初始化构件, 4个数值构件.
 * 这些构件所操作的数据片,
 * 由函数 d_patch_strategy->initializeComponent() 指定.
 *
 *************************************************************************/
void HeatLevelStrategy::initializeLevelIntegrator(
    tbox::Pointer<algs::IntegratorComponentManager<NDIM> > manager) {
  // 初值构件: 管理当前值数据片的内存以及初始化
  d_init_set_value =
      new algs::InitializeIntegratorComponent<NDIM>("INIT", d_patch_strategy, manager);
  d_alloc_multiphysics_data =
      new algs::MemoryIntegratorComponent<NDIM>("ALLOC_MULTIPHYSICS", d_patch_strategy, manager);
  /// 数值构件：组装矩阵.
  d_mat_intc = new algs::NumericalIntegratorComponent<NDIM>("MAT", d_patch_strategy, manager);
  /// 数值构件：组装右端项.
  d_rhs_intc = new algs::NumericalIntegratorComponent<NDIM>("RHS", d_patch_strategy, manager);
  /// 数值构件：后处理计算.
  d_num_intc = new algs::NumericalIntegratorComponent<NDIM>("POST", d_patch_strategy, manager);
  /// 数值构件：计算SG电流密度.
  d_DD_J_SG_intc =
      new algs::NumericalIntegratorComponent<NDIM>("DD_J_SG", d_patch_strategy, manager);

  /// 电学构件：组装矩阵.
  d_E_mat_intc = new algs::NumericalIntegratorComponent<NDIM>("E_MAT", d_patch_strategy, manager);
  /// 电学构件：组装右端项.
  d_E_rhs_intc = new algs::NumericalIntegratorComponent<NDIM>("E_RHS", d_patch_strategy, manager);
  /// 电学构件：加载约束.
  d_E_cons_intc = new algs::NumericalIntegratorComponent<NDIM>("E_CONS", d_patch_strategy, manager);
  /// 电学构件：后处理.
  d_E_post_intc = new algs::NumericalIntegratorComponent<NDIM>("E_POST", d_patch_strategy, manager);
  /// 数值构件：计算电场模值.
  d_Ex_intc =
      new algs::NumericalIntegratorComponent<NDIM>("calculate_Ex", d_patch_strategy, manager);

  /// 温度场构件：组装矩阵.
  d_T_mat_intc = new algs::NumericalIntegratorComponent<NDIM>("T_MAT", d_patch_strategy, manager);
  /// 温度场构件：组装右端项.
  d_T_rhs_intc = new algs::NumericalIntegratorComponent<NDIM>("T_RHS", d_patch_strategy, manager);
  /// 温度场构件：加载约束.
  d_T_cons_intc = new algs::NumericalIntegratorComponent<NDIM>("T_CONS", d_patch_strategy, manager);
  /// 温度场构件：后处理 (计算误差向量并更新温度 plot, 与 3D PossionLevelStrategy 一致)
  d_T_post_intc = new algs::NumericalIntegratorComponent<NDIM>("thermal_POST", d_patch_strategy, manager);
  /// 误差估计构件：多物理场后验误差 (ZZ 梯度恢复, 与 3D ERROR_EST 构件同构)
  d_error_est_intc = new algs::NumericalIntegratorComponent<NDIM>("ERROR_EST", d_patch_strategy, manager);

  /// 归约构件：全局后验误差 (MPI_SUM, 与 3D PossionLevelStrategy 一致)
  d_reduction_intc = new algs::ReductionIntegratorComponent<NDIM>("RED", MPI_SUM, d_patch_strategy, manager);

  d_DD_mat_intc = new algs::NumericalIntegratorComponent<NDIM>("DD_MAT", d_patch_strategy, manager);
  d_DD_rhs_intc = new algs::NumericalIntegratorComponent<NDIM>("DD_RHS", d_patch_strategy, manager);
  d_DD_cons_intc =
      new algs::NumericalIntegratorComponent<NDIM>("DD_CONS", d_patch_strategy, manager);

  /// 数值构件：计算DD迭代误差
  d_DD_iter_post_intc =
      new algs::NumericalIntegratorComponent<NDIM>("DD_ITER_POST", d_patch_strategy, manager);
  /// 数值构件：循环结束后统一更新浓度场
  d_DD_post_intc =
      new algs::NumericalIntegratorComponent<NDIM>("DD_POST", d_patch_strategy, manager);
  /// 步长：计算时间步长.
  d_dt_intc = new algs::DtIntegratorComponent<NDIM>("Dt", d_patch_strategy, manager);
  /// 数值构件：处理约束边界条件.
  d_constraint_intc =
      new algs::NumericalIntegratorComponent<NDIM>("CONS", d_patch_strategy, manager);

  /// 初始化的数值构件
  d_Cell_flag_intc =
      new algs::NumericalIntegratorComponent<NDIM>("CELL_FLAG", d_patch_strategy, manager);
  d_edge_order_intc =
      new algs::EdgeOrderIntegratorComponent<NDIM>("EDGE_ORDER", d_patch_strategy, manager);
  d_com_edge_flag_intc =
      new algs::NumericalIntegratorComponent<NDIM>("COMM_EDGE_FLAG", d_patch_strategy, manager);
  d_int_geom_intc =
      new algs::NumericalIntegratorComponent<NDIM>("INIT_GEOM", d_patch_strategy, manager);
}

void HeatLevelStrategy::initializeLevelData(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                                            const double init_data_time, const bool initial_time) {
  tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;
  d_init_set_value->initializeLevelData(level, init_data_time, initial_time);
  d_Cell_flag_intc->computing(patch_level, init_data_time, 0.);
  d_edge_order_intc->setup(patch_level);
  d_com_edge_flag_intc->computing(patch_level, init_data_time, 0.);
  d_int_geom_intc->computing(patch_level, init_data_time, 0.);
}

/*************************************************************************
 * 获取网格层上的计算时间步长.
 *
 * @note: 求解问题如果是静态的, 则该函数将不被调用, 如果是时间发展的, 则该函数调
 * 用有限元构件的函数getLevelDt(),该构件又进一步自动调用
 * d_fem_strategy->getCellDt(), 完成时间步长的求解.
 *
 ************************************************************************/
double HeatLevelStrategy::getLevelDt(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                                     const double dt_time, const bool initial_time,
                                     const int flag_last_dt, const double last_dt) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!level.isNull());
#endif
  return (d_dt_intc->getLevelDt(level, dt_time, initial_time, flag_last_dt, last_dt, false));
}

/*************************************************************************
 * 向前积分一个时间步.
 *
 * 注解: 该函数调用了有限元构件对象的5个函数，
 * 分别完成计算矩阵右端项, 设置边界条件, 求解线性系统的计算, 误差计算.
 *
 ************************************************************************/
int HeatLevelStrategy::advanceLevel(const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
                                    const double current_time, const double predict_dt,
                                    const double max_dt, const double min_dt, const bool first_step,
                                    const int step_number, double &actual_dt) {
#ifdef DEBUG_CHECK_ASSERTIONS
  assert(!level.isNull());
#endif
  actual_dt = predict_dt;
  /// 从网格片算法类获取参数
  tbox::Pointer<PatchStrategy> p_strategy = d_patch_strategy;
  const tbox::Pointer<hier::PatchLevel<NDIM> > patch_level = level;
  int iter = 0;
  for (int niter = 0; niter < 1000; niter++) {
    d_alloc_multiphysics_data->allocatePatchData(patch_level, current_time + predict_dt);
    /// 浓度场方程求解
    tbox::pout << "开始浓度场的计算......" << endl;
    // S-G方法计算空穴电流密度
    d_DD_J_SG_intc->computing(patch_level, current_time, actual_dt);
    // Drift-Diffusion 矩阵组装
    d_DD_mat_intc->computing(patch_level, current_time, actual_dt);
    // Drift-Diffusion 右端项组装
    d_DD_rhs_intc->computing(patch_level, current_time, actual_dt);
    // Drift-Diffusion 矩阵数据片
    int DD_mat_id = p_strategy->getDDMatrixID();
    // Drift-Diffusion 右端项数据片
    int DD_vec_id = p_strategy->getDDRHSID();
    // Drift-Diffusion 解数据片
    int DD_sol_id = p_strategy->getDDSolutionID();
    // Drift-Diffusion 误差数据片
    int DD_error_id = p_strategy->getDDerrorID();
    d_DD_solver->setMatrix(DD_mat_id);
    d_DD_solver->setRHS(DD_vec_id);
    d_DD_solver->solve(first_step, DD_sol_id, patch_level, d_solver_db->getDatabase("SolverDD"));
    d_DD_iter_post_intc->computing(patch_level, current_time, actual_dt, false);
    d_DD_sol = new JPSOL::JVector<NDIM, double>(patch_level, DD_sol_id);
    d_DD_error = new JPSOL::JVector<NDIM, double>(patch_level, DD_error_id);
    // 计算相对误差收敛判据
    // 计算解向量的L2范数
    double DD_sol_l2norm = d_DD_sol->l2Norm();
    // 计算误差向量的L2范数
    double DD_error_l2norm = d_DD_error->l2Norm();
    double DD_REC = DD_error_l2norm / (EPS + DD_sol_l2norm);
    tbox::pout << "本次浓度场残差: " << DD_REC << endl;
    tbox::pout << "结束浓度场的计算......" << endl;
    /// 电场方程求解
    tbox::pout << "开始电场的计算......" << endl;
    /// Electric 矩阵组装
    d_E_mat_intc->computing(patch_level, current_time, actual_dt);
    /// Electric 右端项组装
    d_E_rhs_intc->computing(patch_level, current_time, actual_dt);
    /// Electric 约束加载
    d_E_cons_intc->computing(patch_level, current_time, actual_dt);
    int E_mat_id = p_strategy->getEMatrixID();
    // Electric 右端项数据片
    int E_vec_id = p_strategy->getERHSID();
    // Electric 解数据片
    int E_sol_id = p_strategy->getESolutionID();
    // Electric 误差数据片
    int E_error_id = p_strategy->getEerrorID();
    d_E_solver->setMatrix(E_mat_id);
    d_E_solver->setRHS(E_vec_id);
    d_E_solver->solve(first_step, E_sol_id, patch_level, d_solver_db->getDatabase("SolverE"));
    // 将解向量写入节点数据片并计算误差
    d_E_post_intc->computing(patch_level, current_time, actual_dt, false);
    // 计算电场模值
    d_Ex_intc->computing(patch_level, current_time, actual_dt);
    d_E_sol = new JPSOL::JVector<NDIM, double>(patch_level, E_sol_id);
    d_E_error = new JPSOL::JVector<NDIM, double>(patch_level, E_error_id);
    // 计算解向量的L2范数
    double E_sol_l2norm = d_E_sol->l2Norm();
    // 计算误差向量的L2范数
    double E_error_l2norm = d_E_error->l2Norm();
    double E_REC = E_error_l2norm / (EPS + E_sol_l2norm);
    tbox::pout << "本次电场残差: " << E_REC << endl;
    /// 温度场求解
    tbox::pout << "结束电场的计算......" << endl;
    tbox::pout << "开始温度场的计算......" << endl;
    /// T 矩阵组装
    d_T_mat_intc->computing(patch_level, current_time, actual_dt);
    /// T 右端项组装
    d_T_rhs_intc->computing(patch_level, current_time, actual_dt);
    /// T 约束加载
    d_T_cons_intc->computing(patch_level, current_time, actual_dt);
    // T 矩阵数据片
    int T_mat_id = p_strategy->getthermalMatrixID();
    // T 右端项数据片
    int T_vec_id = p_strategy->getthermalRHSID();
    // T 解数据片
    int T_sol_id = p_strategy->getthermalSolutionID();
    // T 误差数据片
    int T_error_id = p_strategy->getthermalerrorID();
    d_T_solver->setMatrix(T_mat_id);
    d_T_solver->setRHS(T_vec_id);
    d_T_solver->solve(first_step, T_sol_id, patch_level, d_solver_db->getDatabase("SolverT"));
    // 温度场后处理: 计算误差向量并更新温度 plot (此前缺失, 导致 T_REC 恒为 0)
    d_T_post_intc->computing(patch_level, current_time, actual_dt, false);
    d_T_sol = new JPSOL::JVector<NDIM, double>(patch_level, T_sol_id);
    d_T_error = new JPSOL::JVector<NDIM, double>(patch_level, T_error_id);
    // 计算解向量的L2范数
    double T_sol_l2norm = d_T_sol->l2Norm();
    // 计算误差向量的L2范数
    double T_error_l2norm = d_T_error->l2Norm();
    double T_REC = T_error_l2norm / (EPS + T_sol_l2norm);
    tbox::pout << "结束温度场的计算......" << endl;
    // 收敛判据 
    bool converged = (DD_REC < 1e-5) && (E_REC < 1e-3) && (T_REC < 1e-3);
    tbox::pout << "============================================================" << endl;
    tbox::pout << "  Iteration " << setw(4) << niter
               << " | DD_REC = " << scientific << setprecision(3) << DD_REC
               << " | E_REC  = " << E_REC
               << " | T_REC  = " << T_REC
               << " | " << (converged ? "CONVERGED" : "not converged")
               << fixed << endl;
    tbox::pout << "============================================================" << endl;
    if (converged) {
      tbox::pout << "当前时间步均已收敛" << endl;
      tbox::pout << "DD_REC: " << DD_REC << " E_REC: " << E_REC
                 << " T_REC: " << T_REC << endl;
      iter = niter;
      break;
    }
    // 每轮未收敛时释放多物理场数据片, 下一轮重新分配 (与 3D PossionLevelStrategy::advanceLevel
    // 第337行一致): 否则 DD_RHS 等向量/矩阵数据片跨轮累积不清零, 每轮解 = 上一轮解+固定增量,
    // DD_REC 按 1/(k+2) 调和级数衰减, 外层迭代永不收敛.
    d_alloc_multiphysics_data->deallocatePatchData(patch_level);
  }
  tbox::pout << "外层迭代次数: " << iter << endl;
  if (iter == 1000) tbox::pout << "未收敛(达到最大迭代步数)" << endl;
  /// 循环结束后统一更新浓度场 (DD_temp → DD_plot 滚动到下一时间步)
  d_DD_post_intc->computing(patch_level, current_time, actual_dt, false);

  /// 每个时间步的多物理场后验误差估计 (ZZ 梯度恢复, 逐单元)
  d_error_est_intc->computing(patch_level, current_time, actual_dt, false);

  /// 全局后验误差: ‖e‖ = (Σ_e η²_e)^{1/2}, 跨进程 MPI_SUM 归约
  /// (与 3D PossionLevelStrategy 一致: "RED" 构件 → reduceOnPatch 逐 patch 累加局部 Ση²)
  double glob[3] = {0., 0., 0.};
  d_reduction_intc->reduction(glob, 3, patch_level, current_time, actual_dt);
  d_global_DD_err = sqrt(glob[0]);
  d_global_T_err = sqrt(glob[1]);
  d_global_E_err = sqrt(glob[2]);
  tbox::pout << "全局后验误差 ‖e‖: DD = " << d_global_DD_err
             << "  T = " << d_global_T_err
             << "  E = " << d_global_E_err << endl;

  // 循环结束后释放多物理场数据片内存 (与 3D PossionLevelStrategy::advanceLevel
  // 第347行一致): 否则 RHS 残留会污染下一个时间步的第一轮迭代.
  d_alloc_multiphysics_data->deallocatePatchData(patch_level);

  actual_dt = predict_dt;


  return (0);
}
