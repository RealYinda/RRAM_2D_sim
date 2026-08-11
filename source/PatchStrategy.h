//
// 文件名:     PatchStrategy.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: Tue May 20 08:45:12 2014 $
// 描述　:     网格片策略类派生类, RRAM 多物理场 (电-热-漂移扩散) 2D 版本.
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_PatchStrategy
#define included_PatchStrategy

#include "StandardComponentPatchStrategy.h"
#include "JaVisDataWriter.h"
#include "DOFInfo.h"
#include "RestartManager.h"
#include "UserConstant.h"

using namespace JAUMIN;

class PatchStrategy : public algs::StandardComponentPatchStrategy<NDIM>,
                      tbox::Serializable {
public:
  /*! @brief 构造函数.
   * @param object_name          输入参数, 字符串, 表示对象名称.
   */
  PatchStrategy(const std::string& object_name, bool is_from_restart,
                tbox::Pointer<tbox::Database> input_db);

  /*!
   * @brief 析构函数.
   */
  virtual ~PatchStrategy();

  /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
  // @{

  /*!
   * @brief 初始化指定的积分构件.
   *
   * 注册待填充的数据片或待调度内存空间的数据片到积分构件.
   *
   * @param component 输入参数, 指针, 指向待初始化的积分构件对象.
   */
  void initializeComponent(algs::IntegratorComponent<NDIM>* component) const;

  /**
   * @brief 初始化数据片（支持有限元初值构件）.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示初始化的时刻.
   * @param initial_time   输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
   * @param component_name 输入参数, 字符串, 当前调用该函数的初值构件之名称.
   */
  void initializePatchData(hier::Patch<NDIM>& patch, const double time,
                           const bool initial_time,
                           const std::string& component_name);

  /*************************************************************************
   * 矢量有限元信息生成
   ************************************************************************/
  void initGeometryOnPatch(hier::Patch<NDIM>& patch);

  /*************************************************************************
   * 电处理函数
   ************************************************************************/
  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成矩阵组装（电流连续性）.
   */
  void buildMatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                          const double dt, const std::string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成右端项组装（电流连续性）.
   */
  void buildRHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                       const double dt, const std::string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成右端项组装（电学, E_RHS）.
   */
  void buildERHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const std::string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成约束加载（电压边界）.
   */
  void applyConstraint(hier::Patch<NDIM>& patch, const double time,
                       const double dt, const std::string& component_name);

  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成约束加载（电学, E_CONS）.
   */
  void applyEConstraint(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const std::string& component_name);

  void calculateEx(hier::Patch<NDIM>& patch, const double time, const double dt,
                   const std::string& component_name);

  // 求解时数据片获得
  int getEMatrixID() { return d_E_matrix_id; }
  int getERHSID() { return d_E_rhs_id; }
  int getESolutionID() { return d_E_solution_id; }
  int getEerrorID() { return d_E_error_id; }

  /*************************************************************************
   * 漂移扩散处理函数
   ************************************************************************/
  void buildDDMatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                            const double dt, const std::string& component_name);

  void buildDDRHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                         const double dt, const std::string& component_name);

  void applyDDConstraint(hier::Patch<NDIM>& patch, const double time,
                         const double dt, const std::string& component_name);

  // S-G方法计算载流子电流密度
  void calculateDD_J(hier::Patch<NDIM>& patch, const double time, const double dt,
                     const std::string& component_name);

  // 求解时数据片获得
  int getDDMatrixID() { return d_DD_matrix_id; }
  int getDDRHSID() { return d_DD_rhs_id; }
  int getDDSolutionID() { return d_DD_solution_id; }
  int getDDerrorID() { return d_DD_error_id; }


  /*************************************************************************
   * 电场处理函数
   ************************************************************************/
  void buildEMatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                            const double dt, const std::string& component_name);

  /*************************************************************************
   * 热处理函数
   ************************************************************************/
  void applythermalConstraint(hier::Patch<NDIM>& patch, const double time,
                              const double dt, const std::string& component_name);

  void buildthermalMatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                                 const double dt, const std::string& component_name);

  void buildthermalRHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                              const double dt, const std::string& component_name);

  /*!
   * @brief 温度场矩阵组装 (热传导, T_MAT 构件).
   */
  void buildTMatrixOnPatch(hier::Patch<NDIM>& patch, const double time,
                           const double dt, const std::string& component_name);

  /*!
   * @brief 温度场右端项组装 (焦耳热源, T_RHS 构件).
   */
  void buildTRHSOnPatch(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const std::string& component_name);

  /*!
   * @brief 温度场约束加载 (固定温度, T_CONS 构件).
   */
  void applyTConstraint(hier::Patch<NDIM>& patch, const double time,
                        const double dt, const std::string& component_name);

  // 求解时数据片获得
  int getthermalMatrixID() { return d_thermal_matrix_id; }
  int getthermalRHSID() { return d_thermal_rhs_id; }
  int getthermalSolutionID() { return d_thermal_solution_id; }
  int getthermalerrorID() { return d_thermal_error_id; }

  /*************************************************************************
   * 材料处理函数
   ************************************************************************/
  void calculateSigma(hier::Patch<NDIM>& patch, const double time,
                      const double dt, const std::string& component_name);

  void calculateK(hier::Patch<NDIM>& patch, const double time,
                  const double dt, const std::string& component_name);

  /*************************************************************************
   * 前处理函数
   ************************************************************************/
  // 单元标识
  void computeCellFlagOnPatch(hier::Patch<NDIM>& patch);

  /*************************************************************************
   * 后处理函数
   ************************************************************************/
  /*!
   * @brief 支撑指定名称的数值构件, 在单个网格片上完成后处理计算.
   */
  void postProcess(hier::Patch<NDIM>& patch, const double time, const double dt,
                   const std::string& component_name);

  /*!
   * @brief 电场后处理 (E_POST 构件): 解向量 → 节点数据片 + 误差.
   */
  void postEProcess(hier::Patch<NDIM>& patch, const double time, const double dt,
                    const std::string& component_name);

  void postDDProcess(hier::Patch<NDIM>& patch, const double time, const double dt,
                     const std::string& component_name);

  void postDDIterError(hier::Patch<NDIM>& patch, const double time, const double dt,
                         const std::string& component_name);

  void postthermalProcess(hier::Patch<NDIM>& patch, const double time, const double dt,
                          const std::string& component_name);

  /// 计算各种方法得到的误差
  void calculateErrorOnPatch(hier::Patch<NDIM>& patch, const double time, const double dt,
                             const std::string& component_name);

  /**
   * @brief 注册模型变量, 完成用户变量的注册.
   */
  void registerModelVariable();

  /*!
   * @brief 支撑指定名称的有限元数值构件, 在单个网格片上完成后数值计算.
   *
   * @param patch          输入参数, 网格片类, 表示网格片.
   * @param time           输入参数, 双精度浮点型, 表示当前时刻.
   * @param dt             输入参数, 双精度浮点型, 表示时间步长.
   * @param initial_time   输入参数, BOOL型, 是否为初始时刻.
   * @param component_name 输入参数, 字符串, 表示数值构件的名称.
   */
  void computeOnPatch(hier::Patch<NDIM>& patch, const double time,
                      const double dt, const bool initial_time,
                      const std::string& component_name);

  /*!
   * @brief 支撑指定名称的归约构件, 在单个网格片上执行归约计算.
   */
  void reduceOnPatch(double* vector, int len, hier::Patch<NDIM>& patch,
                     const double time, const double dt,
                     const std::string& component_name);

  /*!
   * @brief 支撑指定名称的有限元步长构件, 在单个网格片上计算时间步长.
   */
  double getPatchDt(hier::Patch<NDIM>& patch, const double time,
                    const bool initial_time, const int flag_last_dt,
                    const double last_dt, const std::string& component_name);

  /**
   * @brief 注册可视化数据.
   */
  void registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer);

  /**
   * @brief 从输入数据库中读取参数, 并设置到计算流程中.
   */
  void setParameter(tbox::Pointer<tbox::Database> input_db);

  /**
   * @brief 从输入数据库中读如参数。
   */
  void getFromInput(tbox::Pointer<tbox::Database> db);

  /**
   * @brief 从重启动数据库中读如参数。
   */
  void getFromRestart(tbox::Pointer<tbox::Database> db);

  /**
   * @brief 将数据写入到重启动数据库。
   */
  void putToDatabase(tbox::Pointer<tbox::Database> db);

  /**
   * @brief 获取自由度信息
   */
  tbox::Pointer<solv::DOFInfo<NDIM> > getDOFInfo() { return d_dof_info; }

  /*************************************************************************
   * 材料数组
   ************************************************************************/
  class Region {
  public:
    int d_id;
    double d_epsilonr;
    double d_mur;
    double d_sigma;
    double d_density;
    double d_K;
    double d_Cp;
  };
  Region d_region_table[6];
  tbox::Array< tbox::Array< int > > material_entity;
  int NumberOfMaterial;

  /*************************************************************************
   * 边界条件
   ************************************************************************/
  // 温度
  class FixTemperature {
  public:
    tbox::Array< int > temp_face;
    double temperature;
  };
  // 固定浓度
  class FixConcentration {
  public:
    tbox::Array< int > con_face;
    double concentration;
  };
  // 固定电压
  class Voltage {
  public:
    tbox::Array< int > vol_face;
    double voltage;
  };
  // 时变电压
  class timeVoltage {
  public:
    tbox::Array< int > vol_face;
  };

  /*************************************************************************
   * 输入参数
   ************************************************************************/
  class BoundaryCondition {
  public:
    // temp boundary
    double T_initial;
    int num_temperature;
    tbox::Array< FixTemperature > FixTemperature_bc;
    // 电压边界
    int num_voltage;
    tbox::Array< Voltage > voltage_bc;
    int num_timevoltage;
    tbox::Array< timeVoltage > timevoltage_bc;
    // 浓度边界
    int num_concentration;
    tbox::Array< FixConcentration > concentration_bc;
  };

  BoundaryCondition BC;
  double pulse[101];

private:
  /*!@brief 对象名.  */
  std::string d_object_name;

  /// 自由度信息
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_info;
  /// 漂移扩散自由度信息
  tbox::Pointer<solv::DOFInfo<NDIM> > d_dof_info_DD;

  /// 有限元计算的形函数类型, 单元类型, 积分器类型.
  std::string d_shape_func_type;
  std::string d_element_type;
  std::string d_integrator_type;
  tbox::Array<std::string> d_constraint_types;
  tbox::Array<int> d_constraint_marks;

  // 时间步长
  double Dt;

  /// 对应于变量的id, 通过id应户可以获取变量数据.
  // 电数据片
  int d_E_solution_id;
  int d_potential_plot_id;
  int d_E_error_id;
  int d_E_matrix_id;
  int d_E_rhs_id;
  int d_Ex_id;
  // 信息生成
  int d_Edge_order_id;
  int d_Edge_flag_id;
  int d_Cell_volume_id;
  int d_Cell_jacobian_id;
  // 漂移扩散方程数据片
  int d_DD_solution_id;
  int d_nD_plot_id;
  int d_DD_nD_temp_id;
  int d_DD_diff_id;
  int d_DD_mobility_id;
  int d_DD_matrix_id;
  int d_DD_rhs_id;
  int d_DD_J_id;
  int d_DD_error_id;
  // 热数据片
  int d_thermal_solution_id;
  int d_temperature_plot_id;
  int d_thermal_old_id;
  int d_thermal_matrix_id;
  int d_thermal_rhs_id;
  int d_thermal_error_id;
  // 材料数据片
  int d_epsilonr_id;
  int d_Sigma_id;
  int d_Density_id;
  int d_Cp_id;
  int d_K_id;
  // 其他
  int d_Cell_flag_id;
  // 误差
  int d_Recovery_basis_id;
  int d_J_error_id;
};
#endif
