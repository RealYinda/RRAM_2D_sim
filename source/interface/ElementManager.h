//
// 文件名:     ElementManager.h
// 软件包:     2D FEM
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     有限元单元管理器
// 类别　:     %Internal File% ( Don't delete this line )
//

#ifndef included_fem_ElementManager
#define included_fem_ElementManager

#include "Array.h"
#include <vector>
using namespace std;

using namespace JAUMIN;
template <int DIM>
class BaseElement;
/**
 * @brief 有限元单元管理类 adpt::ElementManager, 为有限元模块管理单元。
 *
 * 该类维护一个数组, 存储用户注册到JAUMIN框架的单元。
 *
 * 该类提供如下函数, 支持用户从数组中获取指定名字的单元。 - getElement()
 *
 *
 * @see adpt::BaseElement
 */

template <int DIM>
class ElementManager {
public:
  /**
   * @brief 构造函数。
   */
  ElementManager();

  /**
   * @brief 虚析构函数。
   */
  virtual ~ElementManager();

  /**
   * @brief  获取单元理器.
   *
   * @return 指针, 单元管理器.
   */
  static tbox::Pointer<ElementManager<DIM> > getManager();

  /**
   * @brief 将单元注册到单元管理器。
   *
   * @param ele_name 输入参数, 指针, 指向添加的单元.
   */
  void addElement(tbox::Pointer<BaseElement<DIM> > ele);

  /**
   * @brief 获取匹配于指定名称的单元.
   * @param data_name 输入参数, 字符串, 单元名字.
   * @return 指针, 指向匹配的单元.
   */
  tbox::Pointer<BaseElement<DIM> > getElement(const string &element_name) const;

  void setElement(tbox::Array<string> &element_names,
                  tbox::Array<int> &element_marks);

  tbox::Array<string> getElementNames();
  tbox::Array<int> getElementMarks();

private:
  static tbox::Pointer<ElementManager<DIM> > s_element_manager_instance;

  tbox::Array<string> d_element_names;
  tbox::Array<int> d_element_marks;

  /**
   * 存储单元的数组。
   */
  std::vector<tbox::Pointer<BaseElement<DIM> > > d_vec_element;
};

#endif
