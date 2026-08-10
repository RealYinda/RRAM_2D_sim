//
// 文件名:     IntegratorManager.C
// 软件包:     JAUMIN
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     积分器管理类实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "Utilities.h"
#include "BaseIntegrator.h"
#include "IntegratorManager.h"
#ifndef NULL
#define NULL (0)
#endif

template <int DIM>
tbox::Pointer<IntegratorManager<DIM> >
    IntegratorManager<DIM>::s_integrator_manager_instance = 0;

/**************************************************************************
 * 构造函数                                                            *
**************************************************************************/
template <int DIM>
IntegratorManager<DIM>::IntegratorManager() {}

/**************************************************************************
 * 析构函数                                                            *
**************************************************************************/
template <int DIM>
IntegratorManager<DIM>::~IntegratorManager() {}

template <int DIM>
tbox::Pointer<IntegratorManager<DIM> > IntegratorManager<DIM>::getManager() {
  if (s_integrator_manager_instance.isNull()) {
    s_integrator_manager_instance = new IntegratorManager<DIM>();
  }
  return (s_integrator_manager_instance);
}

/************************************************************************
 * 添加积分器.
*************************************************************************/
template <int DIM>
void IntegratorManager<DIM>::addIntegrator(
    tbox::Pointer<BaseIntegrator<DIM> > integrator) {
  const string& name = integrator->getName();
  bool exist = false;
  unsigned int size = d_vec_integrator.size();
  for (unsigned int i = 0; i < size; ++i) {
    exist = d_vec_integrator[i]->isMatching(name);
    if (exist) break;
  }
  if (!exist)
    d_vec_integrator.push_back(integrator);
  else
    TBOX_ERROR("integrator \"" << name << "\" is existed. " << endl);
}

/************************************************************************
* 获取积分器.
*************************************************************************/

template <int DIM>
tbox::Pointer<BaseIntegrator<DIM> > IntegratorManager<DIM>::getIntegrator(
    const string& integrator_name) const {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!integrator_name.empty());
#endif

  tbox::Pointer<BaseIntegrator<DIM> > integrator = NULL;
  bool found = false;
  for (unsigned int i = 0; i < d_vec_integrator.size(); ++i) {
    tbox::Pointer<BaseIntegrator<DIM> > tmp_integrator = d_vec_integrator[i];
    if (tmp_integrator->isMatching(integrator_name)) {
      integrator = tmp_integrator;
      found = true;
    }
  }

  if (!found) {
    TBOX_ERROR("integrator \"" << integrator_name << "\" is not found. "
                               << endl);
  }

  if (integrator.isNull()) {
    TBOX_ERROR("integrator \"" << integrator_name << "\" is NULL. " << endl);
  }

  return (integrator);
}

/************************************************************************
 * 输出类信息到指定流.
*************************************************************************/
template <int DIM>
void IntegratorManager<DIM>::printClassData(ostream& os) const {
  os << "printing IntegratorManager<DIM> data..." << endl;
  os << "IntegratorManager<DIM>: this = " << (IntegratorManager<DIM>*)this
     << endl;

  os << "fem  integrator list: " << endl;
}
template class IntegratorManager<NDIM>;
