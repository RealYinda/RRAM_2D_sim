//
// 文件名:     ShapeFunctionManager.C
// 软件包:     JAUMIN fem
// 版权　:     北京应用物理与计算数学研究所
// 版本号:     $Revision: 0 $
// 修改　:     $Date: 2011-11-01 16:22:08 +0800 (二, 2011-11-01) $
// 描述　:     形函数管理类实现.
// 类别　:     %Internal File% ( Don't delete this line )
//

#include "ShapeFunctionManager.h"
#include "Utilities.h"
#include "BaseShapeFunction.h"
#ifndef NULL
#define NULL (0)
#endif

template <int DIM>
tbox::Pointer<ShapeFunctionManager<DIM> >
    ShapeFunctionManager<DIM>::s_shape_func_manager_instance = 0;

/*
*************************************************************************
*                                                                       *
* Constructor and destructor for ShapeFunctionManager objects.      *
*                                                                       *
*************************************************************************
*/
template <int DIM>
ShapeFunctionManager<DIM>::ShapeFunctionManager() {}

template <int DIM>
ShapeFunctionManager<DIM>::~ShapeFunctionManager() {}

template <int DIM>
tbox::Pointer<ShapeFunctionManager<DIM> >
ShapeFunctionManager<DIM>::getManager() {
  if (s_shape_func_manager_instance.isNull()) {
    s_shape_func_manager_instance = new ShapeFunctionManager<DIM>();
  }
  return (s_shape_func_manager_instance);
}

/************************************************************************
*                                                                       *
* Add shape function to list.                                               *
*                                                                       *
*************************************************************************/

template <int DIM>
void ShapeFunctionManager<DIM>::addShapeFunction(
    tbox::Pointer<BaseShapeFunction<DIM> > function) {
  const string& name = function->getName();
  bool exist = false;
  unsigned int size = d_vec_shape_function.size();
  for (unsigned int i = 0; i < size; ++i) {
    exist = d_vec_shape_function[i]->isMatching(name);
    if (exist) break;
  }
  if (!exist)
    d_vec_shape_function.push_back(function);
  else
    TBOX_ERROR("shape_function \""
               << name << "\" is existed. please choose another name!!"
               << endl);
}

/************************************************************************
*                                                                       *
* Search shape function lists for data matching request.                    *
*                                                                       *
*************************************************************************/

template <int DIM>
tbox::Pointer<BaseShapeFunction<DIM> >
ShapeFunctionManager<DIM>::getShapeFunction(const string& name) const {
#ifdef DEBUG_CHECK_ASSERTIONS
  TBOX_ASSERT(!name.empty());
#endif

  tbox::Pointer<BaseShapeFunction<DIM> > function = NULL;
  bool found = false;
  for (unsigned int i = 0; i < d_vec_shape_function.size(); ++i) {
    tbox::Pointer<BaseShapeFunction<DIM> > tmp_function =
        d_vec_shape_function[i];
    if (tmp_function->isMatching(name)) {
      function = tmp_function;
      found = true;
    }
  }

  if (!found) {
    TBOX_ERROR("shape function \"" << name << "\" is not found. " << endl);
  }

  if (function.isNull()) {
    TBOX_ERROR("shape function \"" << name << "\" is NULL. " << endl);
  }

  return (function);
}

/************************************************************************
*                                                                       *
* Print ShapeFunctionManager class data.                               *
*                                                                       *
*************************************************************************/

template <int DIM>
void ShapeFunctionManager<DIM>::printClassData(ostream& os) const {
  os << "printing ShapeFunctionManager<DIM> data..." << endl;
  os << "ShapeFunctionManager<DIM>: this = " << (ShapeFunctionManager<DIM>*)this
     << endl;

  os << "shape function list: " << endl;
}

template class ShapeFunctionManager<NDIM>;
