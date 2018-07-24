#include "anyode/anyode.hpp"
#include "cvodes_anyode.hpp"
#include "cvodes_cxx.hpp"

BEGIN_NAMESPACE(AnyODE)
struct CvodesPyOdeSys : public AnyODE::PyOdeSys<realtype, indextype> {
    using AnyODE::PyOdeSys<realtype, indextype>::PyOdeSys;
};
END_NAMESPACE(AnyODE)