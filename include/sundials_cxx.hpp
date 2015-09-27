#ifndef SUNDIALS_CXX_HPP_3CSK5Z37F5GSNHG2O23JGOSJWA
#define SUNDIALS_CXX_HPP_3CSK5Z37F5GSNHG2O23JGOSJWA

#include <cstring> // std::memcpy
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */


namespace sundials_cxx {

    namespace nvector_serial {
        class Vector {
        public:
            N_Vector n_vec {nullptr};
            Vector(long int n){
                this->n_vec = N_VNew_Serial(n);
            }
            Vector(long int n, realtype * const data){
                this->n_vec = N_VMake_Serial(n, const_cast<realtype*>(data));
            }
            Vector(std::vector<realtype> v){
                this->n_vec = N_VMake_Serial(v.size(), &v[0]);
            }
            ~Vector(){
                N_VDestroy_Serial(this->n_vec);
            }
            realtype& operator[](long int idx){
                return *(NV_DATA_S(this->n_vec)+idx);
            }
            void dump(realtype * out){
                std::memcpy(out, NV_DATA_S(this->n_vec),
                            NV_LENGTH_S(this->n_vec)*sizeof(realtype));
            }
        };
    }

} // namespace sundials_cxx
#endif /* SUNDIALS_CXX_HPP_3CSK5Z37F5GSNHG2O23JGOSJWA */
