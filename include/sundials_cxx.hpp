#pragma once
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
            Vector(const Vector& v) : n_vec(N_VMake_Serial(v.size(), v.get_data_ptr())) // copy-constructor
            {}
            Vector(long int n, realtype * const data){ // "copy-constructor"
                this->n_vec = N_VMake_Serial(n, const_cast<realtype*>(data));
            }
            Vector(std::vector<realtype> v){ // "copy-constructor"
                this->n_vec = N_VMake_Serial(v.size(), &v[0]);
            }
            Vector& operator=(const Vector& v){ // Copy assignment
                if (v.size() != this->size())
                    throw std::runtime_error("Incompatible sizes");
                for (long int idx=0; idx < v.size(); ++idx)
                    (*this)[idx] = v[idx];
                return *this;
            }
            ~Vector(){
                N_VDestroy_Serial(this->n_vec);
            }
            long int size() const {
                return NV_LENGTH_S(this->n_vec);
            }
            realtype *get_data_ptr() const {
                return NV_DATA_S(this->n_vec);
            }
            realtype& operator[](long int idx) const{
                return *(NV_DATA_S(this->n_vec)+idx);
            }
            void dump(realtype * out) const {
                std::memcpy(out, NV_DATA_S(this->n_vec),
                            NV_LENGTH_S(this->n_vec)*sizeof(realtype));
            }
            void zero_out(){
                for (long int idx=0; idx < this->size(); ++idx)
                    (*this)[idx] = 0;
            }
        };
    }

}
