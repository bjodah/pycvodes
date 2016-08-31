#pragma once
#include <cstring> // std::memcpy
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */


namespace sundials_cxx {

    namespace nvector_serial {
        struct VectorBase_{
            N_Vector n_vec {nullptr};
            VectorBase_(N_Vector v) : n_vec(v) {}
            virtual ~VectorBase_(){
                N_VDestroy_Serial(this->n_vec);
            }
            long int size() const {
                return NV_LENGTH_S(this->n_vec);
            }
            realtype *get_data_ptr() const {
                return NV_DATA_S(this->n_vec);
            }
            void set_data_ptr(realtype * data){
                NV_DATA_S(this->n_vec) = data;
            }
            realtype& operator[](long int idx) const{
                return *(NV_DATA_S(this->n_vec)+idx);
            }
            void dump(realtype * const out) const {
                std::memcpy(out, NV_DATA_S(this->n_vec),
                            NV_LENGTH_S(this->n_vec)*sizeof(realtype));
            }
            void load(const realtype * const in) const {
                std::memcpy(NV_DATA_S(this->n_vec), in,
                            NV_LENGTH_S(this->n_vec)*sizeof(realtype));
            }
            void set_all(realtype value) {
                for (long int idx=0; idx < this->size(); ++idx)
                    (*this)[idx] = value;
            }
            void zero_out(){ set_all(0); }  // deprecated
        };

        struct Vector : public VectorBase_ {
            // Vector owns the memory containing data
            Vector(long int n) :
                VectorBase_(N_VNew_Serial(n)) {}
            Vector(const Vector& v)
                : VectorBase_(N_VNew_Serial(v.size())) { // copy-constructor
                v.dump(get_data_ptr());
            }
            Vector(long int n, const realtype * const data)
                : VectorBase_(N_VNew_Serial(n)) { // "copy-constructor"
                load(data);
            }
            Vector(const std::vector<realtype>& v)
                : VectorBase_(N_VNew_Serial(v.size())) { // "copy-constructor"
                load(&v[0]);
            }
            Vector& operator=(const Vector& v){ // Copy assignment
                if (v.size() != this->size())
                    throw std::runtime_error("Incompatible sizes");
                load(v.get_data_ptr());
                return *this;
            }
        };

        struct VectorView : public VectorBase_ {
            // VectorView DOES NOT own the memory containing data

            VectorView(long int n, realtype * const data)
                : VectorBase_(N_VMake_Serial(n, const_cast<realtype*>(data))) {}

            VectorView(const VectorView& v)
                : VectorBase_(N_VMake_Serial(v.size(), v.get_data_ptr())) {} // copy-constructor

            VectorView(const Vector& v)
                : VectorBase_(N_VMake_Serial(v.size(), v.get_data_ptr())) {}

            VectorView(const std::vector<realtype>& v)
                : VectorBase_(N_VMake_Serial(v.size(), const_cast<realtype*>(&v[0]))) {}

            VectorView& operator=(const Vector& v){ // Copy assignment
                if (v.size() != this->size())
                    throw std::runtime_error("Incompatible sizes");
                set_data_ptr(v.get_data_ptr());
                return *this;
            }
        };
    }
}
