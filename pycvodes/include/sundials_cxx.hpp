#pragma once
#include <cstring> // std::memcpy
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_config.h>
namespace sundials_cxx {
#if SUNDIALS_VERSION_MAJOR >= 3
    const int version_major = SUNDIALS_VERSION_MAJOR;
    const int version_minor = SUNDIALS_VERSION_MINOR;
    const int version_patch = SUNDIALS_VERSION_PATCH;
#else
#  if defined(SUNDIALS_PACKAGE_VERSION)   /* == 2.7.0 */
    const int version_major = 2;
    const int version_minor = 7;
    const int version_patch = 0;
#  else
#    error "Unkown sundials version"
#  endif
#endif
    namespace nvector_serial {
        struct VectorBase_{
            N_Vector n_vec {nullptr};
            VectorBase_(N_Vector v) : n_vec(v) {}
            virtual ~VectorBase_() //= default;
            {
                N_VDestroy_Serial(this->n_vec);
            }

            [[nodiscard]] long int size() const {
                return NV_LENGTH_S(this->n_vec);
            }
            realtype *get_data_ptr() const {
                return NV_DATA_S(this->n_vec);
            }
            void set_data_ptr(realtype * data){
                NV_DATA_S(this->n_vec) = data;
            }
            [[nodiscard]] realtype& operator[](long int idx) const{
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
            Vector(long int n
#if SUNDIALS_VERSION_MAJOR >= 6
                   , SUNContext ctx
#endif
                ) :
                VectorBase_(N_VNew_Serial(n
#if SUNDIALS_VERSION_MAJOR >= 6
                   , ctx
#endif
                                )) {}
            Vector(const Vector& v
#if SUNDIALS_VERSION_MAJOR >= 6
                   , SUNContext ctx
#endif
)
                : VectorBase_(N_VNew_Serial(v.size()
#if SUNDIALS_VERSION_MAJOR >= 6
                   , ctx
#endif
)) { // copy-constructor
                v.dump(get_data_ptr());
            }
            Vector(long int n, const realtype * const data
#if SUNDIALS_VERSION_MAJOR >= 6
                   , SUNContext ctx
#endif
)
                : VectorBase_(N_VNew_Serial(n
#if SUNDIALS_VERSION_MAJOR >= 6
                   , ctx
#endif
)) { // "copy-constructor"
                load(data);
            }
            Vector(const std::vector<realtype>& v
#if SUNDIALS_VERSION_MAJOR >= 6
                   , SUNContext ctx
#endif
)
                : VectorBase_(N_VNew_Serial(v.size()
#if SUNDIALS_VERSION_MAJOR >= 6
                   , ctx
#endif
)) { // "copy-constructor"
                load(&v[0]);
            }
            Vector& operator=(const Vector& v){ // Copy assignment
                if (v.size() != this->size())
                    throw std::runtime_error("Incompatible sizes");
                load(v.get_data_ptr());
                return *this;
            }
            // ~Vector() override {
            // }

        };

        struct VectorView : public VectorBase_ {
            // VectorView DOES NOT own the memory containing data

            VectorView(long int n, realtype * const data
#if SUNDIALS_VERSION_MAJOR >= 6
                   , SUNContext ctx
#endif
)
                : VectorBase_(N_VMake_Serial(n, const_cast<realtype*>(data)
#if SUNDIALS_VERSION_MAJOR >= 6
                   , ctx
#endif
)) {}

            VectorView(const VectorView& v
#if SUNDIALS_VERSION_MAJOR >= 6
                   , SUNContext ctx
#endif
)
                : VectorBase_(N_VMake_Serial(v.size(), v.get_data_ptr()
#if SUNDIALS_VERSION_MAJOR >= 6
                   , ctx
#endif
)) {} // copy-constructor

            VectorView(const Vector& v
#if SUNDIALS_VERSION_MAJOR >= 6
                   , SUNContext ctx
#endif
)
                : VectorBase_(N_VMake_Serial(v.size(), v.get_data_ptr()
#if SUNDIALS_VERSION_MAJOR >= 6
                   , ctx
#endif
)) {}

            VectorView(const std::vector<realtype>& v
#if SUNDIALS_VERSION_MAJOR >= 6
                   , SUNContext ctx
#endif
                )
                : VectorBase_(N_VMake_Serial(v.size(), const_cast<realtype*>(&v[0])
#if SUNDIALS_VERSION_MAJOR >= 6
                   , ctx
#endif
)) {}

            VectorView& operator=(const Vector& v){ // Copy assignment
                if (v.size() != this->size())
                    throw std::runtime_error("Incompatible sizes");
                set_data_ptr(v.get_data_ptr());
                return *this;
            }
        };
    }
}
