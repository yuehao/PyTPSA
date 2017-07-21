cimport cython
cimport libcpp.complex

from libcpp.vector cimport vector
from libcpp.string cimport string

from libcpp cimport bool
cdef extern from "../src/tpsalib.h":
    cdef cppclass CTPS[T]:
        CTPS() except +
        CTPS(const T&) except +
        CTPS(const CTPS[T]&) except +


        @staticmethod
        void Initialize(const int&, const int&)

        @staticmethod
        int Get_Max_Degree()

        @staticmethod
        int Get_TPS_Dim()

        unsigned long findindex(const vector[int]& ) const
        vector[int] findpower(const unsigned long & n) const

        const int get_dim() const
        const unsigned long get_degree() const
        const unsigned long get_terms() const

        const T element(const unsigned long & ) const
        const T element(const vector[int] &) const

        void assign(const T &, const int &)
        void assign(const T &)

        T evaluate(const vector[double] ) const
        CTPS[T] derivative(const int&, const int&) const

        string print_to_string() const

        CTPS[T]& add_to(const CTPS[T] &)
        CTPS[T]& sub_to(const CTPS[T] &)
        CTPS[T]& mul_to(const CTPS[T] &)
        CTPS[T]& div_to(const CTPS[T] &)

        const T cst() const

cdef extern from "../src/tpsalib.h":
    CTPS[double] inv(const CTPS[double] &)
    CTPS[complex] inv(const CTPS[complex] &)

    CTPS[double] exp(const CTPS[double] &)
    CTPS[complex] exp(const CTPS[complex] &)

    CTPS[double] log(const CTPS[double] &)
    CTPS[complex] log(const CTPS[complex] &)

    CTPS[double] sqrt(const CTPS[double] &)
    CTPS[complex] sqrt(const CTPS[complex] &)

    CTPS[double] pow(const CTPS[double] &, const double&)
    CTPS[complex] pow(const CTPS[complex] &, const double&)

    CTPS[double] sin(const CTPS[double] &)
    CTPS[complex] sin(const CTPS[complex] &)

    CTPS[double] cos(const CTPS[double] &)
    CTPS[complex] cos(const CTPS[complex] &)

    CTPS[double] tan(const CTPS[double] &)
    CTPS[complex] tan(const CTPS[complex] &)

    CTPS[double] sinh(const CTPS[double] &)
    CTPS[complex] sinh(const CTPS[complex] &)

    CTPS[double] cosh(const CTPS[double] &)
    CTPS[complex] cosh(const CTPS[complex] &)