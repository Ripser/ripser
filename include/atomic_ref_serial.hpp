#pragma once

#pragma message "Using serial atomic_ref"

// Atuhor: Dmitriy Morozov
// Version 2020-03-15

#include <atomic>
#include <type_traits>
#include <cassert>

namespace mrzv
{

template<class T, class Enable = void>
class atomic_ref;

template<class T>
class atomic_ref_base
{
    public:
        using value_type = T;
        using difference_type = value_type;

    public:

        explicit    atomic_ref_base(T& obj): obj_(&obj)     {}
                    atomic_ref_base(const atomic_ref_base& ref) noexcept =default;

        T           operator=(T desired) const noexcept     { store(desired); return desired; }
        atomic_ref_base&
                    operator=(const atomic_ref_base&)       =delete;

        bool        is_lock_free() const noexcept           { return true; }

        void        store(T desired, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { *obj_ = desired; }
        T           load(std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { return *obj_; }

                    operator T() const noexcept             { return load(); }

        T           exchange(T desired, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { std::swap(*obj_, desired); return desired; }

        // Technically, not quite CAS, but simplified for the only meaningful scenario I can think of
        bool        compare_exchange_weak(T& expected, T desired,
                                          std::memory_order,
                                          std::memory_order) const noexcept
                                                            { assert(*obj_ == expected); *obj_ = desired; return true; }

        bool        compare_exchange_weak(T& expected, T desired,
                                           std::memory_order order = std::memory_order_seq_cst ) const noexcept
                                                            { assert(*obj_ == expected); *obj_ = desired; return true; }

        bool        compare_exchange_strong(T& expected, T desired,
                                            std::memory_order,
                                            std::memory_order) const noexcept
                                                            { assert(*obj_ == expected); *obj_ = desired; return true; }


        bool        compare_exchange_strong(T& expected, T desired,
                                            std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { assert(*obj_ == expected); *obj_ = desired; return true; }

        // would be great to have wait and notify, but unclear how to implement them efficiently with __atomic
        //void wait(T old, std::memory_order order = std::memory_order::seq_cst) const noexcept;
        //void wait(T old, std::memory_order order = std::memory_order::seq_cst) const volatile noexcept;
        //void notify_one() const noexcept;
        //void notify_one() const volatile noexcept;
        //void notify_all() const noexcept;
        //void notify_all() const volatile noexcept;

    protected:
        T* obj_;
};

template<class T>
class atomic_ref<T, typename std::enable_if<std::is_integral<T>::value>::type>: public atomic_ref_base<T>
{
    public:
        using Parent            = atomic_ref_base<T>;
        using value_type        = typename Parent::value_type;
        using difference_type   = typename Parent::difference_type;

    public:
        using Parent::Parent;

        value_type operator++() const noexcept          { return fetch_add(1) + 1; }
        value_type operator++(int) const noexcept       { return fetch_add(1); }
        value_type operator--() const noexcept          { return fetch_sub(1) - 1; }
        value_type operator--(int) const noexcept       { return fetch_sub(1); }

        T operator+=( T arg ) const noexcept            { return fetch_add(arg) + arg; }
        T operator-=( T arg ) const noexcept            { return fetch_sub(arg) - arg; }
        T operator&=( T arg ) const noexcept            { return fetch_and(arg) & arg; }
        T operator|=( T arg ) const noexcept            { return fetch_or(arg) | arg; }
        T operator^=( T arg ) const noexcept            { return fetch_xor(arg) ^ arg; }

        T fetch_add(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { T result = *obj_; *obj_ += arg; return result; }

        T fetch_sub(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { T result = *obj_; *obj_ -= arg; return result; }

        T fetch_and(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { T result = *obj_; *obj_ &= arg; return result; }

        T fetch_or(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { T result = *obj_; *obj_ |= arg; return result; }

        T fetch_xor(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { T result = *obj_; *obj_ ^= arg; return result; }

    protected:
        using Parent::obj_;
};

template<class T>
class atomic_ref<T*>: public atomic_ref_base<T*>
{
    public:
        using Parent            = atomic_ref_base<T*>;
        using value_type        = typename Parent::value_type;
        using difference_type   = typename Parent::difference_type;

    public:
        using Parent::Parent;

        value_type operator++() const noexcept              { return fetch_add(1) + 1; }
        value_type operator++(int) const noexcept           { return fetch_add(1); }
        value_type operator--() const noexcept              { return fetch_sub(1) - 1; }
        value_type operator--(int) const noexcept           { return fetch_sub(1); }

        T* operator+=(std::ptrdiff_t arg) const noexcept    { return fetch_add(arg) + arg; }
        T* operator-=(std::ptrdiff_t arg) const noexcept    { return fetch_sub(arg) - arg; }

        T* fetch_add(std::ptrdiff_t arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { T* result = *obj_; *obj_ += arg; return result; }

        T* fetch_sub(std::ptrdiff_t arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { T* result = *obj_; *obj_ -= arg; return result; }

    protected:
        using Parent::obj_;
};

}
