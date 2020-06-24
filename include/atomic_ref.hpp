#pragma once

// Atuhor: Dmitriy Morozov
// Version 2020-03-15

#if defined(USE_SERIAL_ATOMIC_REF)
#include "atomic_ref_serial.hpp"
#else

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

        bool        is_lock_free() const noexcept           { return __atomic_is_lock_free(sizeof(T), obj_); }

        void        store(T desired, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { __atomic_store_n(obj_, desired, atomic_memory_order(order)); }
        T           load(std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { return __atomic_load_n(obj_, atomic_memory_order(order)); }

                    operator T() const noexcept             { return load(); }

        T           exchange(T desired, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { return __atomic_exchange_n(obj_, desired, atomic_memory_order(order)); }

        bool        compare_exchange_weak(T& expected, T desired,
                                          std::memory_order success,
                                          std::memory_order failure) const noexcept
                                                            { return __atomic_compare_exchange_n(obj_, &expected, desired, true, atomic_memory_order(success), atomic_memory_order(failure)); }

        bool        compare_exchange_weak(T& expected, T desired,
                                           std::memory_order order = std::memory_order_seq_cst ) const noexcept
                                                            { return compare_exchange_weak(expected, desired, order, order); }      // TODO: not quite this simple

        bool        compare_exchange_strong(T& expected, T desired,
                                            std::memory_order success,
                                            std::memory_order failure) const noexcept
                                                            { return __atomic_compare_exchange_n(obj_, &expected, desired, false, atomic_memory_order(success), atomic_memory_order(failure)); }


        bool        compare_exchange_strong(T& expected, T desired,
                                            std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { return compare_exchange_strong(expected, desired, order, order); }    // TODO: not quite this simple

        // would be great to have wait and notify, but unclear how to implement them efficiently with __atomic
        //void wait(T old, std::memory_order order = std::memory_order::seq_cst) const noexcept;
        //void wait(T old, std::memory_order order = std::memory_order::seq_cst) const volatile noexcept;
        //void notify_one() const noexcept;
        //void notify_one() const volatile noexcept;
        //void notify_all() const noexcept;
        //void notify_all() const volatile noexcept;

    protected:
        int         atomic_memory_order(std::memory_order order) const
        {
            if (order == std::memory_order_relaxed)
                return __ATOMIC_RELAXED;
            else if (order == std::memory_order_consume)
                return __ATOMIC_CONSUME;
            else if (order == std::memory_order_acquire)
                return __ATOMIC_ACQUIRE;
            else if (order == std::memory_order_release)
                return __ATOMIC_RELEASE;
            else if (order == std::memory_order_acq_rel)
                return __ATOMIC_ACQ_REL;
            else if (order == std::memory_order_seq_cst)
                return __ATOMIC_SEQ_CST;
            assert(0);
            return __ATOMIC_RELAXED;
        }

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
                                                        { return __atomic_fetch_add(obj_, arg, atomic_memory_order(order)); }

        T fetch_sub(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { return __atomic_fetch_sub(obj_, arg, atomic_memory_order(order)); }

        T fetch_and(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { return __atomic_fetch_and(obj_, arg, atomic_memory_order(order)); }

        T fetch_or(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { return __atomic_fetch_or(obj_, arg, atomic_memory_order(order)); }

        T fetch_xor(T arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                        { return __atomic_fetch_xor(obj_, arg, atomic_memory_order(order)); }

    protected:
        using Parent::obj_;
        using Parent::atomic_memory_order;
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
                                                            { return __atomic_fetch_add(obj_, arg * sizeof(T), atomic_memory_order(order)); }

        T* fetch_sub(std::ptrdiff_t arg, std::memory_order order = std::memory_order_seq_cst) const noexcept
                                                            { return __atomic_fetch_sub(obj_, arg * sizeof(T), atomic_memory_order(order)); }

    protected:
        using Parent::obj_;
        using Parent::atomic_memory_order;
};

}

#endif      // USE_SERIAL_ATOMIC_REF
