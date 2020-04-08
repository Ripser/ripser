#pragma once

// Atuhor: Dmitriy Morozov
// Version 2020-03-12

#include <vector>
#include <utility>

#include <atomic_ref.hpp>

namespace mrzv
{

template <class Key, class T, class H, class E>
class trivial_concurrent_hash_map
{
    public:
        using key_type      = Key;
        using mapped_type   = T;
        using value_type    = std::pair<Key,T>;
        using Storage       = std::vector<value_type>;
        using iterator      = typename Storage::iterator;

        static constexpr size_t storage_multiplier = 8;

    public:
        void                reserve(size_t hint)
        {
            size_t sz = 1;
            while (sz < hint*storage_multiplier)
                sz *= 2;
            storage_.clear();
            storage_.resize(sz, dummy());
        }

        inline iterator     find(Key k);
        inline std::pair<iterator,bool>
                            insert(value_type x);

        // CAS update value
        inline bool         update(iterator it, T& expected, T desired)
                                                    { return atomic_ref<T>(it->second).compare_exchange_weak(expected, desired); }

        iterator            end()                   { return storage_.end(); }

        Key                 key(iterator it)        { return atomic_ref<Key>(it->first).load(); }
        T                   value(iterator it)      { return atomic_ref<T>(it->second).load(); }

        template<class F>
        void foreach(const F& f) const              { for(auto& x : storage_) if (x.first != dummy_key()) f(x); }

    public:
        static Key          dummy_key()             { return static_cast<Key>(-1); }
        static T            dummy_value()           { return static_cast<T>(-1); }
        static value_type   dummy()                 { return value_type { dummy_key(), dummy_value() }; }

    private:
        static size_t       hash(Key k)             { return H()(k); }
        static bool         equal(Key k1, Key k2)   { return E()(k1,k2); }

    private:
        Storage storage_;
};

}


template <class Key, class T, class H, class E>
typename mrzv::trivial_concurrent_hash_map<Key,T,H,E>::iterator
mrzv::trivial_concurrent_hash_map<Key,T,H,E>::
find(Key k)
{
    auto size = storage_.size();
    auto idx = hash(k) % size;

    while(true)
    {
        iterator it = storage_.begin() + idx;

        atomic_ref<Key> ak(it->first);
        Key kk = ak;

        if (equal(kk, k)) {
            while(value(it) == dummy_value());  // make sure a value has been written (possibly should add exponential backoff)
            return it;
        }
        else if (equal(kk, dummy_key()))
            return end();

        idx = (idx + 1) % size;
    }
}

template <class Key, class T, class H, class E>
std::pair<typename mrzv::trivial_concurrent_hash_map<Key,T,H,E>::iterator, bool>
mrzv::trivial_concurrent_hash_map<Key,T,H,E>::
insert(value_type x)
{
    Key k = x.first;

    auto size = storage_.size();
    auto idx = hash(k) % size;

    while(true)
    {
        iterator it = storage_.begin() + idx;

        atomic_ref<Key> ak(it->first);
        Key kk = ak;

        if (equal(kk, k))
            return { it, false };
        else if (equal(kk, dummy_key()))
        {
            Key dk = dummy_key();
            if (ak.compare_exchange_weak(dk, k))
            {
                atomic_ref<T> av(it->second);
                av.store(x.second);
                return { it, true };
            } else if (dk == k)     // somebody overwrote our own key, we lose
                return { it, false };
            // else something is written, but not our key, continue scanning
        }

        idx = (idx + 1) % size;
    }
}
