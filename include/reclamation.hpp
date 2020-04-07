#pragma once

#include <vector>
#include <atomic_ref.hpp>

namespace mrzv
{

template<class T>
struct MemoryManager
{
        std::vector<T*> retired_, to_delete_;
        atomic_ref<int> counter_;
        unsigned        n_threads_;
        bool            even_epoch_ = false;

        MemoryManager(int& epoch_counter, unsigned n_threads):
            counter_(epoch_counter), n_threads_(n_threads)      {}

        ~MemoryManager()
        {
            for(T* p : to_delete_)
                delete p;
            for(T* p : retired_)
                delete p;
        }
        bool is_even_epoch(int counter) const   { return (counter / n_threads_) % 2 == 0; }
        void retire(T* ptr)                     { if(ptr) retired_.push_back(ptr); }
        void quiescent()
        {
            if (even_epoch_ != is_even_epoch(counter_))
            {
                ++counter_;
                even_epoch_ = !even_epoch_;
                for(T* p : to_delete_)
                    delete p;
                retired_.swap(to_delete_);
                retired_.clear();
            }
        }
};

}
