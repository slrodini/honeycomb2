#ifndef THREAD_POOL_HPP
#define THREAD_POOL_HPP

#include <honeycomb2/default.hpp>
#include <mutex>
#include <condition_variable>
#include <thread>

namespace Honeycomb
{

struct ThreadPool {

   ThreadPool(size_t num_threads = 1);
   ~ThreadPool();
   void AddTask(std::function<void()> task);
   void WaitOnJobs();

   std::mutex task_mutex;

private:
   long int _remaining;
   std::mutex _wait_mutex;
   std::condition_variable _done;

   std::vector<std::thread> _threads;
   std::vector<std::function<void()>> _tasks;
   std::mutex _q_mutex;
   std::condition_variable _notifier;
   bool _to_stop = false;
};
} // namespace Honeycomb

#endif // THREAD_POOL_HPP