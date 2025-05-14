#include <honeycomb2/thread_pool.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{
ThreadPool::ThreadPool(size_t num_threads) : _remaining(0)
{
   _am_I_alive = true;
   for (size_t i = 0; i < num_threads; i++) {
      _threads.emplace_back([this] {
         while (true) {

            std::unique_lock<std::mutex> lock(_q_mutex);

            _notifier.wait(lock, [this](void) -> bool {
               return !_tasks.empty() || _to_stop;
            });

            if (_to_stop) {
               if (!_tasks.empty()) {
                  logger(Logger::WARNING, "The threadpool is not empty, but the user requested to stop. "
                                          "Ignoring tasks still in the pool.");
               }
               return;
            }

            std::function<void()> task = std::move(_tasks.back());
            _tasks.pop_back();

            lock.unlock();

            task();
            {
               std::unique_lock<std::mutex> lock(_wait_mutex);
               _remaining--;
            }
            _done.notify_one();
         }
      });
   }
}

ThreadPool::~ThreadPool()
{
   if (_am_I_alive) ShutDown();
}

void ThreadPool::ShutDown()
{
   if (_am_I_alive) {
      std::unique_lock<std::mutex> lock(_q_mutex);
      _to_stop = true;
      lock.unlock();

      _notifier.notify_all();

      for (std::thread &thread : _threads) {
         if (thread.joinable()) thread.join();
      }
   }
}

void ThreadPool::AddTask(std::function<void()> task)
{

   std::unique_lock<std::mutex> lock(_q_mutex);
   _tasks.emplace_back(std::move(task));
   std::unique_lock<std::mutex> lock_done(_wait_mutex);
   _remaining++;
   lock.unlock();
   lock_done.unlock();

   _notifier.notify_one();
}

void ThreadPool::WaitOnJobs()
{
   std::unique_lock<std::mutex> lock(_wait_mutex);
   _done.wait(lock, [this](void) {
      return _remaining == 0;
   });
}

}; // namespace Honeycomb