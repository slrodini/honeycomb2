#include <honeycomb2/thread_pool.hpp>
#include <honeycomb2/utilities.hpp>

namespace Honeycomb
{
ThreadPool::ThreadPool(size_t num_threads)
{
   for (size_t i = 0; i < num_threads; i++) {
      _threads.emplace_back([this] {
         while (true) {

            std::unique_lock<std::mutex> lock(_q_mutex);

            _notifier.wait(lock, [this](void) -> bool {
               return !_tasks.empty() || _to_stop;
            });

            if (_to_stop) {
               if (!_tasks.empty()) {
                  logger(
                      Logger::WARNING,
                      "The threadpool is not empty, but the user requested to stop. Ignoring tasks still in the pool.");
               }
               return;
            }

            std::function<void()> task = std::move(_tasks.back());
            _tasks.pop_back();

            lock.unlock();

            task();
         }
      });
   }
}

ThreadPool::~ThreadPool()
{

   std::unique_lock<std::mutex> lock(_q_mutex);
   _to_stop = true;
   lock.unlock();

   _notifier.notify_all();

   for (std::thread &thread : _threads) {
      thread.join();
   }
}

void ThreadPool::add_task(std::function<void()> task)
{

   std::unique_lock<std::mutex> lock(_q_mutex);
   _tasks.emplace_back(std::move(task));
   lock.unlock();

   _notifier.notify_one();
}

}; // namespace Honeycomb