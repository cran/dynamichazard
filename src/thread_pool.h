/*
 Boost Software License - Version 1.0 - August 17th, 2003

 Permission is hereby granted, free of charge, to any person or organization
 obtaining a copy of the software and accompanying documentation covered by
 this license (the "Software") to use, reproduce, display, distribute,
 execute, and transmit the Software, and to prepare derivative works of the
 Software, and to permit third-parties to whom the Software is furnished to
 do so, all subject to the following:

 The copyright notices in the Software and this entire statement, including
 the above license grant, this restriction and the following disclaimer,
 must be included in all copies of the Software, in whole or in part, and
 all derivative works of the Software, unless such copies or derivative
 works are solely in the form of machine-executable object code generated by
 a source language processor.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.
 */

// Code from:
// Williams, Anthony. C++ concurrency in action. London, 2012

#include <deque>
#include <future>
#include <memory>
#include <functional>
#include <iostream>
#include <iostream>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <type_traits>

// Listing 6.6
template<typename T>
class thread_safe_queue
{
private:
   struct node
   {
      std::shared_ptr<T> data;
      std::unique_ptr<node> next;
   };

   std::mutex head_mutex;
   std::unique_ptr<node> head;
   std::mutex tail_mutex;
   node* tail;

   node* get_tail()
   {
      std::lock_guard<std::mutex> tail_lock(tail_mutex);
      return tail;
   }

   std::unique_ptr<node> pop_head()
   {
      std::lock_guard<std::mutex> head_lock(head_mutex);
      if (head.get() == get_tail())
      {
         return nullptr;
      }
      std::unique_ptr<node> old_head = std::move(head);
      head = std::move(old_head->next);
      return old_head;
   }


public:
   thread_safe_queue() :
      head(new node), tail(head.get())
   {}

   thread_safe_queue(const thread_safe_queue& other) = delete;
   thread_safe_queue& operator=(const thread_safe_queue& other) = delete;

   std::shared_ptr<T> try_pop()
   {
      std::unique_ptr<node> old_head = pop_head();
      return old_head ? old_head->data : std::shared_ptr<T>();
   }

   // Added
   bool try_pop(T& value)
   {
      std::unique_ptr<node> old_head = pop_head();
      if(old_head){
        value = std::move(*old_head->data);
        return true;

      } else{
        return false;

      }
   }

   void push(T new_value)
   {
      std::shared_ptr<T> new_data(
         std::make_shared<T>(std::move(new_value)));
      std::unique_ptr<node> p(new node);
      node* const new_tail = p.get();
      std::lock_guard<std::mutex> tail_lock(tail_mutex);
      tail->data = new_data;
      tail->next = std::move(p);
      tail = new_tail;
   }
};

// Just before listing 8.4
class join_threads
{
  std::vector<std::thread>& threads;
public:
  explicit join_threads(std::vector<std::thread>& threads_);
  ~join_threads();
};

// Listing 9.2
class function_wrapper
{
  struct impl_base {
    virtual void call()=0;
    virtual ~impl_base() {}
  };
  std::unique_ptr<impl_base> impl;
  template<typename F>
  struct impl_type: impl_base
  {
    F f;
    impl_type(F&& f_): f(std::move(f_)) {}
    void call() { f(); }
  };
public:
  template<typename F>
  function_wrapper(F&& f):
    impl(new impl_type<F>(std::move(f)))
  {}

  void operator()() { impl->call(); }
  function_wrapper() = default;
  function_wrapper(function_wrapper&& other):
    impl(std::move(other.impl))
  {}

  function_wrapper& operator=(function_wrapper&& other)
  {
    impl=std::move(other.impl);
    return *this;
  }

  function_wrapper(const function_wrapper&)=delete;
  function_wrapper(function_wrapper&)=delete;
  function_wrapper& operator=(const function_wrapper&)=delete;
};


// Listing 9.7:
class work_stealing_queue
{
private:
  typedef function_wrapper data_type;
  std::deque<data_type> the_queue;
  mutable std::mutex the_mutex;

public:
  work_stealing_queue()
  {}

  work_stealing_queue(const work_stealing_queue& other)=delete;
  work_stealing_queue& operator=(
    const work_stealing_queue& other)=delete;

  void push(data_type data);
  bool empty() const;
  bool try_pop(data_type& res);
  bool try_steal(data_type& res);
};

// Listing 9.8:
class thread_pool
{
  typedef function_wrapper task_type;

  std::atomic_bool done;
  std::atomic_bool start; // added
  thread_safe_queue<task_type> pool_work_queue;
  std::vector<std::unique_ptr<work_stealing_queue> > queues;
  std::vector<std::thread> threads;
  join_threads joiner;

  // A way to get a static __thread member in a header class. See http://stackoverflow.com/a/11711082
  static work_stealing_queue*& local_work_queue(){
    static __thread work_stealing_queue* local_work_queue_;
    return local_work_queue_;
  }
  static unsigned& my_index(){
    static __thread unsigned my_index_;
    return my_index_;
  }

  void worker_thread(unsigned my_index_)
  {
    while(!start){
      std::this_thread::yield();
    }

    my_index()=my_index_;
    local_work_queue()=queues[my_index()].get();

    while(!done)
    {
      run_pending_task();
    }

    local_work_queue() = NULL;
  }

  bool pop_task_from_local_queue(task_type& task)
  {
    return local_work_queue() && local_work_queue()->try_pop(task);
  }

  bool pop_task_from_pool_queue(task_type& task)
  {
    return pool_work_queue.try_pop(task);
  }

  bool pop_task_from_other_thread_queue(task_type& task)
  {
    for(unsigned i=0;i<queues.size();++i)
    {
      unsigned const index=(my_index()+i+1)%queues.size();
      if(queues[index]->try_steal(task))
      {
        return true;
      }
    }

    return false;
  }

public:
  thread_pool(unsigned int n_jobs, unsigned int max_threads = 1):
  done(false), start(false), joiner(threads)
  {
    unsigned const thread_count = std::min(n_jobs, max_threads);

    try
    {
      for(unsigned i=0;i<thread_count;++i)
      {
        queues.emplace_back(new work_stealing_queue);
        threads.push_back(
          std::thread(&thread_pool::worker_thread,this,i));
      }

      start = true;
    }
    catch(...)
    {
      start = true;
      done=true;
      throw;
    }
  }

  ~thread_pool()
  {
    done=true;
  }

  // template<typename ResultType>
  // using task_handle = std::future<ResultType>;
  // future is re-named in the book as far as I gather. using task_handle=std::unique_future<ResultType>;
  template<typename FunctionType>
  std::future<typename std::result_of<FunctionType()>::type> submit(FunctionType f)
  {
    typedef typename std::result_of<FunctionType()>::type result_type;

    std::packaged_task<result_type()> task(f);
    std::future<result_type> res(task.get_future());
    if(local_work_queue())
    {
      local_work_queue()->push(std::move(task));
    }
    else
    {
      pool_work_queue.push(std::move(task));
    }
    return res;
  }

  void run_pending_task()
  {
    task_type task;
    if(pop_task_from_local_queue(task) ||
       pop_task_from_pool_queue(task) ||
       pop_task_from_other_thread_queue(task))
    {
      task();
    }
    else
    {
      std::this_thread::yield();
    }
  }
};
