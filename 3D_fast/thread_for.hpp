#include <algorithm>
#include <thread>
#include <functional>
#include <vector>
using namespace std;

/// @param[in] nb_elements : size of your for loop
/// @param[in] functor(start, end) :
/// your function processing a sub chunk of the for loop.
/// "start" is the first index to process (included) until the index "end"
/// (excluded)
/// @code
///     for(int i = start; i < end; ++i)
///         computation(i);
/// @endcode
/// @param use_threads : enable / disable threads.
///
///
static
void parallel_for(unsigned nb_elements,
				  function<void (int start, int end)> functor,
				  bool use_env = false, bool use_threads = true)
{
	// -------
	unsigned nb_threads_hint = use_env? atoi(getenv("SLURM_CPUS_PER_TASK")): thread::hardware_concurrency();
	
	unsigned nb_threads = nb_threads_hint == 0 ? 8 : (nb_threads_hint);
	
	unsigned batch_size = nb_elements / nb_threads;
	unsigned batch_remainder = nb_elements % nb_threads;
	
	vector<thread> my_threads(nb_threads);
	
	if(use_threads)
	{
		// Multithread execution
		for(unsigned i = 0; i < nb_threads; ++i)
		{
			int start = i * batch_size;
			my_threads[i] = thread(functor, start, start+batch_size);
		}
	}
	else
	{
		// Single thread execution (for easy debugging)
		for(unsigned i = 0; i < nb_threads; ++i){
			int start = i * batch_size;
			functor(start, start + batch_size);
		}
	}
	// Deform the elements left
	int start = nb_threads * batch_size;
	functor(start, start + batch_remainder);
	// Wait for the other thread to finish their task
	if(use_threads)
		for_each(my_threads.begin(), my_threads.end(), mem_fn(&thread::join));
}
