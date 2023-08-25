#ifndef GBWTGRAPH_CONSTRUCTION_H
#define GBWTGRAPH_CONSTRUCTION_H

#include <cstdlib>
#include <functional>

#include <omp.h>

#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/minimizer.h>
#define RYMER

/*
  index.h: Minimizer index construction from GBWTGraph.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

/*
  Index the haplotypes in the graph. Insert the minimizers into the provided index.
  Function argument get_payload is used to generate the payload for each position
  stored in the index.
  The number of threads can be set through OMP.
*/
template<class KeyType>
void
index_haplotypes(const GBWTGraph& graph, MinimizerIndex<KeyType>& index,
                 const std::function<payload_type(const pos_t&)>& get_payload, unsigned int k=21)
{
  typedef typename MinimizerIndex<KeyType>::minimizer_type minimizer_type;

  int threads = omp_get_max_threads();

  // Minimizer caching. We only generate the payloads after we have removed duplicate positions.
  std::vector<std::vector<std::pair<minimizer_type, pos_t>>> cache(threads);
  constexpr size_t MINIMIZER_CACHE_SIZE = 1024;

      // Introduce a shared counter to track the number of indexed rymers.
  std::atomic<size_t> minimizer_count(0);

  auto flush_cache = [&](int thread_id)
  {

    std::vector<std::pair<minimizer_type, pos_t>>& current_cache = cache[thread_id];
    gbwt::removeDuplicates(current_cache, false);
    std::vector<payload_type> payload;
    payload.reserve(current_cache.size());
    for(size_t i = 0; i < current_cache.size(); i++) { payload.push_back(get_payload(current_cache[i].second)); }
    #pragma omp critical (minimizer_index)
    {
      for(size_t i = 0; i < current_cache.size(); i++)
      {
        minimizer_count++;
        index.insert(current_cache[i].first, current_cache[i].second, payload[i]);
      }
    }
    cache[thread_id].clear();
  };

  // Minimizer finding.
  auto find_minimizers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    std::vector<minimizer_type> minimizers = index.minimizers(seq); // Calls syncmers() when appropriate.
    auto iter = traversal.begin();
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();

//        if (std::all_of(minimizers.begin(), minimizers.end(), [](const minimizer_type& minimizer){ return minimizer.empty(); }))
//{
 //   throw std::runtime_error("All minimizers are empty!");
//}


    for(minimizer_type& minimizer : minimizers)
    {
       if(minimizer.empty()) { continue; }

       //cerr << "MINIMIZER KEY: " << minimizer.key << endl;

       std::string minimizer_sequence = minimizer.key.decode(k); // TO CORRECT
       //cerr << "MINIMIZER SEQUENCE (IN MINIMIZER INDEXING)  " << minimizer_sequence << endl;


       auto encoded_key = Key64::encode(minimizer_sequence);

      if (encoded_key.decode(k) != minimizer_sequence){
         throw runtime_error("THERES A PROBLEM WITH THE ENCODING OR DECODING");
      }


      // Find the node covering minimizer starting position.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= minimizer.offset)
      {
        node_start += node_length;
        ++iter;
        node_length = graph.get_length(*iter);
      }
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start };
      if(minimizer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      if(!Position::valid_offset(pos))
      {
        #pragma omp critical (cerr)
        {
          std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large" << std::endl;
        }
        std::exit(EXIT_FAILURE);
      }
      cache[thread_id].emplace_back(minimizer, pos);
    }

    if(cache[thread_id].size() >= MINIMIZER_CACHE_SIZE) { flush_cache(thread_id); }
  };

  /*
    Index the minimizers.
    We do a lot of redundant work by traversing both orientations and finding almost the same minimizers
    in each orientation. If we consider only the windows starting in forward (reverse) orientation,
    we may skip windows that cross from a reverse node to a forward node (from a forward node to a
    reverse node).
  */
  for_each_haplotype_window(graph, index.window_bp(), find_minimizers, (threads > 1), false);
  for(int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }

   // Print or return the counter after the function completes its execution.
  std::cout << "Total minimizers indexed: " << minimizer_count << std::endl;
}


#ifdef RYMER
//-------------------------------------------------------------------------------

/*
  Index the haplotypes in the graph. Insert the minimizers into the provided index.
  Function argument get_payload is used to generate the payload for each position
  stored in the index.
  The number of threads can be set through OMP.
*/
template<class KeyType>
void
index_haplotypes_rymer(const GBWTGraph& graph, MinimizerIndex<KeyType>& index,
                 const std::function<payload_type(const pos_t&)>& get_payload, unsigned int k=21)
{

  typedef typename MinimizerIndex<KeyType>::minimizer_type minimizer_type;

  int threads = omp_get_max_threads();

  // Minimizer caching. We only generate the payloads after we have removed duplicate positions.
  std::vector<std::vector<std::pair<minimizer_type, pos_t>>> cache(threads);
  constexpr size_t MINIMIZER_CACHE_SIZE = 1024;

     // Introduce a shared counter to track the number of indexed rymers.
  std::atomic<size_t> rymer_count(0);

  auto flush_cache = [&](int thread_id)
  {
    std::vector<std::pair<minimizer_type, pos_t>>& current_cache = cache[thread_id];
    gbwt::removeDuplicates(current_cache, false);
    std::vector<payload_type> payload;
    payload.reserve(current_cache.size());
    for(size_t i = 0; i < current_cache.size(); i++) { payload.push_back(get_payload(current_cache[i].second)); }
    #pragma omp critical (minimizer_index)
    {
      for(size_t i = 0; i < current_cache.size(); i++)
      {
        rymer_count++;
        index.insert(current_cache[i].first, current_cache[i].second, payload[i]);
        //cerr << "RYMER INDEX SIZE: " << index.size() << endl;
      }
    }
    cache[thread_id].clear();
  };

  // Minimizer finding.
  auto find_minimizers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {

    //std::cerr << "SEQ SIZE: " << seq.size() << std::endl;

    std::string rymer_seq = gbwtgraph::convertToRymerSpace(seq);
    //std::vector<minimizer_type> minimizers = index.minimizers(seq); // Calls syncmers() when appropriate.
    std::vector<minimizer_type> rymers = index.rymers(rymer_seq);

    for (auto & r : rymers){
      cerr << "RYMER KEY: " << r.key.get_key() << endl;
    }

    auto iter = traversal.begin();
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();

    Key64 thing;
    for(minimizer_type& rymer : rymers)
    {

      if(rymer.empty()){continue;}

      std::string rymer_sequence = rymer.key.decode_rymer(k); // TO CORRECT

      //Key64::key_type original_key = thing.get_original_kmer_key(rymer_sequence);
      //rymer.key = original_key;


      // Find the node covering minimizer starting position.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= rymer.offset)
      {
        node_start += node_length;
        ++iter;
        node_length = graph.get_length(*iter);
      }
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), rymer.offset - node_start };
      if(rymer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      if(!Position::valid_offset(pos))
      {
        #pragma omp critical (cerr)
        {
          std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large" << std::endl;
        }
        std::exit(EXIT_FAILURE);
      }

      cache[thread_id].emplace_back(rymer, pos);
    }

    if(cache[thread_id].size() >= MINIMIZER_CACHE_SIZE) { flush_cache(thread_id); }
  };

  /*
    Index the minimizers.
    We do a lot of redundant work by traversing both orientations and finding almost the same minimizers
    in each orientation. If we consider only the windows starting in forward (reverse) orientation,
    we may skip windows that cross from a reverse node to a forward node (from a forward node to a
    reverse node).
  */
  for_each_haplotype_window(graph, index.window_bp(), find_minimizers, (threads > 1), true);
  for(int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }

  // Print or return the counter after the function completes its execution.
  std::cout << "Total rymers indexed: " << rymer_count << std::endl;

}


} // namespace gbwtgraph
#endif

#endif // GBWTGRAPH_CONSTRUCTION_H
