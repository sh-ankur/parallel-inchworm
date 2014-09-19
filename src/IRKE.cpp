#include <stdlib.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iterator>
#include <math.h>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include "IRKE.hpp"
#include "Fasta_reader.hpp"
#include "sequenceUtil.hpp"
#include "KmerCounter.hpp"
#include "stacktrace.hpp"
#include "irke_common.hpp"
#include "IO_Processor.cpp"
#include "File_Map.cpp"
//#include <scorep/SCOREP_User.h>


IRKE::IRKE () {   // IRKE = Inchworm Recursive Kmer Extension

}

IRKE::IRKE (unsigned int kmer_length, unsigned int max_recursion, float min_seed_entropy, 
			unsigned int min_seed_coverage, float min_any_entropy, 
			bool pacman, bool crawl, unsigned int crawl_length, bool double_stranded) : kcounter(kmer_length, double_stranded) {
	
	MAX_RECURSION = max_recursion;
	MIN_SEED_ENTROPY = min_seed_entropy;
	MIN_SEED_COVERAGE = min_seed_coverage;
	MIN_ANY_ENTROPY = min_any_entropy;
	PACMAN = pacman;
	CRAWL = crawl;
	CRAWL_LENGTH = crawl_length;
	
	INCHWORM_ASSEMBLY_COUNTER = 0;
	
	DOUBLE_STRANDED_MODE = double_stranded;
	PRUNE_SINGLETON_READ_INTERVAL = 0;
	
}


void IRKE::set_prune_singleton_read_interval (unsigned long interval) {
	
	this->PRUNE_SINGLETON_READ_INTERVAL = interval;
	
}


void IRKE::build_graph(const string& fasta_filename, bool reassembleIworm, bool useKmers, float min_ratio_non_error) {
	
	string f = fasta_filename;
	if (useKmers)
		populate_kmers_using_mmap(f, min_ratio_non_error);
	else
		populate_kmers_using_mmap(f, min_ratio_non_error);
	
}

void IRKE::populate_kmers_using_mmap(string filename, float min_ratio_non_error)
{
  cerr << "MMAP: Called" << endl;
  time_t start, end;
  kmer_pair p,q;
  long long count = 0;
  IO_Processor io;
  int PARTITIONS=32;
  int tid;
  start = time(NULL);
  
  File_Map fm(filename, PARTITIONS, 35);
  
  unsigned long long logger[PARTITIONS + 1];
  
#pragma omp parallel for
  for(int i=0; i < (PARTITIONS + 1) ; i++ )
    logger[i] = -1;

int threadcounter = 0;
  
#pragma omp parallel shared(logger, threadcounter) private(io,tid) num_threads(PARTITIONS)
{
	tid = omp_get_thread_num();
    while((tid != (PARTITIONS - 1)) && (logger[tid + 1] == -1));
    io.initialize(fm.pmap[tid],fm.get_limit(tid), logger[tid +1], logger[tid]);
	#pragma omp critical(kmercounter)
	{
	    for(int j = 0 ; j < io.kmer_buffer.size() ; j++)
	    {
	    	try
	    	{
	      		kcounter.add_kmer(io.kmer_buffer[j].first, io.kmer_buffer[j].second);
	      	}
	      	catch(exception e)
	      	{
				cerr << "overlap rejected" << endl;
	      	}
	    }
  	}
}

  end = time(NULL);
  cerr << "kmers parsed:"<< kcounter.size() <<" in " << end - start << "sec" << endl;
  
  ofstream iworm_kmer_count_report_fh;
  iworm_kmer_count_report_fh.open("inchworm.kmer_count");
  iworm_kmer_count_report_fh << kcounter.size() << endl;
  iworm_kmer_count_report_fh.close();
 
}

void IRKE::populate_Kmers_from_kmers(const string& fasta_filename) {
	unsigned int kmer_length = kcounter.get_kmer_length();
	int i, myTid;
	unsigned long sum, 
        *record_counter = new unsigned long[omp_get_max_threads()];
	unsigned long start, end;

  // init record counter
  for (int i = 0; i < omp_get_max_threads(); i++) {
  	record_counter[i] = 0;
  }

	cerr << "-reading Kmer occurences..." << endl;
	start = time(NULL);

	Fasta_reader fasta_reader(fasta_filename); 

  #pragma omp parallel private (myTid)
	{
		myTid = omp_get_thread_num();
		record_counter[myTid] = 0;

		while (true) {
			Fasta_entry fe = fasta_reader.getNext();
			if (fe.get_sequence() == "") break;

			record_counter[myTid]++;

			if (IRKE_COMMON::MONITOR) {
				if (myTid == 0 && record_counter[myTid] % 100000 == 0)
					{
						sum = record_counter[0];
						for (i=1; i<omp_get_num_threads(); i++)
							sum+= record_counter[i];
						cerr << "\r [" << sum/1000000 << "M] Kmers parsed.     ";
					}
			}


			string seq = fe.get_sequence();
			if (seq.length() != kmer_length) {
				continue;
			}

			kmer_int_type_t kmer = kcounter.get_kmer_intval(seq);
			unsigned int count = atoi(fe.get_header().c_str());
			kcounter.add_kmer(kmer, count);
		}
	}
	end = time(NULL);

	sum = record_counter[0];
	for (i=1; i<omp_get_max_threads(); i++)
		sum+= record_counter[i];
	delete [] record_counter;

	cerr << endl << " done parsing " << sum << " Kmers, " << kcounter.size() << " added, taking " << (end-start) << " seconds." << endl;

    ofstream iworm_kmer_count_report_fh;
    iworm_kmer_count_report_fh.open("inchworm.kmer_count");
    iworm_kmer_count_report_fh << kcounter.size() << endl;
    iworm_kmer_count_report_fh.close();
    

	return;
}


void IRKE::populate_Kmers_from_fasta(const string& fasta_filename, bool reassembleIworm) {
	
	unsigned int kmer_length = kcounter.get_kmer_length();
	int i, myTid;
	unsigned long sum, 
                  *record_counter = new unsigned long[omp_get_max_threads()];
	unsigned long start, end;
    
    // init record counter
    for (int i = 0; i < omp_get_max_threads(); i++) {
        record_counter[i] = 0;
    }

	
	cerr << "-storing Kmers..." << endl;
	start = time(NULL);
	
	Fasta_reader fasta_reader(fasta_filename);

    unsigned int entry_num = 0;

    #pragma omp parallel private (myTid)
	{
		myTid = omp_get_thread_num();
		record_counter[myTid] = 0;
		
		while (fasta_reader.hasNext()) {
			Fasta_entry fe = fasta_reader.getNext();
            string accession = fe.get_accession();

            #pragma omp atomic            
            entry_num++;
            record_counter[myTid]++;
			
			if (IRKE_COMMON::MONITOR >= 4) {
				cerr << "[" << entry_num << "] acc: " << accession << ", by thread no: " << myTid << endl;;
			}
			else if (IRKE_COMMON::MONITOR) {
				if (myTid == 0 && record_counter[myTid] % 1000 == 0)
					{
						sum = record_counter[0];
						for (i=1; i<omp_get_num_threads(); i++)
							sum+= record_counter[i];
						cerr << "\r [" << sum << "] sequences parsed.     ";
					}
			}
			
			
			string seq = fe.get_sequence();
			
			if (seq.length() < kmer_length + 1) {
				continue;
			}
			
			if (reassembleIworm) {
				string accession = fe.get_accession();
				string header = fe.get_header();
				// get coverage value from iworm assembly
				vector<string> tokens;
				string_util::tokenize(accession, tokens, ";");
				if (tokens.size() < 2) {
					stringstream err;
					err << "Could not extract coverage value from accession: " << tokens[tokens.size()-1];
					throw(err.str());
				}
				string cov_s = tokens[tokens.size()-1];
				unsigned int cov_val = atoi(cov_s.c_str());
				
				// get Kmer value from header
				vector<string> header_toks;
				string_util::tokenize(header, header_toks, " ");
				if (header_toks.size() < 5) {
					stringstream err;
					err << "Fasta header: " << header << " lacks expected format including Kmer length from previous inchworm assembly run";
					throw(err.str());
				}
				
				unsigned int kmer_val = atoi(header_toks[2].c_str());
				
				unsigned int normalized_coverage_val = static_cast<unsigned int> (cov_val * kmer_val / 25.0 + 0.5);
				
				if (IRKE_COMMON::MONITOR >= 1) {
					cerr << "Adding inchworm assembly " << accession 
						 << " K: " << kmer_val << " Cov: " << cov_val 
						 << " with coverage: " << normalized_coverage_val << endl;
				}
				if (cov_val < 1) {
					stringstream err;
					err << "error parsing coverage value from accession: " << accession;
					throw(err.str());
				}
				kcounter.add_sequence(seq, normalized_coverage_val);
			}
			else {
				kcounter.add_sequence(seq);
			}
			
			// remove singleton kmers at read interval to minimize memory requirements.
			if (PRUNE_SINGLETON_READ_INTERVAL > 0 
				&& 
				myTid == 0
				&&
				record_counter[myTid]/omp_get_num_threads() % PRUNE_SINGLETON_READ_INTERVAL == 0) {
				if (IRKE_COMMON::MONITOR >= 1) {
					cerr << "Reached singleton kmer pruning interval at read count: " << record_counter << endl;
				}
				prune_kmers_min_count(1);
			}
			
			
		}
	}
	end = time(NULL);
	
	sum = record_counter[0];
	for (i=1; i<omp_get_max_threads(); i++)
		sum+= record_counter[i];
    delete [] record_counter;
	
	cerr << endl << " done parsing " << sum << " sequences, extracted " << kcounter.size() << " kmers, taking " << (end-start) << " seconds." << endl;
	
	
	return;
}


void IRKE::prune_kmers_min_count(unsigned int min_count) {
	
	// proxy, sends message to kcounter member
	
	kcounter.prune_kmers_min_count(min_count);
}


void IRKE::prune_kmers_min_entropy(float min_entropy) {
	
	// proxy, send message to kcounter
	
	kcounter.prune_kmers_min_entropy(min_entropy);
	
	
}


bool IRKE::prune_kmer_extensions(float min_ratio_non_error) {
	
	// proxy, send message to kcounter
	
	return(kcounter.prune_kmer_extensions(min_ratio_non_error));
}


bool IRKE::prune_some_kmers(unsigned int min_count, float min_entropy, bool prune_error_kmers, float min_ratio_non_error) {

	return(kcounter.prune_some_kmers(min_count, min_entropy, prune_error_kmers, min_ratio_non_error));
}

bool IRKE::fast_pruning(unsigned int min_count, float min_entropy, bool prune_error_kmers, float min_ratio_non_error) {
	kcounter.fast_pruning(min_count, min_entropy, prune_error_kmers, min_ratio_non_error);
	sorted_kmers.clear();
}


unsigned long IRKE::get_graph_size() {
	
	// proxy call
	
	return(kcounter.size());
	
}


void IRKE::traverse_path(KmerCounter& kcounter, Kmer_Occurence_Pair seed_kmer, Kmer_visitor& visitor,
						 Kmer_visitor& place_holder, float MIN_CONNECTIVITY_RATIO, unsigned int depth) {
	
	if (IRKE_COMMON::MONITOR >= 3) {
		cerr << "traverse_path, depth: " << depth << ", kmer: " << kcounter.get_kmer_string(seed_kmer.first) <<  endl;
	}
	
	
	// check to see if visited already
	if (visitor.exists(seed_kmer.first)) {
		// already visited
		if (IRKE_COMMON::MONITOR >= 3) {
			cout << "\talready visited " << kcounter.get_kmer_string(seed_kmer.first) << endl;
		}
		
		return;
	}
	
	// check if at the end of the max recursion, and if so, don't visit it but set a placeholder.
	if (depth > MAX_RECURSION) {
		place_holder.add(seed_kmer.first);
		return;
	}
	
	visitor.add(seed_kmer.first);
    
	// try each of the forward paths from the kmer:  
	vector<Kmer_Occurence_Pair> forward_candidates = kcounter.get_forward_kmer_candidates(seed_kmer.first);
	
	for(unsigned int i = 0; i < forward_candidates.size(); i++) {
		Kmer_Occurence_Pair kmer = forward_candidates[i];
		
		if (kmer.second	&& exceeds_min_connectivity(kcounter, seed_kmer, kmer, MIN_CONNECTIVITY_RATIO)) {
			
			traverse_path(kcounter, kmer, visitor, place_holder, MIN_CONNECTIVITY_RATIO, depth + 1);
		}
	}
	
	// try each of the reverse paths from the kmer:
	vector<Kmer_Occurence_Pair> reverse_candidates = kcounter.get_reverse_kmer_candidates(seed_kmer.first);
	
	for (unsigned int i = 0; i < reverse_candidates.size(); i++) {
		Kmer_Occurence_Pair kmer = reverse_candidates[i];
		
		if (kmer.second	&& exceeds_min_connectivity(kcounter, seed_kmer, kmer, MIN_CONNECTIVITY_RATIO)) {
			
			traverse_path(kcounter, kmer, visitor, place_holder, MIN_CONNECTIVITY_RATIO, depth + 1);
		}
	}
	
	
	return;
	
}


string add_fasta_seq_line_breaks(string& sequence, int interval) {
    
    stringstream fasta_seq;
    
    int counter = 0;
    for (string::iterator it = sequence.begin(); it != sequence.end(); it++) {
        counter++;
        
        fasta_seq << *it;
        if (counter % interval == 0 && (it + 1) != sequence.end()) {
            fasta_seq << endl;
        }
    }

    return(fasta_seq.str());
}

void IRKE::parallel_compute_sequence_assemblies(float min_connectivity, unsigned int MIN_ASSEMBLY_LENGTH, unsigned int MIN_ASSEMBLY_COVERAGE,
									   bool WRITE_COVERAGE, string COVERAGE_OUTPUT_FILENAME) {

	// use self kcounter
	parallel_compute_sequence_assemblies_same_results(kcounter, min_connectivity, MIN_ASSEMBLY_LENGTH, MIN_ASSEMBLY_COVERAGE, WRITE_COVERAGE, COVERAGE_OUTPUT_FILENAME);
	
	return;
	
}

void IRKE::compute_sequence_assemblies(float min_connectivity, unsigned int MIN_ASSEMBLY_LENGTH, unsigned int MIN_ASSEMBLY_COVERAGE,
									   bool WRITE_COVERAGE, string COVERAGE_OUTPUT_FILENAME) {

	// use self kcounter
	compute_sequence_assemblies(kcounter, min_connectivity, MIN_ASSEMBLY_LENGTH, MIN_ASSEMBLY_COVERAGE, WRITE_COVERAGE, COVERAGE_OUTPUT_FILENAME);
	
	return;
	
}


void IRKE::parallel_compute_sequence_assemblies(KmerCounter& kcounter, float min_connectivity,
									   unsigned int MIN_ASSEMBLY_LENGTH, unsigned int MIN_ASSEMBLY_COVERAGE,
									   bool WRITE_COVERAGE, string COVERAGE_OUTPUT_FILENAME)
{
	int slaves_forward = 2, slaves_reverse = 2;
	
	time_t stt, edd;
	
	stt = time(NULL);

	if (! got_sorted_kmers_flag) 
	{
        stringstream error;
        error << stacktrace() << " Error, must populate_sorted_kmers_list() before computing sequence assemblies" << endl;
        throw(error.str());
	}
    
	vector<Kmer_Occurence_Pair>& kmers = sorted_kmers; 
	
	unsigned long init_size = kcounter.size();
    //SCOREP_USER_REGION_DEFINE(master_handle)
    //SCOREP_USER_REGION_DEFINE(slave_forward_handle)
    //SCOREP_USER_REGION_DEFINE(slave_reverse_handle)
    
    cerr << "Total kcounter hash size: " << init_size << " vs. sorted list size: " << kmers.size() << endl;

	unsigned int kmer_length = kcounter.get_kmer_length();
	ofstream coverage_writer;
	if (WRITE_COVERAGE) 
	{
		coverage_writer.open(COVERAGE_OUTPUT_FILENAME.c_str());
	}
	
	precomputations forward_precomputations, reverse_precomputations;
	typedef precomputations::iterator preiterator;

	int tid;

	#pragma omp parallel num_threads(slaves_forward + slaves_reverse + 1) private(tid) shared(forward_precomputations, reverse_precomputations)
	{	
		//omp_set_nested(1);
		//cerr << omp_get_nested() << endl;
		tid = omp_get_thread_num();
		int mapmissf = 0, mapmissr = 0;
		  
		if(tid == 0)
		{	
			
			//SCOREP_USER_REGION_BEGIN(master_handle, "Master", SCOREP_USER_REGION_TYPE_COMMON)
			int flag_ = 0;
			int contigcounter = 0;
			int fhitcounter = 0, fhitandmiss = 0, rhitcounter = 0, rhitandmiss = 0;
			Path_n_count_pair selected_path_n_pair_forward, selected_path_n_pair_reverse;
			for (unsigned int i = 0; i < kmers.size(); ++i) 
			{
				flag_ = 0;
				unsigned long kmer_counter_size = kcounter.size();
				if (kmer_counter_size > init_size) 
				{
					
					stringstream error;
					error << stacktrace() << "Error, Kcounter size has grown from " << init_size
						  << " to " << kmer_counter_size << endl;
					throw (error.str());
				}
				
				kmer_int_type_t kmer = kmers[i].first;
				unsigned int kmer_count = kcounter.get_kmer_count(kmer);
		        

				if (kmer_count == 0) 
				{
					continue;
				}
				
				if(is_set_left(kcounter.find_kmer(kmer)->second) || is_set_right(kcounter.find_kmer(kmer)->second))
				{
				  continue;
				}
		        
		        
				if (IRKE_COMMON::MONITOR >= 2) 
				{
					cerr << "SEED kmer: " << kcounter.get_kmer_string(kmer) << ", count: " << kmer_count << endl;
				}

		        

				if (kmer == revcomp_val(kmer, kmer_length)) 
				{
					
		            if (IRKE_COMMON::MONITOR >= 2) 
		            {
		                cerr << "SEED kmer: " << kcounter.get_kmer_string(kmer) << " is palidnromic.  Skipping. " << endl;
		            }
		            
		            continue;
				}
				
		        
				if (kmer_count < MIN_SEED_COVERAGE) 
				{
					if (IRKE_COMMON::MONITOR >= 2) 
					{
		                cerr << "-seed has insufficient coverage, skipping" << endl;
		            }
		            
		            continue;
				}
				
				
				float entropy = compute_entropy(kmer, kmer_length);
				
				
				if (entropy < MIN_SEED_ENTROPY) 
				{

				    if (IRKE_COMMON::MONITOR >= 2) 
				    {
					cerr << "-skipping seed due to low entropy: " << entropy << endl;
				    }
		            
				    continue;
				}
				
				int t;

				vector<kmer_int_type_t> joined_path;

				
				int flagf = 0;
				int flagr = 0;
				unsigned int fcount, rcount;
				preiterator itf = forward_precomputations.find(kmer);
				preiterator itr = reverse_precomputations.find(kmer);

				if(itf != forward_precomputations.end())
				{
					fhitcounter += 1;
			    	for(int j = 0 ; j < itf->second.first.size() ; j++)
			      	{
						Kmer_counter_map_iterator itt = kcounter.find_kmer(itf->second.first[j]);
						if(is_set_right(itt->second))
						{
						  flagf = 1;
						  fhitandmiss += 1;
						  break;
						}
			      	}
			  	}
			  	else
			  	{
			  		flagf = 1;
					mapmissf += 1;
			  	}

			  	if(itr != reverse_precomputations.end())
				{
					rhitcounter += 1;
			    	for(int j = 0 ; j < itr->second.first.size() ; j++)
			      	{
						Kmer_counter_map_iterator itt = kcounter.find_kmer(itr->second.first[j]);
						if(is_set_left(itt->second))
						{
						  flagr = 1;
						  rhitandmiss += 1;
						  break;
						}
			      	}
			  	}
			  	else
			  	{
			  		flagr = 1;
					mapmissr += 1;
			  	}

			  	if(flagf)
			  	{
					Kmer_visitor visitorf(kmer_length, DOUBLE_STRANDED_MODE);
					selected_path_n_pair_forward = inchworm(kcounter, 'F', kmer, visitorf, min_connectivity); 
			     	//visitorf.clear();
			  	}

			  	if(flagr)
			  	{
					Kmer_visitor visitorr(kmer_length, DOUBLE_STRANDED_MODE);
					selected_path_n_pair_reverse = inchworm(kcounter, 'R', kmer, visitorr, min_connectivity);
					//visitorr.clear();
			  	}

			  	if(!flagf)
			  	{
					if(!flagr)
					{
						fcount = itf->second.second;
						rcount = itr->second.second;
						joined_path = _join_forward_n_reverse_paths(itr->second.first, kmer, itf->second.first);
					}
					else
					{
						fcount = itf->second.second;
						rcount = selected_path_n_pair_reverse.second;
						joined_path = _join_forward_n_reverse_paths(selected_path_n_pair_reverse.first, kmer, itf->second.first);
					}
			  	}
			  	else
			  	{
					if(!flagr)
					{
						fcount = selected_path_n_pair_forward.second;
						rcount = itr->second.second;
						joined_path = _join_forward_n_reverse_paths(itr->second.first, kmer, selected_path_n_pair_forward.first);
					}
					else
					{
						fcount = selected_path_n_pair_forward.second;
						rcount = selected_path_n_pair_reverse.second;
						joined_path = _join_forward_n_reverse_paths(selected_path_n_pair_reverse.first, kmer, selected_path_n_pair_forward.first);
					}
			  	}

					
				unsigned int total_counts = fcount + rcount + kmer_count; //kcounter.get_kmer_count(kmer); 
				
				vector<unsigned int> assembly_base_coverage;
				string sequence = reconstruct_path_sequence(kcounter, joined_path, assembly_base_coverage);
				
				int avg_cov =  static_cast<int> ( (float)total_counts/(sequence.length()-kcounter.get_kmer_length() +1) + 0.5);
				
				
				if (sequence.length() >= MIN_ASSEMBLY_LENGTH && avg_cov >= MIN_ASSEMBLY_COVERAGE) 
				{
					
					INCHWORM_ASSEMBLY_COUNTER++;
					
					stringstream headerstream;
					
					contigcounter++;
					
					
					headerstream << ">a" << INCHWORM_ASSEMBLY_COUNTER << ";" << avg_cov 
		                         << " total_counts: " <<  total_counts << " Fpath: " << selected_path_n_pair_forward.second << " Rpath: " << selected_path_n_pair_reverse.second << " Seed: " << kmer_count 
								 << " K: " << kmer_length
								 << " length: " << sequence.length();
					
					string header = headerstream.str();
					
					 sequence = add_fasta_seq_line_breaks(sequence, 60);
		            
					cout << header << endl << sequence << endl;
					
					if (WRITE_COVERAGE) 
					{
						
						coverage_writer << header << endl;
						
						for (unsigned int j = 0; j < assembly_base_coverage.size(); ++j) 
						{
							coverage_writer << assembly_base_coverage[j];
							if ( (j+1) % 30 == 0) 
							{
								coverage_writer << endl;
							}
							else 
							{
								coverage_writer << " ";
							}
						}
						coverage_writer << endl;
					}
					
				}
				
				if(flagr)
			  	{
					
			  	}

			  	if(!flagf)
			  	{
			  		for(unsigned int j = 0 ; j < itf->second.first.size() ; ++j)
					{
				    	kmer_int_type_t k = itf->second.first[j];
				    	kcounter.set_forward(k);
					}

					if(!flagr)
					{											
						for(unsigned int j = 0 ; j < itr->second.first.size() ; ++j)
						{
				    		kmer_int_type_t k = itr->second.first[j];
				    		kcounter.set_reverse(k);
						}
					}
					else
					{
						for(unsigned int j = 0 ; j < selected_path_n_pair_reverse.first.size() ; ++j)
						{
				    		kmer_int_type_t k = selected_path_n_pair_reverse.first[j];
				    		kcounter.set_reverse(k);
						}
					}
			  	}
			  	else
			  	{
			  		for(unsigned int j = 0 ; j < selected_path_n_pair_forward.first.size() ; ++j)
					{
				    	kmer_int_type_t k = selected_path_n_pair_forward.first[j];
				   		kcounter.set_forward(k);
					}

					if(!flagr)
					{											
						for(unsigned int j = 0 ; j < itr->second.first.size() ; ++j)
						{
				    		kmer_int_type_t k = itr->second.first[j];
				    		kcounter.set_reverse(k);
						}
					}
					else
					{
						for(unsigned int j = 0 ; j < selected_path_n_pair_reverse.first.size() ; ++j)
						{
				    		kmer_int_type_t k = selected_path_n_pair_reverse.first[j];
				    		kcounter.set_reverse(k);
						}
					}
			  	}
		        
				kcounter.set_seed(kmer);
				
		    }
		    cerr << "Total reported contigs: " << contigcounter << endl;
		    cerr << "forward precompute miss: " << mapmissf << endl << "reverse precompute miss: " << mapmissr << endl;
		    cerr << "forward precomputation hits: " << fhitcounter << " || Hit but recomputed: " << fhitandmiss << endl;
		    cerr << "reverse precomputation hits: " << rhitcounter << " || Hit but recomputed: " << rhitandmiss;
		}
		// slaves
		else if(tid <= slaves_forward)
		{
			//SCOREP_USER_REGION_BEGIN(slave_forward_handle, "Worker-Forward", SCOREP_USER_REGION_TYPE_COMMON)
			int slave_spawn_id;
			int cf = 0;
			Path_n_count_pair selected_path_n_pair_forward;
			for(unsigned int j = tid - 1 ; j < kmers.size() ; j = j + slaves_forward)
			{	
				kmer_int_type_t kmer = kmers[j].first;
				if(is_set_left(kcounter.find_kmer(kmer)->second) || is_set_right(kcounter.find_kmer(kmer)->second))
				{
				  continue;
				}
				Kmer_visitor visitor(kmer_length, DOUBLE_STRANDED_MODE);
				selected_path_n_pair_forward = inchworm(kcounter, 'F', kmer, visitor, min_connectivity);
				visitor.clear();
				#pragma omp critical(forwardmap)
				{
					forward_precomputations.insert(make_pair<kmer_int_type_t, Path_n_count_pair>(kmers[j].first,selected_path_n_pair_forward));
					cf++;
					//cerr << tid << "F ";
				}
			} 
			//#pragma omp critical(print)
			//{
			//cerr << "forward precompute worker " << tid << " computed :" << cf << " contigs" << endl;
			//}//
			//SCOREP_USER_REGION_END(slave_forward_handle)
		}
		else
		{
			//SCOREP_USER_REGION_BEGIN(slave_reverse_handle, "Worker-Reverse", SCOREP_USER_REGION_TYPE_COMMON)
			int slave_spawn_id;
			int cr = 0;
			Path_n_count_pair selected_path_n_pair_reverse;
			for(unsigned int j = tid - slaves_forward - 1 ; j < kmers.size() ; j = j + slaves_reverse)
			{	
				kmer_int_type_t kmer = kmers[j].first;
				if(is_set_left(kcounter.find_kmer(kmer)->second) || is_set_right(kcounter.find_kmer(kmer)->second))
				{
				  continue;
				}
				Kmer_visitor visitor(kmer_length, DOUBLE_STRANDED_MODE);
				selected_path_n_pair_reverse = inchworm(kcounter, 'R', kmer, visitor, min_connectivity);
				visitor.clear();
				#pragma omp critical(reversemap)
				{
					reverse_precomputations.insert(make_pair<kmer_int_type_t, Path_n_count_pair>(kmers[j].first,selected_path_n_pair_reverse));
					cr++;
					//cerr << tid << "F ";
				}
			} 
		}
	
		if (IRKE_COMMON::MONITOR) 
		{
			cerr << endl;
		}
		
		if (WRITE_COVERAGE) 
		{
			coverage_writer.close();
		}
	        
	    // drop sorted kmer list as part of cleanup
	    clear_sorted_kmers_list();
	}
	
	edd = time(NULL);
	
	cerr << "Total time for assembly: " << edd - stt << " seconds" << endl;
	
	return; // end of runIRKE
}

void IRKE::parallel_compute_sequence_assemblies_same_results(KmerCounter& kcounter, float min_connectivity,
									   unsigned int MIN_ASSEMBLY_LENGTH, unsigned int MIN_ASSEMBLY_COVERAGE,
									   bool WRITE_COVERAGE, string COVERAGE_OUTPUT_FILENAME)
{
	int slaves_forward = 2, slaves_reverse = 2;
	
	time_t stt, edd;
	
	stt = time(NULL);

	if (! got_sorted_kmers_flag) 
	{
        stringstream error;
        error << stacktrace() << " Error, must populate_sorted_kmers_list() before computing sequence assemblies" << endl;
        throw(error.str());
	}
    
	vector<Kmer_Occurence_Pair>& kmers = sorted_kmers; 
	
	unsigned long init_size = kcounter.size();
    
    cerr << "Total kcounter hash size: " << init_size << " vs. sorted list size: " << kmers.size() << endl;

	unsigned int kmer_length = kcounter.get_kmer_length();
	ofstream coverage_writer;
	if (WRITE_COVERAGE) 
	{
		coverage_writer.open(COVERAGE_OUTPUT_FILENAME.c_str());
	}
	
	precomputations forward_precomputations, reverse_precomputations;
	typedef precomputations::iterator preiterator;

	int tid;

	#pragma omp parallel num_threads(slaves_forward + slaves_reverse + 1) private(tid) shared(forward_precomputations, reverse_precomputations)
	{	
		tid = omp_get_thread_num();
		int mapmissf = 0, mapmissr = 0;
		  
		if(tid == 0)
		{	
			
			int flag_ = 0;
			int contigcounter = 0;
			int fhitcounter = 0, fhitandmiss = 0, rhitcounter = 0, rhitandmiss = 0;
			Path_n_count_pair selected_path_n_pair_forward, selected_path_n_pair_reverse;
			for (unsigned int i = 0; i < kmers.size(); ++i) 
			{
				flag_ = 0;
				unsigned long kmer_counter_size = kcounter.size();
				if (kmer_counter_size > init_size) 
				{
					
					stringstream error;
					error << stacktrace() << "Error, Kcounter size has grown from " << init_size
						  << " to " << kmer_counter_size << endl;
					throw (error.str());
				}
				
				kmer_int_type_t kmer = kmers[i].first;
				unsigned int kmer_count = kcounter.get_kmer_count(kmer);
		        

				if (kmer_count == 0) 
				{
					continue;
				}
				
				if(is_set_left(kcounter.find_kmer(kmer)->second) || is_set_right(kcounter.find_kmer(kmer)->second))
				{
				  continue;
				}
		        
		        
				if (IRKE_COMMON::MONITOR >= 2) 
				{
					cerr << "SEED kmer: " << kcounter.get_kmer_string(kmer) << ", count: " << kmer_count << endl;
				}

		        

				if (kmer == revcomp_val(kmer, kmer_length)) 
				{
					
		            if (IRKE_COMMON::MONITOR >= 2) 
		            {
		                cerr << "SEED kmer: " << kcounter.get_kmer_string(kmer) << " is palidnromic.  Skipping. " << endl;
		            }
		            
		            continue;
				}
				
		        
				if (kmer_count < MIN_SEED_COVERAGE) 
				{
					if (IRKE_COMMON::MONITOR >= 2) 
					{
		                cerr << "-seed has insufficient coverage, skipping" << endl;
		            }
		            
		            continue;
				}
				
				
				float entropy = compute_entropy(kmer, kmer_length);
				
				
				if (entropy < MIN_SEED_ENTROPY) 
				{

				    if (IRKE_COMMON::MONITOR >= 2) 
				    {
					cerr << "-skipping seed due to low entropy: " << entropy << endl;
				    }
		            
				    continue;
				}
				
				int t;

				vector<kmer_int_type_t> joined_path;

				
				int flagf = 0;
				int flagr = 0;
				unsigned int fcount, rcount;
				preiterator itf = forward_precomputations.find(kmer);
				preiterator itr = reverse_precomputations.find(kmer);

				if(itf != forward_precomputations.end())
				{
					fhitcounter += 1;
			    	for(int j = 0 ; j < itf->second.first.size() ; j++)
			      	{
						Kmer_counter_map_iterator itt = kcounter.find_kmer(itf->second.first[j]);
						if(is_set_right(itt->second) || is_set_right(itt->second))
						{
						  flagf = 1;
						  fhitandmiss += 1;
						  break;
						}
			      	}
			  	}
			  	else
			  	{
			  		flagf = 1;
					mapmissf += 1;
			  	}

			  	if(flagf)
			  	{
					Kmer_visitor visitorf(kmer_length, DOUBLE_STRANDED_MODE);
					selected_path_n_pair_forward = inchworm(kcounter, 'F', kmer, visitorf, min_connectivity); 
			  	}

			  	if(itr != reverse_precomputations.end())
				{
					rhitcounter += 1;
			    	for(int j = 0 ; j < itr->second.first.size() ; j++)
			      	{
						Kmer_counter_map_iterator itt = kcounter.find_kmer(itr->second.first[j]);
						if(is_set_left(itt->second) || is_set_right(itt->second))
						{
						  flagr = 1;
						  rhitandmiss += 1;
						  break;
						}
			      	}
			  	}
			  	else
			  	{
			  		flagr = 1;
					mapmissr += 1;
			  	}


			  	if(flagr)
			  	{
					Kmer_visitor visitorr(kmer_length, DOUBLE_STRANDED_MODE);
					selected_path_n_pair_reverse = inchworm(kcounter, 'R', kmer, visitorr, min_connectivity);
					//visitorr.clear();
			  	}

			  	if(!flagf)
			  	{
					if(!flagr)
					{
						fcount = itf->second.second;
						rcount = itr->second.second;
						joined_path = _join_forward_n_reverse_paths(itr->second.first, kmer, itf->second.first);
					}
					else
					{
						fcount = itf->second.second;
						rcount = selected_path_n_pair_reverse.second;
						joined_path = _join_forward_n_reverse_paths(selected_path_n_pair_reverse.first, kmer, itf->second.first);
					}
			  	}
			  	else
			  	{
					if(!flagr)
					{
						fcount = selected_path_n_pair_forward.second;
						rcount = itr->second.second;
						joined_path = _join_forward_n_reverse_paths(itr->second.first, kmer, selected_path_n_pair_forward.first);
					}
					else
					{
						fcount = selected_path_n_pair_forward.second;
						rcount = selected_path_n_pair_reverse.second;
						joined_path = _join_forward_n_reverse_paths(selected_path_n_pair_reverse.first, kmer, selected_path_n_pair_forward.first);
					}
			  	}

				unsigned int total_counts = fcount + rcount + kmer_count; //kcounter.get_kmer_count(kmer); 
				
				vector<unsigned int> assembly_base_coverage;
				string sequence = reconstruct_path_sequence(kcounter, joined_path, assembly_base_coverage);
				
				int avg_cov =  static_cast<int> ( (float)total_counts/(sequence.length()-kcounter.get_kmer_length() +1) + 0.5);
				
				
				if (sequence.length() >= MIN_ASSEMBLY_LENGTH && avg_cov >= MIN_ASSEMBLY_COVERAGE) 
				{
					
					INCHWORM_ASSEMBLY_COUNTER++;
					
					stringstream headerstream;
					
					contigcounter++;
					
					
					headerstream << ">a" << INCHWORM_ASSEMBLY_COUNTER << ";" << avg_cov 
		                         << " total_counts: " <<  total_counts << " Fpath: " << selected_path_n_pair_forward.second << " Rpath: " << selected_path_n_pair_reverse.second << " Seed: " << kmer_count 
								 << " K: " << kmer_length
								 << " length: " << sequence.length();
					
					string header = headerstream.str();
					
					 sequence = add_fasta_seq_line_breaks(sequence, 60);
		            
					cout << header << endl << sequence << endl;
					
					if (WRITE_COVERAGE) 
					{
						
						coverage_writer << header << endl;
						
						for (unsigned int j = 0; j < assembly_base_coverage.size(); ++j) 
						{
							coverage_writer << assembly_base_coverage[j];
							if ( (j+1) % 30 == 0) 
							{
								coverage_writer << endl;
							}
							else 
							{
								coverage_writer << " ";
							}
						}
						coverage_writer << endl;
					}
					
				}
				
				// remove path

				if(flagr)
			  	{
					
			  	}

			  	if(!flagf)
			  	{
			  		for(unsigned int j = 0 ; j < itf->second.first.size() ; ++j)
					{
				    	kmer_int_type_t k = itf->second.first[j];
				    	kcounter.set_forward(k);
					}

					if(!flagr)
					{											
						for(unsigned int j = 0 ; j < itr->second.first.size() ; ++j)
						{
				    		kmer_int_type_t k = itr->second.first[j];
				    		kcounter.set_reverse(k);
						}
					}
					else
					{
						for(unsigned int j = 0 ; j < selected_path_n_pair_reverse.first.size() ; ++j)
						{
				    		kmer_int_type_t k = selected_path_n_pair_reverse.first[j];
				    		kcounter.set_reverse(k);
						}
					}
			  	}
			  	else
			  	{
			  		for(unsigned int j = 0 ; j < selected_path_n_pair_forward.first.size() ; ++j)
					{
				    	kmer_int_type_t k = selected_path_n_pair_forward.first[j];
				   		kcounter.set_forward(k);
					}

					if(!flagr)
					{											
						for(unsigned int j = 0 ; j < itr->second.first.size() ; ++j)
						{
				    		kmer_int_type_t k = itr->second.first[j];
				    		kcounter.set_reverse(k);
						}
					}
					else
					{
						for(unsigned int j = 0 ; j < selected_path_n_pair_reverse.first.size() ; ++j)
						{
				    		kmer_int_type_t k = selected_path_n_pair_reverse.first[j];
				    		kcounter.set_reverse(k);
						}
					}
			  	}
		        
				kcounter.set_seed(kmer);
				
		    }
		    //SCOREP_USER_REGION_END(master_handle)
		    cerr << "Total reported contigs: " << contigcounter << endl;
		    cerr << "forward precompute miss: " << mapmissf << endl << "reverse precompute miss: " << mapmissr << endl;
		    cerr << "forward precomputation hits: " << fhitcounter << " || Hit but recomputed: " << fhitandmiss << endl;
		    cerr << "reverse precomputation hits: " << rhitcounter << " || Hit but recomputed: " << rhitandmiss;
		}
		// slaves
		else if(tid <= slaves_forward)
		{
			//SCOREP_USER_REGION_BEGIN(slave_forward_handle, "Worker-Forward", SCOREP_USER_REGION_TYPE_COMMON)
			int slave_spawn_id;
			int cf = 0;
			Path_n_count_pair selected_path_n_pair_forward;
			for(unsigned int j = tid - 1 ; j < kmers.size() ; j = j + slaves_forward)
			{	
				kmer_int_type_t kmer = kmers[j].first;
				if(is_set_left(kcounter.find_kmer(kmer)->second) || is_set_right(kcounter.find_kmer(kmer)->second))
				{
				  continue;
				}
				Kmer_visitor visitor(kmer_length, DOUBLE_STRANDED_MODE);
				selected_path_n_pair_forward = inchworm(kcounter, 'F', kmer, visitor, min_connectivity);
				visitor.clear();
				#pragma omp critical(forwardmap)
				{
					forward_precomputations.insert(make_pair<kmer_int_type_t, Path_n_count_pair>(kmers[j].first,selected_path_n_pair_forward));
					cf++;
					//cerr << tid << "F ";
				}
			} 
			//#pragma omp critical(print)
			//{
			//cerr << "forward precompute worker " << tid << " computed :" << cf << " contigs" << endl;
			//}//
			//SCOREP_USER_REGION_END(slave_forward_handle)
		}
		else
		{
			//SCOREP_USER_REGION_BEGIN(slave_reverse_handle, "Worker-Reverse", SCOREP_USER_REGION_TYPE_COMMON)
			int slave_spawn_id;
			int cr = 0;
			Path_n_count_pair selected_path_n_pair_reverse;
			for(unsigned int j = tid - slaves_forward - 1 ; j < kmers.size() ; j = j + slaves_reverse)
			{	
				kmer_int_type_t kmer = kmers[j].first;
				if(is_set_left(kcounter.find_kmer(kmer)->second) || is_set_right(kcounter.find_kmer(kmer)->second))
				{
				  continue;
				}
				Kmer_visitor visitor(kmer_length, DOUBLE_STRANDED_MODE);
				selected_path_n_pair_reverse = inchworm(kcounter, 'R', kmer, visitor, min_connectivity);
				visitor.clear();
				#pragma omp critical(reversemap)
				{
					reverse_precomputations.insert(make_pair<kmer_int_type_t, Path_n_count_pair>(kmers[j].first,selected_path_n_pair_reverse));
					cr++;
					//cerr << tid << "F ";
				}
			} 
		    //#pragma omp critical(print)
		    //{
			//cerr << "reverse precompute worker " << tid - slaves_forward - 1 << " computed :" << cr << " contigs" << endl;
		    //}   //
		    //SCOREP_USER_REGION_END(slave_reverse_handle)
		}
	
		if (IRKE_COMMON::MONITOR) 
		{
			cerr << endl;
		}
		
		if (WRITE_COVERAGE) 
		{
			coverage_writer.close();
		}
	        
	    // drop sorted kmer list as part of cleanup
	    clear_sorted_kmers_list();
	}
	
	edd = time(NULL);
	
	cerr << "Total time for assembly: " << edd - stt << " seconds" << endl;
	
	return; // end of runIRKE
}

void IRKE::compute_sequence_assemblies(KmerCounter& kcounter, float min_connectivity,
									   unsigned int MIN_ASSEMBLY_LENGTH, unsigned int MIN_ASSEMBLY_COVERAGE,
									   bool WRITE_COVERAGE, string COVERAGE_OUTPUT_FILENAME) {
	int contigcounter = 0;
	time_t start, end;
	start = time(NULL);
    if (! got_sorted_kmers_flag) {
        stringstream error;
        error << stacktrace() << " Error, must populate_sorted_kmers_list() before computing sequence assemblies" << endl;
        throw(error.str());
    }
    

	//vector<Kmer_counter_map_iterator>& kmers = sorted_kmers; 
    vector<Kmer_Occurence_Pair>& kmers = sorted_kmers; 
	
	unsigned long init_size = kcounter.size();
    
    cerr << "Total kcounter hash size: " << init_size << " vs. sorted list size: " << kmers.size() << endl;

	unsigned int kmer_length = kcounter.get_kmer_length();
	ofstream coverage_writer;
	if (WRITE_COVERAGE) {
		coverage_writer.open(COVERAGE_OUTPUT_FILENAME.c_str());
	}
	    

	// string s = "before.kmers";
	// kcounter.dump_kmers_to_file(s);
	
	for (unsigned int i = 0; i < kmers.size(); i++) {
		
		// cerr << "round: " << i << endl;
		
		unsigned long kmer_counter_size = kcounter.size();
		if (kmer_counter_size > init_size) {
			
			// string s = "after.kmers";
			// kcounter.dump_kmers_to_file(s);
			
			stringstream error;
			error << stacktrace() << "Error, Kcounter size has grown from " << init_size
				  << " to " << kmer_counter_size << endl;
			throw (error.str());
		}
		
		
		//kmer_int_type_t kmer = kmers[i]->first;
		//unsigned int kmer_count = kmers[i]->second;

        kmer_int_type_t kmer = kmers[i].first;
		// unsigned int kmer_count = kmers[i].second;  // NO!!!  Use for sorting, but likely zeroed out in the hashtable after contig construction
        unsigned int kmer_count = kcounter.get_kmer_count(kmer);
        
        
        // cout << "SEED kmer: " << kcounter.get_kmer_string(kmer) << ", count: " << kmer_count << endl;

        //continue;
        
		if (kmer_count == 0) {
			continue;
		}
        
        
		if (IRKE_COMMON::MONITOR >= 2) {
			cerr << "SEED kmer: " << kcounter.get_kmer_string(kmer) << ", count: " << kmer_count << endl;
		}

        

		if (kmer == revcomp_val(kmer, kmer_length)) {
			// palindromic kmer, avoid palindromes as seeds
			
            if (IRKE_COMMON::MONITOR >= 2) {
                cerr << "SEED kmer: " << kcounter.get_kmer_string(kmer) << " is palidnromic.  Skipping. " << endl;
            }
            
            continue;
		}
		
        
		if (kmer_count < MIN_SEED_COVERAGE) {
			if (IRKE_COMMON::MONITOR >= 2) {
                cerr << "-seed has insufficient coverage, skipping" << endl;
            }
            
            continue;
		}
		
		
		float entropy = compute_entropy(kmer, kmer_length);
		
		
		if (entropy < MIN_SEED_ENTROPY) {

            if (IRKE_COMMON::MONITOR >= 2) {
                cerr << "-skipping seed due to low entropy: " << entropy << endl;
            }
            
            continue;
		}
		
				
		/* Extend to the right */
		
		Kmer_visitor visitor(kmer_length, DOUBLE_STRANDED_MODE);
		Path_n_count_pair selected_path_n_pair_forward = inchworm(kcounter, 'F', kmer, visitor, min_connectivity); 
		
		visitor.clear();
		// add selected path to visitor
		
		vector<kmer_int_type_t>& forward_path = selected_path_n_pair_forward.first;

		if (IRKE_COMMON::MONITOR >= 2) {
            cerr << "Forward path contains: " << forward_path.size() << " kmers. " << endl;
        }


        for (unsigned int i = 0; i < forward_path.size(); i++) {
			kmer_int_type_t kmer = forward_path[i];
			visitor.add(kmer);
            
            if (IRKE_COMMON::MONITOR >= 2) {
                cerr << "\tForward path kmer: " << kcounter.get_kmer_string(kmer) << endl;
            }
            
		}
		
		
		/* Extend to the left */ 
		visitor.erase(kmer); // reset the seed
		
		Path_n_count_pair selected_path_n_pair_reverse = inchworm(kcounter, 'R', kmer, visitor, min_connectivity);
        if (IRKE_COMMON::MONITOR >= 2) {
            vector<kmer_int_type_t>& reverse_path = selected_path_n_pair_reverse.first;
            cerr << "Reverse path contains: " << reverse_path.size() << " kmers. " << endl;
            for (unsigned int i = 0; i < reverse_path.size(); i++) {
                cerr  << "\tReverse path kmer: " << kcounter.get_kmer_string(reverse_path[i]) << endl; 
            }
        }
        
		
		unsigned int total_counts = selected_path_n_pair_forward.second + selected_path_n_pair_reverse.second + kmer_count; //kcounter.get_kmer_count(kmer); 
		
		vector<kmer_int_type_t>& reverse_path = selected_path_n_pair_reverse.first;
		
		vector<kmer_int_type_t> joined_path = _join_forward_n_reverse_paths(reverse_path, kmer, forward_path);
		
		// report sequence reconstructed from path.
		
		vector<unsigned int> assembly_base_coverage;
		string sequence = reconstruct_path_sequence(kcounter, joined_path, assembly_base_coverage);
		
		int avg_cov =  static_cast<int> ( (float)total_counts/(sequence.length()-kcounter.get_kmer_length() +1) + 0.5);
		
		/*
		  cout << "Inchworm-reconstructed sequence, length: " << sequence.length() 
		  << ", avgCov: " << avg_cov
		  << " " << sequence << endl;
		*/
		
		
		
		if (sequence.length() >= MIN_ASSEMBLY_LENGTH && avg_cov >= MIN_ASSEMBLY_COVERAGE) {
			
			INCHWORM_ASSEMBLY_COUNTER++;
			
			stringstream headerstream;
			contigcounter++;
			
			
			headerstream << ">a" << INCHWORM_ASSEMBLY_COUNTER << ";" << avg_cov 
                         << " total_counts: " <<  total_counts << " Fpath: " << selected_path_n_pair_forward.second << " Rpath: " << selected_path_n_pair_reverse.second << " Seed: " << kmer_count 
						 << " K: " << kmer_length
						 << " length: " << sequence.length();
			
			string header = headerstream.str();
			
            sequence = add_fasta_seq_line_breaks(sequence, 60);
            
			cout << header << endl << sequence << endl;
			
			if (WRITE_COVERAGE) {
				
				coverage_writer << header << endl;
				
				for (unsigned int i = 0; i < assembly_base_coverage.size(); i++) {
					coverage_writer << assembly_base_coverage[i];
					if ( (i+1) % 30 == 0) {
						coverage_writer << endl;
					}
					else {
						coverage_writer << " ";
					}
				}
				coverage_writer << endl;
			}
			
		}
		
		// remove path
        
        if (IRKE_COMMON::__DEVEL_zero_kmer_on_use) { 
            
            // dont forget the seed. The forward/reverse path kmers already cleared.
            kcounter.clear_kmer(kmer);
            

        } else {
            
            for (unsigned int i = 0; i < joined_path.size(); i++) {
                
                kmer_int_type_t kmer = joined_path[i];
                
                /*
                  if (DEBUG) {
                  cout << "\tpruning kmer: " << kmer << endl;
                  }
                */
                
                //string kmer_seq = kcounter.get_kmer_string(kmer);
                //cerr << "Purging: " << kcounter.describe_kmer(kmer_seq) << endl;  
    
                kcounter.clear_kmer(kmer);
                

            }
                
        }

        /*
          if (DEBUG) {
		  cout << "done pruning kmers." << endl;
		  }
		*/
		
    }
	
	if (IRKE_COMMON::MONITOR) {
		cerr << endl;
	}
	
	if (WRITE_COVERAGE) {
		coverage_writer.close();
	}
        
    // drop sorted kmer list as part of cleanup
    clear_sorted_kmers_list();
    
	end = time(NULL);
	cerr << "Total contigs generated :" << contigcounter << endl;
	cerr << "Total assembly time by original implementation: " << end - start << " seconds" << endl;
	return; // end of runIRKE
	
}

Path_n_count_pair IRKE::inchworm (KmerCounter& kcounter, char direction, kmer_int_type_t kmer, Kmer_visitor& visitor, float min_connectivity) {
	
	// cout << "inchworm" << endl;
	
	Path_n_count_pair entire_path;
    entire_path.second = 0; // init cumulative path coverage

	unsigned int inchworm_round = 0;
	
	unsigned long num_total_kmers = kcounter.size();
	
	Kmer_visitor eliminator(kcounter.get_kmer_length(), DOUBLE_STRANDED_MODE);
	
	while (true) {


        if (IRKE_COMMON::__DEVEL_rand_fracture) {
            
            // terminate extension with probability of __DEVEL_rand_fracture_prob
            
            float prob_to_fracture = rand() / (float) RAND_MAX;
            //cerr << "prob: " << prob_to_fracture << endl;
            
            if (prob_to_fracture <= IRKE_COMMON::__DEVEL_rand_fracture_prob) {
                
                // cerr << "Fracturing at iworm round: " << inchworm_round << " given P: " << prob_to_fracture << endl;
                
                return(entire_path);
            }
        }
        		
		inchworm_round++;
		eliminator.clear();
		
		if (inchworm_round > num_total_kmers) {
			throw(string ("Error, inchworm rounds have exceeded the number of possible seed kmers"));
		}
		
		if (IRKE_COMMON::MONITOR >= 3) {
			cerr << endl << "Inchworm round(" << string(1,direction) << "): " << inchworm_round << " searching kmer: " << kmer << endl;
			string kmer_str = kcounter.get_kmer_string(kmer);
			cerr << kcounter.describe_kmer(kmer_str) << endl;
		}
		
		visitor.erase(kmer); // seed kmer must be not visited already.
		
		Kmer_Occurence_Pair kmer_pair(kmer, kcounter.get_kmer_count(kmer));
		Path_n_count_pair best_path = inchworm_step(kcounter, direction, kmer_pair, visitor, eliminator, inchworm_round, 0, min_connectivity, MAX_RECURSION);
		

        vector<kmer_int_type_t>& kmer_list = best_path.first;
        unsigned int num_kmers = kmer_list.size();

		if ( (IRKE_COMMON::__DEVEL_zero_kmer_on_use && num_kmers >= 1) || best_path.second > 0) {
			// append info to entire path in reverse order, so starts just after seed kmer
			
			int first_index = num_kmers - 1;
			int last_index = 0;
			if (CRAWL) {
				last_index = first_index - CRAWL_LENGTH + 1;
				if (last_index < 0) {
					last_index = 0;
				}
			}
			
			for (int i = first_index; i >= last_index; i--) {
				kmer_int_type_t kmer_extend = kmer_list[i];
				entire_path.first.push_back(kmer_extend);
				visitor.add(kmer_extend);
				//entire_path.second += kcounter.get_kmer_count(kmer_extend);

                // selected here, zero out:

                
                if (IRKE_COMMON::__DEVEL_zero_kmer_on_use) {
					/*
					#pragma omp critical(HashMap)
						if(direction == 'F')
							set_right_flag(&(kcounter.find_kmer(kmer_extend).second)
						else if(direction == 'R')
							set_left_flag(&(kcounter.find_kmer(kmer_extend).second))
					*/

					//Dont remove the kmer from the hash table

                    kcounter.clear_kmer(kmer_extend);
                }
                
			}
			
			kmer = entire_path.first[ entire_path.first.size() -1 ];
            
            entire_path.second += best_path.second;
            

		}
		else {
			// no extension possible
			break;
		}
	}
	
	if (IRKE_COMMON::MONITOR >= 3) 
		cerr << "No extension possible." << endl << endl;
	
	
	return(entire_path);
}


bool compare (const Path_n_count_pair& valA, const Path_n_count_pair& valB) {
	
#ifdef _DEBUG
    if (valA.second == valB.second)
        return (valA.first > valB.first);
	else
#endif
		return(valA.second > valB.second); // reverse sort.
}

Path_n_count_pair IRKE::inchworm_step (KmerCounter& kcounter, char direction, Kmer_Occurence_Pair kmer, Kmer_visitor& visitor,
									   Kmer_visitor& eliminator, unsigned int inchworm_round, unsigned int depth, 
									   float MIN_CONNECTIVITY_RATIO, unsigned int max_recurse) {
	
	// cout << "inchworm_step" << endl;
	
	if (IRKE_COMMON::MONITOR >= 2) {
		cerr << "\rinchworm: " << string(1,direction) 
			 << " A:" << INCHWORM_ASSEMBLY_COUNTER << " "
			 << " rnd:" << inchworm_round << " D:" << depth << "         "; 
	}
	
	// check to see if kmer exists.  If not, return empty container
	Path_n_count_pair best_path_n_pair;
	best_path_n_pair.second = 0; // init
		
	if ( // !kmer.second || 
        
        visitor.exists(kmer.first) // visited
        || eliminator.exists(kmer.first) // eliminated
		
		 ) {
        
        if (IRKE_COMMON::MONITOR >= 3) {
            cerr << "base case, already visited or kmer doesn't exist." << endl;
            cerr << kmer.first << " already visited or doesn't exist.  ending recursion at depth: " << depth << endl;
        }
		
        return(best_path_n_pair);
		
	}
	
	visitor.add(kmer.first);
	
	if (PACMAN && depth > 0) {
		// cerr << "pacman eliminated kmer: " << kmer << endl;
		eliminator.add(kmer.first);
	}
	
	
	if (depth < max_recurse) {
		
		vector<Kmer_Occurence_Pair> kmer_candidates;
		if (direction == 'F') {
			// forward search
			kmer_candidates = kcounter.get_forward_kmer_candidates(kmer.first);
		}
		else {
			// reverse search
			kmer_candidates = kcounter.get_reverse_kmer_candidates(kmer.first);
		}
		
        if (IRKE_COMMON::MONITOR >= 3) {
            cerr << "Got " << kmer_candidates.size() << " kmer extension candidates." << endl;
        }
        
		bool tie = true;
		unsigned int recurse_cap = max_recurse;
		unsigned int best_path_length = 0;
		while (tie) {

            // keep trying to break ties if ties encountered.
            // this is done by increasing the allowed recursion depth until the tie is broken.
            //  Recursion depth set via: recurse_cap and incremented if tie is found
            

			vector<Path_n_count_pair> paths;
			
			for (unsigned int i = 0; i < kmer_candidates.size(); i++) {
				Kmer_Occurence_Pair kmer_candidate = kmer_candidates[i];
				
				if ( kmer_candidate.second &&
				    
                    !visitor.exists(kmer_candidate.first)  // avoid creating already visited kmers since they're unvisited below...
					&& exceeds_min_connectivity(kcounter, kmer, kmer_candidate, MIN_CONNECTIVITY_RATIO) ) {
					//cout << endl << "\ttrying " << kmer_candidate << endl;
					

                    // recursive call here for extension
					Path_n_count_pair p = inchworm_step(kcounter, direction, kmer_candidate, visitor, eliminator, inchworm_round, depth+1, MIN_CONNECTIVITY_RATIO, recurse_cap);
					
                    if (p.first.size() >= 1) {
                        // only retain paths that include visited nodes.
                        paths.push_back(p);
                    }
					visitor.erase(kmer_candidate.first); // un-visiting
					
                }
				
			} // end for kmer
			
			
			if (paths.size() > 1) {

                sort(paths.begin(), paths.end(), compare);
				
                if (IRKE_COMMON::__DEVEL_no_greedy_extend) {
                    // pick a path at random
                    float p = rand()/(float)((long)RAND_MAX+1);
                    int rand_index = (int) (p * paths.size());
                    tie = false;
                    if (IRKE_COMMON::MONITOR) {
                        cerr << "IRKE_COMMON::__DEVEL_no_greedy_extend -- picking random path index: " << rand_index << " from p: " << p << " and size(): " << paths.size() << endl;
                    }
                    best_path_n_pair = paths[rand_index];
                }
				
                else if (paths[0].second == paths[1].second   // same cumulative coverage values for both paths.
					&&
					// check last kmer to be sure they're different. 
					// Not interested in breaking ties between identically scoring paths that end up at the same kmer.
					paths[0].first[0] != paths[1].first[0]
					
					) {
					
					// got tie, two different paths and two different endpoints:
					if (IRKE_COMMON::MONITOR >= 3) {
						
						
						cerr << "Got tie! " << ", score: " << paths[0].second << ", recurse at: " << recurse_cap << endl;
						vector<unsigned int> v;
						cerr << reconstruct_path_sequence(kcounter, paths[0].first, v) << endl;
						cerr << reconstruct_path_sequence(kcounter, paths[1].first, v) << endl;
												
					}
                    
                    if (IRKE_COMMON::__DEVEL_no_tie_breaking) {
                        tie = false;
                        float p = rand()/(float)((long)RAND_MAX+1);
                        int rand_index = (int) (p * 2); // just consider the first two options for simplicity
                        
                        if (IRKE_COMMON::MONITOR) {
                            cerr << "IRKE_COMMON::__DEVEL_no_tie_breaking, so picking path: " << rand_index << " at random." << endl;
                        }
                        
                        best_path_n_pair = paths[rand_index];
                    }
                    
                    else if (paths[0].first.size() > best_path_length) {
                        // still making progress in extending to try to break the tie.  Keep going.
                        // note, this is the only test that keeps us in this while loop. (tie stays true)
                        recurse_cap++;
						best_path_length = paths[0].first.size();
					}
					else {
						// cerr << "not able to delve further into the graph, though...  Stopping here." << endl;
						tie = false;
                        best_path_n_pair = paths[0]; // pick one
					}
				}
				
				else if ((paths[0].second == paths[1].second   // same cumulative coverage values for both paths.
						  &&
						  paths[0].first[0] == paths[1].first[0] ) // same endpoint
						 ) {
					
					if (IRKE_COMMON::MONITOR >= 3) {
						cerr << "Tied, but two different paths join to the same kmer.  Choosing first path arbitrarily." << endl;
					}
					tie = false;
					best_path_n_pair = paths[0];
				}
				
				else {
					// no tie.
					tie = false;
					best_path_n_pair = paths[0];
				}
				
				
			}
			else if (paths.size() == 1) {
				tie = false;
				best_path_n_pair = paths[0];
			}
			else {
				// no extensions possible.
				tie = false;
			}
			
			
		} // end while tie
	}
	
	// add current kmer to path, as long as not the original seed kmer!
	if (depth > 0) {
		best_path_n_pair.first.push_back(kmer.first);
		best_path_n_pair.second += kmer.second;
        
    }
	
	return(best_path_n_pair);
	
	
}


vector<kmer_int_type_t> IRKE::_join_forward_n_reverse_paths(vector<kmer_int_type_t>& reverse_path, 
															kmer_int_type_t seed_kmer_val, 
															vector<kmer_int_type_t>& forward_path) {
	
	vector<kmer_int_type_t> joined_path;
	
	// want reverse path in reverse order
	
	for (int i = reverse_path.size()-1; i >= 0; i--) {
		joined_path.push_back( reverse_path[i] );
	}
	
	// add seed kmer
	joined_path.push_back(seed_kmer_val);
	
	// tack on the entire forward path.
	
	for (unsigned int i = 0; i < forward_path.size(); i++) {
		joined_path.push_back( forward_path[i] );
	}
	
	return(joined_path);
}


string IRKE::reconstruct_path_sequence(vector<kmer_int_type_t>& path, vector<unsigned int>& cov_counter) {
	// use kcounter member
	return(reconstruct_path_sequence(kcounter, path, cov_counter));
}


string IRKE::reconstruct_path_sequence(KmerCounter& kcounter, vector<kmer_int_type_t>& path, vector<unsigned int>& cov_counter) {
	
	if (path.size() == 0) {
		return("");
	}
	
	string seq = kcounter.get_kmer_string(path[0]);
	cov_counter.push_back( kcounter.get_kmer_count(path[0]) );
	
	for (unsigned int i = 1; i < path.size(); i++) {
		string kmer = kcounter.get_kmer_string(path[i]);
		seq += kmer.substr(kmer.length()-1, 1);
		
		cov_counter.push_back( kcounter.get_kmer_count(path[i]) );
	}
	
	return(seq);
}


bool IRKE::exceeds_min_connectivity (KmerCounter& kcounter, Kmer_Occurence_Pair kmerA, Kmer_Occurence_Pair kmerB, float min_connectivity) {
	
    if (min_connectivity < 1e5) {
        return(true); // consider test off
    }
    
	
	unsigned int kmerA_count = kmerA.second;
	if (kmerA_count == 0)
		return(false);
	unsigned int kmerB_count = kmerB.second;
	if (kmerB_count == 0)
		return(false);
	
	unsigned int minVal;
	unsigned int maxVal;
	
	if (kmerA_count < kmerB_count) {
		minVal = kmerA_count;
		maxVal = kmerB_count;
	}
	else {
		minVal = kmerB_count;
		maxVal = kmerA_count;
	}
	
	float connectivity_ratio = (float) minVal/maxVal;
	
	if (connectivity_ratio >= min_connectivity) {
		return(true);
	}
	else {
		return(false);
	}
}


bool IRKE::exceeds_min_connectivity (KmerCounter& kcounter, string kmerA, string kmerB, float min_connectivity) {

	kmer_int_type_t valA = kmer_to_intval(kmerA);
	kmer_int_type_t valB = kmer_to_intval(kmerB);

	Kmer_Occurence_Pair pairA(valA, kcounter.get_kmer_count(valA));
	Kmer_Occurence_Pair pairB(valB, kcounter.get_kmer_count(valB));

	return exceeds_min_connectivity(kcounter, pairA, pairB, min_connectivity);

}



string IRKE::thread_sequence_through_graph(string& sequence) {
	
	// describe each of the ordered kmers in the input sequence as they exist in the graph.
	
	unsigned int kmer_length = kcounter.get_kmer_length();
	
	if (sequence.length() < kmer_length) {
		cerr << "Sequence length: " << sequence.length() << " is too short to contain any kmers." << endl;
		return("");
	}
	
	stringstream s;
	
	for (unsigned int i=0; i <= sequence.length() - kmer_length; i++) {
		
		string kmer = sequence.substr(i, kmer_length);
		
		s << kcounter.describe_kmer(kmer) << endl;
	}
	
	return(s.str());
	
	
}



bool IRKE::sequence_path_exists(string& sequence, unsigned int min_coverage, float min_entropy, float min_connectivity, 
								vector<unsigned int>& coverage_counter) {
	
	unsigned int kmer_length = kcounter.get_kmer_length();
	
	if (sequence.length() < kmer_length) {
		return(false);
	}
	
	bool path_exists = true;
	
	string prev_kmer = sequence.substr(0, kmer_length);
	if (contains_non_gatc(prev_kmer) || ! kcounter.kmer_exists(prev_kmer)) {
		path_exists = false;
		coverage_counter.push_back(0);
	}
	else {
		unsigned int kmer_count = kcounter.get_kmer_count(prev_kmer);
		coverage_counter.push_back(kmer_count);

		float entropy = compute_entropy(prev_kmer);

		if (kmer_count < min_coverage || entropy < min_entropy) {
			path_exists = false;
		}
	}
	
	
	for (unsigned int i=1; i <= sequence.length() - kmer_length; i++) {
		
		string kmer = sequence.substr(i, kmer_length);
		
		if (contains_non_gatc(kmer) || ! kcounter.kmer_exists(kmer)) {
			path_exists = false;
			coverage_counter.push_back(0);
		}
		else {
			unsigned int kmer_count = kcounter.get_kmer_count(kmer);
			coverage_counter.push_back(kmer_count);

			float entropy = compute_entropy(kmer);

			if (kmer_count < min_coverage || entropy < min_entropy) {
				path_exists = false;
			}
		}

		
		if (path_exists && ! exceeds_min_connectivity(kcounter, prev_kmer, kmer, min_connectivity)) {
			path_exists = false;
		}
		
		prev_kmer = kmer;
		
	}
	
	return(path_exists);
}



void IRKE::describe_kmers() {
	
	return(describe_kmers(kcounter));
}

void IRKE::describe_kmers(KmerCounter& kcounter) {
	
	// proxy call to Kcounter method.
	
	kcounter.describe_kmers();
	
	return;
}

void IRKE::populate_sorted_kmers_list() {
    
    sorted_kmers = kcounter.get_kmers_sort_descending_counts();
    got_sorted_kmers_flag = true;

    return;
}

void IRKE::clear_sorted_kmers_list() {

    sorted_kmers.clear();
    got_sorted_kmers_flag = false;
    
    return;
}
