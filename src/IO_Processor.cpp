# include "IO_Processor.hpp"
int KMER_SIZE = 25;
void IO_Processor::initialize(char* mapping, long length, unsigned long long bound, unsigned long long &log)
{
	this->map = mapping;
	this->length = length;
	this->bound_kmer = bound;
	this->generate_kmer_buffer(log);
}

void IO_Processor::generate_kmer_buffer(unsigned long long &log) {
	long i = 0, j, counter = 0;
	kmer_pair temp;
	string ripped_str("");

	while (this->length > 0 && this->map[0] != '>') {
		this->map++;
		this->length--;
	}

	while (i < this->length) {
		
		if(kmer_buffer.size() >= 1)
		  log = this->kmer_buffer[0].first;
	  
		if (this->map[i] == '>') {

			i++;
			j = 0;
			while (j != KMER_SIZE && i < this->length) {
				ripped_str += this->map[i];
				if (this->map[i] == 'A' || this->map[i] == 'T'
						|| this->map[i] == 'G' || this->map[i] == 'C')
					j++;
				i++;
			}
			if (j == KMER_SIZE) {
				temp = extract_kmer(ripped_str);
				if(temp.first == bound_kmer)
					break;
				if(temp.second >= 1)
					this->kmer_buffer.push_back(temp);
				else 
					counter++;
				ripped_str.clear();
			}
		} else
			i++;
	}

	//cerr << "Parsed :" << this->kmer_buffer.size() << "kmers, Pruned:"  << counter << "kmers" <<endl; 

}

