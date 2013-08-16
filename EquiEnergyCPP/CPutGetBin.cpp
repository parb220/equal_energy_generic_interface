#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <glob.h>
#include <algorithm>
#include "dw_rand.h"

#include "CSampleIDWeight.h"
#include "CPutGetBin.h"
string CPutGetBin::GetFileNameForDump() const
{
	stringstream convert;
        convert << id << "." << nDumpFile << "." << suffix;
        string file_name = filename_prefix + convert.str();
	return file_name; 
}

vector <string> CPutGetBin::GetFileNameForFetch() const
{
	stringstream convert;
        convert.str(string());
        convert << id << ".*.*"; // << ".*" // all cluster node; 

        string filename_pattern = filename_prefix + convert.str();

        glob_t glob_result;
        glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	vector <string> filename_fetch; 
	if (glob_result.gl_pathc > 0)
	{
		for (unsigned int i=0; i<glob_result.gl_pathc; i++)
			if (NumberRecord(string(glob_result.gl_pathv[i])) >= capacity)
				filename_fetch.push_back(string(glob_result.gl_pathv[i])); 
	}
	else 
		filename_fetch.clear(); 
        globfree(&glob_result);
        return filename_fetch; 
}

vector <string> CPutGetBin::GetFileNameForConsolidate() const
{
	stringstream convert;
        convert.str(string());
        convert << id << ".*.*"; // << ".*" // all cluster node; 

        string filename_pattern = filename_prefix + convert.str();

        glob_t glob_result;
        glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	vector <string> filename_consolidate; 
	if (glob_result.gl_pathc > 0)
	{
		for (unsigned int i=0; i<glob_result.gl_pathc; i++)
			if (NumberRecord(string(glob_result.gl_pathv[i])) < capacity)
				filename_consolidate.push_back(string(glob_result.gl_pathv[i])); 
	}
	else 
		filename_consolidate.clear(); 
        globfree(&glob_result);
        return filename_consolidate; 
}

bool CPutGetBin::Dump(const string &_filename)
{
	string file_name; 
	if (_filename.empty())
		file_name = GetFileNameForDump(); 
	else 
		file_name = _filename; 
	fstream oFile(file_name.c_str(), ios::out | ios::binary); 
	if (!oFile)
	{
		cerr << "Error in dumping samples to " << file_name << endl; 
		return false; 
	}
	for (unsigned int i=0; i<nPutUsed; i++)
		write(oFile, &(dataPut[i])); 
	oFile.flush(); 
	oFile.close();  
	return true; 
}

bool CPutGetBin::Fetch(const vector<string> &filename_fetch)
{
	if (filename_fetch.empty() && capacity > nPutUsed)
	{
		// cerr << "Error in Fetch(): files are empty and bins are not full yet.\n"; 
		return false;
	}

	size_t nFetchFile = filename_fetch.size(); 
	if (dataGet.size() < capacity)
		dataGet.resize(capacity);  
	
	unsigned int select; 
	vector < vector <unsigned int > > select_per_file(nFetchFile+1); 
	// First nFetchFile(): files
	// Last: cache
	for (unsigned int i=0; i<capacity; i++)
	{
		select=dw_uniform_int(nFetchFile*capacity+nPutUsed); 
		select_per_file[select/capacity].push_back(select%capacity); 
	}
	int counter =0; 
	for (unsigned int i=0; i<nFetchFile; i++)	// Read data from file
	{
		if (!select_per_file[i].empty())
		{
			if (!ReadFromOneFile(filename_fetch[i], counter, select_per_file[i]) )
			{
				cerr << "Error in reading data from " << filename_fetch[i] << endl; 
				return false; 
			}
		}
	}
	if (!select_per_file[nFetchFile].empty() ) // Read data from cache
	{
		for (unsigned int n=0; n<select_per_file[nFetchFile].size(); n++)
		{
			dataGet[counter] = dataPut[select_per_file[nFetchFile][n]]; 
			counter ++;
		}
	}
	return true; 
}

bool CPutGetBin::ReadFromOneFile(const string &file_name, int &counter, const vector <unsigned int> &index) 
{	
	fstream iFile; 
        iFile.open(file_name.c_str(), ios::in|ios::binary);
        if (!iFile)
        	return false;
	// determine the dimension of CSampleIDWeight
	CSampleIDWeight temp_data; 
	read(iFile, &temp_data); 

        for (unsigned int n=0; n<index.size(); n++)
        {
        	iFile.seekg(temp_data.GetSize_Data()*index[n], ios_base::beg);
                read(iFile, &(dataGet[counter]));
                counter ++;
        }
      	iFile.close();
	return true; 
}

vector <CSampleIDWeight> CPutGetBin::ReadSampleFromFile(const string &file_name) const
{
	unsigned int nRecord = NumberRecord(file_name); 
	if (nRecord <= 0)
		return vector<CSampleIDWeight>(0);

	fstream iFile;
        iFile.open(file_name.c_str(), ios::in|ios::binary);
        if (!iFile)
        	return vector<CSampleIDWeight>(0); 

        vector <CSampleIDWeight> sample(nRecord); 
        for(unsigned int n=0; n<nRecord; n++)
		read(iFile, &(sample[n]));
        iFile.close();
	return sample; 
}

unsigned int CPutGetBin::NumberRecord(const string &file_name) const
{
	ifstream iFile; 
	iFile.open(file_name.c_str(), ios::in|ios::binary);
        if (!iFile)
		return 0;
	
	// To determine size of each record
	CSampleIDWeight temp; 
	read(iFile, &temp); 

	iFile.seekg(0, ios::beg); 
	iFile.seekg(0, ios::end); 
	unsigned int lenFile = iFile.tellg(); 
	unsigned int number_record = lenFile/temp.GetSize_Data(); 
	iFile.close(); 
	return number_record; 
}

bool compare_CSampleIDWeight(const CSampleIDWeight &i, const CSampleIDWeight &j) 
{
	return i.weight < j.weight; 
}

bool CPutGetBin::Load_K_MostWeightSample(size_t K, const string &file_name, vector<CSampleIDWeight> &best_samples) const
{
	vector <CSampleIDWeight> samples = ReadSampleFromFile(file_name); 
	if (samples.empty())
		return false; 
	samples.insert(samples.end(), best_samples.begin(), best_samples.end()); // merge sampels and best_samples together
	sort(samples.begin(), samples.end(), compare_CSampleIDWeight); // ascending based on weight
	if (samples.size() <=K)
		best_samples = samples; 
	else 
	{
		best_samples.resize(K); 
		copy(samples.end()-K, samples.end(), best_samples.begin()); 
	}
	reverse(best_samples.begin(), best_samples.end()); 
	return true; 
}

bool CPutGetBin::Load_K_LeastWeightSample(size_t K, const string &file_name, vector<CSampleIDWeight> &best_samples) const
{
	vector <CSampleIDWeight > samples = ReadSampleFromFile(file_name); 
	if (samples.empty()) 
		return false; 
	samples.insert(samples.end(), best_samples.begin(), best_samples.end()); // merge sampels and best_samples together
	sort(samples.begin(), samples.end(), compare_CSampleIDWeight); // ascending based on weight
	if (samples.size() <= K)
		best_samples = samples; 
	else 
	{
		best_samples.resize(K); 
		copy(samples.begin(), samples.begin()+K, best_samples.begin()); 
	}
	return true; 
}

bool CPutGetBin::LoadLeastWeightSample(const string &file_name, CSampleIDWeight &best) const
{
	vector <CSampleIDWeight > samples = ReadSampleFromFile(file_name); 
	if ( samples.empty() )
		return false; 
	best = samples[0]; 
	for (unsigned int i=1; i<samples.size(); i++)
	{
		if (samples[i].weight < best.weight )
			best = samples[i]; 
	}
	return true; 
}

bool CPutGetBin::LoadMostWeightSample(const string &file_name, CSampleIDWeight &best) const
{
	vector <CSampleIDWeight> samples = ReadSampleFromFile(file_name); 
	if ( samples.empty() )
		return false; 
	best = samples[0]; 
	for (unsigned int i=1; i<samples.size(); i++)
	{
		if (samples[i].weight > best.weight )
			best = samples[i]; 
	}
	return true; 
}

CPutGetBin::CPutGetBin(int _id, unsigned int _nDumpFile, size_t _capacity, string _grandPrefix, int _suffix ) : 
suffix(_suffix), id(_id),  nDumpFile(_nDumpFile), capacity(_capacity), 
nPutUsed(0), nGetUsed(_capacity), filename_prefix(_grandPrefix), 
dataPut(vector<CSampleIDWeight>(0)), dataGet(vector<CSampleIDWeight>(0)) 
{
}

CPutGetBin::CPutGetBin(const CPutGetBin &right) : 
suffix(right.suffix), id(right.id), nDumpFile(right.nDumpFile), capacity(right.capacity),
nPutUsed(right.nPutUsed), nGetUsed(right.nGetUsed), filename_prefix(right.filename_prefix),
dataPut(right.dataPut), dataGet(right.dataGet)
{
}
 
CPutGetBin & CPutGetBin::operator=(const CPutGetBin &right)
{
	suffix = right.suffix; 
	id = right.id; 
	nDumpFile = right.nDumpFile; 
	capacity = right.capacity; 
	nPutUsed = right.nPutUsed; 
	nGetUsed = right.nGetUsed; 
	filename_prefix = right.filename_prefix; 
	dataPut = right.dataPut; 
	dataGet = right.dataGet; 
	return *this; 
}

CPutGetBin::~CPutGetBin()
{
}

size_t CPutGetBin::GetNumberFileForFetch() const
{
	vector <string> file_fetch = GetFileNameForFetch();
	return file_fetch.size(); 
}


size_t CPutGetBin::GetNumberFileForDump() const
{
	stringstream convert; 
	convert.str(string()); 
	convert << id << ".*." << suffix; // << "." << cluster_node; 

	string filename_pattern = filename_prefix + convert.str(); 
	
	glob_t glob_result; 
	glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result); 
	size_t final_result = glob_result.gl_pathc; 
	globfree(&glob_result); 
	return final_result; 
}

size_t CPutGetBin::GetNumberFileForConsolidate() const
{
	vector <string> file_consolidate = GetFileNameForConsolidate(); 
	return file_consolidate.size(); 
}


unsigned int CPutGetBin::DepositSample(const CSampleIDWeight &sample)
{
	unsigned int index =  nPutUsed;
	if (dataPut.size() <= nPutUsed)
		dataPut.push_back(sample); 
	else
		dataPut[index] = sample; 
	nPutUsed ++; 
	if (nPutUsed == capacity)
	{
		Dump(); 
		nDumpFile ++; 
		nPutUsed = 0; 
	}
	
	return nDumpFile*capacity+nPutUsed; 
}

bool CPutGetBin::Draw_K_MostWeightSample(size_t K, vector<CSampleIDWeight> &sample) const
{
	vector <string> file = GetFileNameForConsolidate();
        vector <string> fetch_file = GetFileNameForFetch();
        file.insert(file.end(), fetch_file.begin(), fetch_file.end());

	if (file.empty())
		return false; 
        
        for (unsigned int i=0; i<file.size(); i++)
        {
                if (!Load_K_MostWeightSample(K, file[i], sample) )
                {
                        cerr << "Draw_K_MostWeightSample : Error in opening " << file[i] << " for reading data.\n";
                        return false;
                }
        }
        return true;

}

bool CPutGetBin::Draw_K_LeastWeightSample(size_t K, vector<CSampleIDWeight> &sample) const
{
	vector <string> file = GetFileNameForConsolidate();
        vector <string> fetch_file = GetFileNameForFetch();
	file.insert(file.end(), fetch_file.begin(), fetch_file.end());

	if (file.empty())
		return false; 	

	for (unsigned int i=0; i<file.size(); i++)
	{
		if (!Load_K_LeastWeightSample(K, file[i], sample) )
		{
			cerr << "Draw_K_LeastWeightSample : Error in opening " << file[i] << " for reading data.\n"; 
			return false; 
		}
	}
	return true; 
}

bool CPutGetBin::DrawLeastWeightSample(CSampleIDWeight &sample)  const
{
	CSampleIDWeight temp_sample; 

	vector <string> file = GetFileNameForConsolidate();
	vector <string> fetch_file = GetFileNameForFetch(); 
	file.insert(file.end(), fetch_file.begin(), fetch_file.end()); 
	if (file.empty())
		return false; 	 

	for (unsigned int i=0; i<file.size(); i++)
	{
		if (!LoadLeastWeightSample(file[i], temp_sample)) 
		{
			cerr << "DrawLeastWeightSample : Error in opening " << file[i] << " for reading data.\n"; 
			return false; 
		}
		if (i==0 || temp_sample.weight < sample.weight)
			sample = temp_sample; 
	}
	return true; 	
}

bool CPutGetBin::DrawMostWeightSample(CSampleIDWeight &sample)  const
{
        CSampleIDWeight temp_sample;

        vector <string> file = GetFileNameForConsolidate();
        vector <string> fetch_file = GetFileNameForFetch();
        file.insert(file.end(), fetch_file.begin(), fetch_file.end());
	
	if (file.empty())
		return false; 

        for (unsigned int i=0; i<file.size(); i++)
        {
                if (!LoadMostWeightSample(file[i], temp_sample))
                {
                        cerr << "Error in opening " << file[i] << " for reading data.\n";
                        return false;
                }
                if (i==0 || temp_sample.weight > sample.weight)
                        sample = temp_sample;
        }
        return true;
}

bool CPutGetBin::DrawSample(CSampleIDWeight &sample)
{
	unsigned int index; 
	if (nGetUsed < capacity)	// if data in the memory has not been exhaustively used 
	{
		index = dw_uniform_int(capacity);
                nGetUsed ++;
                sample = dataGet[index];
		return true; 
	}
	else	// replenish and get data
	{
		vector <string> filename_fetch = GetFileNameForFetch(); 
		size_t nFetchFile = filename_fetch.size(); 
		if (nFetchFile == 0 && nPutUsed == 0 )
		{
			// cerr << "Error in DrawSample() : files and bins are empty.\n"; 
			return false; 
		}
		else if (nFetchFile == 0 ) // && nPutUsed > 0) 
		/* when data have not been dumped to files,
 		will directly get a data from dataPut
 		*/
		{
			/*
 			nGetUsed is initialized as capacity and is reset to 0 after each Fetch. 
			This part is executated when there is no files for fetch (and therefore 
			nGetUsed cannot be 0). So we let nGetUsed continue to increase until the 
			net increment, nGetUsed - capacity exceeds the number of samples in 
			nPutUsed.
 			*/
			/*if (nGetUsed >= capacity + nPutUsed)
			{
				cerr << "Error in DrawSample(): files are empty while samples in bins have been exhausitively used.\n"; 
				return false; 
			}*/

			/*
 			If DrawSample is attempted even after all data in dataPut have been exhausted,
			we will let dataPut continue to be used
 			*/
			index = dw_uniform_int(nPutUsed); 
			sample = dataPut[index];
			nGetUsed ++;  
			return true; 
		}
		else	// nFetchFile > 0 and/or nPutUsed > 0, get data from files
		{
			Fetch(filename_fetch);
                        nGetUsed = 0;
			index = dw_uniform_int(capacity);
                	nGetUsed ++;
                	sample = dataGet[index];
                	return true;
		}
	}
}

void CPutGetBin::finalize()
{
	if (nPutUsed > 0)
	{
		Dump();
		nDumpFile++;
		nPutUsed =0; 
	}  
}

void CPutGetBin::consolidate()
{
	vector < string> file_consolidate = GetFileNameForConsolidate(); 
	if (file_consolidate.size() > 1)
	{
		vector <CSampleIDWeight> sample_consolidate; 
		vector <CSampleIDWeight> temp_sample; 
		for (unsigned int i=0; i<file_consolidate.size(); i++)
		{
			temp_sample = ReadSampleFromFile(file_consolidate[i]); 
			sample_consolidate.insert(sample_consolidate.end(), temp_sample.begin(), temp_sample.end()); 
			remove(file_consolidate[i].c_str()); 
		}
		size_t nComplete = sample_consolidate.size()/capacity; 
		size_t nRemaining = sample_consolidate.size()%capacity; 
		for (unsigned int iComplete=0; iComplete<nComplete; iComplete ++)
		{
			dataPut.resize(capacity); 
			for (unsigned int j=0; j<capacity; j++)
				dataPut[j] = sample_consolidate[j+iComplete*capacity]; 
			nPutUsed = capacity; 
			Dump(file_consolidate[iComplete]); 
		}
		if (nRemaining > 0)
		{
			if (dataPut.size() < nRemaining)
				dataPut.resize(nRemaining); 
			for (unsigned int j=0; j<nRemaining; j++)
				dataPut[j] = sample_consolidate[j+nComplete*capacity]; 
			nPutUsed = nRemaining; 
			Dump(file_consolidate[nComplete]); 
		}
	}
}

void CPutGetBin::restore()
{
	nDumpFile = GetNumberFileForDump(); 
	if (nDumpFile > 0)
	{
       	 	string file_name;
        	stringstream convert;

        	convert.str(string());
        	convert << id << "." << nDumpFile-1 << "." << suffix; 
        	file_name = filename_prefix + convert.str();
		dataPut = ReadSampleFromFile(file_name); 
		nPutUsed = dataPut.size(); 
 
		if (nPutUsed >0 && nPutUsed < capacity)
			nDumpFile --; 
		else if (nPutUsed >= capacity)
			nPutUsed =0; 
	}
}

void CPutGetBin::RestoreForFetch()
{
	vector <string> filename_partial = GetFileNameForConsolidate(); 
	// After consolidation, there should be at most 1 partial file
	if (!filename_partial.empty())
	{
		vector<CSampleIDWeight> tempSample = ReadSampleFromFile(filename_partial[0]); 
		if (!tempSample.empty())
		{
			nPutUsed = tempSample.size(); 
			if (dataPut.size() < tempSample.size())
				dataPut.resize(capacity); 
			for (unsigned int j=0; j<tempSample.size(); j++)
				dataPut[j] = tempSample[j]; 
		}
	}
}
 
bool CPutGetBin::empty() const
{
        stringstream convert;
        convert.str(string());
        convert << id << ".*.*"; // << ".*" // all cluster node; 

        string filename_pattern = filename_prefix + convert.str();

        glob_t glob_result;
        glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	bool return_result; 
        if (glob_result.gl_pathc > 0)
		return_result = false; 
	else 
		return_result = true; 

        globfree(&glob_result);
	return return_result; 
}

void CPutGetBin::DisregardHistorySamples()
{
	vector <string> filename_fetch = GetFileNameForFetch(); 
	size_t nFetchFile = filename_fetch.size(); 
		
	for (unsigned int iFile =0; iFile<nFetchFile; iFile++)
		remove(filename_fetch[iFile].c_str());
}

void CPutGetBin::ClearDepositDrawHistory()
{
	nDumpFile = 0; 
 	nGetUsed = capacity;
	nPutUsed = 0; 
}

size_t CPutGetBin::GetTotalNumberRecord() const
{
	stringstream convert;
        convert.str(string());
        convert << id << ".*.*"; 
	string filename_pattern = filename_prefix + convert.str();

        glob_t glob_result;
        glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	int nRecord=0; 
        for (unsigned int i=0; i<glob_result.gl_pathc; i++)
        	nRecord += NumberRecord(string(glob_result.gl_pathv[i]));
	return nRecord; 
}



