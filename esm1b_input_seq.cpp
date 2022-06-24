/*
  parsing the protein sequences to format input for ESM-1b
  divide sequences into segments if len>=1022
  ESM-1b adds an tok in the front and at the end of the sequence, max limit is 1024
  output fasta
*/
#include "class_list.h"

int main(int argc, char *argv[]) {
	//arguments check
	if (argc<2) {
		cout << "Please provide input file name!\n";
		return 1;
	}

	size_t found;
	unsigned int seg, len, pos[2];
	unsigned long long i;
	string temp, path, file, seq;
	ifstream input, search;
	ofstream output;
	
	path = argv[1]; // input is a list of transcript IDs, this is for deletions, so wild type sequences
        found = path.find_last_of('.');
        if (found==string::npos) {
                temp = path;
                path += ".txt";
                temp += ".fasta";
        }
        else {
                temp = path.substr(0,found);
                temp += ".fasta";
        }

	input.open(path.c_str());
        if (input.is_open()==0) {
                cout << path << " variant file is not open\n";
                return 1;
        }
	
	output.open(temp.c_str());
	if (output.is_open()==0) {
		cout << temp << " sequence file is not open\n";
		return 1;
	}

	path = "~/resource/ensembl/MANE_ensembl_protein/";

	while (input.good()&&!input.eof()) {
		getline(input,temp); // transcript
		if (temp=="\0") continue;
	
		file = path;
		file += temp;
		file += ".fasta";
		search.open(file.c_str());
		if (search.is_open()==0) {
			cout << file << " is not open\n";
			continue;
		}

		getline(search,seq); // header
                getline(search,seq); // sequence
		search.close();

		len = seq.length();
		if (len<1023)
			output << '>' << temp << endl << seq << endl;
		else {
			seg = len / 1022; // number of segments - 1 (always round down)
			pos[0] = 0;
			for (i=0;i<=seg;i++) {
				pos[1] = len * (i+1) / (seg+1);
				output << '>' << temp << '_' << i << endl;
				output << seq.substr(pos[0],pos[1]-pos[0]) << endl;
				pos[0] = pos[1];
			}
		}
	}

	input.close();
	output.close();
	return 0;
}

