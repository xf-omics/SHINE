/*
 * input includes common_name, scientific_name, gene accession, protein accession, sequence and barcode (level of tree)
 * 
 */

#include "/share/terra/Users/xf2193/resource/CPP/class_list.h"

int main(int argc, char *argv[]) {
	//arguments check
	if (argc<2) {
		cout << "There is no input file basename!\n";
		return 1;
	}

	unsigned int i=0, depth=0;
//	float len[2];
	size_t found;
	string temp, path, key, common, scientific, gene, protein, seq, barcode;
	ifstream input;
	ofstream output[2],list;
	
	path = "/share/terra/Users/xf2193/resource/ensembl/msa/";
	path += argv[1];	// ***.txt
	path += ".txt";
	input.open(path.c_str());

        path = "/share/terra/Users/xf2193/resource/ensembl/a3m/";
        path += argv[1];        // ***.a3m
        path += ".a3m";	
	output[0].open(path.c_str());
	
	path = "/share/terra/Users/xf2193/resource/ensembl/script/temp";
        output[1].open(path.c_str());

	path = "/share/terra/Users/xf2193/resource/ensembl/ensembl_genetree_gene_protein_id.txt";
	list.open(path.c_str(),ios_base::app);

	if (input.is_open()==0||list.is_open()==0||output[0].is_open()==0||output[1].is_open()==0) {
		cout << "input or output file is not open\n";
		return 1;
	}
	
	
	while (input.good()&&!input.eof()) {
		getline(input,temp);
		if (temp=="\0") continue;
		
		found = temp.find(' ');	// key
        	if (found==string::npos) {
			cout << "wrong format!\n";
            		cout << temp << endl;
	        	break;
        	}

		i++;
        	key = temp.substr(0,found);
        	temp = temp.substr(found+1);
                found = temp.find_last_of(' '); // content
                if (found==string::npos) {
                        cout << "wrong format!\n";
                        cout << temp << endl;
                        break;
                }
		
		if (key == "common_name") {
			common = temp.substr(0,found);
			barcode = temp.substr(found+1);
		}
		else {
			temp = temp.substr(0,found);
			if (key == "scientific_name")
				scientific = temp;
			else if (key == "seq") {
				seq = temp;
				found = seq.find('*');
				while (found!=string::npos) {
					seq[found] = '-';
					found = seq.find('*',found);
				}
		     	}
			else if (key == "accession") {
				found = temp.find("P0");
				if (found!=string::npos)
					protein = temp;
				else {
					found = temp.find("G00");
					if (found!=string::npos)
						gene = temp;
					else if (gene=="")
						gene = temp;
					else
						protein = temp;
				}
			}
			else cout << "wrong key: " << key << endl;
		}

		if (i==5) {
			if (scientific==""||common==""||gene==""||protein==""||seq=="") {
				cout << "One of keys is missing!\n";
				cout << argv[1] << ": " << scientific << ", " << common << ", " << gene << ", " << protein << endl;
				break;
			}
//			if (common=="Human") { // there may be more than 1 human sequences
			if (gene==argv[1]) {
				output[0] << ">" << scientific << ", " << common << ", " << gene << ", " << protein << ", " << barcode << endl << seq << endl;
				path = gene;
				path += "\t";
				path += protein;
				output[0].close();
			} else {
//				len[1] = 0;
//				len[0] = seq.length();
//				for (i=0;i<len[0];i++) {
//					if (seq[i]!='-')
//					len[1]++;
//				}
//				if (len[1]/len[0]>0.1) {
					output[1] << ">" << scientific << ", " << common << ", " << gene << ", " << protein << ", " << barcode << endl << seq << endl;
					depth++;
//				}
			}
			scientific = "";
			common = "";
			gene = "";
			protein = "";
			seq = "";
			i = 0;
		}
	}
	
	list << path << "\t" << depth << endl;

	list.close();
	input.close();
	output[1].close();
	return 0;
}

