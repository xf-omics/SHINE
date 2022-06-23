/*
 * first input a3m with original species from genetree
 * second input number of original species
 * trim the phylogenetic trees to 300 species
 */

#include "/share/terra/Users/xf2193/resource/CPP/class_list.h"

int main(int argc, char *argv[]) {
	//arguments check
	if (argc<3) {
		cout << "There is no input file basename or number of desirable species!\n";
		return 1;
	}

	int th;
	unsigned int flag=0, depth=1024, level=2;
//	float len[2];
	size_t found;
	string temp, path, key, barcode[2];
	ifstream input;
	ofstream output,list;
	
	th = atoi(argv[2]);

        path = "/share/terra/Users/xf2193/resource/ensembl/a3m_trim/";
        path += argv[1];        // ***.a3m
        path += ".a3m";	
	output.open(path.c_str());
	
	path = "/share/terra/Users/xf2193/resource/ensembl/ensembl_genetree_gene_protein_depth.txt";
	list.open(path.c_str(),ios_base::app);

	if (list.is_open()==0||output.is_open()==0) {
		cout << "input or output file is not open\n";
		return 1;
	}

	while (depth>th) {
		path = "/share/terra/Users/xf2193/resource/ensembl/a3m/";
	        path += argv[1];        // ***.txt
        	path += ".a3m";
	        input.open(path.c_str());
		
		if (input.is_open()==0) {
			cout << path << " is not open\n";
			break;
		}
		
		level++;
		flag = 0;
		depth = 0;

		while (input.good()&&!input.eof()) {
			getline(input,temp);
			if (temp=="\0") continue;
		
			found = temp.find_last_of(' ');	// barcode
        		if (found==string::npos) {
				cout << "wrong format!\n";
            			cout << temp << endl;
	        		break;
        		}

                	if (flag==0) {
				barcode[0] = temp.substr(found+1,level);
                	        flag = 1;
        	                getline(input,temp);
				continue;
	                }
			else {
				barcode[1] = temp.substr(found+1);
				if (barcode[1].length()>=level&&barcode[1].substr(0,level)==barcode[0])
					depth++;
			}
			getline(input,temp);
		}
		input.close();
	}

	path = "/share/terra/Users/xf2193/resource/ensembl/a3m/";
        path += argv[1];        // ***.txt
        path += ".a3m";
        input.open(path.c_str());

        if (input.is_open()==0) {
               cout << path << " is not open\n";
               return 1;
        }

	while (input.good()&&!input.eof()) {
        	getline(input,temp);
        	if (temp=="\0") continue;
		getline(input,path);

        	found = temp.find_last_of(' '); // barcode
        	if (found==string::npos) {
                       cout << "wrong format!\n";
		       cout << temp << endl;
                       break;
        	}
		
		barcode[1] = temp.substr(found+1);
                if (barcode[1].length()>=level&&barcode[1].substr(0,level)==barcode[0])
                	output << temp << endl << path << endl;
	}

	list << argv[1] << "\t" << depth << endl;

	list.close();
	output.close();
	return 0;
}

