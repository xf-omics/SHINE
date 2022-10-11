// divide MSA into pieces if length > 1023

#include "class_list.h"

int main(int argc, char *argv[]) {
	//arguments check
	if (argc<2) {
		cout << "There is no input file basename!\n";
		return 1;
	}

	unsigned int i, j, flag=1, seg=0, pos[4];
	size_t found[2];
	float len;
	string temp, path, header, seq;
	ifstream input;
	ofstream output,list;
	
	path = "~/resource/ensembl/a3m_trim/";
	path += argv[1];	// ***.txt
	path += ".a3m";
	input.open(path.c_str());

	path = "~/resource/ensembl/ensembl_genetree_trim_gene_protein_index.txt";
	list.open(path.c_str(),ios_base::app);

        path = "~/resource/ensembl/a3m_trim_seg/";

	if (input.is_open()==0||list.is_open()==0) {
		cout << "input or output file is not open\n";
		return 1;
	}
	
	
	while (input.good()&&!input.eof()) {
		getline(input,header);
		if (header=="\0") continue;
                getline(input,seq);
                if (seq=="\0") {
			cout << "wrong format!\n";
            		cout << header << endl;
	        	break;
        	}

		if (flag) { // human
			found[1] = header.find(',');
			found[0] = header.find(',',found[1]+1);
                        found[1] = header.find(',',found[0]+1);
			list << header.substr(found[0]+2,found[1]-found[0]-2); // gene stable id
                        found[0] = header.find(',',found[1]+1);
			list << "\t" << header.substr(found[1]+2,found[0]-found[1]-2); // protein stable id

			pos[1] = seq.length();

			seg = ceil((pos[1])/1023.0);
			len = float(pos[1]) / seg; // length of each segment
			
			flag = 0;
	                for (i=0;i<seg;i++) {
        	                temp = path;
                	        temp += argv[1];        // ***.a3m
                        	temp += "_";
	                        temp += to_string(i);
        	                temp += ".a3m";
                	        output.open(temp.c_str());
                        	if (output.is_open()==0) {
                                	cout << temp << " output file is not open\n";
                                	cout << header << endl;
					break;
                        	}

				output << header << endl;
                        	pos[2] = round(i*len);
                        	pos[3] = round((i+1)*len);
               	        
				temp = seq.substr(pos[2],pos[3]-pos[2]);
				for (j=0;j<temp.length();j++) {
					if (temp[j]!='-')
						flag++;
				}
				list << "\t" << flag;

				output << temp << endl;
        	                output.close();
	                }

			list << endl;
			list.close();
			flag = 0;
			continue;
		}

		for (i=0;i<seg;i++) {
			temp = path;
			temp += argv[1];        // ***.a3m
        		temp += "_";
			temp += to_string(i);
			temp += ".a3m";
        		output.open(temp.c_str(),ios::app);
			if (output.is_open()==0) {
		                cout << temp << " output file is not open\n";
                		break;
       			}

			output << header << endl;
			pos[2] = round(i*len);
			pos[3] = round((i+1)*len);
			output << seq.substr(pos[2],pos[3]-pos[2]) << endl;
			output.close();
		}
	}
	
	input.close();
	return 0;
}

