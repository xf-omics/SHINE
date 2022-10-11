/*
  parsing the variant file to format input for ESM1b
  check if the reference amino acids match the protein sequences
  segment protein sequences if longer than 1022
  output fasta
  1/31/2022
  add -std=c++0x when compiling for to_string
*/
#include "/share/terra/Users/xf2193/resource/CPP/class_list.h"

int main(int argc, char *argv[]) {
	//arguments check
	if (argc<2) {
		cout << "Please provide input file name!\n";
		return 1;
	}

	size_t found[2];
	unsigned int seg, len, pos[2];
        unsigned long long i, j;
	string temp, path, file, id, info, gene, transcript, sequence, amino[2];
	ifstream input, search;
	ofstream output, list;
	
	path = argv[1];
	found[0] = path.find_last_of('.');
        if (found[0]==string::npos) {
                temp = path;
		file = path;
                path += ".var";
                temp += "_1b.var";
		file += "_1b.fasta";
        }
        else {
                temp = path.substr(0,found[0]);
        	file = temp;
		temp += "_1b.var";
		file += "_1b.fasta";
        }

	input.open(path.c_str());
        if (input.is_open()==0) {
                cout << path << " variant file is not open\n";
                return 1;
        }
	
	output.open(file.c_str());
	if (output.is_open()==0) {
		cout << file << " sequence file is not open\n";
		return 1;
	}

        list.open(temp.c_str());
        if (list.is_open()==0) {
                cout << temp << " list file is not open\n";
                return 1;
        }

	path = "/share/terra/Users/xf2193/resource/ensembl/MANE_ensembl_protein/";

// ID      CHROM   POS     REF     ALT     Pathogenicity   indel   Gene    Feature Protein.stable.ID       Symbol  Func    Protein_position        Amino_acids     gnomAD  gnomAD_popmax   gnomAD3 UKBB
	list << "VarID\tCHROM\tPOS\tREF\tALT\tLabel\tindel\tGene\tProtein\tSymbol\tFunc\tgnomAD2\tgnomAD2_max\tgnomAD3\tUKBB\tTranscript\tAA_pos\tAA\tindex\n";
	getline(input,temp);	// header
	while (input.good()&&!input.eof()) {
		getline(input,temp);
		if (temp=="\0") continue;
	
		found[0] = temp.find("\t"); // ID
		for (i=0;i<4;i++) {
			found[1] = found[0];
			found[0] = temp.find("\t",found[0]+1);
			if (i>0) {
				id += "_";
				id += temp.substr(found[1]+1,found[0]-found[1]-1);
			}
			else
				id = temp.substr(found[1]+1,found[0]-found[1]-1);
		}
		found[1] = temp.find("\t");
		info = id;
		info += temp.substr(found[1],found[0]-found[1]+1);
		found[1] = temp.find("Benign",found[0]+1);	// Pathogenicity

                if (found[1]!=string::npos)
                        info += "-1";
                else {
                        found[1] = temp.find("Unaffected",found[0]+1);
                        if (found[1]!=string::npos)
                                info += "-1";
                        else
                                info += "1";
                }

		found[0] = temp.find("\t",found[0]+1);
		found[1] = found[0];
                for (i=0;i<2;i++)
                        found[1] = temp.find("\t",found[1]+1);
		info += temp.substr(found[0],found[1]-found[0]);

                found[0] = temp.find("\t",found[1]+1);
		gene = temp.substr(found[1]+1,found[0]-found[1]-1);
		if (transcript!=gene) {
			transcript = gene;
			file = path + transcript;
			file += ".fasta";
			search.open(file.c_str());
		        if (search.is_open()==0) {
                		cout << file << " fasta file is not available\n";
                		continue;
        		}
			getline(search,sequence); // transcript id
			getline(search,sequence); // sequence
			len = sequence.length();       	
			search.close();
		}

		found[1] = found[0];
                for (i=0;i<3;i++)
                        found[1] = temp.find("\t",found[1]+1);
                info += temp.substr(found[0],found[1]-found[0]);

		found[0] = temp.find_first_of("-	",found[1]+1);
		pos[0] = atoi(temp.substr(found[1]+1,found[0]-found[1]-1).c_str());
		if (temp[found[0]]=='-') {
			found[1] = temp.find("\t",found[0]+1);
			pos[1] = atoi(temp.substr(found[0]+1,found[1]-found[0]-1).c_str());
			found[0] = found[1];
		} else {
			pos[1] = pos[0];
		}
		
		found[1] = temp.find("/",found[0]+1);
		amino[0] = temp.substr(found[0]+1,found[1]-found[0]-1);
		found[0] = temp.find("\t",found[1]+1);
		amino[1] = temp.substr(found[1]+1,found[0]-found[1]-1);
		info += temp.substr(found[0]);

		found[1] = amino[1].find('*');
		if (found[1]!=string::npos) continue; // stop-gain

		if (amino[0]!="-"&&sequence.substr(pos[0],amino[0].length())!=amino[0]) {
                	pos[0]--;
			pos[1]--;
			if (sequence.substr(pos[0],amino[0].length())!=amino[0]) {
				cout << temp << " amino acid does not match the reference " << sequence.substr(pos[0]-1,amino[0].length()+3) << endl;
				continue;
			}
                }
		
		if (amino[0]!="-"&&amino[1]!="-") {
			if (amino[0].length()>amino[1].length()) { // deletion
				while (amino[0][0]==amino[1][0]) {
					pos[0]++;
					amino[0].erase(0,1);
					if (amino[1].length()==1) {
						amino[1] = '-';
						break;
					} else
						amino[1].erase(0,1);
				} 
				while (amino[0][amino[0].length()-1]==amino[1][amino[1].length()-1]) {
					pos[1]--;
					amino[0].erase(amino[0].length()-1,1);
					if (amino[1].length()==1) {
						amino[1] = '-';
						break;
					} else
						amino[1].erase(amino[1].length()-1,1);
				}
				if (amino[1]!="-") {
				//	cout << temp << " deletion across codons\n";
					continue;
				}
			} else {
			if (amino[0].length()<amino[1].length()) { // insertion
				while (amino[0][0]==amino[1][0]) {
					pos[1]++;
					amino[1].erase(0,1);
				       	if (amino[0].length()==1) {
						amino[0] = '-';
						break;
					} else
						amino[0].erase(0,1);
				}
				while (amino[0][amino[0].length()-1]==amino[1][amino[1].length()-1]) {
					pos[0]--;
					amino[1].erase(amino[1].length()-1,1);
				       	if (amino[0].length()==1) {
						amino[0] = '-';
						break;
					} else
						amino[0].erase(amino[0].length()-1,1);
				}
                                if (amino[0]!="-") {
                                //        cout << temp << " insertion across codons\n";
                                        continue;
                                }
			} else
				cout << temp << " not an INDEL\n";
			}
		}

		temp = sequence;
                if (amino[1]!="-") { // insertion
			if (amino[0]!="-") // indel
				continue;
			temp.insert(pos[0],amino[1]);
			len = temp.length();
                        if (len>=500) {
				i = (pos[0] + pos[1])/2;
				if (i>250) i -= 250;
				else i = 0;	// sequence start position
				temp = temp.substr(i,500);
				pos[0] -= i;
				pos[1] -= i;
			}
			output << '>' << id << endl;
			output << temp << endl;
			list << info << "\t" << transcript << "\t" <<  pos[1] << "\t" << amino[1] << "\t0\n"; // position 1-based
                } else { // deletioni
			j = 0;
			if (len>=1023) {
				seg = len / 1022;
				for (i=0;i<=seg;i++) {
					if (pos[0] <= len * (i+1) / (seg+1)) {
						if (pos[1] <= len * (i+1) / (seg+1)) {
							transcript += '_';
							transcript += to_string(i);
							pos[0] -= len * i / (seg+1);
                                                        pos[1] -= len * i / (seg+1);
						} else { // deletions in two separate segments
                                                        transcript += '_';
							for (pos[0]=pos[0]-(len*i/(seg+1));pos[0]<=len*(i+1)/(seg+1);pos[0]++) {
								list << info << "\t" << transcript << i << "\t" << pos[0] << "\t" << amino[0][j] << "\t" << j << endl;
								j++;
							}
                                                        i++;
							transcript += to_string(i);
                                                        pos[0] = 1;
                                                        pos[1] -= len * i / (seg+1);
						}
						break;
					}
				}
			}
			for (;pos[0]<=pos[1];pos[0]++) {
				list << info << "\t" << transcript << "\t" << pos[0] << "\t" << amino[0][j] << "\t" << j << endl;
                                j++;
                        }
		}
// ID     CHROM   POS     REF     ALT     Label   indel   Gene Protein.stable.ID     Symbol  Func   gnomAD  gnomAD_popmax   gnomAD3 UKBB transcript AA_pos AA_ref/AA_alt
	}

	input.close();
	output.close();
	list.close();
	return 0;
}

