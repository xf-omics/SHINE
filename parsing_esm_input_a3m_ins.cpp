/*
  parsing the insertion variant file, use together with the msa segments
  1/3/2022
  input: variant file, index file, and msa directory
  remove deletion
  update amino acid positions to map msa positions

  add -std=c++0x when compiling for to_string
*/
// #include "/share/mendel/home/xf2193/Script/CPP/class_list.h"
#include "/share/terra/Users/xf2193/resource/CPP/class_list.h"

int main(int argc, char *argv[]) {
	//arguments check
	if (argc<4) {
		cout << "Please provide input file name, index file and directory of msa sequences!\n";
		return 1;
	}

	size_t found[3];
	
	unsigned int i, j, len, pos[5], num;
	unsigned long long seg;
	string temp, index, info, path, segment, file, id, transcript, gene, sequence, amino[2];
	ifstream input, search;
	ofstream list;
	
	path = argv[1];
	found[0] = path.find_last_of('.');
	if (found[0]==string::npos) {
		temp = path;
                path += ".var";
		temp += "_msa.var";
	}
	else {
		temp = path.substr(0,found[0]);
                temp += "_msa.var";
	}

	input.open(path.c_str());
        if (input.is_open()==0) {
                cout << path << " variant input file is not open\n";
                return 1;
        }
	
        list.open(temp.c_str());
        if (list.is_open()==0) {
                cout << temp << " variant output file is not open\n";
                return 1;
        }

//        index = "/share/terra/Users/xf2193/resource/ensembl/ensembl_genetree_gene_protein_index.txt";
//	path = "/share/terra/Users/xf2193/resource/ensembl/a3m_seg/";
	index = argv[2];
	path = argv[3];
	
// ID	CHROM	POS	REF	ALT	Pathogenicity	indel	Gene	Feature	Protein.stable.ID	Symbol	Func	Protein_position	Amino_acids	gnomAD	gnomAD_popmax	gnomAD3	UKBB
	getline(input,temp); // header
	list << "ID\tCHROM\tPOS\tREF\tALT\tLabel\tindel\tTranscript\tProtein\tSymbol\tFunc\tgnomAD2\tgnomAD2_max\tgnomAD3\tUKBB\tGene\tAA_msa_pos\tAA_alt\n";
//	list << temp << endl;
	while (input.good()&&!input.eof()) {
		getline(input,temp);
		if (temp=="\0") continue;

		found[0] = 0;

                for (i=0;i<5;i++)
                        found[0] = temp.find("\t",found[0]+1); // before Pathogenicity
                info = temp.substr(0,found[0]+1);
		found[1] = temp.find("Benign",found[0]+1);

                if (found[1]!=string::npos)
                        info += "-1";
                else {
                        found[1] = temp.find("Unaffected",found[0]+1);
                        if (found[1]!=string::npos)
                                info += "-1";
			else
				info += "1";
                }

                found[1] = temp.find("\t",found[0]+1); // after pathogenicity
		found[0] = temp.find("\t",found[1]+1); // before gene stable ID
		info += temp.substr(found[1],found[0]-found[1]); // indel
		
		found[1] = temp.find("\t",found[0]+1);	// gene ensembl ID
		transcript = temp.substr(found[0]+1,found[1]-found[0]-1);

		found[0] = found[1];
                for (i=0;i<4;i++)
                        found[1] = temp.find("\t",found[1]+1); // after Func
		info += temp.substr(found[0],found[1]-found[0]);

		found[0] = temp.find_first_of("-	",found[1]+1);  // Protein_position
                pos[0] = atoi(temp.substr(found[1]+1,found[0]-found[1]-1).c_str());
                if (temp[found[0]]=='-') {
                        found[1] = temp.find("\t",found[0]+1);
                        pos[1] = atoi(temp.substr(found[0]+1,found[1]-found[0]-1).c_str());
                } else {
                        found[1] = found[0];
                        pos[1] = pos[0];
                }

                found[0] = temp.find("/",found[1]+1);   // Amino_acids
                amino[0] = temp.substr(found[1]+1,found[0]-found[1]-1);	// reference
                found[1] = temp.find("\t",found[0]+1);
                amino[1] = temp.substr(found[0]+1,found[1]-found[0]-1); // alternate
		
		info += temp.substr(found[1]); // AF
		info += "\t";
		info += transcript;
		info += "_";
            
		if (amino[0]!="-") {
                        if (amino[0][0]==amino[1][0]) {
                                if (amino[0][0]=='*') {
					amino[0]="-";
					pos[0]--;
		//	insertion is at the end of protein
				} else
					pos[1]++;	// did not modify amino[0] in order to check the reference allele
				amino[1].erase(0,1);
				pos[4] = pos[0];	// position of the reference allele
                        } else {
                                i = amino[0].length()-1;
                                j = amino[1].length()-1;
                                if (amino[0][i]==amino[1][j]) {
                                //        amino[0].erase(i,1); save to check reference allele
					amino[1].erase(j,1);
					pos[0]--;
				//	if (amino[0]=="\0")
				//		amino[0] = '-';
					pos[4] = pos[1];
                                } else {
					cout << temp << " insertion in the middle of a codon\n";
					continue;
				}
                        }
			if (amino[0].length()>1) {
				cout << "complicated insertion\n";
				cout << temp << endl;
				continue;
			}
                } else
			pos[4] = pos[0];
		
		if (gene!=transcript) {
			gene = transcript;
			search.open(index.c_str());
			if (search.is_open()==0) {
                		cout << index << " file is not open\n";
                		break;
        		}

			while (search.good()&&!search.eof()) {
				getline(search,segment);
				if (segment=="\0") continue;
				found[0] = segment.find("\t"); // after gene id
				if (segment.substr(0,found[0])==gene)
					break;
			}
			search.close();

                	if (segment.substr(0,found[0])!=gene) {
                        	cout << gene << " is not found in the protein segment index file\n";
                        	gene = "";
				continue;
                	}
		} else {
			found[0] = segment.find("\t"); // after gene id
		}

		found[2] = segment.find("\t",found[0]+1); // after protein id
		seg=0;
		pos[2] = 0;
		while (1) {
			found[1] = segment.find("\t",found[2]+1); // amino acid for each segment
			if (found[1]!=string::npos)
				pos[3] = atoi(segment.substr(found[2]+1,found[1]-found[2]-1).c_str());
			else
				pos[3] = atoi(segment.substr(found[2]+1).c_str());
			
			if (pos[0]<=pos[3])
				break;
                        if (found[1]==string::npos) {
				cout << "insertion is beyond the size of the protein\n";
				cout << segment << ", " << seg << endl;
				cout << temp << endl;
				break;
			}
			found[2] = found[1];
			pos[2] = pos[3];
			seg++;
		}

		file = path + transcript;
		file += "_";
		file += to_string(seg);
		file += ".a3m";
		search.open(file.c_str());
	        if (search.is_open()==0) {
                	cout << file << " msa file is not available\n";
                	continue;
        	}
		getline(search,sequence); // transcript id
	        getline(search,sequence); // sequence
			
		search.close();
		len = sequence.length();
		for (i=0;i<len;i++) {
			if (sequence[i]!='-') {
				pos[2]++;
				if (pos[4]-pos[2]==1) j = i; // record the last position
                                if (pos[4]==pos[2]) { // break after reference check
                                        if(amino[0]!="-") {
                                               if (sequence[i]!=amino[0][0]) {
						       if (sequence[j]!=amino[0][0]) {
                                               			cout << "amino acid reference does not match\n";
                                                	        cout << temp << "\t" << pos[4] << "\t" << sequence.substr(i-1,3) << endl;
                                                        	j = -1; // fail flag
                                                        	break;
                                                	} else {
								if (pos[4]>pos[0]) {
									cout << "amino acid reference check\n";
                                                                	cout << temp << "\t" << pos[4] << "\t" << sequence.substr(j-1,3) << endl;
                                                                	j = -1; // fail flag
                                                                	break;
								}
								i = j;
							}
					       }
					}
					if (pos[4]>pos[0]) i = j;
                                        break;
                                }
			}
		}

		if (j==-1) continue;
		if (++i<len) // i is 0-based
			list << info << seg << "\t" << i << "\t" << amino[1] << endl;
		else {
//			cout << "insertion falls in-between two segments\n" << temp << endl;
        		file = path + transcript;
                	file += "_";
	                file += to_string(seg+1);
	                file += ".a3m";
       			search.open(file.c_str());
	                if (search.is_open()==1) {
       			        list << info << seg+1 << "\t" << 0 << "\t" << amino[1] << endl;;
 				search.close();
	                } else
				list << info << seg << "\t" << i-1 << "\t" << amino[1] << endl;
                }
	}

	list.close();
	input.close();
	return 0;
}

