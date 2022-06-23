#include "/share/terra/Users/xf2193/resource/CPP/class_list.h"

int main(int argc, char *argv[]) {
	//arguments check
	if (argc<2) {
		cout << "There is no input file basename!\n";
		return 1;
	}

	unsigned int i, j, depth=0, len, flag;
	string path, header[301], seq[301];
	ifstream input;
	ofstream output;
	
	path = "/share/terra/Users/xf2193/resource/ensembl/a3m_trim/";
	path += argv[1];	// ***.a3m
	path += ".a3m";
	input.open(path.c_str());

	path = "/share/terra/Users/xf2193/resource/ensembl/script/temp";
	output.open(path.c_str());

	if (input.is_open()==0||output.is_open()==0) {
		cout << "input or output file is not open\n";
		return 1;
	}
	
	while (input.good()&&!input.eof()) {
                getline(input,path);
                if (path=="\0") continue;
		if (depth>300) {
			cout << argv[1] << " has more than 301 species\n";
			return 1;
		}
		header[depth] = path;
		getline(input,seq[depth]);
		if (seq[depth++]=="\0") {
			cout << argv[1] << ' '<< depth << " wrong format!\n";
			return 1;
		}
	}

	len = seq[0].length();
	for (i=0;i<len;i++) {
		flag = 0;
		for (j=0;j<depth;j++)
			if (seq[j][i]!='-') {
				flag = 1;
				break;
			}
		if (flag==0) {
			for (j=0;j<depth;j++)
				seq[j].erase(i,1);
			i--;
			len--;
		}
	}

	for (i=0;i<depth;i++) {
		output << header[i] << endl;
		output << seq[i] << endl;
	}

	output.close();
	input.close();
	return 0;
}

