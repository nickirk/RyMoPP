#ifndef _LINE_
#define _LINE_
istream& ignoreline(ifstream& in, ifstream::pos_type& pos){
	pos = in.tellg();
	return in.ignore(numeric_limits<streamsize>::max(), '\n');
}

string getLastLine(ifstream& in){
	ifstream::pos_type pos = in.tellg();
	ifstream::pos_type lastPos;
	while (in >> ws && ignoreline(in, lastPos))
	    pos = lastPos;
	in.clear();
	in.seekg(pos);
	string line;
	getline(in, line);
	return line;
} // get the last line of the file where the energy is.

#endif
