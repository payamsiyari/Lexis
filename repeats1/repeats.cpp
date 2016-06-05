/**
 *	@author Matthias Gall\'e, Centre de Recherche INRIA Rennes - Bretagne Atlantique
 *	@date	Jan 2009
 *	@brief Computing of Largest Maximal Repeats (a.k.a. Near Supermaximal Repeats), Maximal Repeats and Simple Repeats														
 * 	@update: Feb 2012: optionally, take as input a list of integers
 * 	@update: Aug 2012: combines integer input with multilines (separators are negative numbers)
 * 	@update: Aug 2012: add functionality of outputing only largest occurrences of LMR
 * 	@update: Jan 2013: add functions to compute xkcd-repeats
 */

#include<stack>
#include<queue>
#include<iostream>
#include<iterator>
#include<assert.h>		// assert()
#include<math.h>

#include<set>
#include<algorithm>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include <unistd.h>
#include"sarray.h"
#include"sais.h"

using namespace std;

static const char* pygramFormat = "p";

static int DEBUG_LEVEL = 10;

#define MAX(x,y) ( (x)>(y) ? (x) : (y) )

/**
	prints a message according to the current debug level
	@param str message to print
	@param level message is print if this is greater then @DEBUG_LEVEL
*/
//static void debug(char* str, int level){
static void debug(string str, int level){
	if (level > DEBUG_LEVEL){
		cerr << "DEBUG: " << str << endl;
	}
}

class repeat{


};

class node_xkcd{
public:
	node_xkcd(int start, int length, int rightctxt, set<int>* leftctxt){
		st=start;
		len=length;
		rc=rightctxt;
		lc=leftctxt;
	}
	int st;
	int len;
	int rc;
	set<int>* lc;
	
	~node_xkcd(){
		delete lc;
	}
};
	

/**
 A tuple of three elements. Ugly. 
 */
template <class A, class B, class C>
class three{
public:
	three(A f, B s, C t){
		first = f;
		second = s;
		third = t;
	}
	A first;
	B second;
	C third;
};



/** Enhanced Suffix Array.
 */
class esarray{
private:
	int* sarray, *inverse, *lcp, *bwt;
	int* sequence;
	int len, sigmasize;
	char* printFormat;					// specifies print format
	bool integerfile;					// flag: True if input file contains integer IDs
	void* gsequence;
	
	void print(){ //ostream out){
		for (int i=0;i<=len;i++){
			cerr << lcp[i] << "\t";
			cerr << sarray[i] << "\t";
			copy(sequence+sarray[i],sequence+len,ostream_iterator<int>(cerr," "));
			cerr << endl;
		}

	}
	int* translate(char* csequence){
		/** 
			char* -> int* 
		*/
		for(int i=0; i<len-1;i++){
			sequence[i] = (int) (csequence[i] & 0x00FF);
			/*			switch (csequence[i]){
			 case 'A':
			 case 'a':
			 sequence[i] = 1;
			 break;
			 case 'C':
			 case 'c':
			 sequence[i] = 2;
			 break;
			 case 'G':
			 case 'g':
			 sequence[i] = 3;
			 break;
			 case 'T':
			 case 't':
			 sequence[i] = 4;
			 break;
			 }*/
		}
		return sequence;
	}
	
public:

	/**
		set print Format, containing:
			'l': length
			'w': word
			'o': number of occcurrences
			' ': skip
	*/
	void setPrintFormat(char* s){
		printFormat = s;
	}	
	
	/**
	 creates an enhancend suffix array (sarray, inverse, lcp, bwt) and saves the sequence in an integer sequence.
	 inverse and integer could be deleted after the creation of lcp and sarray
	 Execution time O(n) (? depends on which algorithm I'm using)
	 @param multiline: if True, then all occurrences of @param separator are replaced by a unique symbol
	 @param gsequence is either a char* or int* sequence (depending on the flag on @param integer)
	 if @param integer, then @param multiline (and therefore @param separator) are ignored)
	 */
	esarray(void* gseq, int l,bool multiline=false, char separator='\n', int integer=false){
		printFormat = (char*) "w  ";
		len = l; //strlen(csequence) + 1;
		integerfile=integer;
		gsequence=gseq;
		sequence = new int[len+5];		
		if (!integer){
			char* csequence=(char*)gsequence;
			translate((char*) csequence);
			if (multiline){
				int isep = 256;
				for(int i=0;i <len; i++){
					if (csequence[i] == separator){
						sequence[i] = isep++;
					}
				}
			}
		} else{
			copy((int*)gsequence,(int*)gsequence+len,sequence);
		}
		sequence[len-1] = 0;									// sentinel
		sequence[len] = sequence[len+1] = sequence[len+2] = 0;	// requirement of algorithm
		sigmasize = *max_element(sequence, sequence+len);
		/* sarray, inverse, lcp, bwt */
		sarray = new int[len+5];
//		suffixArray(sequence, sarray, len, sigmasize);
		sais_int(sequence, sarray, len, sigmasize+1);
		sarray[len]=0;
		inverse = new int[len];
		for(int i=0; i<len; i++)
			inverse[sarray[i]] = i;
		lcp = new int[len+1];
		lcpa(sarray, sequence, lcp, len, inverse);
		lcp[len] = 0;
		bwt = new int[len+1]; //TODO Can't this be replace by a constant function, reducing storage (or is this better for cache reasons?)
		for(int i=0; i<len; i++)
			bwt[i] = (sarray[i] == 0 ? 0 : sequence[sarray[i]-1]);
		bwt[len]=0;
//		delete[] sequence; sequence = NULL;
		delete[] inverse; inverse = NULL;
	}

	inline int LEletter(int l1, int l2){
		return (l1 == 0 or l1 != l2 ? 0 : l1);
	}



	inline int gbwt(int pos){
		/** as pred_char of sarray[pos], but returns 0 instead of -1 if pos==0*/
		int pred = pred_char(sequence, sarray[pos]);
		return (pred == -1 ? 0 : pred);
	}

	void smr(int minsize=1){
		/** 
			outputs all super-maximal repeats
		*/
		int* alph = new int[sigmasize+1]; /* LAST array in the trick of Puglisi et al 2008 */
		fill(alph, alph+sigmasize+1, 0);
		int i = 1;
		while(i <= len){
			if (lcp[i] > lcp[i-1]){
				int start = i;
				bool all_diff = true;
				if (sarray[i-1] != 0){
					alph[gbwt(i-1)] = start;
				}
				while(i <= len && lcp[i] == lcp[start]){
					if (sarray[i] != 0){
						all_diff &= alph[gbwt(i)] != start;
						alph[gbwt(i)] = start;
					}
					i++;
				}
				if (i > len || (lcp[i] < lcp[start] && all_diff && lcp[start] >= minsize)){ // is maximal over lcp (or, equivalently, arrived the end) and previous characters are different
					output(lcp[start], i-start+1 ,start-1);
					if (i <= len) i++;
				}
			} else i++;

		}
		delete[] alph;
	}



	bool enoughCtxt  (int lc, int rc, int o, float threshold){
		/* 	@param lc and @param rc are the number of left and right ctxt, respectively
			@param o is the number of occurrences
			@param threshold is the threshold line
		*/
//		return lc>=1 && rc>=1 && (lc+rc)/(2.0*o) > 1.0/3.0;
		threshold=MAX(2,(float)o/3.0); //sqrt((float)o);
		return 	lc >= threshold && rc >= threshold; //threshold;
	}

	void xkcd(int x, int minLength=1){
		/** extracts all k--context-diverse repeats 
 *			runs in O(log(|\Sigma|)*|\Sigma|*n) time and outputs (lctxt size, rctxt size, #occs, length)
 * 		*/

		stack<node_xkcd*> st;	
		set<int>* t=new set<int>(); t->insert(gbwt(0));
		st.push(new node_xkcd(0,0,1,t)); // asserts that the stack will never be empty, and avoids cumbersome extreme cases and ifs
		node_xkcd *r;
		int stpos;
		set<int>* savet; savet=NULL;
		for (int i = 1; i<=len; i++){
			stpos=i-1;
			delete savet;
			savet=new set<int>(); savet->insert(gbwt(stpos));
			/* the usual: first process all the repeats that finish here. Then check if another continues or a new one starts */
			while (/*!st.empty() && */ st.top()->len > lcp[i]){

				r=st.top(); st.pop();
				stpos=r->st; /* update start position */
				delete savet;
				savet=new set<int>(); savet->insert(r->lc->begin(),r->lc->end());
//				cout << r->len << endl;
//				if (r->lc->size() > 1 && r->rc > 1)
					cout << r->lc->size() << " " << r->rc << " " << i-r->st << " " << r-> len << endl; //--- "; // << endl;
//				if (enoughCtxt(r->lc->size(),r->rc,i-r->st, (float)x) && r->len >= minLength){
//					output(r->len,i - r->st , r->st);
//				}
				assert(!st.empty()) ;
//				if (!st.empty()) 
				st.top()->lc->insert(r->lc->begin(),r->lc->end()); // inherits all left contexts from child
				delete r;
			}
				assert(!st.empty()) ;
			if (/*!st.empty() && */st.top()->len==lcp[i]){
				st.top()->rc++;
				st.top()->lc->insert(gbwt(i));
			} else{ // isempty or st.top()->len < lcp[i+1]
				// the new repeat already has stpos occurrences
/*				set<int>* t=new set<int>(); 
				for(int j=stpos;j<=i;j++)
					t->insert(gbwt(j)); */
				savet->insert(gbwt(i));
				// Cool fact: there are only 2 diff right context: the one of the childs, and the current
				set<int>* tmp=new set<int>(); tmp->insert(savet->begin(),savet->end());
				st.push(new node_xkcd(stpos,lcp[i],2,tmp));
			}

		}
		if (savet!=NULL) delete savet;
		// there is now only one remaining tuple in the stack
		delete st.top(); 
		st.pop(); 

	}

	/**
	 * outputs all maximal repeats, by traversing the lcp-interval tree
	 * based on "Fast Optimal Algorithms for Computing All the Repeats in a String", Puglisi, et al
	 */
	void mr(int minsize=1){
		int stpos, ibwt, bwt1, bwt2;
		three<int,int,int>* t; /* length, startpos, ? */
		stpos = 0;
//		ilcp = 0; //lcp[0];
		bwt1 = bwt[0];
		stack<three<int,int,int>* > st; // = new stack<three<int,int,int>* >();
		st.push(new three<int,int,int>(0, stpos, bwt1)); //0 == lcp[0]
		for (int i = 1; i <= len; i++){
			stpos = i-1;
			bwt2 = bwt[i];
			ibwt = LEletter(bwt1, bwt2);
			bwt1 = bwt2;
			while(st.top()->first > lcp[i]){
				t = st.top(); st.pop();
				int prevbwt = t->third;
				if (prevbwt == 0 && t->first >= minsize){
					/* between t->second and i+1 (non inclusive) over the sarray there is a maximal repeat of size t->first */
					output(t->first, i - t->second, t->second);
				}
				stpos = t->second;
				if (!st.empty()){
					st.top()->third = LEletter(prevbwt, st.top()->third);
				}
				ibwt = LEletter(prevbwt, ibwt);
				delete t;
			}
			if (st.top()->first == lcp[i]){
				st.top()->third = LEletter(st.top()->third, ibwt);
			} else{ // top()->length < lcp[i]
				three<int,int,int> * t = new three<int,int,int>(lcp[i],stpos,ibwt);
				if (t==NULL){
					cerr << "Not enough space for a three" << endl;
					exit(-1);
				}
				st.push(t);
			}
		}

		//cerr << st.top()->first << " " << st.top()->second << " " << st.top()->third << endl;

		delete st.top();
		st.pop();
//		cerr << st.size() << endl;
//		cerr << st.top()->first << " " << st.top()->second << " " << st.top()->third << endl;
		assert (st.empty()); 
	}



	void r(int minsize=1){
		int stpos;
		pair<int,int>* t;
		stack<pair<int,int>* > s; // elements in s correspond to <len,startpos> of repetitions which end position is not yet founded
		s.push(new pair<int,int>(0,0));
		for(int i=1; i<=len; i++){
			stpos = i-1;
			while(s.top()->first > lcp[i]){
				/* repeat over [s.top()->second, i] of length s.top()->first */
				t=s.top();
				stpos = t->second;
				output(t->first, i - t->second, t->second);
				delete t;
				s.pop();
			}
			if (s.top()->first < lcp[i]){
				int limit = s.top()->first + 1;
				for(int len = limit; len <= lcp[i]; len++)
					s.push(new pair<int,int>(len,stpos));
			} // don't do nothing if equal
		}
		while(s.size()>1){ // last one is the empty one
			t=s.top();
			output(t->first, len - t->second + 1, t->second);
			delete t;
			s.pop();
		}
		delete s.top(); 	s.pop();

		assert(s.empty());
	}

/*	void r(int minsize=1){
		stack<pair<int,int>* > s; // elements in s correspond to <len,startpos> of repetitions which end position is not yet founded
		for(int i=0; i<=len; i++){
			int limit = minsize; //s.empty() ? minsize : s.top()->first;
			if(lcp[i] >= limit){
				for(int len = limit; len <= lcp[i]; len++)
					s.push(new pair<int,int>(len,i-1));
			} else{
				while(s.size()>1 && lcp[i] < s.top()->first){
					// repeat over [s.top()->second, i] of length s.top()->first 
					output(s.top()->first, i - s.top()->second, s.top()->second);
					delete s.top();
					s.pop();
				}
			}
		}
//		while(!s.empty()){
		while(s.size()>1){ // last one is the empty one
			output(s.top()->first, len - s.top()->second + 1, s.top()->second);
			delete s.top();
			s.pop();
		}
		delete s.top(); 	s.pop();

		assert(s.empty());
	}*/

	
	
	
	/*
	 Computes all the largest maximal repeats.
	 @pre calculation of lcp, bwt (and sarray if the output shows information over the sequence).
	 @pre lcp must be of length len+1 (lcp[len] should be 0).
	 
	 O(n) execution time, uses 4*n of auxiliary storage + n*sizeof(three) during computation
			(2012, Aug 27th update: at last term there is an additional n term (I also store which occurrences are lcunique, but they can be at most n of these)
	 
	 */
	void lmr(int minLength = 1,bool onlylargestoccurrences=false){
		#define isLCUnique(i) (prev[(i)] < trees2.top()->first && next[(i)] > *(endTree[(i)]))
		int i;
		/* prev and next */
		int* next = new int[len];
		int* prev = new int[len];
		int lastchar[sigmasize+1]; for(i=0;i<sigmasize+1;i++) lastchar[i] = -1;
		for(i=0;i<len;i++){
			prev[i] = lastchar[bwt[i]];
			lastchar[bwt[i]] = i;
		}
		for(i=0;i<sigmasize+1;i++) lastchar[i]=len+1;

		for(i=len-1;i>=0;i--){
			next[i] = lastchar[bwt[i]];
			lastchar[bwt[i]] = i;
		}
		delete[] bwt; bwt = NULL;
		
		debug("auxiliary structures created",1);
		
		/* end of tree: computing of the final position over the suffix array of each tree */
		stack<three<int,int,int*>*> trees;
		int** endTree = new int*[len]; for(int i=0; i<len;i++) endTree[i] = NULL;
		int* final = new int[1]; *final = len;
		trees.push(new three<int,int,int*>(0,0,final));								/* the $ interval (lcp = 0) */
		int** toClean = new int*[len]; for(int i=0; i<len;i++) toClean[i] = NULL;	/* to keep track of allocated memory */
		int lb;
		debug("starting first traversal",1);
		/* traversing of the lcp-interval tree */
		for(i=0;i<len;i++){
			lb =i ;
			if ( trees.top()->second > lcp[i+1]){
				endTree[i] = trees.top()->third;
				while( trees.top()->second > lcp[i+1]){
					lb = trees.top()->first;
					*(trees.top()->third) = i;										/* i is the end position of the tree */
					delete trees.top();	trees.pop();
				}
			}
			if (trees.top()->second < lcp[i+1]){
				int* future = new int[1];
				toClean[i] = future;
				if (endTree[i] == NULL)	endTree[i] = future;
				trees.push(new three<int,int,int*>(lb,lcp[i+1],future));
			} else{ // ==
				if (endTree[i] == NULL) endTree[i] = trees.top()->third;
			}
		}
		delete trees.top(); trees.pop();											/* the $ interval (lcp = 0) */
		
		debug("end of inmediate tree found",1);
				
		/* computation of LMR: second traversing of the lcp-interval tree */
		stack<three<int,int,vector<int>* >* > trees2;
		trees2.push(new three<int,int,vector<int>* >(0,0,new vector<int>()));
//		bool b;
		for(i=0;i<len;i++){
			lb =i ;

			if (lcp[i+1] < trees2.top()->second){
				if (isLCUnique(i)){ //(prev[i] < trees2.top()->first && next[i] > *(endTree[i])); 
					trees2.top()->third->push_back(sarray[i]);
				}
//				if (b) cerr << sarray[trees2.top()->first] << ' ' << trees2.top()->second << " -> " << sarray[i] << endl;
				while(lcp[i+1] < trees2.top()->second){
					/* a candidate to be LMR */
					lb = trees2.top()->first;
//					cout << "not a lmr " << i << ' '; output(trees2.top()->second, i - trees2.top()->first + 1, trees2.top()->first);
					if (trees2.top()->third->size()>0){	

						/* Bingo!: a LMR of length trees2.top()->second, 
							where all the occurrences are in the interval [trees2.top()->first, i] over the suffix array
						 */
						if (trees2.top()->second >= minLength){
							output(trees2.top()->second, i - trees2.top()->first + 1, trees2.top()->first, onlylargestoccurrences ? trees2.top()->third : NULL);
						}
					} 
					delete trees2.top();	trees2.pop();
					assert(i == *(endTree[i])); 
				}
				
			}
			if (lcp[i+1] > trees2.top()->second){
				trees2.push(new three<int,int,vector<int>* >(lb,lcp[i+1],new vector<int>()));

			}

			if (i == lb){
				if (prev[i] < trees2.top()->first && next[i] > *(endTree[i])){
					trees2.top()->third->push_back(sarray[i]);
				}
//				if (b) cerr << sarray[trees2.top()->first] << " "  << trees2.top()->second << " -> " << sarray[i]<< endl;
			}

		}
		delete trees2.top()->third; delete trees2.top(); trees2.pop(); delete[] final;
		
		for(int i=0;i<len;i++) if (toClean[i] != NULL) delete[] toClean[i];
		delete[] toClean;
		
		delete[] endTree;
		delete[] next;
		delete[] prev;
	}
	


	void printsubseq(int startpos, int length){
		if (integerfile)
			copy((int*)gsequence+sarray[startpos], (int*)gsequence+sarray[startpos]+length, ostream_iterator<int>(cout," "));
		else
			copy((char*)gsequence+sarray[startpos], (char*)gsequence+sarray[startpos]+length, ostream_iterator<char>(cout,""));

	}
	/**
		Execution time depends of the selected option. If it is size (s) or numb. of occurrences (o) it is constant. 
	*/
	void output(int length, int occurrences, int startpos, vector<int>* soccs=NULL){
		if (strcmp(printFormat, pygramFormat) == 0){
			/* pygram format */
			cout << length << ":m:"; /* mark all as maximal repeat */
			printsubseq(startpos,length);
			cout << ":-:0(" << length << '#';
			vector<int>* occs; 
			if (soccs==NULL){
				occs=new vector<int>(occurrences);
				for(int i=0; i<occurrences; i++){
					occs->at(i)=(sarray[startpos + i]);
				}
			} else{
				occs=soccs;
			}
			sort(occs->begin(), occs->end());
			for(vector<int>::iterator i=occs->begin(); i != occs->end() - 1; i++)
				cout << *i << ',';
			cout << *(occs->end()-1) << ")." << endl;
			if (soccs==NULL) delete occs;
		}
		else{
			for (int i = 0; i<(int) strlen(printFormat); i++){
				switch(printFormat[i]){
					case 's':
						cout << length << ' ';
						break;
					case 'w':
						printsubseq(startpos,length);
						cout << ' ';
						break;
					case 'o':
						cout << occurrences << ' ';
						break;
					case 'l':
						vector<int>* occs; 
						if (soccs==NULL){
							occs=new vector<int>(occurrences);
							for(int i=0; i<occurrences; i++){
								occs->at(i)=(sarray[startpos + i]);
							}
						} else{
							occs=soccs;
						}
						sort(occs->begin(), occs->end());
						cout << "(";
						for(vector<int>::iterator i=occs->begin(); i != occs->end() - 1; i++)
							cout << *i << ',';
						cout << *(occs->end()-1) << ")";
						if (soccs==NULL) delete occs;
						break;
				}
			}
			cout << endl;
		}
		return;
	}
	
	~esarray(){
		/* release memory */
		if (sequence != NULL) delete[] sequence;
		delete[] sarray;
		if (inverse != NULL) delete[] inverse;
		delete[] lcp;
		if (bwt != NULL) delete[] bwt;
	}
	
	int getSequenceSize(){
		return len-1;
	}
	int getLengthLongestRepeat(){
		return *(max_element(lcp, lcp+len));
	}
};

/**
	read file character by character
*/
char* readCFile(char* path, int& len){
	FILE* f = fopen(path, "r");
	if (f== NULL){
		perror("filename");
		exit(-1);
	}

	/* get length */
	fseek (f , 0 , SEEK_END);
	len = ftell (f);
	rewind (f);

	/* read it */
	char* csequence = new char[len+4];
	fread(csequence, len, 1, f);
	csequence[len] = '\0';
	fclose(f);
	return csequence;
}

/**
	read file integer by integer. 
	Have to worry about separators (in char case it is handled by the esarray creator).
	if @param separator is true, then any negative symbol will be considered a separator
*/
int* readIFile(char* path, int& len, bool separator=false){
	vector<int> v;
	int i;
	int* isequence;
	int maxsym=-1;
	//FILE* f = fopen(path, "r");
	//if (f== NULL){
	//	perror("filename");
	//	exit(-1);
	//}
	const int BUF_SIZE = 1024;
	char buffer[BUF_SIZE];
	size_t contentSize = 1; // includes NULL
	/* Preallocate space.  We could just allocate one char here,
	 but that wouldn't be efficient. */
	char *content = (char*)malloc(sizeof(char) * BUF_SIZE);
	if(content == NULL)
	{
		perror("Failed to allocate content");
		exit(1);
	}
	content[0] = '\0'; // make null-terminated
	while(fgets(buffer, BUF_SIZE, stdin))
	{
		char *old = content;
		contentSize += strlen(buffer);
		content = (char*)realloc(content, contentSize);
		if(content == NULL)
		{
			perror("Failed to reallocate content");
			free(old);
			exit(2);
		}
		strcat(content, buffer);
	}
	
	if(ferror(stdin))
	{
		free(content);
		perror("Error reading from stdin.");
		exit(3);
	}
	/* save result in vector (no idea how big it is */
	int offset;
	while( sscanf(content,"%i%n",&i,&offset)==1){
	//while( fscanf(f,"%i",&i)==1){
		//printf("%d",i);
		content += offset;
		v.push_back(i);
		if (i>maxsym) maxsym=i;
	}
	//fclose(f);
	len=v.size();
	isequence=new int[len+4];
	for(int i=0;i<len;i++){
		if (v[i]<0) v[i]=++maxsym;
		isequence[i]=v[i];
	}
	return isequence;
}

typedef enum{REPEAT, MAXREPEAT, LMAXREPEAT, LMAXOCCUR, SMAXREPEAT, XKCDREPEAT} repeats_t;
int main(int argc, char** argv){
	char c;
	int minsize = 1;
	char* printFormat = (char*)"w  ";
	repeats_t wordtype = MAXREPEAT;
	bool multiline=false;
	bool isinteger=false;
	char separator = '\n';
	int* isequence;
	char* csequence; /* only one of both will be used */

	while((c = getopt(argc, argv, "ihmn:p:r:vx:")) != -1){
		switch(c){
		case 'i':
			isinteger=true;
			break;
		case 'm':
			multiline = true;
			break;
		case 'n':
			minsize = atoi(optarg);
			if (minsize == 0){
				cerr << "Argument \"" << optarg << "\" must be a positive integer" << endl;
				exit(-1);
			}
			break;
		case 'p':
			printFormat = optarg;
			break;
		case 'r':
			if (strcmp(optarg, "r") == 0) wordtype = REPEAT;
			else if (strcmp(optarg, "mr") == 0) wordtype = MAXREPEAT;
			else if (strcmp(optarg, "lmr") == 0) wordtype = LMAXREPEAT;
			else if (strcmp(optarg, "Lmr") == 0) wordtype = LMAXOCCUR;
			else if (strcmp(optarg, "smr") == 0) wordtype = SMAXREPEAT;
			else if (strcmp(optarg, "xkcd") ==0) wordtype = XKCDREPEAT;


			else cerr << "Repeat type " << optarg << " not recognized. Using maximal repetas instead" << endl;
			break;
		case 'v':
			DEBUG_LEVEL = -1;
			debug("verbose mode: ON",0);
			break;
		case 'x':
			separator = optarg[0];
			break;

		case 'h':
help:			cerr << "Usage: " << argv[0] << " [-v | -h | -n <minsize> | -v | -r (r | mr | lmr | Lmr | smr | xkcd) | -p slowp -x sep] <filename>" << endl;
				cerr	<< "\t -v: verbose mode" << endl
						<< "\t -h: help" << endl
						<< "\t -i: integer sequences" << endl
						<< "\t -n: minimal size of repeat" << endl
						<< "\t -m: multiline: permits several sequences. Uses as separator @separator (-x) or any negative integer (if -i is specified)" << endl
						<< "\t -r kind of repeat (normal repeat, maximal repeat, largest maximal repeat (Lmr is for only largest-occurrences of LMR's), x-kcd" << endl
						<< "\t -p what to print:" << endl
						<< "\t\t s: size of repeat" << endl
						<< "\t\t l: list of occurrences" << endl
						<< "\t\t o: number of occurrences" << endl
						<< "\t\t w: the word itself" << endl
						<< "\t\t p: pygram format" << endl
						<< "\t -x <separator>: indicates the separator (newline by default) " << endl
				;
				
						
			exit(0);
		}
	}
	if (argc < 2) goto help;
	
	/* read file */
	char* filePath = argv[argc-1];
	int len=1;
	if (isinteger)
		isequence =  readIFile(filePath,len,multiline);
	else
		csequence = readCFile(filePath,len);
	debug("File read",1);
	esarray * esa = new esarray((isinteger?(void*) isequence:(void*) csequence),len+1, multiline, separator,isinteger);
	esa->setPrintFormat(printFormat);
	debug("ESA creataed",1);
	
	if (strcmp(printFormat, pygramFormat) == 0){
		cout << ">00>sequences:1" << endl << ">01>0:";
		cout << esa->getSequenceSize() << ':' << filePath << endl << ">02>data:" << esa->getLengthLongestRepeat() << endl;
	}
	cout << "#>length= " << len << endl;
	switch(wordtype){
		case REPEAT:
			esa->r(minsize);
			break;
		case MAXREPEAT:
			esa->mr(minsize);
			break;
		case LMAXREPEAT:
		case LMAXOCCUR:
			esa->lmr(minsize,wordtype==LMAXOCCUR);
			break;
		case XKCDREPEAT:
			esa->xkcd(2,minsize);
			break;
		case SMAXREPEAT:
			esa->smr(minsize);
			break;
	}
	
	debug("all found, now clean up",1);
	delete esa;
	if (isinteger) delete[] isequence; else delete[] csequence;
}

