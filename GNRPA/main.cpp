#include <iostream>
#include <fstream>
#include <vector>

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <list>
#include <iterator>
#include <sstream>
#include <math.h>
#include <time.h>

using namespace std;


const int MaxPlayoutLength = 5000;
const int MaxLegalMoves = 2;

vector<vector<int> > neighbours;
vector<float> proba;

int maxNeighbours = 0;

bool useFree = true;

class Move
{
public:
    int next_vertex;
    bool p;//play or not
};

class Board
{
public:
  set<int> current_set;
  Move rollout [MaxPlayoutLength];
  int l, length;
  bool check [MaxPlayoutLength];
  int free [MaxPlayoutLength];
  int support [MaxPlayoutLength];
  bool tried [MaxPlayoutLength];

  Board()
    {
      init ();
    }

  void init ()
    {
      l = 0;
      length = 0;
      for (int i=0; i<neighbours.size (); i++) {
	check [i] = false;
	tried [i] = false;
	free [i] = neighbours [i].size ();
	support [i] = 0;
	/*
	for (int j =0 ; j < neighbours [i].size (); j++)
	  support [i] += neighbours [neighbours [i] [j]].size ();
	*/
      }
    }

  int code (Move m,  float & b) {
    b = 0.0;
    if (m.p)
      b = 2*(2 * proba[m.next_vertex]-1) ;
    return m.next_vertex * 2 + m.p;

  }

/*

  int code (Move m,  float & b) {
    b = 0.0;
    if (m.p)
        b = 1 * proba[m.next_vertex];
    return m.next_vertex * 2 + m.p;
  }


  int code (Move m,  float & b) {
    b = -2 * log(1+(2*proba[m.next_vertex] - 1));
    if (m.p)
        b = -b;
    return m.next_vertex * 2 + m.p;
  }
*/
  int legalMoves (Move moves [MaxLegalMoves]) {
    int next_vertex = 0, best = MaxPlayoutLength + 1;

    for (int i = 0; i < neighbours.size (); i++) {
      if (!check [i] && !tried [i])
	if (free [i] < best) {
	  next_vertex = i;
	  best = free [i];
	}
    }

    if (best == MaxPlayoutLength + 1)
      return 0;

    if (check[next_vertex]==1)
      {
          moves [0].next_vertex=next_vertex;
          moves [0].p = 0;
          return 1;
      }
      else
      {
        moves [0].next_vertex=next_vertex;
        moves [1].next_vertex=next_vertex;

        moves[0].p=0;
        moves[1].p=1;
        return 2;
      }
  }

  void play(Move move) {
    tried [move.next_vertex] = 1;
    if (move.p == 1) {
      current_set.insert (move.next_vertex);
      check [move.next_vertex] = 1;
      for (int i=0 ; i<neighbours[move.next_vertex].size() ; i++) {
	int m=neighbours[move.next_vertex][i];
	if (check [m] == 0) {
	  check [m] = 1;
	  if (useFree)
	    for (int j = 0; j < neighbours [m].size (); j++)
	      free [neighbours [m] [j]]--;
	}
      }
      l++;
    }
    rollout [length] = move;
    length++;
  }

  bool terminal ()
    {
      Move moves [MaxLegalMoves];
      int nb = legalMoves (moves);
      return nb == 0;
    }

  int score ()
  {
    /**/
    bool seen = true;
    while (seen) {
      int mini = 1000000;
      int best;
      for (int i = 0; i < neighbours.size (); i++) {
	if (!check [i])
	  if (free [i] < mini) {
	    mini = free [i];
	    best = i;
	  }
      }
      seen = false;
      if (mini < 1000000) {
	seen = true;
	Move m;
	m.next_vertex = best;
	m.p = 1;
	play (m);
      }
    }
    /**/
    /*
    for (int i = 0; i < neighbours.size (); i++) {
      if (!check [i]) {
	Move m;
	m.next_vertex = i;
	m.p = 1;
	play (m);
      }
    }
    /**/
    return l;
  }

  void print_vertices()
  {
    set <int>::iterator itr;
    for (itr = current_set.begin(); itr != current_set.end(); ++itr)
      {
	cout<<*itr<<endl;
      }
    cout<<endl;
  }

  void print (FILE * fp)
  {
    set <int>::iterator itr;
    for (itr = current_set.begin(); itr != current_set.end(); ++itr)
      {
	fprintf (fp, "%d,", *itr);
      }
    fprintf (fp, "\n");
  }


};


//#include "nestedSimple.c"
#include "gnrpa.pol.opt.cpp"
//#include "gnrpa.pol.opt.beam.c"
//#include "gnrpa.pol.opt.diversity.c"
//#include "gnrpa.pol.opt.stuck.c"
//#include "gnrpa.set.pol.opt.beam.stuck.parallel.c"


int main(int argc, char *argv [])
{
  //graph file
  char * nom1 = argv [1];
  //probability file
  char * nom2 = argv [2];

  ifstream infile(nom1);

  int n ,m , NbVertices,NbEdges;

  infile >> NbVertices;
  cout<<"Graph has "<<NbVertices<<" vertices."<<endl;
  infile >> NbEdges;
  cout<<"Graph has "<<NbEdges<<" edges."<<endl;

  for (int i=0; i<NbVertices; i++) {
    vector<int> neighbor;
    neighbours.push_back(neighbor);
  }

  //Find neighbours
  for(int i=0; i<NbEdges; i++)
    {
      infile >> n;
      infile >> m;

      neighbours[n].push_back(m);
      neighbours[m].push_back(n);
    }

  maxNeighbours = 0;
  for (int i=0; i<NbVertices; i++) {
    if (neighbours [i].size () > maxNeighbours)
      maxNeighbours = neighbours [i].size ();
  }


  ifstream fichier(nom2, ios::in);

  for (int i=0 ; i<NbVertices ; i++)
   {

    float p;
    fichier >> p;
    //std::cout << "Mon entier vaut : " << p << std::endl;
    proba.push_back(p);

    }



  /**/
  printLevel = 1;
  P = 4;
  Policy pol;
  nrpa (3, pol);

  /**/

  float score;
  bestBoard.print_vertices();

 

  return 0;
}


