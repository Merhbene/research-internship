#include <float.h>

float ALPHA = 1.0;
int N = 100;
//int N [10] = {20, 20, 20, 20, 20, 20, 20, 20, 50, 50};

bool reinit = false;
float regul = 0.0;
float r = 0.0;
float dropout = 0.0;
bool readapt = false;
float k = 0.0;
int P = 1;
float tau = 1.0;
int limitRestart = 100;

const int MaxLevel = 6;


float scoreBestRollout [10];
int lengthBestRollout [10];
Move bestRollout [10] [MaxPlayoutLength];
int bestCodeBestRollout [10] [MaxPlayoutLength];
int nbMovesBestRollout [10] [MaxPlayoutLength];
int codeBestRollout [10] [MaxPlayoutLength] [MaxLegalMoves];
int betaBestRollout [10] [MaxPlayoutLength] [MaxLegalMoves];

class Sequence {
 public:
  float scoreBestRollout;
  int lengthBestRollout;
  Move bestRollout [MaxPlayoutLength];
  int bestCodeBestRollout [MaxPlayoutLength];
  int nbMovesBestRollout [MaxPlayoutLength];
  int codeBestRollout [MaxPlayoutLength] [MaxLegalMoves];
  int betaBestRollout [MaxPlayoutLength] [MaxLegalMoves];

  bool operator== (const Sequence & seq) {
    if (lengthBestRollout != seq.lengthBestRollout)
      return false;
    for (int i = 0; i < lengthBestRollout; i++)
      if (bestCodeBestRollout [i] != seq.bestCodeBestRollout [i])
	return false;
    return true;
  }

  float score () {
    return scoreBestRollout;
  }
};

Board bestBoard;

clock_t startClockNRPA, stopClockNRPA;
float bestScoreNRPA = -DBL_MAX;
float nextTimeNRPA = 0.01;
int indexTimeNRPA;
int indexSearch;
float valueAfterTimeNRPA [10000] [100];
float sumValueAfterTimeNRPA [100];
int nbSearchTimeNRPA [100];

int nbTimesNRPA = 14;
bool stopOnTime = false;
float firstTimeNRPA = 0.01;
int nbSearchesNRPA = 200;

int SizeBeam [MaxLevel] = {1, 1, 1, 1, 1, 1};
int startLearning = 0;
int printLevel = 2;
float constante = 0.4;
float kAMAF = 1.0;
float minNorm = -0.1, maxNorm = 2.0;

#include <vector>

//static union{float d; struct{int j, i;} n;} e;
//#define A (1048576/M_LN2)
//#define C 60801
//#define expon(y) (e.n.i=A*(y)+(1072693248-C),e.d)

float expon (float x) {
  return exp (x);
  //union { float f; int i; } y;
  //y.i = (int)(x * 0xB5645F + 0x3F7893F5);
  //return (y.f);
}

using namespace std;

//const int SizeTablePolicy = 65535;
const int SizeTablePolicy = 511;

class MoveStat {
 public:
  int code;
  float sumScores;
  int nbPlayouts;
};

class PolicyAMAF {
 public:
  vector<MoveStat> table [SizeTablePolicy + 1];
  float maxScore, minScore;

  PolicyAMAF () {
    maxScore = 0.0;
    minScore = DBL_MAX;
  }

  void add (int code, float score) {
    if (score > maxScore)
      maxScore = score;
    if (score < minScore)
      minScore = score;
    int index = code & SizeTablePolicy;
    for (int i = 0; i < table [index].size (); i++) {
      if (table [index] [i].code == code) {
	table [index] [i].sumScores += score;
	//if (score > table [index] [i].sumScores)
	//table [index] [i].sumScores = score;
	table [index] [i].nbPlayouts++;
	return;
      }
    }
    MoveStat p;
    p.code = code;
    p.sumScores = score;
    p.nbPlayouts = 1;
    table [index].push_back (p);
  }

  float get (int code) {
    int index = code & SizeTablePolicy;
    for (int i = 0; i < table [index].size (); i++)
      if (table [index] [i].code == code) {
	return (table [index] [i].sumScores / table [index] [i].nbPlayouts - minScore) / (maxScore - minScore);
	//return table [index] [i].sumScores / maxScore;
      }
    return 0.0;
  }

};

PolicyAMAF policyAMAF;

class ProbabilityCode {
 public:
  int code;
  float proba;
  bool zero;
};

class Policy {
 public:
  vector<ProbabilityCode> table [SizeTablePolicy + 1];

  void set (int code, float proba) {
    int index = code & SizeTablePolicy;
    for (int i = 0; i < table [index].size (); i++) {
      if (table [index] [i].code == code) {
	table [index] [i].proba = proba;
	return;
      }
    }
    ProbabilityCode p;
    p.code = code;
    p.proba = proba;
    p.zero = false;
    table [index].push_back (p);
  }

  float get (int code) {
    int index = code & SizeTablePolicy;
    for (int i = 0; i < table [index].size (); i++)
      if (table [index] [i].code == code) {
	return table [index] [i].proba;
      }
    //return minNorm + (maxNorm - minNorm) * policyAMAF.get (code);
    return 0.0;
  }

  void setZero (int code, bool z) {
    int index = code & SizeTablePolicy;
    for (int i = 0; i < table [index].size (); i++) {
      if (table [index] [i].code == code) {
	table [index] [i].zero = z;
	return;
      }
    }
  }

  bool zero (int code) {
    int index = code & SizeTablePolicy;
    for (int i = 0; i < table [index].size (); i++)
      if (table [index] [i].code == code) {
	return table [index] [i].zero;
      }
    return false;
  }

  bool element (int code) {
    int index = code & SizeTablePolicy;
    for (int i = 0; i < table [index].size (); i++)
      if (table [index] [i].code == code) {
	return true;
      }
    return false;
  }

  void add (Policy & pol) {
    for (int index = 0; index < SizeTablePolicy + 1; index++)
      for (int i = 0; i < pol.table [index].size (); i++) {
	int code = pol.table [index] [i].code;
	bool present = false;
	for (int j = 0; j < table [index].size (); j++) {
	  if (table [index] [j].code == code) {
	    present = true;
	    table [index] [j].proba = pol.table [index] [i].proba;
	    break;
	  }
	}
	if (!present)
	  table [index].push_back (pol.table [index] [i]);
      }
  }

  float maximum () {
    float m = -DBL_MAX;
    bool empty = true;
    for (int index = 0; index < SizeTablePolicy + 1; index++)
      for (int i = 0; i < table [index].size (); i++)
	if (table [index] [i].proba > m) {
	  empty = false;
	  m = table [index] [i].proba;
	}
    if (empty)
      return 0.0;
    return m;
  }

  float minimum () {
    float m = DBL_MAX;
    bool empty = true;
    for (int index = 0; index < SizeTablePolicy + 1; index++)
      for (int i = 0; i < table [index].size (); i++)
	if (table [index] [i].proba < m) {
	  empty = false;
	  m = table [index] [i].proba;
	}
    if (empty)
      return 0.0;
    return m;
  }

  void normalize (float mi, float ma) {
    float maxi = maximum (), mini = minimum ();
    for (int index = 0; index < SizeTablePolicy + 1; index++)
      for (int i = 0; i < table [index].size (); i++)
	table [index] [i].proba =
	  mi + (ma - mi) * (table [index] [i].proba - mini) / (maxi - mini);
  }

  void regularize () {
    for (int index = 0; index < SizeTablePolicy + 1; index++)
      for (int i = 0; i < table [index].size (); i++)
	table [index] [i].proba -= r * table [index] [i].proba;
  }
};


float playoutNRPA (Board & board, Policy & pol) {
  int nbMoves = 0;
  Move moves [MaxLegalMoves];
  float probaMove [MaxLegalMoves];

  while (true) {
    if (board.terminal ()) {
      float score = board.score ();
      scoreBestRollout [0] = score;
      lengthBestRollout [0] = board.length;
      for (int k = 0; k < board.length; k++)
	bestRollout [0] [k] = board.rollout [k];
      //for (int k = 0; k < board.length; k++)
      //policyAMAF.add (board.code (board.rollout [k]), score);
      if (score > bestScoreNRPA) {
	bestScoreNRPA = score;
	bestBoard = board;
	//board.print (stderr);
      }
      stopClockNRPA = clock ();
      float time = ((float)(stopClockNRPA - startClockNRPA)) / CLOCKS_PER_SEC;
      if (time > nextTimeNRPA) {
	while (time > 2 * nextTimeNRPA) {
	  indexTimeNRPA++;
	  nextTimeNRPA *= 2;
	}
	valueAfterTimeNRPA [indexSearch] [indexTimeNRPA] = bestScoreNRPA;
	sumValueAfterTimeNRPA [indexTimeNRPA] += bestScoreNRPA;
	nbSearchTimeNRPA [indexTimeNRPA]++;
	indexTimeNRPA++;
	nextTimeNRPA *= 2;
      }
      return score;
    }
    nbMoves = board.legalMoves (moves);
    nbMovesBestRollout [0] [board.length] = nbMoves;
    for (int i = 0; i < nbMoves; i++) {
      float b = 0.0;
      int c = board.code (moves [i], b);
      /**/
      //float p = board.penalty (moves [i]);
      //probaMove [i] = expon (pol.get (c) + moves [i].penalty);
      //if (moves [i].penalty == -1000.0)
      //fprintf (stderr, "tabu move, probaMove [%d/%d] = %lf, ", i, nbMoves, probaMove [i]);
      //if (!pol.element (c))
      //probaMove [i] = expon (-p);
      //else
      //probaMove [i] = expon (pol.get (c));
      /**/
      probaMove [i] = expon (pol.get (c) / tau + b);
      codeBestRollout [0] [board.length] [i] = c;
      betaBestRollout [0] [board.length] [i] = b;
    }

    float sum = probaMove [0];
    for (int i = 1; i < nbMoves; i++)
      sum += probaMove [i];
    float r = (rand () / (RAND_MAX + 1.0)) * sum;
    int j = 0;
    float s = probaMove [0];
    while ((s < r) && (j < nbMoves - 1)) {
      j++;
      s += probaMove [j];
    }
    bestCodeBestRollout [0] [board.length] = codeBestRollout [0] [board.length] [j];
    board.play (moves [j]);
  }
  return 0.0;
}

float probaMove [10] [MaxPlayoutLength] [MaxLegalMoves];
float z [MaxPlayoutLength];

void adaptLevel (int length, int level, Policy & pol) {
  for (int i = 0; i < length; i++) {
    z [i] = 0.0;
    for (int j = 0; j < nbMovesBestRollout [level] [i]; j++) {
      float b = betaBestRollout [level] [i] [j];
      int c = codeBestRollout [level] [i] [j];
      probaMove [level] [i] [j] = expon (pol.get (c) / tau + b);
      z [i] += probaMove [level] [i] [j];
    }
  }
  for (int i = 0; i < length; i++) {
    pol.set (bestCodeBestRollout [level] [i], pol.get (bestCodeBestRollout [level] [i]) + ALPHA / tau);
    for  (int j = 0; j < nbMovesBestRollout [level] [i]; j++)
      pol.set (codeBestRollout [level] [i] [j], pol.get (codeBestRollout [level] [i] [j]) - (ALPHA / tau) * probaMove [level] [i] [j] / z [i]);
    //polAdapt.regularize ();
  }
  //pol.regularize ();
}
/**/

const int MaxBeam = 10;

class Beam {
 public:
  int nb, maxi;
  Sequence seq [MaxBeam];
  Policy pol [MaxBeam];

  Beam (int m = 1) {
    nb = 0;
    maxi = m;
  }

  void init (int m = 1) {
    nb = 0;
    maxi = m;
  }

  void add (Sequence & s, Policy & p) {
    for (int k = 0; k < nb; k++)
      if (seq [k] == s)
	return;
    if (nb < maxi) {
      if (nb == 0) {
	seq [0] = s;
	pol [0] = p;
	nb++;
      }
      else {
	for (int i = nb - 1; i >=0; i--)
	  if (s.score () > seq [i].score ()) {
	    seq [i + 1] = seq [i];
	    pol [i + 1] = pol [i];
	    if (i == 0) {
	      seq [0] = s;
	      pol [0] = p;
	    }
	  }
	  else {
	    seq [i + 1] = s;
	    pol [i + 1] = p;
	    break;
	  }
	nb++;
      }
    }
    else {
      for (int i = maxi - 1; i >=0; i--)
	if (s.score () > seq [i].score ()) {
	  if (i < maxi - 1) {
	    seq [i + 1] = seq [i];
	    pol [i + 1] = pol [i];
	  }
	  if (i == 0) {
	    seq [0] = s;
	    pol [0] = p;
	  }
	}
	else {
	  if (i < maxi - 1) {
	    seq [i + 1] = s;
	    pol [i + 1] = p;
	  }
	  break;
	}
    }
  }
};

static Policy polLevel [10];

float getAlpha (int last, int current) {
  // return 1.0; // nrpa
  // if (last == current) // edelkamp
  // return 1.0;
  // else
  // return 0.0;
  return ALPHA * expon (-k * (current - last));
}

int stepCopyPolicy = 0;

float nrpa(int level, Policy & pol) {
  scoreBestRollout [level] = -DBL_MAX;
  if (level == 0) {
    //polLevel [level].add (pol);
    //return playoutNRPA (polLevel [level]);
    Board b;
    float s = playoutNRPA (b, pol);
    //fprintf (stderr, "%2.3f ", s);
    //b.print_vertices();
    return s;
  }
  else {
    polLevel [level] = pol;
    //polLevel [level].normalize (pol.minimum (), pol.maximum ());
    //polLevel [level].normalize (pol.minimum () / minNorm, pol.maximum () / minNorm);
    //polLevel [level].normalize (polLevel [level].minimum () / minNorm, polLevel [level].maximum () / minNorm);

    //polLevel [level].normalize (minNorm, maxNorm);
    //polLevel [level].add (pol);
    int last = 0;
    int n = N;
    int Level = 1;
    if (level <= Level)
      n = P * N;
    for (int i = 0; i < n; i++) {
      float score = nrpa (level - 1, polLevel [level]);
      //if ((i == stepCopyPolicy) && (level == 3))
      //polLevel [level] = polLevel [level - 1];
      if (score >= scoreBestRollout [level]) {

	if (readapt) {
	  if (score >= scoreBestRollout [level]) {
	    last = i;
	    scoreBestRollout [level] = score;
	    lengthBestRollout [level] = lengthBestRollout [level - 1];
	    for (int k = 0; k < lengthBestRollout [level]; k++)
	      bestRollout [level] [k] = bestRollout [level - 1] [k];
	    for (int k = 0; k < lengthBestRollout [level]; k++) {
	      bestCodeBestRollout [level] [k] = bestCodeBestRollout [level - 1] [k];
	      nbMovesBestRollout [level] [k] = nbMovesBestRollout [level - 1] [k];
	      for (int l = 0; l < nbMovesBestRollout [level - 1] [k]; l++) {
		codeBestRollout [level] [k] [l] = codeBestRollout [level - 1] [k] [l];
		betaBestRollout [level] [k] [l] = betaBestRollout [level - 1] [k] [l];
	      }
	    }
	    if (score > scoreBestRollout [level]) {
	      if (level > 1) {
		last = i;
		int delta = 1;
		if (level == 1)
		  delta = 1;
		polLevel [level] = pol;
		ALPHA = 5.0;
		adaptLevel (lengthBestRollout [level], level, polLevel [level]);
		if (false)
		  for (int k = 0; k < i ; k += delta) {
		    ALPHA = 1.0;
		    adaptLevel (lengthBestRollout [level], level, polLevel [level]);
		  }
		ALPHA = 1.0;
	      }
	    }
	  }
	}
	else {
	  last = i;
	  scoreBestRollout [level] = score;
	  lengthBestRollout [level] = lengthBestRollout [level - 1];
	  for (int k = 0; k < lengthBestRollout [level]; k++)
	    bestRollout [level] [k] = bestRollout [level - 1] [k];
	  for (int k = 0; k < lengthBestRollout [level]; k++) {
	    bestCodeBestRollout [level] [k] = bestCodeBestRollout [level - 1] [k];
	    nbMovesBestRollout [level] [k] = nbMovesBestRollout [level - 1] [k];
	    for (int l = 0; l < nbMovesBestRollout [level - 1] [k]; l++) {
	      codeBestRollout [level] [k] [l] = codeBestRollout [level - 1] [k] [l];
	      betaBestRollout [level] [k] [l] = betaBestRollout [level - 1] [k] [l];
	    }
	  }
	}
	if (level > printLevel) {
	  for (int t = 0; t < level - 1; t++)
	    fprintf (stderr, "\t");
	  fprintf(stderr,"Level : %d, N:%d, score : %f\n", level, i, score);
	}
      }

      if ((level > Level) || ((level <= Level) && (i % P == P - 1))) {
	adaptLevel (lengthBestRollout [level], level, polLevel [level]);
      }
      //if (level > 2)
      //polLevel [level].add (polLevel [level - 1]);
	//polLevel [level] = polLevel [level - 1];
      //alpha *= 0.998;
      if (stopOnTime && (indexTimeNRPA > nbTimesNRPA))
	   {   //b.print_vertices();
	       return scoreBestRollout [level];}
    }
    //b.print_vertices();
    return scoreBestRollout [level];

  }
}

void writeValues (int nbThreads, const char * prefix) {
  char s [1000];
  //sprintf (s, "NRPA.time.nbThreads=%d.nbSearches=%d.k=%2.3f.plot", nbThreads, nbSearchesNRPA, k);
  sprintf (s, "%s.nbThreads=%d.nbSearches=%d.plot", prefix, nbThreads, nbSearchesNRPA);
  FILE * fp = fopen (s, "w");
  if (fp != NULL) {
    fprintf (fp, "# %d searches\n", indexSearch + 1);
    float t = firstTimeNRPA;
    for (int j = 0; j <= nbTimesNRPA; j++) {
      float sum = 0.0;
      int nbValues = 0;
      for (int l = 0; l < indexSearch + 1; l = l + nbThreads)
	if (l + nbThreads - 1 < indexSearch + 1) {
	  float maxi = -DBL_MAX;
	  for (int m = 0; m < nbThreads; m++)
	    if (valueAfterTimeNRPA [l + m] [j] > maxi)
	      maxi = valueAfterTimeNRPA [l + m] [j];
	  sum += maxi;
	  nbValues++;
	}
      fprintf (fp, "%f %f\n", t, sum / nbValues);
      t *= 2;
    }
  }
  fclose (fp);
}

void writeAllValues (int index, const char * prefix) {
  char s [1000];
  //sprintf (s, "NRPA.time.nbThreads=%d.nbSearches=%d.k=%2.3f.plot", nbThreads, nbSearchesNRPA, k);
  sprintf (s, "%s.AllValues.nbTimes.%d.nbSearches.%d.plot", prefix, nbTimesNRPA, nbSearchesNRPA);
  FILE * fp = fopen (s, "w");
  if (fp != NULL) {
    fprintf (fp, "%d\n", index + 1);
    fprintf (fp, "%d\n", nbTimesNRPA);
    for (int l = 0; l < index + 1; l++) {
      float t = firstTimeNRPA;
      for (int j = 0; j <= nbTimesNRPA; j++) {
	fprintf (fp, "%f %f\n", t, valueAfterTimeNRPA [l] [j]);
	t *= 2;
      }
    }
  }
  fclose (fp);
}

void testTimeNRPA (int level, const char * prefix) {
  char s [1000];
  for (int i = 0; i < 100; i++) {
    for (int j = 0; j < 10000; j++)
      valueAfterTimeNRPA [j] [i] = 0.0;
    sumValueAfterTimeNRPA [i] = 0.0;
    nbSearchTimeNRPA [i] = 0;
  }
  stopOnTime = true;
  for (indexSearch = 0; indexSearch < nbSearchesNRPA; indexSearch++) {
    bestScoreNRPA = -DBL_MAX;
    nextTimeNRPA = firstTimeNRPA;
    indexTimeNRPA = 0;
    startClockNRPA = clock ();
    while (true) {
      Policy pol;
      for (int i = 0; i <= level; i++)
	polLevel [i] = pol;
      nrpa (level, pol);
      if (indexTimeNRPA > nbTimesNRPA)
	break;
    }
    float t = firstTimeNRPA;
    fprintf (stderr, "level %d, iteration %d\n", level, indexSearch + 1);
    for (int j = 0; j <= nbTimesNRPA; j++) {
      if (nbSearchTimeNRPA [j] >= (7 * (indexSearch + 1)) / 10)
	fprintf (stderr, "%f %f\n", t, sumValueAfterTimeNRPA [j] / nbSearchTimeNRPA [j]);
      t *= 2;
    }
    writeAllValues (indexSearch, prefix);
    //writeValues (1, prefix);
    //writeValues (2);
    //writeValues (4);
    //writeValues (8);
    //writeValues (16);
  }
  /*
  sprintf (s, "NRPA.time.Level=%d.nbSearches=%d.k=%2.3f.alpha=%2.2f.N=%d.plot", level, nbSearchesNRPA, k, alpha, N);
  FILE * fp = fopen (s, "w");
  if (fp != NULL) {
    fprintf (fp, "# %d searches\n", nbSearchesNRPA);
    float t = firstTimeNRPA;
    for (int j = 0; j <= nbTimesNRPA; j++) {
      if (nbSearchTimeNRPA [j] >= (7 * nbSearchesNRPA) / 10)
	fprintf (fp, "%f %f\n", t, sumValueAfterTimeNRPA [j] / nbSearchTimeNRPA [j]);
      t *= 2;
    }
  }
  fclose (fp);
  */
}

double nestedPolicy (Board & board, int n, Policy & pol) {
  int nbMoves = 0;
  Move moves [MaxLegalMoves];

  lengthBestRollout [n] = -1;
  scoreBestRollout [n] = -DBL_MAX;
  float res;
  while (true) {
    if (board.terminal ())
      return 0.0;
      //return board.score ();
    nbMoves = board.legalMoves (moves);
    for (int i = 0; i < nbMoves; i++) {
      Board b = board;
      double score;
      if (n == 1) {
        b.play (moves [i]);
        score = playoutNRPA (b, pol);
      }
      else {
        b.play (moves [i]);
        nestedPolicy (b, n - 1, pol);
        score = b.score ();
      }
      if (score > scoreBestRollout [n]) {
        scoreBestRollout [n] = score;
        lengthBestRollout [n] = b.length;
        for (int k = 0; k < b.length; k++)
          bestRollout [n] [k] = b.rollout [k];
        if (n > printLevel) {
          for (int t = 0; t < n - 1; t++)
            fprintf (stderr, "\t");
          fprintf (stderr, "n = %d, progress = %d, score = %f\n", n, board.length, scoreBestRollout [n]);
          int depth = 0;
          //b.print (stderr);
          //fprintf (stderr, "\n");
        }
      }
      if (stopOnTime && (indexTimeNRPA > nbTimesNRPA))
	return scoreBestRollout [n];
    }
    board.play (bestRollout [n] [board.length]);
  }
  return 0.0;
}
