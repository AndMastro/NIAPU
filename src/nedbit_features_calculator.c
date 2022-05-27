#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

/* .....................................
   .....GENERATORE DI NUMERI RANDOM.....
   .(da usare solo con interi a 64 bit).
   ..................................... */
/* ..................................... */
#ifdef UNIX
#define geneRand drand48()
int seed = 34652, startseed;
#else
#define Max_Rand 2147483647L
#define F_Rand 4.656612875e-10
unsigned long long seed = 839264236, startseed;
#define geneRand ((double)(seed = (16807L * seed) % Max_Rand) * F_Rand)
#endif

typedef struct Node{
  int id;
  char *name;
  int class;
  double score;
} Node;

typedef struct Link{
  int node1;
  int node2;
  double weight;
  struct Link *next;
} Link;

void PutInLinkList(Link **testa, Link **coda, Link *elemento) {
  elemento->next = NULL;
  if (!(*testa)) {
    *testa = elemento;
    *coda = elemento;
  }
  else {
    (*coda)->next = elemento;
    *coda = elemento;
  }
}

int NotExist(Link *testa, int n1, int n2) {
  int ret;
  Link *app;

  ret = 1;
  app = testa;
  while (app && ret) {
    if (app->node1 == n1) {
      if (app->node2 == n2) {
	ret = 0;
      }
    }
    else if (app->node1 == n2) {
      if (app->node2 == n1) {
	ret = 0;
      }
    }
    app = app->next;
  }
  return(ret);
} 

#define MAXINPUTLINE 256

int nnodi;
int nlink;
Link *ReadRegularLink(char *filename) {
  int l1, l2;
  Link *newlinkhead, *newlinktail, *elem;
  char stbuf[MAXINPUTLINE];
  FILE *fp;

  if ((fp = fopen(filename, "r")) == NULL) {
    fprintf(stderr,"[%s]: can't open %s\n",__FUNCTION__,filename);
    exit(EXIT_FAILURE);
  }

  nlink = 0;
  nnodi = 0;
  newlinkhead = newlinktail = NULL;
  while (fgets(stbuf, MAXINPUTLINE, fp)) {
    if (sscanf(stbuf, "%d %d", &l1, &l2) == 2) {
      if (nnodi < l1) nnodi = l1;
      if (nnodi < l2) nnodi = l2;
      elem = (Link *)malloc(sizeof(struct Link));
      elem->node1 = l1;
      elem->node2 = l2;
      PutInLinkList(&newlinkhead, &newlinktail, elem);
      nlink++;
    }
    else {
      fprintf(stderr, "[%s]: ill-formed link %s\n",__FUNCTION__, stbuf); 
      exit(EXIT_FAILURE);
    }
  }
  fclose(fp);
  nnodi++;
  return(newlinkhead);
}

int *hash_idgene;
Link *createPositiveList(Link *link, Node *genes) {
  Link *newlinkhead, *newlinktail, *elem;

  newlinkhead = newlinktail = NULL;
  // crea field per la connessione
  while (link) {
    if ((genes[hash_idgene[link->node1]].class == 1) && 
	(genes[hash_idgene[link->node2]].class == 1)) {
      elem = (Link *)malloc(sizeof(struct Link));
      elem->node1 = link->node1;
      elem->node2 = link->node2;
      PutInLinkList(&newlinkhead, &newlinktail, elem);
    }
    link = link->next;
  }
  return(newlinkhead);
}

int ngenes;
int nseedgenes;
double totscore;
Node *ReadRegularGenes(char *filename) {
  int count, id;
  double score;
  Node *node;
  char stbuf[MAXINPUTLINE], gene[MAXINPUTLINE];
  FILE *fp;

  /* open file for reading */
  if( (fp = fopen(filename,"r")) == NULL ) {
    fprintf(stderr,"[%s]: can't open %s\n",__FUNCTION__,filename);
    exit(EXIT_FAILURE);
  }

  /* count the number of elements */
  count = 0;
  while (fgets(stbuf, MAXINPUTLINE, fp)) {
    if (sscanf(stbuf, "%d %s %lf", &id, gene, &score) != 3) {
      fprintf(stderr, "[%s]: ill-formed gene field %s\n",__FUNCTION__, stbuf); 
      exit(EXIT_FAILURE);
    }
    else {
      count++;
    }
  }
  fprintf(stderr, "%d sono le linee di %s\n", count, filename);

  node = (Node *)malloc(count * sizeof(Node));
  memset(node, 0, count * sizeof(Node));
  hash_idgene = (int *)malloc(count * sizeof(int));
  memset(hash_idgene, 0, count * sizeof(int));
  nseedgenes = 0;
  totscore = 0.0;
  count = 0;
  rewind(fp);
  while (fgets(stbuf, MAXINPUTLINE, fp)) {
    sscanf(stbuf, "%d %s %lf", &id, gene, &score);
    hash_idgene[id] = count;
    node[count].name = strdup(gene);
    node[count].id = id;
    node[count].score = score;
    if (node[count].score > 0.0) {
      node[count].class = 1;
      totscore += score;
      nseedgenes++;
    }
    else node[count].class = 0;
    count++;
  }
  fclose(fp);
  ngenes = count;
  return(node);
}

void cutSeedGenes(Node *genes, double percentage) {
  int i, elemtocut, icut, elem, *tocut;

  elemtocut = nseedgenes - (int)(round(percentage * nseedgenes));
  if (elemtocut > 0) {
    tocut = (int *)malloc(nseedgenes * sizeof(int));
    memset(tocut, 0, nseedgenes * sizeof(int));
    for (i = 0; i < elemtocut; i++) {
      do {
	elem = (int)(geneRand * nseedgenes);
      } while (tocut[elem]);
      tocut[elem] = 1;
    }
    icut = 0;
    for (i = 0; i < ngenes; i++) {
      if (genes[i].class == 1) {
	if (tocut[icut]) genes[i].class = -1;
	icut++;
      }
    }
    free(tocut);
  }
}

void PrintLink(Link *present) {
  printf("%d--%d\n", present->node1, present->node2);
}

int Connected(Link *link) {
  int *field, dimclust, change;
  Link *elem;

  // crea field per la connessione
  field = (int *)malloc(nnodi * sizeof(int));
  memset(field, 0, nnodi * sizeof(int));
  field[0] = 1;
  dimclust = 1;
  do {
    elem = link;
    change = 0;
    while (elem) {
      if (field[elem->node1] ^ field[elem->node2]) {
	field[elem->node1] = 1;
	field[elem->node2] = 1;
	dimclust++;
	change = 1;
      }
      elem = elem->next;
    }
  } while (change);
  free(field);
  return(dimclust);
}

int *Clusters(Link *link, int *nc) {
  int i, *field, nclust, dimclust, sumclust, change;
  Link *elem;

  // crea field per la connessione
  field = (int *)malloc(nnodi * sizeof(int));
  memset(field, 0, nnodi * sizeof(int));
  nclust = 1;
  field[0] = nclust;
  dimclust = 1;
  sumclust = 0;
  do {
    do {
      elem = link;
      change = 0;
      while (elem) {
	if (((field[elem->node1] == nclust) && (!field[elem->node2])) || 
	    ((field[elem->node2] == nclust) && (!field[elem->node1]))) {
	  field[elem->node1] = nclust;
	  field[elem->node2] = nclust;
	  dimclust++;
	  change = 1;
	}
	elem = elem->next;
      }
    } while (change);
    sumclust += dimclust;
    if (sumclust < nnodi) {
      nclust++;
      dimclust = 1;
      i = 0;
      while (field[i]) i++;
      field[i] = nclust;
    }
  } while (sumclust < nnodi);
  *nc = nclust;
  return(field);
}

int *ClustersSeedGenes(Link *link, Node *genes, int *nc) {
  int i, *field, nclust, dimclust, sumclust, change;
  Link *elem;

  // crea field per la connessione
  field = (int *)malloc(nnodi * sizeof(int));
  memset(field, 0, nnodi * sizeof(int));
  nclust = 1;
  i = 0;
  while (genes[hash_idgene[i]].class != 1) i++;
  field[i] = nclust;
  dimclust = 1;
  sumclust = 0;
  do {
    do {
      elem = link;
      change = 0;
      while (elem) {
	if (((field[elem->node1] == nclust) && (!field[elem->node2])) || 
	    ((field[elem->node2] == nclust) && (!field[elem->node1]))) {
	  field[elem->node1] = nclust;
	  field[elem->node2] = nclust;
	  dimclust++;
	  change = 1;
	}
	elem = elem->next;
      }
    } while (change);
    sumclust += dimclust;
    nclust++;
    dimclust = 1;
    i = 0;
    while ((i < nnodi)&&((field[i]) ||
			 (genes[hash_idgene[i]].class != 1))) i++;
    if (i < nnodi) field[i] = nclust;
  } while (i < nnodi);
  *nc = nclust;
  return(field);
}


#define score2rank(A) (1.0 - (A) / maxscore)
int *ring, *count;
double *rank, minscore = 0.01, maxscore = 0.33;
double alpha = 0.5;

void netRank(Link *link, Node *genes, int *degree) {
  int i, change, nring, stat;
  double sum;
  Link *elem;

  // crea field per la connessione e il rank
  memset(ring, 0, nnodi * sizeof(int));
  memset(count, 0, nnodi * sizeof(int));
  memset(rank, 0, nnodi * sizeof(double));

  // assign a rank to the seed genes
  elem = link;
  change = 0;
  while (elem) {
    if (genes[hash_idgene[elem->node1]].class == 1) {
      ring[elem->node1] = 1;
      count[elem->node1]++;
      rank[elem->node1] += score2rank(genes[hash_idgene[elem->node2]].score);
    }
    if (genes[hash_idgene[elem->node2]].class == 1) {
      ring[elem->node2] = 1;
      count[elem->node2]++;
      rank[elem->node2] += score2rank(genes[hash_idgene[elem->node1]].score);
    }
    elem = elem->next;
  }
  for (i = 0; i < nnodi; i++) {
    if (genes[hash_idgene[i]].class == 1) {
      fprintf(stderr, "%s %f %f %d = ", genes[hash_idgene[i]].name, 
	      genes[hash_idgene[i]].score, rank[i], count[i]);
      rank[i] = alpha * score2rank(genes[hash_idgene[i]].score) + 
	         (1 - alpha) * rank[i] / count[i];
      fprintf(stderr, "%f\n", rank[i]);
    }
  }

  nring = 1;
  do {
    elem = link;
    change = 0;
    while (elem) {
      if (ring[elem->node1] == nring) {
	if ((ring[elem->node2] == nring + 1) || 
	    (!ring[elem->node2])) {
	  ring[elem->node2] = nring + 1;
	  count[elem->node2]++;
	  rank[elem->node2] += rank[elem->node1] - (nring - 1);
	  change = 1;
	}
      }
      else if (ring[elem->node2] == nring) {
	if ((ring[elem->node1] == nring + 1) || 
	    (!ring[elem->node1])) {
	  ring[elem->node1] = nring + 1;
	  count[elem->node1]++;
	  rank[elem->node1] += rank[elem->node2] - (nring - 1);
	  change = 1;
	}
      }
      elem = elem->next;
    }
    nring++;
    sum = 0;
    stat = 0;
    for (i = 0; i < nnodi; i++) {
      if (ring[i] == nring) {
	rank[i] = (nring - 1) + (rank[i] + degree[i] - count[i]) / degree[i];
	sum += rank[i];
	stat++;
      }
    }
    fprintf(stdout, "%d %f %d\n", nring, sum/stat, stat);
  } while (change);
}

double *weightNS;
void netShort(Link *link, Node *genes) {
  int i, j, ij, ji, k;
  double *d;
  Link *elem;

  fprintf(stderr, "computing initial genes weight\n");
  memset(weightNS, 0, nnodi * sizeof(double));
  for (i = 0; i < nnodi; i++) {
    if (genes[hash_idgene[i]].class == 1) weightNS[i] = genes[hash_idgene[i]].score / maxscore;
    else if (genes[hash_idgene[i]].class == -99) weightNS[i] =  0.0;
    else weightNS[i] =  alpha * minscore / maxscore;
  }
  d = (double *)malloc(nnodi * nnodi * sizeof(double));
  memset(d, 0, nnodi * nnodi * sizeof(double));
  fprintf(stderr, "computing initial distance\n");
  elem = link;
  while (elem) {
//    fprintf(stderr, "%d (%f) %d (%f)\n", elem->node1, weightNS[elem->node1], elem->node2, weightNS[elem->node2]);
    if (weightNS[elem->node1] + weightNS[elem->node2] > 0)
      elem->weight = 2.0 / (weightNS[elem->node1]+weightNS[elem->node2]);
    else elem->weight = DBL_MAX;
    ij = elem->node1 * nnodi + elem->node2;
    ji = elem->node2 * nnodi + elem->node1;
    if ((ij > nnodi * nnodi) || (ji > nnodi * nnodi)) {
      fprintf(stderr, "problema %d o %d maggiore di %d\n", ij, ji, nnodi * nnodi);
      exit(-1);
    }
    d[ij] = elem->weight;
    d[ji] = elem->weight;
//    fprintf(stderr, "(%d) %d, %d = %f\n", nnodi * nnodi, ij, ji, elem->weight);
    elem = elem->next;
  }
  fprintf(stderr, "reset two step distance\n");
  for (i = 0; i < nnodi; i++) {
    for (j = 0; j < nnodi; j++) {
      if ((i != j) && (d[i * nnodi + j] == 0)) d[i * nnodi + j] = DBL_MAX;
    }
  }

  fprintf(stderr, "computing shortest distance\n");
  for (k = 0; k < nnodi; k++) {
    fprintf(stderr, "%d su %d\n", k, nnodi);
    for (i = 0; i < nnodi; i++) {
      for (j = 0; j < nnodi; j++) {
	if (d[i * nnodi + j] > d[i * nnodi + k] + d[k * nnodi + j])
	  d[i * nnodi + j] = d[i * nnodi + k] + d[k * nnodi + j];
      }
    }
  }
  /*
  fprintf(stderr, "print distance\n");
  for (i = 0; i < nnodi; i++) {
    for (j = 0; j < nnodi; j++) {
      fprintf(stdout, "%s %s %f\n", genes[hash_idgene[i]].name, genes[hash_idgene[j]].name, d[i * nnodi + j]);
    }
  }
  */
  memset(weightNS, 0, nnodi * sizeof(double));
  for (i = 0; i < nnodi; i++) {
    for (j = 0; j < nnodi; j++) {
      if (i != j) weightNS[i] += 1.0 / d[i * nnodi + j];
    }
  }
}



double *weight;
void myNetShort(Link *link, Node *genes) {
  int i, j, change, stat;
  double *weightns, appw;
  Link *elem;

  memset(weightNS, 0, nnodi * sizeof(double));
  for (i = 0; i < nnodi; i++) {
    if (genes[hash_idgene[i]].class == 1) weightNS[i] = genes[hash_idgene[i]].score / maxscore;
    else if (genes[hash_idgene[i]].class == -99) weightNS[i] =  0.0;
    else weightNS[i] =  alpha * minscore / maxscore;
  }
  elem = link;
  while (elem) {
    elem->weight = 2.0 / (weightNS[elem->node1]+weightNS[elem->node2]);
    elem = elem->next;
  }
  memset(weightNS, 0, nnodi * sizeof(double));

  weightns = (double *)malloc(nnodi * sizeof(double));
  for (i = 0; i < nnodi; i++) {
    memset(weightns, 0, nnodi * sizeof(double));
    // first step: assign weights to nearest neighbors of i
    elem = link;
    change = 0;
    while (elem) {
      if (elem->node1 == i) {
	weightns[elem->node2] = elem->weight;
	change++;
      }
      else if (elem->node2 == i) {
	weightns[elem->node1] = elem->weight;
	change++;
      }
      elem = elem->next;
    }
    if (change == 0) {
      fprintf(stdout, "isolated node %d\n", i);
    }
    else {
      stat = 0;
      do {
	elem = link;
	change = 0;
	while (elem) {
	  if (weightns[elem->node1] > 0.0) {
	    appw = weightns[elem->node1] + elem->weight;
	    if (((weightns[elem->node2] == 0) || (weightns[elem->node2] > appw))&&(elem->node2 != i)) {
	      weightns[elem->node2] = appw;
	      change++;
	    }
	  }
	  else if (weightns[elem->node2] > 0.0) {
	    appw = weightns[elem->node2] + elem->weight;
	    if (((weightns[elem->node1] == 0) || (weightns[elem->node1] > appw))&&(elem->node1 != i)) {
	      weightns[elem->node1] = appw;
	      change++;
	    }
	  }
	  elem = elem->next;
	}
	stat++;
	fprintf(stdout, "node %d stat %d change %d\n", i, stat, change);
      } while (change > 0);
    }
    for (j = 0; j < nnodi; j++) {
      if (weightns[j] > 0) weightNS[i] += 1.0 / weightns[j];
    }
  }
}


double *diffus;
double wdt;
void diffusionHeat(Link *link, int *degree, double **fieldnew, double **fieldold) {
  int i;
  double *appfield;

  memset(diffus, 0, nnodi * sizeof(double));
  while (link) {
    if (degree[link->node2] != 0)
      diffus[link->node1] += (*fieldold)[link->node2] / degree[link->node2];
    if (degree[link->node1] != 0)
    diffus[link->node2] += (*fieldold)[link->node1] / degree[link->node1];
    link = link->next;
  }

  for (i = 0; i < nnodi; i++) {
    (*fieldnew)[i] = (1 - wdt) * (*fieldold)[i] + diffus[i] * wdt;
  }

  appfield = *fieldnew;
  *fieldnew  = *fieldold;
  *fieldold = appfield;
}

void diffusionInfo(Link *link, int *degree, double **fieldnew, double **fieldold) {
  int i;
  double *appfield;

  memset(diffus, 0, nnodi * sizeof(double));
  while (link) {
    diffus[link->node1] += (*fieldold)[link->node2];
    diffus[link->node2] += (*fieldold)[link->node1];
    link = link->next;
  }

  for (i = 0; i < nnodi; i++) {
    (*fieldnew)[i] = (1 - wdt * degree[i])*(*fieldold)[i] + diffus[i] * wdt;
  }

  appfield = *fieldnew;
  *fieldnew  = *fieldold;
  *fieldold = appfield;
}

int oneRing(Link *present, int *ringnew, int *ringold, int step) {
  int nchanged;
  
  fprintf(stdout, "nuova versione step %d\n", step);
  nchanged = 0;
  while (present) {
    if (ringold[present->node1] == step) {
      ringnew[present->node1] = step;
      if (ringold[present->node2] == 0) {
	ringnew[present->node2] = step + 1;
	nchanged++;
      }
    }
    if (ringold[present->node2] == step) {
      ringnew[present->node2] = step;
      if (ringold[present->node1] == 0) {
	ringnew[present->node1] = step + 1;
	nchanged++;
      }
    }
    present = present->next;
  }
  memcpy(ringold, ringnew, nnodi * sizeof(int));
  return(nchanged);
}

void covDegree(Link *present, int *covdeg, int *ringGene) {
  while (present) {
    if (ringGene[present->node1] == 1) {
      covdeg[present->node2]++;
    }
    else if (ringGene[present->node2] == 1) {
      covdeg[present->node1]++;
    }
    present = present->next;
  }
}

int *computeDegree(Link *present) {
  int *grado;

  grado = (int *)malloc(nnodi * sizeof(int));
  memset(grado, 0, nnodi * sizeof(int));
  while (present) {
    grado[present->node1]++;
    grado[present->node2]++;
    present = present->next;
  }
  return(grado);
}


#define nARG 3
#define ARGfileLink 1
#define ARGfileGene 2
#define ARGfileOut 3

int main (int argc, char *argv[]) {
  Link *linklista, *linkseed;
  Node *genes;
  int *clusters, *elemClus, *degree, i, step, nc, stat;
  double sumdeg, vardeg, coeff;
  double *fieldHeat, *fieldInfo, *app;
  FILE *fw;

  if (argc != nARG + 1) {
    fprintf(stderr, "[%s]: Uso: %s filelink filegene fileout\n",
	    __FUNCTION__, argv[0]);
    exit(1);
  }

  // Seed generatori random
#ifndef FIXEDSEED
#ifdef UNIX
  seed=time(NULL);
  srand48(seed);
  startseed = seed;
  fprintf(stdout, "UNIX-SEME = %d\n", startseed);
#else
  seed=(unsigned long long)time(NULL);
  startseed = seed;
  __mingw_fprintf(stdout, "WINDOWS-SEME = %llu\n", startseed);
#endif
#endif

  linklista = ReadRegularLink(argv[ARGfileLink]);
  genes = ReadRegularGenes(argv[ARGfileGene]);

  //
  // Controllo connessione e stampa geni isolati
  // ...........................................
  fprintf(stdout, "major cluster %d elements over %d\n", 
	  Connected(linklista), nnodi);
  clusters = Clusters(linklista, &nc);
  elemClus = (int *)malloc((nc + 1) * sizeof(int));
  memset(elemClus, 0, (nc + 1) * sizeof(int));
  for (i = 1; i < nnodi; i++) {
    elemClus[clusters[i]]++;
  }
  fprintf(stdout, "clusters:\n");
  for (i = 1; i < nc + 1; i++) {
    fprintf(stdout, "%d %d\n", i, elemClus[i]);
  }

  fprintf(stdout, "singoli:\n");
  for (i = 0; i < nnodi; i++) {
    if (clusters[i] != 1) fprintf(stdout, "%s %d\n", genes[i].name, genes[i].class);
  }
  free(clusters);
  free(elemClus);

  // Blocco per il calcolo del grado e della clusterizzazione 
  // del sottografo dei seed genes
  fprintf(stdout, "compute seed subgraph:\n");
  linkseed = createPositiveList(linklista, genes);
  fprintf(stdout, "clusterization of seed subgraph:\n");
  clusters = ClustersSeedGenes(linkseed, genes, &nc);
  elemClus = (int *)malloc((nc + 1) * sizeof(int));
  memset(elemClus, 0, (nc + 1) * sizeof(int));
  for (i = 1; i < nnodi; i++) {
    elemClus[clusters[i]]++;
  }
  for (i = 1; i < nc + 1; i++) {
    fprintf(stdout, "%d %d\n", i, elemClus[i]);
  }
  fprintf(stdout, "compute degree of seed subgraph:\n");
  stat = 0;
  sumdeg = 0;
  vardeg = 0;
  degree = computeDegree(linkseed);
  for (i = 0; i < nnodi; i++) {
    if (genes[hash_idgene[i]].class == 1) {
      fprintf(stdout, "%d %s %d %d\n", i, genes[hash_idgene[i]].name, 
	      degree[i], clusters[i]); 
    }
  }

  free(degree);

  printf("Parametri del sistema:\n");
  printf("          file link:   %s\n", argv[ARGfileLink]);
  printf("             n link:   %d\n", nlink);
  printf("          file geni:   %s\n", argv[ARGfileGene]);
  printf("             n geni:   %d\n", ngenes);
  printf("        n geni seed:   %d\n", nseedgenes);
  fflush(stdout);

  // Controllo connessione
  fprintf(stdout, "major cluster %d elements over %d\n", 
	  Connected(linklista), nnodi);

  clusters = Clusters(linklista, &nc);
  elemClus = (int *)malloc((nc + 1) * sizeof(int));
  memset(elemClus, 0, (nc + 1) * sizeof(int));

  for (i = 1; i < nnodi; i++) {
    elemClus[clusters[i]]++;
  }

  fprintf(stdout, "clusters:\n");
  for (i = 1; i < nc + 1; i++) {
    fprintf(stdout, "%d %d\n", i, elemClus[i]);
  }

  /*
  fprintf(stdout, "singoli:\n");
  for (i = 0; i < nnodi; i++) {
    if ((clusters[i] != 1)&&(genes[i].class == 1)) {
      fprintf(stdout, "(%d %s %d) %d\n", i + 1, 
	      genes[i].name, genes[i].class, 
	      clusters[i]);
    }
    if (clusters[i] != 1) fprintf(stdout, "%s %d\n", genes[i].name, genes[i].class);
  }
  */


  fprintf(stdout, "compute degree\n");
  stat = 0;
  sumdeg = 0;
  vardeg = 0;
  degree = computeDegree(linklista);
  for (i = 0; i < nnodi; i++) {
//    fprintf(stdout, "%d %s %d %d\n", i, genes[hash_idgene[i]].name, degree[i], clusters[i]); 
    if (clusters[i] == 1) {
      sumdeg += degree[i];
      vardeg += degree[i] * degree[i];
      stat++;
    }
    else {
      genes[hash_idgene[i]].class = -99;
      fprintf(stdout, "%d %s %d\n", i, genes[hash_idgene[i]].name, 
	                               genes[hash_idgene[i]].class); 
    }
  }
  sumdeg /= (double)stat;
  vardeg = sqrt((vardeg - sumdeg * sumdeg * nnodi)/(nnodi - 1.0));
  fprintf(stdout, "nodes %d, average degree: %f, standard deviation: %f\n", stat, sumdeg, vardeg);

  // allocate memory
  ring = (int *)malloc(nnodi * sizeof(int));
  count = (int *)malloc(nnodi * sizeof(int));
  rank = (double *)malloc(nnodi * sizeof(double));
  weightNS = (double *)malloc(nnodi * sizeof(double));
  fieldHeat = (double *)malloc(nnodi * sizeof(double));
  fieldInfo = (double *)malloc(nnodi * sizeof(double));
  diffus = (double *)malloc(nnodi * sizeof(double));
  app = (double *)malloc(nnodi * sizeof(double));

  fprintf(stdout, "net short:\n");
  netShort(linklista, genes);

  fprintf(stdout, "net rank:\n");
  netRank(linklista, genes, degree);

  coeff = (double)nseedgenes / (double)totscore;
  wdt = 0.001;
  step = 263;
  fprintf(stdout, "heat diffusion %f %d\n", wdt, step); 
  memset(fieldHeat, 0, nnodi * sizeof(double));
  memset(app, 0, nnodi * sizeof(double));
  for (i = 0; i < nnodi; i++) {
    if (genes[hash_idgene[i]].class == 1) {
      //fieldHeat[i] = 1.0;
      fieldHeat[i] = coeff * genes[hash_idgene[i]].score;
    }
  }
  for (i = 0; i < step; i++) {
    diffusionHeat(linklista, degree, &app, &fieldHeat);
  }

  wdt = 0.0001;
  step = 50;
  fprintf(stdout, "info diffusion %f %d\n", wdt, step); 
  memset(fieldInfo, 0, nnodi * sizeof(double));
  memset(app, 0, nnodi * sizeof(double));
  for (i = 0; i < nnodi; i++) {
    if (genes[hash_idgene[i]].class == 1) {
//      fieldInfo[i] = 1.0;
      fieldInfo[i] = coeff * genes[hash_idgene[i]].score;
    }
  }
  for (i = 0; i < step; i++) {
    diffusionInfo(linklista, degree, &app, &fieldInfo);
  }

  fw = fopen(argv[ARGfileOut], "w");
  fprintf(fw, "name,class,degree,ring,NetRank,NetShort,HeatDiff,InfoDiff\n");
  for (i = 0; i < nnodi; i++) {
    fprintf(fw, "%s,%d,%d,%d,%e,%e,%e,%e\n", genes[hash_idgene[i]].name, 
	    genes[hash_idgene[i]].class, degree[i], ring[i], 
	    rank[i], weightNS[i], fieldHeat[i], fieldInfo[i]);
  }
  fclose(fw);

  free(ring);
  free(count);
  free(rank);
  free(diffus);

  return(0);
}
