#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define MAXINPUTLINE 1024

void uccomma2space(char *inputstr) {
  while(*inputstr) {
    if ((*inputstr >= 97) && (*inputstr <= 122)) {
      *inputstr -= 32;
    }
    else if (*inputstr == ',') *inputstr = 32;
    inputstr++;
  }
}

char *StrSep(char **inputstr, int ch) {
  char *app, *retstr;

  if ((*inputstr == NULL)||(**inputstr == '\0')) return(NULL);
  while(**inputstr && (**inputstr==' ' || **inputstr=='\t')) (*inputstr)++;
  app = strchr(*inputstr, ch);
  retstr = *inputstr;
  if (app != NULL) {
    *inputstr = app + 1;
    while(**inputstr && (**inputstr==' ' || **inputstr=='\t')) (*inputstr)++;
    *app = '\0';
  }
  else {
    *inputstr = NULL;
  }
  return(retstr);
}

double **createMatrix(int nRow, int nCol) {
  int i;
  double **mat;

  mat = (double **)malloc(nRow * sizeof(double *));
  for (i = 0; i < nRow; i++) {
    mat[i] = (double *)malloc(nCol * sizeof(double));
    memset(mat[i], 0, nCol * sizeof(double));
  }
  return(mat);
}

int **createMatrixInt(int nRow, int nCol) {
  int i, **mat;

  mat = (int **)malloc(nRow * sizeof(int *));
  for (i = 0; i < nRow; i++) {
    mat[i] = (int *)malloc(nCol * sizeof(int));
    memset(mat[i], 0, nCol * sizeof(int));
  }
  return(mat);
}

void freeMatrix(double **mat, int nRow) {
  int i;

  for (i = 0; i < nRow; i++) {
    free(mat[i]);
  }
  free(mat);
}

void freeMatrixInt(int **mat, int nRow) {
  int i;

  for (i = 0; i < nRow; i++) {
    free(mat[i]);
  }
  free(mat);
}

double *createVect(int nElem) {
  double *vec;

  vec = (double *)malloc(nElem * sizeof(double));
  memset(vec, 0, nElem * sizeof(double));
  return(vec);
}

int partition(double *vec, int low, int high) {
  int i, j;
  double mid, app;
  
  mid = vec[high];
  i = low - 1;
  for (j = low; j <= high -1; j++) {
    if (vec[j] < mid) {
      i++;
      if (i != j) {
	app = vec[i];
	vec[i] = vec[j];
	vec[j] = app;
      }
    }
  }
  if ((i + 1) != (high)) {
    app = vec[i + 1];
    vec[i + 1] = vec[high];
    vec[high] = app;
  }
  return(i + 1);
}

void qSort(double *vec, int low, int high) {
  int mid;

  if (low < high) {
    mid = partition(vec, low, high);
    qSort(vec, low, mid - 1);
    qSort(vec, mid + 1, high);
  }
}

void setFeature(double **features, int col, double *x, int nRow) {
  int i;
  double min, max, cof;

  min = x[0];
  max = x[0];
  for (i = 1; i < nRow; i++) {
    if (min > x[i]) min = x[i];
    if (max < x[i]) max = x[i];
  }
  cof = 1.0 / (max - min);
  for (i = 0; i < nRow; i++) features[i][col] = cof * (x[i] - min);
}


#define alpha 0.8
#define threshold 0.000001

#define nARG 5
#define _ARGfileIn_ 1
#define _ARGheader_ 2
#define _ARGfileOut_ 3
#define _sQuantile_ 4
#define _rnQuantile_ 5

int main (int argc, char *argv[]) {
  int i, j, k, nnodi, nfeature, nseed, n1, n2;
  int ie, nRelNeg, nLabel, nRejected, header;
  int *class;
  double *aveFeat, *Wvec, *D, *relNeg, *relNegSort, *label, *appLabel;
  double *G0, *Gt_1, *Gt;
  double sumG0, sumGt;
  char **gene;
  double **feature, **normFeatures, **W, **Wr, dist, min, max, cof;
  double sQuantile, rnQuantile, renorm;
  double soglia, sogliaLNWN, sogliaWNLP;
  char stbuf[MAXINPUTLINE], *app, sep = ',';
  FILE *fr, *fw;

  fprintf(stderr, "Number of arguments: %d\n", argc-1); 
  for (k = 1; k < argc; k++) {
     fprintf(stderr, "%d --> <%s>\n", k  , argv[k]); 
  }

  if (argc != nARG + 1) {
    fprintf(stderr, "Uso: %s fileIn flagHeader fileOut sQuantile rnQuantile\n",
	    argv[0]);
    exit(1);
  }

  /* open file for reading */
  if( (fr = fopen(argv[_ARGfileIn_], "r")) == NULL ) {
    fprintf(stderr,"[%s]: can't open %s\n",__FUNCTION__, argv[_ARGfileIn_]);
    exit(EXIT_FAILURE);
  }

  /* reads header flag */
  header = atoi(argv[_ARGheader_]);

  /* reads quantile for link */
  sQuantile = atof(argv[_sQuantile_]);

  /* reads quantile for Reliable Negative */
  rnQuantile = 1.0 - atof(argv[_rnQuantile_]);

  /* skip the header, if any */
  if (header) fgets(stbuf, MAXINPUTLINE, fr);

  /* count the number of genes and the number of feature */
  nfeature = 0;
  fgets(stbuf, MAXINPUTLINE, fr);
  app = stbuf;
  while (StrSep(&app, sep)) nfeature++;
  nfeature -= 2;

  nnodi = 1;
  while (fgets(stbuf, MAXINPUTLINE, fr)) {
    nnodi++;
  }
  rewind(fr);

  fprintf(stderr, "           file in:   %s", argv[_ARGfileIn_]);
  if (header) fprintf(stderr, " with header\n");
  else fprintf(stderr, " without header\n");
  fprintf(stderr, "          file out:   %s\n", argv[_ARGfileOut_]);
  fprintf(stderr, "        quantile S:   %s\n", argv[_sQuantile_]);
  fprintf(stderr, "       quantile RN:   %s\n", argv[_rnQuantile_]);
  fprintf(stderr, "            ngenes:   %d\n", nnodi);

  // allocate memory
  gene = (char **)malloc(nnodi * sizeof(char *));
  class = (int *)malloc(nnodi * sizeof(int));
  feature = createMatrix(nfeature, nnodi);
  nseed = 0;

  /* ignore the header */
  if (header) fgets(stbuf, MAXINPUTLINE, fr);

  /* reads the features */
  for (i = 0; i < nnodi; i++) {
    fgets(stbuf, MAXINPUTLINE, fr);
    app = stbuf;
    gene[i] = strdup(StrSep(&app, sep));
    class[i] = atoi(StrSep(&app, sep));
    for (j = 0; j < nfeature; j++) feature[j][i] = atof(StrSep(&app, sep));
    nseed += class[i];
  }
  fclose(fr);
  fprintf(stderr, "            nseeds:   %d\n\n", nseed);

  /* debug print*/
  // for (i = 0; i < 10; i++) {
  //   fprintf(stderr, "%s %d", gene[i], class[i]);
  //   for (j = 0; j < nfeature; j++) {
  //     fprintf(stderr, " %f", feature[j][i]);
  //   }
  //   fprintf(stderr, "\n");
  // }

  // Caricamento features
  normFeatures = createMatrix(nnodi, nfeature);
  for (j = 0; j < nfeature; j++) {
    setFeature(normFeatures, j, feature[j], nnodi);
  }

  // Step 1
  fprintf(stderr, "\n---step 1: creation of distance matrix W\n");
  W = createMatrix(nnodi, nnodi);
  aveFeat = createVect(nfeature);
  min = 2.0;
  max = 0.0;
  for (n1 = 0; n1 < nnodi; n1++) {
    if (class[n1]) {
      for (k = 0; k < nfeature; k++) {
	aveFeat[k] += normFeatures[n1][k] / (double)nseed;
      }
    }
    for (n2 = n1 + 1; n2 < nnodi; n2++) {
      dist = 0;
      for (k = 0; k < 4; k++) {
	dist += (normFeatures[n1][k] - normFeatures[n2][k]) *
	        (normFeatures[n1][k] - normFeatures[n2][k]);
      }
      W[n1][n2] = sqrt(dist);
      if (min > W[n1][n2]) min = W[n1][n2];
      else if (max < W[n1][n2]) max = W[n1][n2];
    }
  }

  i = 0;
  ie = nnodi * (nnodi-1) / 2;
  Wvec = createVect(ie);
  cof = 1.0 / (max - min);
  for (n1 = 0; n1 < nnodi; n1++) {
    W[n1][n1] = 1;
    for (n2 = n1 + 1; n2 < nnodi; n2++) {
      W[n1][n2] = 1.0 - cof * (W[n1][n2] - min);
      W[n2][n1] = W[n1][n2];
      Wvec[i++] = W[n1][n2];
    }
  }

  // Step 2
  fprintf(stderr, "\n---step 2: computes reduced normalized matrix Wr\n");
  qSort(Wvec, 0, ie - 1);
  soglia = Wvec[(int)(ie * sQuantile)];
  fprintf(stderr, "   threshold value %f at %d su %d\n", 
	  soglia, (int)(ie * sQuantile), ie);
  nRejected = 0;
  D = createVect(nnodi);
  for (n1 = 0; n1 < nnodi; n1++) {
    D[n1] += 1;
    for (n2 = n1 + 1; n2 < nnodi; n2++) {
      if (W[n1][n2] < soglia) {
	W[n1][n2] = 0;
	W[n2][n1] = 0;
	nRejected += 2;
      }
      else {
	D[n1] += W[n1][n2];
	D[n2] += W[n2][n1];
      }
    }
  }
  fprintf(stderr, "   rejected %d values in the matrix W\n", nRejected);
  free(Wvec);

  Wr = createMatrix(nnodi, nnodi);
  for (n1 = 0; n1 < nnodi; n1++) {
    for (n2 = 0; n2 < nnodi; n2++) {
      Wr[n1][n2] = W[n1][n2] / D[n1];
    }
  }
  freeMatrix(W, nnodi);

  // Step 3
  fprintf(stderr, "\n---step 3: compute reliable negative genes and initial condition\n");
  relNeg = createVect(nnodi);
  relNegSort = createVect(nnodi);
  for (n1 = 0; n1 < nnodi; n1++) {
    if (!class[n1]) {
      for (k = 0; k < 4; k++) {
	relNeg[n1] += (aveFeat[k] - normFeatures[n1][k]) *
	              (aveFeat[k] - normFeatures[n1][k]);
      }
      relNeg[n1] = sqrt(relNeg[n1]);
      relNegSort[n1] = relNeg[n1];
    }
  }
  qSort(relNegSort, 0, nnodi - 1);
  soglia = relNegSort[(int)(nnodi * rnQuantile)];

  nRelNeg = 0;
  G0 = createVect(nnodi);
  for (n1 = 0; n1 < nnodi; n1++) {
    if (class[n1]) {
      G0[n1] = 1;
    }
    else {
      if (relNeg[n1] > soglia) {
	G0[n1] = -1;
	nRelNeg++;
      }
    }
  }
  fprintf(stderr, "   reliable negative gene %d over %d\n", nRelNeg, nnodi);
  renorm = (double)nseed / (double)nRelNeg;
  sumG0 = 0;
  for (n1 = 0; n1 < nnodi; n1++) {
    if (G0[n1] == -1) G0[n1] *= renorm;
    sumG0 += G0[n1];
  }
  fprintf(stderr, "   sum G0 = %f\n", sumG0);
  free(relNeg);
  free(relNegSort);

  // Step 4
  fprintf(stderr, "\n---step 4: markov dynamic with restart\n");
  Gt = createVect(nnodi);
  Gt_1 = createVect(nnodi);
  memcpy(Gt_1, G0, nnodi * sizeof(double));
  do {
    sumGt = 0;
    for (n1 = 0; n1 < nnodi; n1++) {
      Gt[n1] = alpha * G0[n1];
      for (n2 = 0; n2 < nnodi; n2++) {
	Gt[n1] += (1 - alpha) * Wr[n2][n1] * Gt_1[n2];
      }
      sumGt += Gt[n1];
    }
    dist = 0;
    for (n1 = 0; n1 < nnodi; n1++) {
      dist += fabs(Gt[n1] - Gt_1[n1]);
    }
    fprintf(stderr, "   sum Gt = %f,\tdistance |Gt - Gt_1|=%f\n", 
	    sumGt, dist);
    memcpy(Gt_1, Gt, nnodi * sizeof(double));
  } while (dist > threshold);

  // Step 5
  fprintf(stderr, "\n---step 5: labelling procedure\n");
  i = 0;
  nLabel = nnodi - nseed - nRelNeg;
  appLabel = createVect(nLabel);
  for (n1 = 0; n1 < nnodi; n1++) {
    if (G0[n1] == 0) {
      appLabel[i] = Gt[n1];
      i++;
    }
  }
  qSort(appLabel, 0, nLabel - 1);
  sogliaLNWN = appLabel[(int)(0.33333 * nLabel)];
  sogliaWNLP = appLabel[(int)(0.66666 * nLabel)];
  free(appLabel);

  label = createVect(nnodi);
  for (n1 = 0; n1 < nnodi; n1++) {
    if (G0[n1] == 1) label[n1] = 1;
    else if (G0[n1] < 0) label[n1] = 5;
    else if (Gt[n1] < sogliaLNWN) label[n1] = 4;
    else if (Gt[n1] < sogliaWNLP) label[n1] = 3;
    else label[n1] = 2;
  }

  /* open file for writing */
  if( (fw = fopen(argv[_ARGfileOut_], "w")) == NULL ) {
    fprintf(stderr,"[%s]: can't open %s\n",__FUNCTION__, argv[_ARGfileOut_]);
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < nnodi; i++) {
    fprintf(fw, "%s %f %d\n", gene[i], Gt[i], (int)(label[i]));
  }
  fclose(fw);

  return(0);
}
