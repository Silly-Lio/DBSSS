/***********************************************************************/
/************************* DISCLAIMER **********************************/
/***********************************************************************/
/** This UCR Suite software is copyright protected ?2012 by          **/
/** Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,            **/
/** Gustavo Batista and Eamonn Keogh.                                 **/
/**                                                                   **/
/** Unless stated otherwise, all software is provided free of charge. **/
/** As well, all software is provided on an "as is" basis without     **/
/** warranty of any kind, express or implied. Under no circumstances  **/
/** and under no legal theory, whether in tort, contract,or otherwise,**/
/** shall Thanawin Rakthanmanon, Bilson Campana, Abdullah Mueen,      **/
/** Gustavo Batista, or Eamonn Keogh be liable to you or to any other **/
/** person for any indirect, special, incidental, or consequential    **/
/** damages of any character including, without limitation, damages   **/
/** for loss of goodwill, work stoppage, computer failure or          **/
/** malfunction, or for any and all other damages or losses.          **/
/**                                                                   **/
/** If you do not agree with these terms, then you you are advised to **/
/** not use this software.                                            **/
/***********************************************************************/
/***********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

#define INF 1e20       //Pseudo Infitinte number for this code
#define N_MAX_ROWS 200000

using namespace std;

/// Data structure for sorting the query
typedef struct Index
{
    double value;
    int    index;
} Index;

/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void* a, const void* b)
{
    Index* x = (Index*)a;
    Index* y = (Index*)b;
    return (x->value) - (y->value);   // low to high
}

typedef struct IndexDTW {
    int index;
    double dist;
    double dist_1d;
    IndexDTW *next;
}IndexDTW;

struct deque
{
    int* dq;
    int size, capacity;
    int f, r;
};

/// Initial the queue at the begining step of envelop calculation
void init(deque* d, int capacity)
{
    d->capacity = capacity;
    d->size = 0;
    d->dq = (int*)malloc(sizeof(int) * d->capacity);
    d->f = 0;
    d->r = d->capacity - 1;
}

/// Destroy the queue
void destroy(deque* d)
{
    free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct deque* d, int v)
{
    d->dq[d->r] = v;
    d->r--;
    if (d->r < 0)
        d->r = d->capacity - 1;
    d->size++;
}

/// Delete the current (front) element from queue
void pop_front(struct deque* d)
{
    d->f--;
    if (d->f < 0)
        d->f = d->capacity - 1;
    d->size--;
}

/// Delete the last element from queue
void pop_back(struct deque* d)
{
    d->r = (d->r + 1) % d->capacity;
    d->size--;
}

/// Get the value at the current position of the circular queue
int front(struct deque* d)
{
    int aux = d->f - 1;

    if (aux < 0)
        aux = d->capacity - 1;
    return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct deque *d)
int back(struct deque* d)
{
    int aux = (d->r + 1) % d->capacity;
    return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct deque* d)
{
    return d->size == 0;
}
void lower_upper_lemire(double* t, int len, int r, double* l, double* u)
{
    struct deque du, dl;

    init(&du, 2 * r + 2);
    init(&dl, 2 * r + 2);

    push_back(&du, 0);
    push_back(&dl, 0);

    for (int i = 1; i < len; i++)
    {
        if (i > r)
        {
            u[i - r - 1] = t[front(&du)];
            l[i - r - 1] = t[front(&dl)];
        }
        if (t[i] > t[i - 1])
        {
            pop_back(&du);
            while (!empty(&du) && t[i] > t[back(&du)])
                pop_back(&du);
        }
        else
        {
            pop_back(&dl);
            while (!empty(&dl) && t[i] < t[back(&dl)])
                pop_back(&dl);
        }
        push_back(&du, i);
        push_back(&dl, i);
        if (i == 2 * r + 1 + front(&du))
            pop_front(&du);
        else if (i == 2 * r + 1 + front(&dl))
            pop_front(&dl);
    }
    for (int i = len; i < len + r + 1; i++)
    {
        u[i - r - 1] = t[front(&du)];
        l[i - r - 1] = t[front(&dl)];
        if (i - front(&du) >= 2 * r + 1)
            pop_front(&du);
        if (i - front(&dl) >= 2 * r + 1)
            pop_front(&dl);
    }
    destroy(&du);
    destroy(&dl);
}

double dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
}
double dist(double x1, double y1, double x2, double y2)
{
    return (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
}
double dist(double x, double y)
{
    return (x - y) * (x - y);
}

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
double lb_keogh_cumulative(double* t, double* uo, double* lo, int len, double best_so_far = INF)
{
    double lb = 0;
    double x, d;

    for (int i = 0; i < len && lb < best_so_far; i++)
    {
        x = t[i];
        d = 0;
        if (x > uo[i])
            d = dist(x, uo[i]);
        else if (x < lo[i])
            d = dist(x, lo[i]);
        lb += d;
    }
    return lb;
}


/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double* X_t, double* Y_t, double* X_q, double* Y_q,  int m, int r, double bsf = INF)
{
    // cost[k] - cummulative cost
    double *cost;
    double *cost_prev;
    double *cost_tmp;
    int i,j,k;
    double x,y,z,min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost[k]=INF;

    cost_prev = (double*)malloc(sizeof(double)*(2*r+1));
    for(k=0; k<2*r+1; k++)    cost_prev[k]=INF;

    for (i=0; i<m; i++)
    {
        k = max(0,r-i);
        min_cost = INF;

        for(j=max(0,i-r); j<=min(m-1,i+r); j++, k++)
        {
            /// Initialize all row and column
            if ((i==0)&&(j==0))
            {
                cost[k]=dist(X_t[0],Y_t[0], X_q[0], Y_q[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j-1<0)||(k-1<0))     y = INF;
            else                      y = cost[k-1];
            if ((i-1<0)||(k+1>2*r))   x = INF;
            else                      x = cost_prev[k+1];
            if ((i-1<0)||(j-1<0))     z = INF;
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min( min( x, y) , z) + dist(X_t[i], Y_t[i], X_q[j], Y_q[j]); //last loop: end at k = r
            
            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {
                min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (min_cost >= bsf)
        {
            free(cost);
            free(cost_prev);
            return min_cost;
        }

        /// Move current array to previous array. ( and  cost_prev point to cost, cost point to cost_prev)
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--; // k = (r+1) - 1 = r 

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}

double dtw(double* A, double* B, int m, int r, double bsf = INF)
{

    double* cost;
    double* cost_prev;
    double* cost_tmp;
    int i, j, k;
    double x, y, z, min_cost;

    /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array of size O(r).
    cost = (double*)malloc(sizeof(double) * (2 * r + 1));
    for (k = 0; k < 2 * r + 1; k++)    cost[k] = INF;

    cost_prev = (double*)malloc(sizeof(double) * (2 * r + 1));
    for (k = 0; k < 2 * r + 1; k++)    cost_prev[k] = INF;

    for (i = 0; i < m; i++)
    {
        k = max(0, r - i);
        min_cost = INF;

        for (j = max(0, i - r); j <= min(m - 1, i + r); j++, k++)
        {
            /// Initialize all row and column
            if ((i == 0) && (j == 0))
            {
                cost[k] = dist(A[0], B[0]);
                min_cost = cost[k];
                continue;
            }

            if ((j - 1 < 0) || (k - 1 < 0))     y = INF;
            else                      y = cost[k - 1];
            if ((i - 1 < 0) || (k + 1 > 2 * r))   x = INF;
            else                      x = cost_prev[k + 1];
            if ((i - 1 < 0) || (j - 1 < 0))     z = INF;
            else                      z = cost_prev[k];

            /// Classic DTW calculation
            cost[k] = min(min(x, y), z) + dist(A[i], B[j]);

            /// Find minimum cost in row for early abandoning (possibly to use column instead of row).
            if (cost[k] < min_cost)
            {
                min_cost = cost[k];
            }
        }

        /// We can abandon early if the current cummulative distace with lower bound together are larger than bsf
        if (min_cost >= bsf)
        {
            free(cost);
            free(cost_prev);
            return min_cost;
        }

        /// Move current array to previous array.
        cost_tmp = cost;
        cost = cost_prev;
        cost_prev = cost_tmp;
    }
    k--;

    /// the DTW distance is in the last cell in the matrix of size O(m^2) or at the middle of our array.
    double final_dtw = cost_prev[k];
    free(cost);
    free(cost_prev);
    return final_dtw;
}

void AddData(IndexDTW** indexdtw, int index, double dist, double dist_1d, bool replace) {
    IndexDTW* p, * new_indexdtw;
    new_indexdtw = (IndexDTW*)malloc(sizeof(IndexDTW));

    p = *indexdtw;
    if (dist > p->dist) {
        new_indexdtw->index = index;
        new_indexdtw->dist = dist;
        new_indexdtw->dist_1d = dist_1d;

        new_indexdtw->next = *indexdtw;
        *indexdtw = new_indexdtw;
        return;
    }

    while (p->next && (dist <= p->next->dist)) {
        p = p->next;
    }

    new_indexdtw->index = index;
    new_indexdtw->dist = dist;
    new_indexdtw->dist_1d = dist_1d;

    new_indexdtw->next = p->next;
    p->next = new_indexdtw;

    if (replace) {
        p = *indexdtw; 
        *indexdtw = p->next;
        free(p);
    }
}

/// Print function for debugging
void printArray(double *x, int len)
{   for(int i=0; i<len; i++)
        printf(" %6.2lf",x[i]);
    printf("\n");
}

/// If expected error happens, teminated the program.
void error(int id)
{
    if(id==1)
        printf("ERROR : Memory can't be allocated!!!\n\n");
    else if ( id == 2 )
        printf("ERROR : File not Found!!!\n\n");
    else if ( id == 3 )
        printf("ERROR : Can't create Output File!!!\n\n");
    else if ( id == 4 )
    {
        printf("ERROR : Invalid Number of Arguments!!!\n");
        printf("Command Usage:  UCR_DTW.exe  file-path  data-file  query-file   m   R       topK  std (tag)\n\n");
        printf("For example  :  UCR_DTW.exe  D:/        data.txt   query.txt   128  0.05    20    T/F (tag)\n");
    }
    exit(1);
}

void standardize(double *&ary, double mean, double std, int m) {
    for (int i = 0; i < m; i++) ary[i] = (ary[i] - mean) / std;
}

/// Main Function
int main(  int argc , char *argv[] )
{
    FILE *fp;            /// data file pointer
    FILE *qp;            /// query file pointer

    double *X_t, *Y_t, *X_q, *Y_q; 
    double dist = 0, dist_1d = 0; 
    //double bsf = INF; 
    double lb_keogh1, lb_keogh2;

    long long i , j, index;
    int m=-1, r=-1, topK=-1;
    double t1,t2;
    bool isStd;
    Index* Q_tmp;
    
    IndexDTW* topkdtw = (IndexDTW*)malloc(sizeof(IndexDTW));
    
    
    /// If not enough input, display an error.
    if (argc<=4)
        error(4);

    /// read size of the query
    if (argc>4)
        m = atol(argv[4]);

    /// read warping windows
    if (argc>4)
    {   double R = atof(argv[5]);
        if (R<=1)
            r = floor(R*m);
        else
            r = floor(R);
    }
    
    if (argc > 6) 
        topK = atoi(argv[6]);
    
    if (topK <= 0) {
        free(topkdtw);
    }
    
    switch (argv[7][0])
    {
    case 'T': 
    case 't': isStd = true; break;
    case 'F': 
    case 'f': isStd = false; break;
    default:
        cout << "illegal option." << endl;
        error(4);
    }
    
    //string tag(argv[8]);

    string dataFilePath(argv[1]);
    string queryFilePath(argv[1]);
    dataFilePath += string(argv[2]);
    queryFilePath += string(argv[3]);

    fp = fopen(dataFilePath.c_str(),"r");
    if( fp == NULL )
        error(2);

    qp = fopen(queryFilePath.c_str(),"r");
    if( qp == NULL )
        error(2);

    /// malloc everything here
    
    X_q = (double*)malloc(sizeof(double) * m);
    if (X_q == NULL)
        error(1);
    Y_q = (double*)malloc(sizeof(double) * m);
    if (Y_q == NULL)
        error(1);
    X_t = (double*)malloc(sizeof(double) * m);
    if (X_t == NULL)
        error(1);
    Y_t = (double*)malloc(sizeof(double) * m);
    if (Y_t == NULL)
        error(1);

    double *Merged_q, *Merged_t;

    Merged_q = (double*)malloc(sizeof(double) * m * 2);
    if (Merged_q == NULL)
        error(1);
    Merged_t = (double*)malloc(sizeof(double) * m * 2);
    if (Merged_t == NULL)
        error(1);
    
    Q_tmp = (Index*)malloc(sizeof(Index) * N_MAX_ROWS);
    if (Q_tmp == NULL)
        error(1);

    // fscanf -> sscanf 
    /// Read query file
    double x, y;
    char *stringLine, *loc;
    char row[4*5000*30]; //一行所包含字符最多 个
    int offset;
    double xmean = 0, xstd = 0, ymean = 0, ystd = 0;
    
    fscanf(qp, "%s\n", &row);
    stringLine = (char*)malloc(sizeof(char) * strlen(row));
    strncpy(stringLine, row, strlen(row));
    // printf("%s", row);
    // printf("%s", stringLine);

    i = 0;
    loc = stringLine; //save original pointer
    while (sscanf(stringLine, "(%lf,%lf)%n", &x, &y, &offset)) {
        if ((offset < strlen(stringLine) - 1) && (i <= m - 1)) stringLine = stringLine + offset;
        else break;

        X_q[i] = x;
        Y_q[i] = y;
        
        if(isStd){
            xmean += x;
            ymean += y;
            xstd += x * x;
            ystd += y * y;
        }
        
        i++;
    }
    
    if (isStd) {
        xmean /= i;
        ymean /= i;
        xstd /= i;
        ystd /= i;
        xstd = sqrt(xstd - xmean * xmean);
        ystd = sqrt(ystd - ymean * ymean);
        standardize(X_q, xmean, xstd, i);
        standardize(Y_q, ymean, ystd, i);
    }
    //free(loc);
    
    if (argc > 4)
    {
        double R = atof(argv[4]);
        if (R <= 1)
            r = floor(R * i);
        else
            r = floor(R);
    }
    ///拼接 X_q 与 Y_q 
    for (j = 0; j < i; j++) {
        Merged_q[j] = X_q[j]; 
        Merged_q[j + i] = Y_q[j];
    }

    ///计算拼接序列的
    ////上界和下界
    double* lb_q, * ub_q, * lb_t, * ub_t;

    lb_q = (double*)malloc(sizeof(double) * m * 2);
    if (lb_q == NULL)
        error(1);
    ub_q = (double*)malloc(sizeof(double) * m * 2);
    if (ub_q == NULL)
        error(1);
    lb_t = (double*)malloc(sizeof(double) * m * 2);
    if (lb_t == NULL)
        error(1);
    ub_t = (double*)malloc(sizeof(double) * m * 2);
    if (ub_t == NULL)
        error(1);
    

    lower_upper_lemire(Merged_q, i*2, r, lb_q, ub_q);  


    /// Read data file
    int scanned = 0; 
    int keogh1c = 0, keogh2c = 0, dtwc = 0;
    
    while (fscanf(fp, "%s\n", &row) != EOF) {
        
        stringLine = (char*)malloc(sizeof(char) * strlen(row));
        
        strncpy(stringLine, row, strlen(row));
        i = 0;
        loc = stringLine; 
        
        sscanf(stringLine, "%lld:%n", &index, &offset);
        stringLine = stringLine + offset;
        xmean = 0, xstd = 0, ymean = 0, ystd = 0;
        while (sscanf(stringLine, "(%lf,%lf)%n", &x, &y, &offset)) {
            if ((offset < strlen(stringLine) - 1) && (i<= m-1)) stringLine = stringLine + offset;
            else break;

            X_t[i] = x;
            Y_t[i] = y;

            if (isStd) {
                xmean += x;
                ymean += y;
                xstd += x * x;
                ystd += y * y;
            }

            i++;
        }

        if (isStd) {
            xmean /= i;
            ymean /= i;
            xstd /= i;
            ystd /= i;
            xstd = sqrt(xstd - xmean * xmean);
            ystd = sqrt(ystd - ymean * ymean);
            standardize(X_t, xmean, xstd, i);
            standardize(Y_t, ymean, ystd, i);
        }

        if (topK < 1) {
            Q_tmp[scanned].index = index;
            dist = dtw(X_t, Y_t, X_q, Y_q, i, r); 
            Q_tmp[scanned].value = sqrt(dist);
            
        } else {
            //索引条件
            ///拼接 X_t 与 Y_t 
            
            for (j = 0; j < i; j++) {
                Merged_t[j] = X_t[j];
                Merged_t[j + i] = Y_t[j];
            }
            
            if (scanned < 0) {
                printf("variable <scanned> ERROR.");
                return 1;
            }
            else if (scanned == 0) {
                topkdtw->index = index;
                topkdtw->dist = dtw(X_t, Y_t, X_q, Y_q, i, r);
                topkdtw->dist_1d = dtw(Merged_q, Merged_t, 2*i, r, INF);
                topkdtw->next = NULL;
            }
            else if (scanned < topK) {
                dist = dtw(X_t, Y_t, X_q, Y_q, i, r);
                dist_1d = dtw(Merged_q, Merged_t, 2*i, r, INF); 
                AddData(&topkdtw, index, dist, dist_1d, false);
            } 
            else {
                lb_keogh1 = lb_keogh_cumulative(Merged_t, ub_q, lb_q, i*2, topkdtw->dist_1d);
                if (lb_keogh1 < topkdtw->dist_1d) {
                    ///计算拼接序列的
                    ////上界和下界
                    lower_upper_lemire(Merged_t, i * 2, r, lb_t, ub_t);
                    lb_keogh2 = lb_keogh_cumulative(Merged_q, ub_t, lb_t, i * 2, topkdtw->dist_1d);

                    if (lb_keogh2 < topkdtw->dist_1d) {
                        dist = dtw(X_t, Y_t, X_q, Y_q, i, r, topkdtw->dist); 
                        if (dist < topkdtw->dist)
                        {
                            dist_1d = dtw(Merged_q, Merged_t, 2*i, r, INF); 
                            AddData(&topkdtw, index, dist, dist_1d, true);
                        }
                        else dtwc++;
                    }
                    else keogh2c++;
                }
                else keogh1c++;
            } 

        }

        //free(loc);
        
        scanned++; 
    }
    /*free(lb_q);
    free(lb_t);
    free(ub_q);
    free(ub_t);*/
    
    /*
    string outFilePath("./CdtwOutput/");
    string outFileK(argv[6]);
    string outFileId(argv[3]);
    //string outFileLeng(argv[4]);
    string underline("_");
    int pos = outFileId.rfind('.');

    outFilePath = outFilePath + tag + underline + outFileId.substr(0,pos) + underline + outFileK ;
    if (isStd) {
        string string_std("std");
        outFilePath = outFilePath + underline + string_std + outFileId.substr(pos);
    }
    else {
        outFilePath += outFileId.substr(pos);
    }
    */
    string outFilePath("./CdtwOutput/");
    string outFileK(argv[6]);
    string outFileId(argv[3]);
    string underline("_");
    int pos = outFileId.rfind('.');

    outFilePath = outFilePath + outFileId.substr(0, pos) + underline + outFileK;
    if (isStd) {
        string string_std("std");
        outFilePath = outFilePath + underline + string_std + outFileId.substr(pos);
    }
    else {
        outFilePath += outFileId.substr(pos);
    }
    
    FILE* op = fopen(outFilePath.c_str(), "w");
    if (op == NULL)
        error(2);

    if (topK < 1){
        // Sort Part
        qsort(Q_tmp, scanned, sizeof(Index), comp);
        
        // Save Part
        fprintf(op, "tazid,dist\n");
        for (i = 0; i < scanned; i++) {
            fprintf(op, "%d,%lf\n", Q_tmp[i].index, Q_tmp[i].value);
        }
    }
    else {
        Q_tmp = (Index*)malloc(sizeof(Index) * topK);

        IndexDTW* temp = topkdtw;
        i = 0; 
        while (topkdtw) {
            Q_tmp[i].index = topkdtw->index; 
            Q_tmp[i].value = sqrt(topkdtw->dist);
            topkdtw = topkdtw->next;
            i++; 
        }
        
        fprintf(op, "tazid,dist\n");
        while (i > 0) {
            i--; 
            fprintf(op, "%d,%lf\n", Q_tmp[i].index, Q_tmp[i].value);
        }
        free(temp);
    }
    fclose(op);

    fclose(qp);
    fclose(fp);

    free(Q_tmp);
    
    free(X_q);
    free(Y_q);
    free(X_t);
    free(Y_t);
    free(Merged_q);
    free(Merged_t);

    return 0;
}
