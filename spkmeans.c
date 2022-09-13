#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int num_of_cord(char filename[]);
int num_of_lines (char filename[]);
void init_the_mat_of_vec(char filename[], double *mat, int size_of_mat);
double calculate_the_norm (double *arr1, double *arr2, int len);
double **wam(double **input_matrix, int row, int col);
double **ddg(double **input_matrix, int row, int col);
double sum_of_row(double *row, int len);
double **lnorm(double **input_matrix, int row, int col);
double **jacobi(double **input_matrix, int size);
double **multi_mat(double **left_mat, double **right_mat, int size);
int *pivot(double **mat, int size);
double off(double **mat, int size);
void mult_the_matrix_with_the_rotate(double **mat, int size, double c, double s);
double **create_rotate_matrix(double **mat, int size);
void swap(double* xp, double* yp);
void selection_sort(double arr[], int n);
double **spk(int k, int rows,double **lnorm_mat, double **jacobi_mat);
int eigengap_heuristic(double **input_matrix, int size);
double **renormalize(double** mat, int row, int column);
double unit_length(double *vector, int size);
void spkmeans(char *input_file, char goal[]);
void kmeans(int k, int max_iter, double eps, double **mat_of_all_vec, double **arr_of_centroid, int vector_size, int num_of_line);
void create_mat_of_new_centroid (double **arr_of_sum, int arr_of_counters[], double **mat_of_vec,
                                 double **arr_of_centroid, int num_of_lines, int K, int size_of_vector);
int find_the_closest_centroid(double vector[], double **arr_of_centroid, int K, int size_of_line);
double find_the_distance (double *arr1, double *arr2, int len);
int find_the_biggest_norma(double **arr1, double **arr2, int K, int len_of_vec, double eps);
double **create_double_mat(int rows, int cols);

int main(int argc,char *argv[]) {

    if (argc != 3){
        printf("%s", "Invalid Input!");
        printf("%c", '\n');
        return 1;
    }
    spkmeans(argv[2], argv[1]);
    return 0;
}

void spkmeans(char *input_file, char goal[]){ /* The main function that call to the fourth first goals */
    int i, j, rows, cols;
    double **output_mat, **input_mat;

    rows = num_of_lines(input_file);
    cols = num_of_cord(input_file);

    input_mat = create_double_mat(rows, cols);

    init_the_mat_of_vec(input_file, input_mat[0], rows * cols);

    if (strcmp(goal, "wam") == 0){
        output_mat = wam(input_mat, rows, cols);
    }

    else if (strcmp(goal, "ddg") == 0){
        output_mat = ddg(input_mat, rows, cols);
    }

    else if(strcmp(goal, "lnorm") == 0){
        output_mat = lnorm(input_mat, rows, cols);
    }

    else if(strcmp(goal, "jacobi")==0){
        output_mat = jacobi(input_mat, rows);
        for (i = 0; i < rows; ++i) {
            printf("%.4f", input_mat[i][i]);
            if (i<rows-1){
                printf("%c", ',');
            }
        }
        printf("%c", '\n');
    }
    else{
        printf("%s", "Invalid Input!");
        printf("%c", '\n');
        exit(1);
    }

    for (i = 0; i < rows ; i++) {
        for  (j = 0; j < rows; j++) {
            printf("%.4f", output_mat[i][j]);
            if (j<rows-1){
                printf( "%c", ',');
            }
        }
        printf("%c", '\n');
    }

    free(input_mat[0]);
    free(input_mat);
    free(output_mat[0]);
    free(output_mat);

}

double **spk(int k, int rows,double **lnorm_mat, double **jacobi_mat){ /*The function get the jacobi mat (after the third part of apk algorithm amd compute the mat we return to python*/
    int i, j;
    double *arr_of_eigenvalues;
    double **mat_of_k_biggest_eigenvectors, **renormalizing_mat;

    arr_of_eigenvalues = calloc(2*rows, sizeof(double));

    mat_of_k_biggest_eigenvectors = create_double_mat(rows, k);


    for (i = 0; i < rows; ++i) {
        arr_of_eigenvalues[i] = lnorm_mat[i][i];
    }

    selection_sort(arr_of_eigenvalues, rows); /*sort the eigenvalues and save the fit indexes*/

    for (i = 0; i < rows; ++i) {  /*Init the matrix of the k biggest eigenvectors of lnorm*/
        for (j = 0; j < k; j++) {
            mat_of_k_biggest_eigenvectors[i][j] = jacobi_mat[i][(int)arr_of_eigenvalues[rows+j]];
        }
    }

    renormalizing_mat = renormalize(mat_of_k_biggest_eigenvectors, rows, k);

    free(lnorm_mat[0]);
    free(lnorm_mat);
    free(arr_of_eigenvalues);
    free(jacobi_mat[0]);
    free(jacobi_mat);
    free(mat_of_k_biggest_eigenvectors[0]);
    free(mat_of_k_biggest_eigenvectors);

    return renormalizing_mat;
}

double **wam(double **input_matrix, int row, int col){  /*Calculate the Weighted Adjacency Matrix and return it in a new matrix*/
    int i, j;
    double **wam_mat;

    wam_mat = create_double_mat(row, row);

    for (i = 0; i<row; i++){  /*Calculate the return matrix*/
        for(j = 0; j<row; j++){
            if(i==j){
                wam_mat[i][j] = 0;
            }
            else{wam_mat[i][j] = calculate_the_norm(input_matrix[i], input_matrix[j], col);}
        }
    }
    return wam_mat;
}

double **ddg(double **input_matrix, int row, int col){
    int i;
    double **ddg_mat;
    double **wam_matrix = wam(input_matrix, row, col); /*Calculate the weighted adjacency matrix*/

    ddg_mat = create_double_mat(row, row);

    for (i=0; i<row; i++){  /*Calculate the diagonal degree matrix*/
        ddg_mat[i][i] = sum_of_row(wam_matrix[i], row);
    }

    free(wam_matrix[0]);
    free(wam_matrix);

    return ddg_mat;

}

double **lnorm(double **input_matrix, int row, int col){ /* res = I - ddg_mat^-0.5*wam_mat*ddg_mat^-0.5 */
    int i, j;
    double **wam_mat = wam(input_matrix, row, col);
    double **ddg_mat = ddg(input_matrix, row, col);

    for (i=0; i<row; i++){
        ddg_mat[i][i] = pow(ddg_mat[i][i], -0.5);
    }

    for (i=0; i<row; i++){
        for (j = 0; j < row; j++) {
            if (i==j){
                wam_mat[i][j] = 1 - (wam_mat[i][j] * ddg_mat[i][i] * ddg_mat[j][j]);
            }
            else{
                wam_mat[i][j] = -wam_mat[i][j] * ddg_mat[i][i] * ddg_mat[j][j];
            }
        }
    }

    free(ddg_mat[0]);
    free(ddg_mat);

    return wam_mat;
}

double **jacobi(double **input_matrix, int size){
    double epsilon = 0.00001;
    double **p, **new_p;
    double *id;
    double prev_off = off(input_matrix, size);
    double current_off;
    int counter = 1;
    int i;

    if(prev_off==0){
        id = calloc(size*size, sizeof(double));
        p = calloc(size, sizeof(double *));
        if (p== NULL||id == NULL){
            printf("%s", "An Error Has Occurred");
            printf("%c", '\n');
            exit(1);
        }
        for (i=0; i<size; i++){
            p[i] = id + i*size;
        }
        for (i = 0; i < size; i++) {
            p[i][i] = 1;
        }
        return p;
    }

    p = create_rotate_matrix(input_matrix, size);  /*create the first rotate matrix*/
    current_off = off(input_matrix, size);

    while (counter < 100 && prev_off-current_off>epsilon){
        prev_off = current_off;
        new_p = create_rotate_matrix(input_matrix, size); /*calculate the new rotate matrix*/
        p = multi_mat(p, new_p, size); /*create the new p matrix and free the memory of the prev p and new_p*/
        current_off = off(input_matrix, size);
        counter++;
        }


    return p;  /*the input matrix is diagonal matrix*/
}

int eigengap_heuristic(double **matrix, int size){ /*choose k in case of k=0*/
    int k = 0;
    int i;
    double diff = 0;
    double *arr;

    arr = calloc(size*2, sizeof(double));

    if (arr == NULL){
        printf("%s", "An Error Has Occurred");
        printf("%c", '\n');
        exit(1);
    }

    for (i = 0; i < size; ++i) {
        arr[i] = matrix[i][i];
    }
    selection_sort(arr, size); /*sort the eigenvalues*/
    for (i = 0; i < size/2; ++i) {
        if (fabs(arr[i]-arr[i+1])>diff){
            k = i+1;
            diff = fabs(arr[i]-arr[i+1]);
        }
    }

    free(arr);

    return k;
}

void swap(double* xp, double* yp){
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void selection_sort(double arr[], int n){ /*sort the given array and save the original indexes before the sorting*/
    int i, j, max_idx;
    double *indexes_arr;

    indexes_arr = calloc(n, sizeof(double));

    if (indexes_arr == NULL){
        printf("%s", "An Error Has Occurred");
        printf("%c", '\n');
        exit(1);
    }

    for (i = 0; i < n; ++i) {
        indexes_arr[i] = i;
    }

    for (i = 0; i < n - 1; i++) {
        max_idx = i;
        for (j = i + 1; j < n; j++)
            if (arr[j] > arr[max_idx])
                max_idx = j;

        swap(&arr[max_idx], &arr[i]);
        swap(&indexes_arr[max_idx], &indexes_arr[i]);
    }

    for (i = n; i < 2*n; i++) {
        arr[i] = indexes_arr[i-n];
    }
    free(indexes_arr);
}


int num_of_cord(char filename[]){
    FILE *file;
    int num_of_cord = 1;
    int c;
    file = fopen(filename, "r");
    if (file == NULL){
        printf("%s", "Invalid Input!");
        printf("%c", '\n');
        exit(1);
    }

    while (c != '\n'){
        c = fgetc(file);
        if (c == 44){
            num_of_cord += 1;
        }
    }
    fclose(file);
    return num_of_cord;
}

int num_of_lines (char filename[]){
    FILE *file;
    int num_of_lines = 0;
    char c;
    file = fopen(filename, "r");
    if (file == NULL){
        printf("%s", "Invalid Input!");
        printf("%c", '\n');
        exit(1);
    }

    while (c != EOF){
        c = fgetc(file);
        if (c == '\n'){
            num_of_lines += 1;
        }
    }
    fclose(file);
    return num_of_lines;
}

void init_the_mat_of_vec(char filename[], double *mat, int size_of_mat){ /* the function get matrix and file and fill the mat with the number of the file*/
    int i=0;
    FILE *fp;
    double cord;
    int read;
    fp = fopen(filename, "r");

    while (i < size_of_mat) {
        read = fscanf(fp, "%lf", &cord);
        if (read == 0){
            read = fscanf(fp, "%*c");
        }
        else{
            mat[i] = cord;
            i++;
        }
    }
    fclose(fp);
}

double calculate_the_norm (double *arr1, double *arr2, int len) {
    double distance = 0;
    int i;
    for (i = 0; i < len; i++) {
        distance += pow((arr1[i] - arr2[i]), 2);
    }

    return exp(-(pow(distance, 0.5)/2));
}

double sum_of_row(double *row, int len){  /*calculate the sum of the coordinates in the given vector*/
    int i;
    double sum = 0;
    for (i=0; i<len; i++){
       sum += row[i];
    }
    return sum;
}

double **multi_mat(double **left_mat, double **right_mat, int size){ /*The function multiple 2 matrix and return the result in new matrix*/
    int i, j, k;
    double **mat_to_return;

    mat_to_return = create_double_mat(size, size);

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            for (k=0; k<size; k++){
                mat_to_return[i][j] += left_mat[i][k]*right_mat[k][j];
            }
        }
    }
    free(left_mat[0]);
    free(left_mat);
    free(right_mat[0]);
    free(right_mat);
    return mat_to_return;
}

int *pivot(double **mat, int size){ /*find the off-diagonal element with the largest absolute value*/
    int *indexes;
    int i, j;
    indexes = calloc(2, sizeof(int));
    if (indexes == NULL){
        printf("%s", "An Error Has Occurred");
        printf("%c", '\n');
        exit(1);
    }
    indexes[0] = 0;
    indexes[1] = 1;
    for (i = 0; i<size; i++){
        for (j=i; j<size; j++){
            if (i != j){
                if (fabs(mat[i][j])> fabs(mat[indexes[0]][indexes[1]])){
                    indexes[0] = i;
                    indexes[1] = j;
                }
            }
        }
    }
    return indexes;
}

double off(double **mat, int size){
    int i, j;
    double sum = 0;
    for (i = 0; i < size; i++) {
        for (j=0; j<size; j++){
            if (i != j){
                sum += pow(mat[i][j], 2);
            }
        }
    }
    return sum;
}

void mult_the_matrix_with_the_rotate(double **mat, int size, double c, double s){ /* new_mat = rotate * mat * transpose(rotate) */
    int* indexes = pivot(mat, size);
    int i = indexes[0];
    int j = indexes[1];
    int k;
    double aki, akj;
    double aii = mat[i][i];
    double ajj = mat[j][j];
    double aij = mat[i][j];
    for (k=0; k<size; k++){
        aki = mat[k][i];
        akj = mat[k][j];
        if (k == i){
            mat[k][i] = c*c*aii + s*s*ajj - 2*s*c*aij;
        }

        else if (k == j){
            mat[k][j] = s*s*aii + c*c*ajj + 2*s*c*aij;
        }

        else{
            mat[k][i] = c*aki - s*akj;
            mat[i][k] = mat[k][i];
            mat[k][j] = c*akj +s*aki;
            mat[j][k] = mat[k][j];
        }

    }
    mat[i][j] = 0;
    mat[j][i] = 0;

    free(indexes);

}

double **create_rotate_matrix(double **mat, int size){ /*create the rotated matrix of the current pivot and call another function to create A'*/
    int *indexes = pivot(mat, size);
    int i, j, sign;
    double **mat_to_return;
    double c, t, s, theta;

    mat_to_return = create_double_mat(size, size);

    for(i=0; i<size; i++){
        mat_to_return[i][i] = 1;
    }
    i = indexes[0];
    j = indexes[1];
    theta = (mat[j][j] - mat[i][i])/(2*mat[i][j]);
    if (theta == 0){
        sign = 1;
    }
    else{
        sign = (int) (theta/ fabs(theta));
    }
    t = sign/(fabs(theta) + pow(theta*theta+1, 0.5));
    c = 1/pow(t*t+1, 0.5);
    s = t*c;

    mat_to_return[i][i] = c;
    mat_to_return[j][j] = c;
    mat_to_return[i][j] = s;
    mat_to_return[j][i] = -s;

    mult_the_matrix_with_the_rotate(mat, size, c, s);

    free(indexes);

    return mat_to_return;
}

double **renormalize(double** mat, int row, int column){ /*renormalize each row of the given mat and return new mat*/
    double **renormalize_mat;
    int i, j;
    double unit;

    renormalize_mat = create_double_mat(row, column);

    for (i=0; i<row; i++) {
        for (j = 0; j < column; j++) {
            unit = unit_length(mat[i], column);
            if (unit == 0){
                renormalize_mat[i][j] = mat[i][j];
            }
            else{
                renormalize_mat[i][j] = mat[i][j] / (pow(unit, 0.5));
            }
        }
    }

    return renormalize_mat;
}

double unit_length(double *vector, int size){
    int i;
    double num = 0;

    for (i = 0; i < size; i++) {
        num += pow(vector[i], 2);
    }
    return num;

}

/*Here start the part of the project that use the algorithm from HW2*/

void kmeans(int k, int max_iter, double eps, double **mat_of_all_vec, double **arr_of_centroid, int vector_size, int num_of_line){
    int num_of_cord = vector_size;
    int num_of_lines = num_of_line;
    int i, j;
    double **arr_of_sum, **copy_of_centroid;
    int *arr_of_counters;

    if (k<1 || k>num_of_lines || max_iter<0){
        printf("%s", "Invalid Input!");
        printf("%c", '\n');
        exit(1);
    }

    copy_of_centroid = create_double_mat(num_of_cord, k);

    arr_of_sum = create_double_mat(num_of_cord, k);

    arr_of_counters = calloc(k, sizeof(int));
    if (arr_of_counters == NULL){
        printf("An Error Has Occurred");
        printf("%c", '\n');
        exit(1);
    }


    for (i=0; i<k*num_of_cord; i++){
        copy_of_centroid[0][i] = arr_of_centroid[0][i];
    }

    create_mat_of_new_centroid(arr_of_sum, arr_of_counters, mat_of_all_vec,
                               arr_of_centroid,num_of_lines, k, num_of_cord);

    j=1;

    while (find_the_biggest_norma(arr_of_centroid, copy_of_centroid, k, num_of_cord, eps) && j < max_iter){
        for (i=0; i<k*num_of_cord; i++){
            copy_of_centroid[0][i] = arr_of_centroid[0][i];
        }
        create_mat_of_new_centroid(arr_of_sum, arr_of_counters, mat_of_all_vec,
                                   arr_of_centroid, num_of_lines,k, num_of_cord);

        j++;
    }

    free(arr_of_sum[0]);
    free(arr_of_sum);
    free(mat_of_all_vec[0]);
    free(mat_of_all_vec);
    free(arr_of_counters);
    free(copy_of_centroid[0]);
    free(copy_of_centroid);

}

void create_mat_of_new_centroid (double **arr_of_sum, int arr_of_counters[], double **mat_of_vec,
                                 double **arr_of_centroid, int num_of_lines, int K, int size_of_vector) {
    int i, index, j;

    for (i = 0; i < num_of_lines; i++) {
        index = find_the_closest_centroid(mat_of_vec[i] , arr_of_centroid, K, size_of_vector);
        for (j = 0; j < size_of_vector; j++) {
            arr_of_sum[index][j] += mat_of_vec[i][j];
        }
        arr_of_counters[index]++;
    }
    for (i = 0; i < K ; i++) {
        for (j = 0; j < size_of_vector ; j++) {
            arr_of_centroid[i][j] = arr_of_sum[i][j] / arr_of_counters[i];
            arr_of_sum[i][j] = 0;
        }
        arr_of_counters[i] = 0;
    }
}


int find_the_closest_centroid(double vector[], double **arr_of_centroid, int K, int size_of_line){
    int i;
    double distance = 1000000;
    double res;
    int index = 0;
    for (i = 0;  i<K ; i++) {
        res = find_the_distance(vector, arr_of_centroid[i], size_of_line);
        if (res < distance){
            distance = res;
            index = i;
        }
    }
    return index;
}

double find_the_distance (double *arr1, double *arr2, int len) {
    double distance = 0;
    int i;
    for (i = 0; i < len; i++) {
        distance += pow((arr1[i] - arr2[i]), 2);
    }

    return distance;
}


int find_the_biggest_norma(double **arr1, double **arr2, int K, int len_of_vec, double eps){
    int i;
    double res;
    for (i=0; i<K; i++){
        res = pow(find_the_distance(arr1[i], arr2[i], len_of_vec), 0.5);
        if (res >= eps){
            return 1;
        }
    }
    return 0;
}

double **create_double_mat(int rows, int cols){
    double *mat_to_return;
    double **pointers_to_return;
    int i;

    mat_to_return = calloc(rows*cols, sizeof(double)); /*Init the return matrix*/
    pointers_to_return = calloc(rows, sizeof(double *));
    if (mat_to_return == NULL || pointers_to_return == NULL){
        printf("%s", "An Error Has Occurred");
        printf("%c", '\n');
        exit(1);
    }
    for (i=0; i<rows; i++){
        pointers_to_return[i] = mat_to_return + i*cols;
    }

    return pointers_to_return;
}