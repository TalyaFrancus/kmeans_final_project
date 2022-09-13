#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include "spkmeans.h"

PyObject *c_mat_to_list_py(double **c_mat, int row, int col){ /*convert c array to python list*/
    PyObject *py_list, *py_vector;
    int i, j;

    py_list = PyList_New(row);
    for (i=0; i<row; i++){
        py_vector = PyList_New(col);
        for(j=0; j<col; j++){
            PyList_SET_ITEM(py_vector, j, PyFloat_FromDouble(c_mat[i][j]));
        }
        PyList_SET_ITEM(py_list, i, py_vector);
    }

    return py_list;
}

double ** py_list_to_c_mat(PyObject *py_list, int row, int col){ /*convert python list to c array*/
    double **mat_to_return;
    int i, j;

    mat_to_return = create_double_mat(row, col);

    for (i = 0; i < row; i++) {
        for (j = 0; j < col; j++) {
            mat_to_return[i][j] = PyFloat_AS_DOUBLE(PyList_GetItem(PyList_GetItem(py_list, i), j));
            if (PyErr_Occurred()){
                puts("An Error Has Occurred");
                exit(1);
            }
        }
    }

    return mat_to_return;
}

static PyObject* wam_ddg_lnorm_jacobi (PyObject* self, PyObject* args) {
    char *file_name, *goal;

    if (!PyArg_ParseTuple(args, "ss", &file_name, &goal)){
        return NULL;
    }

    spkmeans(file_name, goal);

    Py_RETURN_NONE;
}

static PyObject* spk_module (PyObject* self, PyObject* args){
    char *file_name;
    int k, rows, cols;
    double **mat_to_return, **input_mat, **lnorm_mat, **jacobi_mat;
    PyObject *py_list_to_return;

    if (!PyArg_ParseTuple(args, "si", &file_name, &k)){
        return NULL;
    }

    cols = num_of_cord(file_name);
    rows = num_of_lines(file_name);

    input_mat = create_double_mat(rows, cols);

    init_the_mat_of_vec(file_name, input_mat[0], rows*cols); /*init the input mat*/

    lnorm_mat = lnorm(input_mat, rows, cols); /*calculate the lnorm of the input*/

    jacobi_mat = jacobi(lnorm_mat, rows); /*calculate the eigenvectors of the lnorm matrix*/

    if (k == 0){
        k = eigengap_heuristic(lnorm_mat, rows);
    }

    if(k>rows){
        printf("%s%c", "Invalid Input!", '\n');
        exit(1);
    }

    mat_to_return = spk(k, rows, lnorm_mat, jacobi_mat);

    py_list_to_return = c_mat_to_list_py(mat_to_return, rows, k);

    free(input_mat[0]);
    free(input_mat);
    free(mat_to_return[0]);
    free(mat_to_return);

    return py_list_to_return;
}

static PyObject* fit (PyObject* self, PyObject* args){
    int num_of_vectors, vector_size, k, max_iter;
    double eps;
    PyObject *list_of_centroids, *list_of_all_vectors, *list_to_return;
    double **arr_of_centroids, **input_mat;

    if (!PyArg_ParseTuple(args, "iiiidOO", &num_of_vectors, &vector_size, &k, &max_iter,
                          &eps, &list_of_centroids, &list_of_all_vectors)){
        return NULL;
    }

    arr_of_centroids = py_list_to_c_mat(list_of_centroids, k, vector_size);

    input_mat = py_list_to_c_mat(list_of_all_vectors, num_of_vectors, vector_size);

    kmeans(k, max_iter, eps, input_mat, arr_of_centroids, vector_size, num_of_vectors);

    list_to_return = c_mat_to_list_py(arr_of_centroids, k, vector_size);

    free(arr_of_centroids[0]);
    free(arr_of_centroids);

    return list_to_return;
}

static PyMethodDef spkmeansMethods[] = {
        {"wam_ddg_lnorm_jacobi",
                (PyCFunction) wam_ddg_lnorm_jacobi,
                     METH_VARARGS,
                        PyDoc_STR("wam_ddg_lnorm_jacobi")},
        {"spk_module",
                (PyCFunction) spk_module,
                     METH_VARARGS,
                        PyDoc_STR("spkmeans method")},
        { "fit",
                (PyCFunction) fit,
                     METH_VARARGS,
                        PyDoc_STR("find K clusters")},

        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeansmodule",
        NULL,
        -1,
        spkmeansMethods
};


PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}