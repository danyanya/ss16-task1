
/*
 * Superprac 1. Algowiki Lanczos.
 * Code rewrited by Sliusar Daniil and Grigoriev Mihail.
 *
 */

#include <mpi.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>

#define elem_t int
#define MPI_ELEM MPI_INTEGER

using namespace std;

const double pi = 3.141592653589793;
static unsigned long x = 123456789, y = 362436069, z = 521288629;
unsigned long xorshf96(void) {
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
}

elem_t *new_vec(int rows)
{
    return  new elem_t[rows]();
}
void del_vec(elem_t *vec)
{
    delete[] vec;
}
void gen_matrix(elem_t **matrix, int rows, int cols)
{
    for (auto i1 = 0; i1 < rows; ++i1)
        for (auto j1 = 0; j1 < cols; ++j1)
            matrix[i1][j1] = xorshf96();
}
void gen_vector(elem_t *vec, int cols)
{
    for (auto i = 0; i < cols; ++i)
        vec[i] = sqrt(2.0 / cols + 1) * sin(pi * i / cols + 1);
}
elem_t * normalize_vector(elem_t *in_vec, const int size)
{
    elem_t divider;
    int i = 0;
    for (i = 0; i < size; ++i)
        divider += (elem_t) (in_vec[i] * in_vec[i]);

    divider = sqrt(divider);

    for (i = 0; i < size; ++i)
        in_vec[i] /= divider < 1e-16 ? 1 : divider;

    return in_vec;
}

/*
void count_part_sum(elem_t **matrix, int rows, int cols, )
{
    int i = 0;
    for (i = 0; i < rows_per_master; ++i) {
        z_part[i] = 0;
        for (j = 0; j < index; ++j) {
            z_part[i] += a_rows[i][j] * q[j];
        }
    }
    alpha_part_sum = 0.0;
    for (i = 0; i < rows_per_master; ++i) {
        alpha_part_sum += z_part[i] * q[index_start + i];
    }
}
*/

int main (int argc, char **argv)
{
    int ierr, proc_id, proc_num, mpi_tag;
    MPI_Status status;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &proc_num);

    if (proc_id == 0) {
        cout << "MPI inited.\n The number of processes is " << proc_num << endl;
    }

    elem_t *z_part;
    elem_t *q, *q_prev;
    elem_t *alpha, *beta;
    elem_t *b_init;
    elem_t *buf;

    elem_t beta_prev = 0.0;
    elem_t alpha_cur;
    elem_t alpha_part_sum, beta_part_sum;

    int index_start;
    int i, j;
    int index;
    int workers_num = proc_num - 1;
    int rows_per_worker, rows_per_master;

    const int iter_count = 25;
    const int one_step = 10000;
    const int start_index = 50000;
    const int stop_index = 400000;

    if (proc_id == 0) {
    //        printf("Parameters: shift %d; limit %d/%d; k %d\n", one_step, stop_index, max_lim, iter_count);

    }

    for (index = start_index; index <= stop_index; index += one_step) {

        rows_per_worker = index / workers_num;
        rows_per_master = index % workers_num;

        elem_t **a_rows;
        q_prev = new_vec(index);

        int seed = 0;
        if (proc_id == 0) seed = xorshf96();

        ierr = MPI_Bcast(&seed, 1, MPI_ELEM, 0, MPI_COMM_WORLD);

        /* Initialisation of new metrix and new vector for produce test */
        if (proc_id == 0) {
            alpha = new_vec(index);
            beta = new_vec(index);
            b_init = new_vec(index);
            a_rows = new elem_t*[rows_per_master];
            for (auto i = 0; i < rows_per_master; ++i)
                a_rows[i] = new_vec(index);
            z_part = new_vec(rows_per_master + 1);
            buf = new_vec(rows_per_worker + 1);

            index_start = 0;
            gen_matrix(a_rows, rows_per_master, index);
            gen_vector(b_init, index);
            q = normalize_vector(b_init, index);

        } else {
            q = new_vec(index);
            a_rows = new elem_t*[rows_per_worker];
            for (auto i = 0; i < rows_per_worker; ++i)
                a_rows[i] = new_vec(index);
            z_part = new_vec(rows_per_worker);
            index_start = rows_per_master + (proc_id - 1) * rows_per_worker;
            gen_matrix(a_rows, rows_per_worker, index);
        }

        /* Produce iter_count number of iterations */
        double start_time = MPI_Wtime();
        for (int iter = 0; iter < iter_count; ++iter) {
            ierr = MPI_Bcast(q, index, MPI_ELEM, 0, MPI_COMM_WORLD);
            if (proc_id == 0) {
                for (i = 0; i < rows_per_master; ++i) {
                    z_part[i] = 0;
                    for (j = 0; j < index; ++j) {
                        z_part[i] += a_rows[i][j] * q[j];
                    }
                }
                alpha_part_sum = 0.0;
                for (i = 0; i < rows_per_master; ++i) {
                    alpha_part_sum += z_part[i] * q[index_start + i];
                }
            } else {
                for (i = 0; i < rows_per_worker; ++i) {
                    z_part[i] = 0;
                    for (j = 0; j < index; ++j) {
                        z_part[i] += a_rows[i][j] * q[j];
                    }
                }
                alpha_part_sum = 0.0;
                for (i = 0; i < rows_per_worker; ++i) {
                    alpha_part_sum += z_part[i] * q[index_start + i];
                }

                mpi_tag = proc_id;
                ierr = MPI_Send(&alpha_part_sum, 1, MPI_ELEM, 0,
                                mpi_tag, MPI_COMM_WORLD);
            }

            if (proc_id == 0) {
                alpha_cur = alpha_part_sum;
                for (i = 0; i < workers_num; ++i) {
                    ierr = MPI_Recv(&alpha_part_sum, 1, MPI_ELEM,
                        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    alpha_cur += alpha_part_sum;
                }
                alpha[iter] = alpha_cur;
            }

            ierr = MPI_Bcast(&alpha_cur, 1, MPI_ELEM, 0, MPI_COMM_WORLD);

            if (proc_id == 0) {
                for (i = 0; i < rows_per_master; ++i) {
                    z_part[i] -= alpha_cur * q[index_start + i] +
                            beta_prev * q_prev[index_start + i];
                }
                beta_part_sum = 0;
                for (i = 0; i < rows_per_master; ++i) {
                    beta_part_sum += z_part[i] * z_part[i];
                }

            } else {
                for (i = 0; i < rows_per_worker; ++i) {
                    z_part[i] -= alpha_cur * q[index_start + i] +
                            beta_prev * q_prev[index_start + i];
                }
                beta_part_sum = 0;
                for (i = 0; i < rows_per_worker; ++i) {
                    beta_part_sum += z_part[i] * z_part[i];
                }

                mpi_tag = proc_id;
                ierr = MPI_Send(&beta_part_sum, 1, MPI_ELEM, 0, mpi_tag,
                                MPI_COMM_WORLD);
            }

            if (proc_id == 0) {
                beta_prev = beta_part_sum;
                for (i = 0; i < workers_num; ++i) {
                    ierr = MPI_Recv(&beta_part_sum, 1, MPI_ELEM,
                        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    beta_prev += beta_part_sum;
                }
                beta_prev = sqrt(beta_prev);
                beta[iter] = beta_prev;
            }

            ierr = MPI_Bcast(&beta_prev, 1, MPI_ELEM, 0, MPI_COMM_WORLD);

            if (beta_prev < 1e-16) break;

            if (proc_id == 0) {
                for (i = 0; i < rows_per_master; ++i) z_part[i] /= beta_prev;
                std::swap(q_prev, q);
            } else {
                for (i = 0; i < rows_per_worker; ++i) z_part[i] /= beta_prev;
                mpi_tag = index_start;
                ierr = MPI_Send(z_part, rows_per_worker, MPI_ELEM,
                                0, mpi_tag, MPI_COMM_WORLD);
                std::swap(q_prev, q);
            }

            if (proc_id == 0) {
                for (i = 0; i < workers_num; ++i) {
                    ierr = MPI_Recv(buf, rows_per_worker, MPI_ELEM,
                        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    mpi_tag = status.MPI_TAG;
                    memcpy(q + mpi_tag, buf, rows_per_worker * sizeof *q);
                }
            }
        }

        double uptime = MPI_Wtime() - start_time;
        if (proc_id == 0) {
            cout << index << " " << proc_num << " " << uptime << endl;

            del_vec(alpha);
            del_vec(beta);
            del_vec(buf);
            for (auto i = 0; i < rows_per_master; ++i) del_vec(a_rows[i]);
        } else {
            for (auto i = 0; i < rows_per_worker; ++i) del_vec(a_rows[i]);
        }

        del_vec((elem_t*)a_rows);
        del_vec(q);
        del_vec(q_prev);
        del_vec(z_part);
    }

    ierr = MPI_Finalize();

    return 0;
}
