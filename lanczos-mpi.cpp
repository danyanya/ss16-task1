
/*
 * Superprac 1. Algowiki Lanczos.
 * Code rewrited by Sliusar Daniil and Grigoryev Mikhail.
 * First realisation was presented by MATVEC
 * https://people.sc.fsu.edu/~jburkardt/cpp_src/mpi/matvec_mpi.cpp
 *
 */

#include <mpi.h>
#include <string>
#include <cmath>
#include <ctime>
#include <iostream>

#define elem_t int
#define MPI_ELEM MPI_INTEGER

using namespace std;

const double pi = 3.141592653589793;

/* For generating values was found linear method
 * https://en.wikipedia.org/wiki/Linear_congruential_generator
 */
const long long int a = 84589;
const long long int b = 45989;

unsigned long random(int row, int col, int seed)
{
    return a * (row ^ seed) + a * b * (row ^ col) + seed;
}

elem_t *new_vec(int rows)
{
    return  new elem_t[rows]();
}
void del_vec(elem_t *vec)
{
    delete[] vec;
}
void gen_matrix(elem_t **matrix, int rows, int cols, int seed,
                int shift_x = 0, int shift_y = 0)
{
    for (auto i = 0; i < rows; ++i)
        for (auto j = 0; j < cols; ++j)
            matrix[i][j] = random(i + shift_x, j + shift_y, seed);
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
        in_vec[i] /= divider < 1e-4 ? 1 : divider;

    return in_vec;
}

int main (int argc, char **argv)
{
    int ierr, proc_id, cpu_num, mpi_tag;
    elem_t *z_part, *q, *q_prev, *alpha, *beta, *b_init, *result_buffer;
    elem_t beta_prev, alpha_cur, alpha_part_sum, beta_part_sum;

    int shift_x, i, j, index, slaves_num, rows_worker, rows_master;

    const int iter_count = 25;
    const int one_step = 10000;
    const int start_index = 50000;
    const int stop_index = 400000;

    MPI_Status status;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &cpu_num);

    slaves_num = cpu_num - 1;

    if (proc_id == 0) {
        cout << "Lanczos algorithm start to work.." << endl;
    }
    for (index = start_index; index <= stop_index; index += one_step) {

        rows_worker = index / slaves_num;
        rows_master = index % slaves_num;

        elem_t **a_rows;
        q_prev = new_vec(index);

        int seed = 0;
        if (proc_id == 0) seed = time(0);

        ierr = MPI_Bcast(&seed, 1, MPI_ELEM, 0, MPI_COMM_WORLD);

        /* Initialisation of new matrix and new vector for produce test */
        if (proc_id == 0) {
            alpha = new_vec(index);
            beta = new_vec(index);
            b_init = new_vec(index);
            a_rows = new elem_t*[rows_master];
            for (auto i = 0; i < rows_master; ++i)  a_rows[i] = new_vec(index);
            z_part = new_vec(rows_master + 1);
            result_buffer = new_vec(rows_worker + 1);

            shift_x = 0;
            gen_matrix(a_rows, rows_master, index, seed);
            gen_vector(b_init, index);
            q = normalize_vector(b_init, index);

        } else {
            q = new_vec(index);
            a_rows = new elem_t*[rows_worker];
            for (auto i = 0; i < rows_worker; ++i) a_rows[i] = new_vec(index);
            z_part = new_vec(rows_worker);
            shift_x = rows_master + (proc_id - 1) * rows_worker;
            gen_matrix(a_rows, rows_worker, index, seed, shift_x);
        }

        /* Produce iter_count number of iterations */
        double start_time = MPI_Wtime();
        for (int iter = 0; iter < iter_count; ++iter) {
            ierr = MPI_Bcast(q, index, MPI_ELEM, 0, MPI_COMM_WORLD);

            /* Compute part summ for vector alpha */
            if (proc_id == 0) {
                for (i = 0; i < rows_master; ++i) {
                    z_part[i] = 0;
                    for (j = 0; j < index; ++j) {
                        z_part[i] += a_rows[i][j] * q[j];
                    }
                }
                alpha_part_sum = 0;
                for (i = 0; i < rows_master; ++i)
                    alpha_part_sum += z_part[i] * q[shift_x + i];
            } else {
                for (i = 0; i < rows_worker; ++i) {
                    z_part[i] = 0;
                    for (j = 0; j < index; ++j) {
                        z_part[i] += a_rows[i][j] * q[j];
                    }
                }
                alpha_part_sum = 0;
                for (i = 0; i < rows_worker; ++i) {
                    alpha_part_sum += z_part[i] * q[shift_x + i];
                }

                mpi_tag = proc_id;
                ierr = MPI_Send(&alpha_part_sum, 1, MPI_ELEM, 0,
                                mpi_tag, MPI_COMM_WORLD);
            }

            /* Send alpha vector to slaves */
            if (proc_id == 0) {
                alpha_cur = alpha_part_sum;
                for (i = 0; i < slaves_num; ++i) {
                    ierr = MPI_Recv(&alpha_part_sum, 1, MPI_ELEM,
                        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    alpha_cur += alpha_part_sum;
                }
                alpha[iter] = alpha_cur;
            }

            ierr = MPI_Bcast(&alpha_cur, 1, MPI_ELEM, 0, MPI_COMM_WORLD);

            /* Compute beta vector part summ */
            if (proc_id == 0) {
                for (i = 0; i < rows_master; ++i) {
                    z_part[i] -= alpha_cur * q[shift_x + i] +
                            beta_prev * q_prev[shift_x + i];
                }
                beta_part_sum = 0;
                for (i = 0; i < rows_master; ++i) {
                    beta_part_sum += z_part[i] * z_part[i];
                }

            } else {
                for (i = 0; i < rows_worker; ++i) {
                    z_part[i] -= alpha_cur * q[shift_x + i] +
                            beta_prev * q_prev[shift_x + i];
                }
                beta_part_sum = 0;
                for (i = 0; i < rows_worker; ++i) {
                    beta_part_sum += z_part[i] * z_part[i];
                }

                mpi_tag = proc_id;
                ierr = MPI_Send(&beta_part_sum, 1, MPI_ELEM, 0, mpi_tag,
                                MPI_COMM_WORLD);
            }

            /* Send beta vector to slaves*/
            if (proc_id == 0) {
                beta_prev = beta_part_sum;
                for (i = 0; i < slaves_num; ++i) {
                    ierr = MPI_Recv(&beta_part_sum, 1, MPI_ELEM,
                        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    beta_prev += beta_part_sum;
                }
                beta_prev = sqrt(beta_prev);
                beta[iter] = beta_prev;
            }

            ierr = MPI_Bcast(&beta_prev, 1, MPI_ELEM, 0, MPI_COMM_WORLD);

            if (beta_prev < 1e-4) break;

            /* Normalize vector z and saving next part result */
            if (proc_id == 0) {
                for (i = 0; i < rows_master; ++i) z_part[i] /= beta_prev;
                std::swap(q_prev, q);
            } else {
                for (i = 0; i < rows_worker; ++i) z_part[i] /= beta_prev;
                mpi_tag = shift_x;
                ierr = MPI_Send(z_part, rows_worker, MPI_ELEM,
                                0, mpi_tag, MPI_COMM_WORLD);
                std::swap(q_prev, q);
            }

            /* Save what slaves computed */
            if (proc_id == 0) {
                for (i = 0; i < slaves_num; ++i) {
                    ierr = MPI_Recv(result_buffer, rows_worker, MPI_ELEM,
                        MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    mpi_tag = status.MPI_TAG;
                    memcpy(q + mpi_tag, result_buffer, rows_worker * sizeof *q);
                }
            }
        }

        /* Write answer to stdout and free all memory */
        double uptime = MPI_Wtime() - start_time;
        if (proc_id == 0) {
            cout << index << " " << cpu_num << " " << uptime << endl;

            del_vec(alpha);
            del_vec(beta);
            del_vec(result_buffer);
            for (auto i = 0; i < rows_master; ++i) del_vec(a_rows[i]);
        } else {
            for (auto i = 0; i < rows_worker; ++i) del_vec(a_rows[i]);
        }

        del_vec((elem_t*)a_rows);
        del_vec(q);
        del_vec(q_prev);
        del_vec(z_part);
    }

    ierr = MPI_Finalize();

    return 0;
}
