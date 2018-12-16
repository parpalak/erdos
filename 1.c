#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>

#define N 20
#define size 5
#define delta 1e-4
#define step 0.2e-4
#define max_grad 10
#define rand_amount 28.0e-6

double sqr(double x) {
    return x * x;
}

double pair_energy(double dist2) {
    double d4 = sqr(1 - dist2);
    return exp(- d4 * 11) + 0.002/(0.01 + d4) - 100*exp(- dist2 * 60);
}

struct vector
{
     double x;
     double y;
};

struct vector data[N];
struct vector grad[N];

double my_rand() {
    return (double)rand()/(double)(RAND_MAX/size);
}

void init () {
    int seed = time(NULL);
    srand(seed);

    for (int i = N; i-- ;) {
        data[i].x = my_rand();
        data[i].y = my_rand();
    }
}

double distance2(struct vector a, struct vector b) {
    return sqr(a.x - b.x) + sqr(a.y - b.y);
}

double get_energy(struct vector a[]) {
    double w = 0;

    for (int i = N; i--; ) {
        for (int j = i; j--; ) {
            double dist2 = distance2(a[i], a[j]);
            w += pair_energy(dist2);
        }
    }

    return w;
}

double get_part_energy(struct vector a[], int current) {
    double w = 0;

    for (int i = N; i--; ) {
        if (i != current) {
            double dist2 = distance2(a[i], a[current]);
            w += pair_energy(dist2);
        }
    }

    return w;
}

double get_narrow_energy(struct vector a[]) {
    double w = 0;

    for (int i = N; i--; ) {
        for (int j = i; j--; ) {
            double dist2 = distance2(a[i], a[j]);
            w += exp(- sqr(1 - dist2) * 1000);
//            printf("%f w = %f\n", 1 - dist2,  w);
        }
    }

    return w;
}

struct vector cut (struct vector a) {
    double l2 = sqr(a.x) + sqr(a.y);

    if (l2 > sqr(max_grad)) {
        double l = sqrt(l2);
        a.x = a.x * max_grad / l;
        a.y = a.y * max_grad / l;
    }

    return a;
}

struct vector b[N];

void do_step(struct vector a[], struct vector grad[]) {
    double w0 = get_energy(a);
    double w2, w1;

    for (int i = N; i-- ;) {
        b[i].x = a[i].x;
        b[i].y = a[i].y;
    }

    for (int i = N; i-- ;) {
        b[i].x = a[i].x + delta;
        w2 = get_part_energy(b, i);

        b[i].x = a[i].x - delta;
        w1 = get_part_energy(b, i);

        grad[i].x = (w2 - w1) / delta;
        b[i].x = a[i].x;


        b[i].y = a[i].y + delta;
        w2 = get_part_energy(b, i);

        b[i].y = a[i].y - delta;
        w1 = get_part_energy(b, i);

        grad[i].y = (w2 - w1) / delta;
        b[i].y = a[i].y;

        grad[i] = cut(grad[i]);
    }

    for (int i = N; i-- ;) {
        a[i].x += grad[i].x * step + rand_amount * (my_rand() - 0.5*size);
        a[i].y += grad[i].y * step + rand_amount * (my_rand() - 0.5*size);
    }
}

void dump () {
    for (int i = N; i-- ;) {
        printf("%f, %f\n", data[i].x, data[i].y);
    }
}

void dump_points_to_file (char* name) {
    FILE *f;
    f = fopen(name, "w");
    for (int i = N; i-- ;) {
        fprintf(f, "%f\t%f\t%f\t%f\n", data[i].x, data[i].y, 0.1 * grad[i].x, 0.1 * grad[i].y);
    }
    fclose(f);
}

int dump_lines_to_file (char* name, struct vector a[]) {
    FILE *f;
    int n = 0;

    f = fopen(name, "w");
    for (int i = N; i-- ;) {
        for (int j = i; j--; ) {
            double dist2 = distance2(a[i], a[j]);
            if (fabs(1 - dist2) < 1.0e-3) {
                n++;
                fprintf(f, "%f\t%f\t%f\t%f\t%f\n", a[i].x, a[i].y, a[j].x, a[j].y, dist2);
            }
        }
    }
    fclose(f);

    return n;
}

static volatile int keepRunning = 1;

void int_handler(int dummy) {
    keepRunning = 0;
}

int main() {
    signal(SIGINT, int_handler);
    init();

    for (int j = 0; j < 10000000; j++) {
        do_step(data, grad);
        if (!(j % 201)) {
            dump_points_to_file("plane.dat");
            int lines_num = dump_lines_to_file("lines.dat", data);
            printf("%d) w = %f, n = %d\n", j, get_energy(data), lines_num);
            if (!keepRunning) {
                break;
            }
        }
    }

//    dump();
    dump_points_to_file("plane.dat");
}
