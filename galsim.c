#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "graphics/graphics.h"  
#include <X11/Xlib.h>           
#include <string.h>
#include <sys/time.h>
#include <unistd.h>

//Reference: Lab5, Lab6, Lab7, Assignment2, chatgpt, deepseek

// Struct of Arrays (SOA)
typedef struct Particles {
    double *x;      
    double *y;      
    double *mass;   
    double *v_x;    
    double *v_y;   
    double *bright; 
} Particles;

//measuring time, learn from Lab5_Task4
//Here we use inline to reduce functions calls,learn from Lab5_Task5
static inline double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec / 1000000;
}

//Read particles from a binary file, learn from Assignment2_part2
int readfile(const char *filename, Particles *particles, const int N) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        fprintf(stderr, "Cannot open input file: %s\n", filename);
        return 1;
    }
    for (int i = 0; i < N; i++) {
        if (fread(&particles->x[i], sizeof(double), 1, file) != 1 ||
            fread(&particles->y[i], sizeof(double), 1, file) != 1 ||
            fread(&particles->mass[i], sizeof(double), 1, file) != 1 ||
            fread(&particles->v_x[i], sizeof(double), 1, file) != 1 ||
            fread(&particles->v_y[i], sizeof(double), 1, file) != 1 ||
            fread(&particles->bright[i], sizeof(double), 1, file) != 1) {
            fprintf(stderr, "Error reading particle data at index %d\n", i);
            fclose(file);
            return 1;
        }
    }
    fclose(file);
    return 0;
}
    


// // accept a_x and a_y as arguments, which avoids repeated memory allocation, improving perfomance, it comes from Deepseek.
void simulationStep(Particles *restrict particles, const int N, const double dt, double *restrict a_x, double *restrict a_y) {
    const double G = 100.0 / N;
    const double eps = 1e-3;

     // Reset accelerations to zero, learn from Lab5_Task4
     memset(a_x, 0, N * sizeof(double));
     memset(a_y, 0, N * sizeof(double));

    // Compute forces and update accelerations
    // we use temporal locality, reuse data, improve the speed
    for (int i = 0; i < N; i++) {
        double tmpax = a_x[i];
        double tmpay = a_y[i];
        double tmpx = particles->x[i];
        double tmpy = particles->y[i];
        double tempmass=particles->mass[i];
        for (int j = i + 1; j < N; j++) { 
            double r_x = tmpx - particles->x[j];
            double r_y = tmpy - particles->y[j];
            double r   = sqrt(r_x * r_x + r_y * r_y);
            double r_e = (r + eps) * (r + eps) * (r + eps);
            double temp = -G  / r_e;
            tmpax += temp *particles->mass[j]*r_x;
            tmpay += temp *particles->mass[j]*r_y;
            a_x[j] -= temp*tempmass *r_x;
            a_y[j] -= temp*tempmass *r_y;
        }
        a_x[i] = tmpax;
        a_y[i] = tmpay;
        
    }

    // Update velocities and positions
    for (int i = 0; i < N; i++) {
        particles->v_x[i] += dt * a_x[i];
        particles->v_y[i] += dt * a_y[i];
        particles->x[i]   += dt * particles->v_x[i];
        particles->y[i]   += dt * particles->v_y[i];
    

      
    }

}

//Write particles to a binary file, which we use chatgpt to modify our code
int binary_particle(const char *filename, const Particles *particles, const int N) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        fprintf(stderr, "Cannot open output file: %s\n", filename);
        return 1;
    }
    for (int i = 0; i < N; i++) {
        if (fwrite(&particles->x[i], sizeof(double), 1, file) != 1 ||
            fwrite(&particles->y[i], sizeof(double), 1, file) != 1 ||
            fwrite(&particles->mass[i], sizeof(double), 1, file) != 1 ||
            fwrite(&particles->v_x[i], sizeof(double), 1, file) != 1 ||
            fwrite(&particles->v_y[i], sizeof(double), 1, file) != 1 ||
            fwrite(&particles->bright[i], sizeof(double), 1, file) != 1) {
            fprintf(stderr, "Error writing particle data at index %d\n", i);
            fclose(file);
            return 1;
        }
    }
    fclose(file);
    return 0;
}

int main(int argc, char *argv[]) {
    double startTime1 = get_wall_seconds();
    if (argc != 6) {
        printf("Usage: galsim N filename nsteps delta_t graphics\n");
        return 1;
    }

    const int N = atoi(argv[1]);
    const char *filename = argv[2];
    const int nsteps = atoi(argv[3]);
    const double dt = atof(argv[4]);
    const int graphicsEnabled = atoi(argv[5]);
    int step = 0;
    
    // Allocate memory for the Particles struct (SOA)
    Particles particles;
    particles.x = (double*) malloc(N * sizeof(double));
    particles.y = (double*) malloc(N * sizeof(double));
    particles.v_x = (double*) malloc(N * sizeof(double));
    particles.v_y = (double*) malloc(N * sizeof(double));
    particles.mass = (double*) malloc(N * sizeof(double));
    particles.bright = (double*) malloc(N * sizeof(double));

    if (particles.x == NULL || particles.y == NULL || particles.v_x == NULL || 
        particles.v_y == NULL || particles.mass == NULL || particles.bright == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    readfile(filename, &particles, N);

    

    // Set window dimensions for graphics
    const int windowWidth = 800;
    const int windowHeight = 600;
    // Initialize graphics if enabled
    if (graphicsEnabled == 1) {
        InitializeGraphics(argv[0], windowWidth, windowHeight); 
    }

    //using buffer for stroring a_x and a_y contiguously in memory, learn from Lab6_Task1
    double *buffer = (double *)malloc(2 * N * sizeof(double));
    if (buffer == NULL) {
        fprintf(stderr, "Memory allocation failed for buffer\n");
        return 1;
    }
    double *restrict a_x = buffer;
    double *restrict a_y = buffer + N;

    


    // Main simulation loop
    while (step < nsteps) {
        simulationStep(&particles, N, dt, a_x, a_y);
        step++;

        if (graphicsEnabled == 1) {
            ClearScreen();
            float min_bright = 0.0f;
            float max_bright = 1.0f;

            // when I use the input data's brightness, which is out of range, so we 
            // need to normalize brightness to [0, 1]
            for (int i = 0; i < N; i++) {
                if (particles.bright[i] < min_bright) min_bright = particles.bright[i];
                if (particles.bright[i] > max_bright) max_bright = particles.bright[i];
            }
            for (int i = 0; i < N; i++) {
                particles.bright[i] = (particles.bright[i] - min_bright) / (max_bright - min_bright);
            }

            // Draw each particle as a circle.
            // Map particle coordinates (assumed in [0,1]) to screen coordinates, we learn it from chatgpt
            for (int i = 0; i < N; i++) {
                float screenX = (float)(particles.x[i] * windowWidth);
                float screenY = (float)(particles.y[i] * windowHeight);
                //making the radius proportional to the particle's mass,  which comes from deepseek
                float radius = 5.0f * (float)particles.mass[i];
                //use particle's brightness field for rendering, which comes from deepseek
                float color = (float)particles.bright[i];
                DrawCircle(screenX, screenY, (float)windowWidth, (float)windowHeight, radius, color);
            }
            Refresh();
            //control the simulation speed, add a small delay, which comes from chatgpt
            usleep(20000);
            if (CheckForQuit()) break;
        }
    }

    // Clean up
    if (graphicsEnabled == 1) {
        FlushDisplay();
        CloseDisplay();
    } else {
        printf("Simulation complete.\n");
    }

    // Write results to file
    binary_particle("result.gal", &particles, N);

    // Free memory
    free(particles.x);
    free(particles.y);
    free(particles.v_x);
    free(particles.v_y);
    free(particles.mass);
    free(particles.bright);
    free(buffer);

    // Print total time taken
    double totalTimeTaken = get_wall_seconds() - startTime1;
    printf("totalTimeTaken = %f\n", totalTimeTaken);
    return 0;
}