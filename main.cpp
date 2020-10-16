// Derivative work of Farouk Ounanes' original work
// found here: https://github.com/ChinksofLight/mandelbrot_cpp?source=post_page-----c7ad6a1bf2d9----------------------
//
// Revised by: Daniel Carns
// This program has been modified to run with OpenMP and
// includes multiple colors instead of only black and red

#include <iostream>
#include <fstream>
#include <complex>
#include <omp.h>

using namespace std;

// threads to request
const int REQUESTED_THREADS = 50;

// iterations to perform on a number to see
// if it is in the mandelbrot set
const int MAX_ITERATIONS = 100;

// dimensions of the image
const int WIDTH = 1000;
const int HEIGHT = 1000;


// calculate the color value to return for a given pixel location in the image
int colorizePixel (int x, int y)  {
    complex<float> point((float)x/WIDTH-1.5, (float)y/HEIGHT-0.5);
    complex<float> z(0, 0);
    int nb_iter = 0;
    
    
    while (abs (z) < 2 && nb_iter <= MAX_ITERATIONS) {
        z = z * z + point;
        nb_iter++;
    }
    
    if (nb_iter < MAX_ITERATIONS)
        return nb_iter;
    else
       return 0;
}



// write the PPM image file using our value function
int main()  {
    
    // used to normalize our histogram data
    double u[MAX_ITERATIONS];

    // used to time our parallel region
    double start1 = 0.0, start2 = 0.0, stop1 = 0.0, stop2 = 0.0;
    
    // set up the file out stream to ppm file
    ofstream my_Image("mandelbrot.ppm");
    

    // our image buffer, stores strings which contain
    // rows of RGB values
    string image_buffer[REQUESTED_THREADS];

    // used as a histogram to calculate u values
    double count[MAX_ITERATIONS] = { 0.0 };
    int *values[HEIGHT];
    for (int row = 0; row < HEIGHT; row++)
        values[row] = (int *)malloc(sizeof(int) * WIDTH);

    // start timing the first parallel region
    start1 = omp_get_wtime();

#   pragma omp parallel num_threads(REQUESTED_THREADS)
    {
        // the current thread's ID
        int threadNum = omp_get_thread_num();
        int val = 0;
        // total threads given
        int given_thread_ct = omp_get_num_threads();
        
        int th_start = (long)(threadNum * WIDTH)/given_thread_ct;
        int th_stop =(long)((threadNum+1) * WIDTH)/given_thread_ct;
        
        // we'll build the file by dividing it into multiple chunks
        // of rows and write the rows to an array of strings
        for (int row = th_start; row < th_stop; row++) {
            // we do WIDTH / 2 because the image is symmetrical about the Y axis, so we only need to
            // count one half the number of total points
            for (int col = 0; col < (WIDTH/2); col++) {
                val = colorizePixel(row, col);
                values[row][col] = val;
                values[row][WIDTH-col] = val;
#               pragma omp atomic
                count[val]+=2;
            }
        }
    }
    stop1 = omp_get_wtime();
    

    
    // create the histogram
    for (int i = 2; i < MAX_ITERATIONS; i++)
        count[i] += count[i-1];
    
    
    // get data normalized
    for (int i = 1; i < MAX_ITERATIONS; i++)
        u[i] = (float)(count[i] / count[MAX_ITERATIONS-1]);
    

    
    start2 = omp_get_wtime();
    
    
    if (my_Image.is_open ()) {

        my_Image << "P3\n" << WIDTH << " " << HEIGHT << " 255\n";
        my_Image << "# Derivative work of Farouk Ouname's Mandelbrot program\n";
        my_Image << "# Revised by: Daniel Carns\n";

    
#       pragma omp parallel num_threads(REQUESTED_THREADS)
        {
            // we'll use a string to write to the file later because it's very fast
            // and we have the RAM available to do so
            string rowString = "";
            int rgb = 0;
            int val = 0;
            
            // to be used to extract rgb values from our calculated value
            int r = 0, g = 0, b = 0;
            
            
            // the current thread's ID
            int threadNum = omp_get_thread_num();
            
            // total threads given
            int given_thread_ct = omp_get_num_threads();
            
            int th_start = (long)(threadNum * WIDTH)/given_thread_ct;
            int th_stop =(long)((threadNum+1) * WIDTH)/given_thread_ct;
        
            
            // we'll build the file by dividing it into multiple chunks
            // of rows and write the rows to an array of strings
            for (int row = th_start; row < th_stop; row++) {
                rowString = "";
                for (int col = 0; col < WIDTH; col++) {
                    
                    // lets write 5 sets of RGB values per line
                    // because that's what is excected in each row
                    val = values[row][col];
                    
                    rgb = 16777215 * u[val];
                    
                    r = ((rgb & 0x00FF0000) >> 16);
                    g = 0;
                    b = 0;
                    
                    if ((col+1) % 5 == 0)
                        rowString += to_string(r) + " " + to_string(g) + " " + to_string(b) + "\n";
                    else
                        rowString += to_string(r) + " " + to_string(g) + " " + to_string(b) + " ";

                }
                image_buffer[threadNum] += rowString;
            }

            
        }
        // our parallel region is done so get the stop time
        stop2 = omp_get_wtime();
      

        // write the file serially
        for (int i = 0; i < REQUESTED_THREADS; i++)
            my_Image << image_buffer[i];
        
        my_Image.close();

    
        for (int row = 0; row < HEIGHT; row++)
            free(values[row]);
    
   
    } else {
        // if we cant open the file, let the user know
        cerr << "Could not open the file, failed execution." << endl;
        exit(1);
    }
    
    cout << "Execution time 1st region: " << (stop1 - start1) << endl;
    cout << "Execution time 2nd region: " << (stop2 - start2) << endl;
    cout << "Execution time minus file writing: " << (stop2 - start1) << endl;
    cout << "Program executed normally.\n";
    return 0;
}
