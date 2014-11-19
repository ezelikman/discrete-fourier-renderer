#include <complex>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <OpenGL/GL.h>
#include <SDL2/SDL.h>
#include <stdio.h>

using namespace std;
typedef complex<double> compDouble;

vector<compDouble> slicer(vector<compDouble>& x, bool parity) // Return every other odd or even value in a vector
{
    vector<compDouble> slicedval(x.size()/2); // Creates an array with half the size of the original
    int b=0;
    for (int ii = (int)parity; ii<x.size();ii+=2) slicedval[b++]=x[ii]; // Defines the new array as the set of values indexed by whether it's odd or even.
    return slicedval;
}

double extremaArray (vector<compDouble> datain, double inverseTransformn, string extremaType)
{
    double maxivalf=0;
    double minivalf=0;
    for (int i = 0; i < datain.size(); i++) // Defines maxvalf and minvalf as the minimum values of the array
    {
        double realmaxx=real(datain[i])/(datain.size());
        double imagmaxy=inverseTransformn*imag(datain[i])/(datain.size());
        maxivalf=max(maxivalf,realmaxx);
        maxivalf=max(maxivalf,imagmaxy);
        minivalf=min(minivalf,realmaxx);
        minivalf=min(minivalf,imagmaxy);
    }
    if (extremaType=="min") return minivalf;
    else if (extremaType=="max") return maxivalf;
    throw;
}

vector<compDouble> loadArray(string filename) // Load the complex vector array's real part from a file
{
    
    ifstream stream(filename, ios::in | ios::binary);
    vector<double> realData;
    if (stream.is_open())
    {
        string line;
        while (getline(stream,line, ' ')) // Delimiter definition
        {
            double d = atof(line.c_str());
            realData.push_back(d);
        }
    }
    else // Alert if unable to open file
    {
        cout << "Input file not found" << endl;
        throw;
    }
    vector<compDouble> x(realData.size());
    for (int i = 0; i<x.size();i++) x[i]=compDouble(realData[i],0);
    return x;
}

void FT(vector<compDouble>& x)
{
    const double N = x.size();
    if (fmod(log(N)/log(2.0),1)>.00001) // Abandon if not a power of 2,
    {
        cout << N << ": Invalid size, must be a power of 2" << endl;
        throw;
    }
    if (N <= 1) return; // The function starts at a function of two values
    vector<compDouble> even = slicer(x,false); // Defines the even-numbered components of x
    vector<compDouble> odd  = slicer(x,true); //  Defines the odd -numbered components of x
    
    FT(even); // Performs transform on the even subsection of the function
    FT(odd); //  Performs transform on the odd  subsection of the function
    for (double k = 0; k < N/2; k++)
    {
        complex<double> t = polar(1.0,-6.28318530717958*k/N)*odd[k]; // Rectangular form of the addition of the polar angles of odd-numbered values by corresponding circle fractions.
        x[k]     = even[k] + t; // Replaces the first  half of the values with the even values in order, with the addition    of the polar modified odd values.
        x[k+N/2] = even[k] - t; // Replaces the second half of the values with the even values in order, with the subtraction of the polar modified odd values.
    }
}

int main(int argc, const char* argv[])
{
    double precOut=10; // Output precision
    string fileLoc="";
    if (argc>=2) fileLoc = string(argv[1]); // Tests for the passing of argument corresponding to non-default datafile.
    vector<compDouble> data = loadArray((fileLoc=="") ? "data.txt" : fileLoc); // Loads the datafile to data (Possible TODO: Add support for complex input)
    double asize=data.size(); // Record array size
    string conNext=""; // Next parameters
    if (argc>=3) conNext=       string(argv[2]); // Defines inverse parameter if a third argument exists
    double inverseTransform =   (conNext=="ifft" || conNext=="inverse" || conNext=="ift" || conNext=="-i") ? -1 : 1; // Tests if there is a flag for inverse transform
    bool approximation =        !(conNext=="-p" || conNext=="precis" || conNext=="precise" || conNext=="p"); // Tests if there is a flag for approximation
    if (argc>=4)
    {
        conNext=string(argv[3]);
        if (inverseTransform!=-1) inverseTransform = (conNext=="ifft" || conNext=="inverse" || conNext=="ift" || conNext=="-i") ? -1 : 1;
        approximation |= (conNext=="-a" || conNext=="approx" || conNext=="approximate" || conNext=="a");
    }
    
    double winsizex=1368;
    double winsizey=725;
    if (argc>=5) // Fetches dimensions if four parameters or more
    {
        cout << "Window width ";
        cin >> winsizex;
        cout << "Window height ";
        cin >> winsizey;
    }
    ofstream fout("foutn.txt");
    SDL_Window *window;
    SDL_Init (SDL_INIT_VIDEO); // Initialize SDL
    SDL_Renderer *renderer;
    
    SDL_CreateWindowAndRenderer(winsizex, winsizey, 0, &window, &renderer); // Create window with renderer given size (Infinitely scalable)
    double maxvalf=extremaArray(data,inverseTransform, "max"); // Array extrema defined
    double minvalf=extremaArray(data,inverseTransform, "min");
    
    SDL_SetRenderDrawColor(renderer,255,0,0,255); // Defines drawing color of original points
    for (int i = 0; i < asize; i++) // Draws the points with a line connecting them
    {
        if (i>0) SDL_RenderDrawLine(renderer,
                                    (i-1)*winsizex/asize, winsizey*(1-(real(data[i-1])/asize-minvalf)/(maxvalf-minvalf)),
                                    (i  )*winsizex/asize, winsizey*(1-(real(data[i  ])/asize-minvalf)/(maxvalf-minvalf)));
    }
    
    SDL_SetRenderDrawColor(renderer,0,100,100,122); // x-axis color
    SDL_RenderDrawLine(renderer, 0, winsizey-winsizey*(-minvalf)/(maxvalf-minvalf), winsizex, winsizey-winsizey*(-minvalf)/(maxvalf-minvalf)); //x-axis

    SDL_SetRenderDrawColor(renderer,30,30,0,20); // Data marker line color
    if (asize<winsizex/2) for (int i = 0; i < asize; i++) SDL_RenderDrawLine(renderer, i*winsizex/asize, 0, i*winsizex/asize, winsizey); // Marks datapoints with lines
    
    FT(data); // The transform is performed
    
    double maxval=extremaArray(data,inverseTransform, "max"); // New extreme found
    double minval=extremaArray(data,inverseTransform, "min");
    
    for (int i = 0; i < asize; i++)
    {
        fout<<floor(.5+precOut*real(data[i]))/precOut<<"Cos[-160*Pi*"<<i<<"*t]+"<<floor(.5+inverseTransform*precOut*imag(data[i]))/precOut<<"Sin[-160*Pi*"<<i<<"*t]+"; // Writes values to file with given precision
    }
    
    SDL_SetRenderDrawColor(renderer,50,100,200,60); // Transform axis color
    SDL_RenderDrawLine(renderer, 0, winsizey-1+winsizey*minval/(maxval-minval), winsizex, winsizey-1+winsizey*minval/(maxval-minval)); // Draws the transform axis

    for (int i = 0; i < asize; i++) // Draws the transform
    {
        SDL_SetRenderDrawColor(renderer,200,100,50,60); // Color of the real or cosine part
        SDL_RenderDrawLine(renderer,
            i*winsizex/asize+1,     winsizey*(1+minval/(maxval-minval))-1,
            i*winsizex/asize+1,     winsizey*(1-(                 real(data[i])/asize-minval)/(maxval-minval))-1
                           );
        SDL_SetRenderDrawColor(renderer,100,200,0,60);  // Color of the imaginary or sine part
        SDL_RenderDrawLine(renderer,
            i*winsizex/asize+3,     winsizey*(1+minval/(maxval-minval))-1,
            i*winsizex/asize+3,     winsizey*(1-(inverseTransform*imag(data[i])/asize-minval)/(maxval-minval))-1
                           );
    }
    
    double sumxnp = 0;
    double modn=1+approximation; // Notes the modification depending on approximation
    SDL_SetRenderDrawColor(renderer,0,0,255,20); // Sets Fourier series curve color to blue
    for (int i = 0; i < winsizex;i++) // Draws a line for every pixel
    {
        double sumxn=0;
        for (int l = 0; l<asize/modn;l++) // Evaluation of the series formula at a point
        {
            sumxn+=real(data[l])*cos(-2*3.14159265*i*l/(winsizex));
            sumxn+=imag(data[l])*sin(-2*3.14159265*i*l/(winsizex));
        }
        if (i>0) SDL_RenderDrawLine(renderer, // Drawing of the evaluated formula
            i-1,    winsizey*((modn+1)/2-modn*((sumxnp/(asize*asize)-minvalf)/(maxvalf-minvalf))), // Scaling by window and min and max of original
            i,      winsizey*((modn+1)/2-modn*((sumxn /(asize*asize)-minvalf)/(maxvalf-minvalf)))
                                    );
        sumxnp=sumxn;
    }
        SDL_RenderPresent(renderer);
                if (window == NULL) {
                printf ("Could not create window:%s \n", SDL_GetError ());
                return 1; 
            }
    bool open=true;
    while (open)
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
            if (event.type==SDL_QUIT)
            {
                open=!open;
            }
    }
        return 0;
    }
