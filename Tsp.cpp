// Student: Thy Vu
// CSS 534 Autumn 2016
// Professor Munehiro Fukuda

#include <iostream>  // cout
#include <fstream>   // ifstream
#include <string.h>  // strncpy
#include <stdlib.h>  // rand
#include <math.h>    // sqrt, pow
#include <omp.h>     // OpenMP
#include "Timer.h"
#include "Trip.h"

using namespace std;

// Already implemented. see the actual implementations below
void initialize( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] );
void select( Trip trip[CHROMOSOMES], Trip parents[TOP_X] );
void populate( Trip trip[CHROMOSOMES], Trip offsprings[TOP_X] );
extern int compareTrip(const void* _tripA, const void* _tripB);

// Implemented by Thy Vu
extern void evaluate( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] );
extern void crossover( Trip parents[TOP_X], Trip offsprings[TOP_X], int coordinates[CITIES][2] );
extern double getDistance( char cityA, char cityB, int coordinates[][2]);
extern char getUnvisitedCity(bool visitedCities[CITIES]);
extern void mutate( Trip offsprings[TOP_X] );

/*
 * MAIN: usage: Tsp #threads
 */
int main( int argc, char* argv[] ) {
  Trip trip[CHROMOSOMES];       // all 50000 different trips (or chromosomes)
  Trip shortest;                // the shortest path so far
  int coordinates[CITIES][2];   // (x, y) coordinates of all 36 cities:
  int nThreads = 1;
  
  // verify the arguments
  if ( argc == 2 )
    nThreads = atoi( argv[1] );
  else {
    cout << "usage: Tsp #threads" << endl;
    if ( argc != 1 )
      return -1; // wrong arguments
  }
  cout << "# threads = " << nThreads << endl;

  // shortest path not yet initialized
  shortest.itinerary[CITIES] = 0;  // null path
  shortest.fitness = -1.0;         // invalid distance

  // initialize 5000 trips and 36 cities' coordinates
  initialize( trip, coordinates );

  // start a timer 
  Timer timer;
  timer.start( );

  // change # of threads
  omp_set_num_threads( nThreads );

  // find the shortest path in each generation
  for ( int generation = 0; generation < MAX_GENERATION; generation++ ) {

    // evaluate the distance of all 50000 trips
    evaluate( trip, coordinates );

    // just print out the progress
    if ( generation % 20 == 0 )
      cout << "generation: " << generation << endl;

    // whenever a shorter path was found, update the shortest path
    if ( shortest.fitness < 0 || shortest.fitness > trip[0].fitness ) {

      strncpy( shortest.itinerary, trip[0].itinerary, CITIES );
      shortest.fitness = trip[0].fitness;

      cout << "generation: " << generation 
       << " shortest distance = " << shortest.fitness
       << "\t itinerary = " << shortest.itinerary << endl;
    }




    //define TOP_X parents and offsprings.
    Trip parents[TOP_X], offsprings[TOP_X];

    // choose TOP_X parents from trip
    select( trip, parents );

    // generates TOP_X offsprings from TOP_X parenets
    crossover( parents, offsprings, coordinates );
    

    // mutate offsprings
    mutate( offsprings );

    // populate the next generation.
    populate( trip, offsprings );
     
    
  }

  // stop a timer
  cout << "elapsed time = " << timer.lap( ) << endl;
  return 0;
}



/*
 * Initializes trip[CHROMOSOMES] with chromosome.txt and coordiantes[CITIES][2] with cities.txt
 *
 * @param trip[CHROMOSOMES]:      50000 different trips
 * @param coordinates[CITIES][2]: (x, y) coordinates of 36 different cities: ABCDEFGHIJKLMNOPQRSTUVWXYZ
 */
void initialize( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] ) {
  // open two files to read chromosomes (i.e., trips)  and cities
  ifstream chromosome_file( "chromosome.txt" );
  ifstream cities_file( "cities.txt" );
  
  // read data from the files
  // chromosome.txt:                                                                                           
  //   T8JHFKM7BO5XWYSQ29IP04DL6NU3ERVA1CZG                                                                    
  //   FWLXU2DRSAQEVYOBCPNI608194ZHJM73GK5T                                                                    
  //   HU93YL0MWAQFIZGNJCRV12TO75BPE84S6KXD
  for ( int i = 0; i < CHROMOSOMES; i++ ) {
    chromosome_file >> trip[i].itinerary;
    trip[i].fitness = 0.0;
  }

  // cities.txt:                                                                                               
  // name    x       y                                                                                         
  // A       83      99                                                                                        
  // B       77      35                                                                                        
  // C       14      64                                                                                        
  for ( int i = 0; i < CITIES; i++ ) {
    char city;
    cities_file >> city;
    int index = ( city >= 'A' ) ? city - 'A' : city - '0' + 26;
    cities_file >> coordinates[index][0] >> coordinates[index][1];
  }

  // close the files.
  chromosome_file.close( );
  cities_file.close( );

  // just for debugging
  if ( DEBUG ) {
    for ( int i = 0; i < CHROMOSOMES; i++ )
      cout << trip[i].itinerary << endl;
    for ( int i = 0; i < CITIES; i++ )
      cout << coordinates[i][0] << "\t" << coordinates[i][1] << endl;
  }
}


/*
 * Based on coordinates of each city, evaluate total distance for each trip's itinerary
 * (Trip starts from coordinate (0,0)) .
 * and assign that value as the trip's fitness value.
 * Sort the trips from shortest to longest.
 * 
 * @param trip[CHROMOSOMES]:       50000 different trips
 * @param coordinates[CITIES][2]:  (x, y) coordinates of 36 different cities
 */
extern void evaluate( Trip trip[CHROMOSOMES], int coordinates[CITIES][2] ){
  
    #pragma omp parallel for shared(trip)
    // calculate total distance for each trip, starting from (0,0)
    for (int i = 0; i < CHROMOSOMES; i++){
        double totalDistance = 0.0;
        
        // trip starts at (0,0)
        int prevX = 0;  
        int prevY = 0; 

        for (int j = 0; j < CITIES; j++){
            // calculate distance to current city based on its and last city's coordinates
            char currCity = trip[i].itinerary[j];            
            int cityIndex = ( currCity >= 'A' ) ? currCity - 'A' : currCity - '0' + 26;
            int deltaX = coordinates[cityIndex][0] - prevX;
            int deltaY = coordinates[cityIndex][1] - prevY; 
            prevX = coordinates[cityIndex][0];
            prevY = coordinates[cityIndex][1];
            totalDistance += sqrt((double)(deltaX * deltaX + deltaY * deltaY));
        }
        
        // set this trip's fitness
        trip[i].fitness = (float) totalDistance;
             
    }
    
    
    // sort trips from shortest to longest
    qsort((void*) trip, CHROMOSOMES, sizeof(trip[0]), compareTrip);
    
}



/*
 * Function to use for qsort() used for Trips : Return 1 if tripA is longer than tripB, 
 * or -1 if it is shorter than tripB, or 0 if they are equal.
 * 
 * @param _tripA:   pointer to a Trip object with assigned total distance (fitness)
 * @param _tripB:   pointer to a Trip object with assigned total distance (fitness)
 * Credit: this function was created with callback tutorial at www.newty.de/fpt/callback.html
 */
extern int compareTrip(const void* _tripA, const void* _tripB){
    const Trip* tripA = (const Trip*) _tripA;
    const Trip* tripB = (const Trip*) _tripB;
    
    if (tripA->fitness > tripB->fitness) {
        return 1;
    }else{
        if(tripA->fitness == tripB->fitness)
            return 0;
        else
            return -1;
    }
}


/*
 * Select the first TOP_X parents from trip[CHROMOSOMES]
 *
 * @param trip[CHROMOSOMES]: all trips
 * @param parents[TOP_X]:    the firt TOP_X parents
 */
void select( Trip trip[CHROMOSOMES], Trip parents[TOP_X] ) {
  // just copy TOP_X trips to parents
  for ( int i = 0; i < TOP_X; i++ )
    strncpy( parents[i].itinerary, trip[i].itinerary, CITIES + 1 );
}



/*
 * Replace the bottom TOP_X trips with the TOP_X offsprings
 */
void populate( Trip trip[CHROMOSOMES], Trip offsprings[TOP_X] ) {
  // just copy TOP_X offsprings to the bottom TOP_X trips.
  for ( int i = 0; i < TOP_X; i++ )
    strncpy( trip[ CHROMOSOMES - TOP_X + i ].itinerary, offsprings[i].itinerary, CITIES + 1 );

  // for debugging
  if ( DEBUG ) {
    for ( int chrom = 0; chrom < CHROMOSOMES; chrom++ ) 
      cout << "chrom[" << chrom << "] = " << trip[chrom].itinerary 
       << ", trip distance = " << trip[chrom].fitness << endl;
  }
}



/*
 * Generate 25000  offsprings from 25000 parents. For a pair of parents,
 * their 2 offsprings are complementary of each other.
 *
 * @param parents[TOP_X]:          25000 different parent trips
 * @param offsprings[TOP_X]:       (to be populated with) 25000 offspring trips
 * @param coordinates[CITIES][2]:  (x, y) coordinates of 36 different cities
 */
extern void crossover( Trip parents[TOP_X], Trip offsprings[TOP_X], int coordinates[CITIES][2] ){

    // complementary strands: e.g. A complements 9, K complements Z
    string cities(    "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789");
    string compCities("9876543210ZYXWVUTSRQPONMLKJIHGFEDCBA");

    #pragma omp parallel for firstprivate(cities,compCities) shared(offsprings)
 
    // generate 2 offsprings for each parent pair
    for (int i = 0; i < TOP_X; i+=2){ 
        
        string parentAItin, parentBItin; // parent[i]'s itinerary, parent[i+1]'s itinerary
        parentAItin.assign(parents[i].itinerary, CITIES+1);
        parentBItin.assign(parents[i+1].itinerary, CITIES+1);

        bool visitedCities[CITIES];   // for picking a random unvisited city
        
        // offspring[i] has not visited any city, so set false for all cities
        for (int k = 0; k < CITIES; k++){
            visitedCities[k] = false;
        }
        
        
        // first city offspring[i] visits
        char currCity;
        currCity = parents[i].itinerary[0];
        offsprings[i].itinerary[0] = currCity;
        
        // offspring[i+1] visits the complementary city of offspring[i] visited
        offsprings[i+1].itinerary[0] = compCities.at(cities.find(currCity));

        
        string childItin;      // itinerary of offspring[i] to keep track of visited cities
        childItin.append(1,currCity); 
        int currCityIndex = ( currCity >= 'A' ) ? currCity - 'A' : currCity - '0' + 26;
        visitedCities[currCityIndex] = true;    
        
        int currX= coordinates[currCityIndex][0];
        int currY= coordinates[currCityIndex][1];
        

        // check parents' next cities
        for (int j = 1; j < CITIES; j++){
        
            
            char parentANext,parentBNext;   // next city according to each parent's itinerary
            bool parentAItinEnded, parentBItinEnded;  // true if current city is parent's last
            
            // Find parent A's next city from current city 
            // If current city is the last city, set bool flag
            if (parentAItin.find(currCity) < CITIES-1) { 
                parentANext = parentAItin.at(parentAItin.find(currCity)+1);
                parentAItinEnded = false;
            }else{  
                parentAItinEnded = true;
            }
            
            // Find parent B's next city from current city 
            // If current city is the last city, set bool flag
            if (parentBItin.find(currCity) < CITIES-1) {    // Parent B
                parentBNext = parentBItin.at(parentBItin.find(currCity)+1);
                parentBItinEnded = false;
            }else{  
                parentBItinEnded = true;
            }
            
            char nextCity;      // next city offspring[i] must visit
            int nextCityIndex;  // index on coordinates[CITIES][2] of next city
            
            // Edge cases: either or both parents' itinerary ended
            if (parentAItinEnded){
                if(parentBItinEnded ||  childItin.find(parentBNext) !=  string::npos){
                    // Parent B ended or child has visited B's next.
                    nextCity = getUnvisitedCity(visitedCities);
                }else{
                    // child has not visited Parent B's next 
                    nextCity = parentBNext;
                }
            } else if (parentBItinEnded){
                if (parentAItinEnded ||  childItin.find(parentANext) !=  string::npos){
                    // Parent A ended or child has visited A's next.
                    nextCity = getUnvisitedCity(visitedCities);
                }else{
                    // Child has not visit Parent A's next
                    nextCity = parentANext;
                }
                
            // Non-edge cases: neither parents' itinerary ended
            }else if (parentANext == parentBNext && childItin.find(parentANext) == string::npos){
                // if parent A's and parent B's next city is same and unvisited, use that city 
                nextCity = parentANext;
                
            }else if( childItin.find(parentANext) !=  string::npos) { // child has visited A's next
                
                if( childItin.find(parentBNext) != string::npos){ // child has visited B's & A's nexts
                    nextCity = getUnvisitedCity(visitedCities);
                    
                }else{ // child didn't visit B's next (but visited A's next)
                    nextCity = parentBNext;
                }
            } else if ( childItin.find(parentBNext) != string::npos) {
                // child visited B's next but did not visit A's next
                nextCity = parentANext;
            
            } else { // child did not visit B's next nor A's next
                // choose the city closer to current city
                if (getDistance (currCity, parentANext, coordinates) < getDistance (currCity, parentBNext, coordinates)){
                    nextCity = parentANext;
                } else { 
                    // distance to either cities are equal or currCity to B's next is shorter
                    nextCity = parentBNext;
                }
            
            }
            
            // set next city for offspring[i]'s and [i+1]'s itineraries
            offsprings[i].itinerary[j] = nextCity;
            offsprings[i+1].itinerary[j] = compCities.at(cities.find(nextCity));
            
            
            // on visitedCities: set true (which means offspring[i] visited) for next city
            nextCityIndex = ( nextCity >= 'A' ) ? nextCity - 'A' : nextCity - '0' + 26;
            visitedCities[nextCityIndex] = true;
            
            // update string of visited cities
            childItin.append(1,nextCity);

            // update currCity for next iteration
            currCity = nextCity;
            currCityIndex = nextCityIndex;

        }
  
    }

}


/*
 * Return distance value between cityA and cityB according to their coordinates.
 *
 * @param cityA:   a char with a value from 'A' to 'Z' or '0' to '9'
 * @param cityB:   a char with a value from 'A' to 'Z' or '0' to '9'
 * @param coordinates[CITIES][2]: (x, y) coordinates of 36 different cities
 */
extern double getDistance( char cityA, char cityB, int coordinates[][2]){

    // find corresponding index on coordinates array for each city
    int cityAIndex = ( cityA >= 'A' ) ? cityA - 'A' : cityA - '0' + 26; 
    int cityBIndex = ( cityB >= 'A' ) ? cityB - 'A' : cityB - '0' + 26;

    // find delta X and delta Y for 2 cities
    int deltaX = coordinates[cityAIndex][0] - coordinates[cityBIndex][0];
    int deltaY = coordinates[cityAIndex][1] - coordinates[cityBIndex][1]; 

    // return distance between 2 cities
    return sqrt((double)(deltaX * deltaX + deltaY * deltaY));

}

/*
 * Return a random (char) unvisited city from visitedCities bool array
 *
 * @param visitedCities[CITIES]:  36 cities, false value = unvisited, true value = visited
 */
extern char getUnvisitedCity(bool visitedCities[CITIES]){
    string allUnvisited;   // string contains unvisited cities only
    char unvisited;        // one unvisited city
    
    // concatenate all unvisited city to allUnvisited string
    for (int i = 0; i < CITIES; i++){
        if (!visitedCities[i]){
            unvisited = ( i < 26? 'A'+ i : i-26+'0');
            allUnvisited.append(1, unvisited);
        }
    }
    
    // get random index on allUnvisited string
    int randomIndex = rand() % (allUnvisited.length());
    
    // return a random unvisited city
    return allUnvisited.at(randomIndex);
}

/*
 * For each offspring's itinerary, randomly choose 2 cities with a given probability
 * and swap them.
 * 
 * @param offsprings[TOP_X]:      25000 different trips
 */
extern void mutate( Trip offsprings[TOP_X] ){
    for (int i = 0; i < TOP_X; i++){
        if( rand()%100 < MUTATE_RATE) {
            // get 2 random city indexes
            int city1Index = rand()%36;
            int city2Index = rand()%36;
            
            // swap cities
            char city1 = offsprings[i].itinerary[city1Index];
            offsprings[i].itinerary[city1Index] = offsprings[i].itinerary[city2Index];
            offsprings[i].itinerary[city2Index] = city1;
        }
    
    }
}

