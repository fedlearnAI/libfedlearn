#include "../include/logger.h"

void log(){};

template<typename First, typename ...Rest>
 void log(First && first, Rest && ...rest)
{   
        std::cout << std::forward<First>(first);
        log(std::forward<Rest>(rest)...);

}

template<typename First, typename ...Rest>
 void log(First const && first, Rest const && ...rest)
{   
        std::cout << std::forward<First>(first);
        log(std::forward<Rest>(rest)...);

}

template<typename First, typename ...Rest> 
 void log(bool verbose, First && first, Rest && ...rest)
{   if(verbose) {
        std::cout << std::forward<First>(first);
        log(std::forward<Rest>(rest)...);
    }
}

template<typename First, typename ...Rest> 
 void log(bool const verbose, First const && first, Rest const && ...rest) {
     if(verbose) {
        std::cout << std::forward<First>(first);
        log(std::forward<Rest>(rest)...);
    }
 }