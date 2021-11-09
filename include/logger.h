#ifndef LOG_HPP
#define LOG_HPP

#include <iostream>
#include <utility>

void log();

template<typename First, typename ...Rest> 
 void log(bool verbose, First && first, Rest && ...rest);

template<typename First, typename ...Rest>
 void log(First const && first, Rest const && ...rest);

 template<typename First, typename ...Rest>
 void log(First && first, Rest && ...rest);

template<typename First, typename ...Rest> 
 void log(bool const verbose, First const && first, Rest const && ...rest); 

#endif //LOG_HPP