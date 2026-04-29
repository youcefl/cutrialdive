/*
* Youcef Lemsafer
* 2026.04.27
*/
#pragma once

#include <ostream>
#include <string>
#include <chrono>

class timer
{
public:
    timer(std::string prefix, std::ostream & out);
    ~timer();
private:
    std::string prefix_;
    std::ostream & out_;
    decltype(std::chrono::high_resolution_clock::now()) start_time_;
};
