/*
* Youcef Lemsafer
* 2026.04.27
*/
#include "timer.hpp"

timer::timer(std::string prefix, std::ostream & out)
    : prefix_(prefix)
    , out_(out)
    , start_time_(std::chrono::high_resolution_clock::now())
{
}

timer::~timer()
{
    out_ << prefix_ << std::chrono::duration<double>(
                std::chrono::high_resolution_clock::now() - start_time_) << std::endl;
}
