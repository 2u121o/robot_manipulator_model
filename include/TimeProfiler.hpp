#ifndef TIME_PROFILER_HPP
#define TIME_PROFILER_HPP

#include <iostream>
#include <chrono>

namespace profiler
{

class TimeProfiler
{
    public:
        TimeProfiler(const std::string &function_name);
        TimeProfiler(const TimeProfiler&) = delete;
        TimeProfiler(TimeProfiler&&) = delete;
        TimeProfiler& operator=(const TimeProfiler&) = delete;
        TimeProfiler& operator=(TimeProfiler&&) = delete;

        void printCurrentTime();

        ~TimeProfiler();
    private:
        const std::string FUNCTION_NAME;
        const std::chrono::steady_clock::time_point START;

};

}//profiler

#if USE_PROFILER
#define MEASURE_FUNCTION() profiler::TimeProfiler timer{__func__};
#else
#define MEASURE_FUNCTION()
#endif

#endif