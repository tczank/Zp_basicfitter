#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include <cmath>

void show_progress_bar(int time, const std::string &message, char symbol)
{
    std::string progress_bar;
    const double progress_level = 1.42;

    std::cout << message << "\n\n";

    for (double percentage = 0; percentage <= 100; percentage += progress_level)
    {
        progress_bar.insert(0, 1, symbol);
        std::cout << "\r [" << std::ceil(percentage) << '%' << "] " << progress_bar;
        std::this_thread::sleep_for(std::chrono::milliseconds(time));
    }
    std::cout << "\n\n";
}
