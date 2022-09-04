#include <iostream>
#include <thread>

class State{
    public:
        State(int x, int y): x(x), y(y){}
        int x;
        int y;
};

void start_thread(int x, int y){
    State state(x, y);
    int i = 0;
    while(i < 10){
        std::cout << "x: " << state.x << " y: " << state.y << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(1));
        i++;
    }
}

void out_func(int x, int y){
    std::thread t0(start_thread, x, y);
    std::cout << "Thread started" << std::endl;
    t0.join();
}

