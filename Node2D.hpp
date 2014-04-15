#ifndef NODE2DHEADERDEF
#define NODE2DHEADERDEF

#include "Node.hpp"

class Node2D : public Node {
public:
    struct coordinates {
        double x;
        double y;
    };
    coordinates N, S, W, E, C; // North, south, west, east, center
};

#endif

