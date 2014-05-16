#ifndef NODEHEADERDEF
#define NODEHEADERDEF

class Node {
public:
    struct coordinates {
        double x;
        double y;
        double z;
    };
    coordinates N, S, W, E, C; // North, south, west, east, center
    int num; /*# of the node*/
    double value; // value of u for boundary
};

#endif