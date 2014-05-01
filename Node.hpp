#ifndef NODEHEADERDEF
#define NODEHEADERDEF

class Node {
public:
    struct coordinates {
        double x;
        double y;
        double z;
        double N; /*# of the node*/
    };
    coordinates N, S, W, E, C; // North, south, west, east, center
};

#endif