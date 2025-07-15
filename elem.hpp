#ifndef ELEM_HPP
class Elem {
    public:
        Elem(int i, int j) : x_left(i*h), x_right((i+1)*h), y_down(j*h), y_up((j+1)*h) {
            x_center = (x_left + x_right) / 2;
            y_center = (y_down + y_up) / 2;
        }
        void ref2phy(const Coord ref_coord, Coord &phy_coord) {
            phy_coord[0] = ref_coord[0]*h + x_left; 
            phy_coord[1] = ref_coord[1]*h + y_down; 
        }

        int edge_sign[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
        double x_left;
        double x_right;
        double x_center;
        double y_down;
        double y_up;
        double y_center;
};

#define ELEM_HPP
#endif

