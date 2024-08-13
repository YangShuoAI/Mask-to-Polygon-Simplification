#include "apsc.h"
#include <vector>
#include <iostream>
#include "OpencvLibrary.h"

int initVertexTree(const std::vector<cv::Point2f>& pts, VertexRelation *vertex_tree)
{
    //std::vector<int> ID;
    //std::vector<cv::Point2f> XY;
    //std::vector<double> error;
    //std::vector<int> parent;
    //std::vector<int> LC;
    //std::vector<int> RC;
    //std::vector<int> LS;
    //std::vector<int> RS;
    //std::vector<bool> current;

    int n = (int)pts.size();

    // Initialize default values for each vertex
    for (int i = 0; i < n; i++)
    {
        vertex_tree->ID.push_back(i);
        vertex_tree->XY.push_back(pts[i]);
        vertex_tree->error.push_back(0);
        vertex_tree->parent.push_back(-1);
        vertex_tree->LC.push_back(-1);
        vertex_tree->RC.push_back(-1);
        vertex_tree->LS.push_back(i - 1);
        vertex_tree->RS.push_back(i + 1);
        vertex_tree->current.push_back(true);
    }

    // If not a closed loop, set start/end left/right siblings to -1
    vertex_tree->RS[n - 1] = -1;
    vertex_tree->LS[0] = -1;

    return 0;
}

int clockwise(cv::Point2f A, cv::Point2f B, cv::Point2f C)
{
    // Determines if points A,B & C are sequenced clockwise around a triangle.
    // Args :
    // A, B, C : (x, y) tuples
    // Returns :
    // True if points are sequenced clockwise.Result is unpredictable if points
    // are collinear.
    //std::cout << ((C.y - A.y) * (B.x - A.x)) << " " << ((B.y - A.y) * (C.x - A.x)) << std::endl;
    return ((C.y - A.y) * (B.x - A.x)) > ((B.y - A.y) * (C.x - A.x));

}

int intersect(cv::Point2f A, cv::Point2f B, cv::Point2f C, cv::Point2f D)
{
    // Quickly determines if segment AB intersects segment CD.
    // Args :
    // A, B, C, D : (x, y) tuples
    // Returns:
    // True if segments intersect.Result is unpredictable if
    // lines are parallel or segments meet at endpoint.

    //std::cout << clockwise(A, C, D) << std::endl;
    //std::cout << clockwise(B, C, D) << std::endl;
    //std::cout << clockwise(A, B, C) << std::endl;
    //std::cout << clockwise(A, B, D) << std::endl;
    //std::cout << (clockwise(A, C, D) != clockwise(B, C, D)) && (clockwise(A, B, C) != clockwise(A, B, D));


    return (clockwise(A, C, D) != clockwise(B, C, D)) && (clockwise(A, B, C) != clockwise(A, B, D));
}

cv::Point2f intersection(cv::Point2f A, cv::Point2f B, cv::Point2f C, cv::Point2f D, bool infinite = true)
{
    // Returns the intersection point of two lines AB & CD.
    // A, B, C, D and return value are all lists of two coordinates each.
    // If lines are parallel or do not intersect, returns a pair of Nones.
    // Code taken from Stephen Wise, pp. 48 - 9
    // Args :
    // A, B, C, D : (x, y) tuples
    // infinite : If False, will return (None, None)
    // Returns :
    // (x, y) tuple of coordinates of intersection point
    // if ultra_precise :
    //     A = (Decimal(A[0]), Decimal(A[1])) # dec_pt(A);
    //     B = (Decimal(B[0]), Decimal(B[1])) # dec_pt(B);
    //     C = (Decimal(C[0]), Decimal(C[1])) # dec_pt(C);
    //     D = (Decimal(D[0]), Decimal(D[1])) # dec_pt(D);
    float xp = 0.0f;
    float yp = 0.0f;

    if (A.x == B.x)
    {
        if (C.x == D.x)
        {
            xp = -1; // lines are parallel
            yp = -1;
        }
        else // first line vertical
        {
            float b2 = (D.y - C.y) / (D.x - C.x);
            float a2 = C.y - b2 * C.x;
            xp = A.x;
            yp = a2 + b2 * xp;
        }

    }
    else
    {
        if (C.x == D.x) // second line vertical
        {
            if ((B.y == -1) || (A.y == -1) || (B.x == -1) || (A.x == -1))
            {
                std::cout << "yeah!";
            }

            float b1 = (B.y - A.y) / (B.x - A.x);
            float a1 = A.y - b1 * A.x;
            xp = C.x;
            yp = a1 + b1 * xp;
        }
        else // neither line vertical
        {
            float b1 = (B.y - A.y) / (B.x - A.x);
            float b2 = (D.y - C.y) / (D.x - C.x);
            float a1 = A.y - b1 * A.x;
            float a2 = C.y - b2 * C.x;
            if (b1 == b2)
            {
                xp = -1;
                yp = -1; // lines are parallel
            }
            else
            {
                xp = -(a1 - a2) / (b1 - b2);
                yp = a1 + b1 * xp;
            }
        }
        


    }

    // test whether intersection point falls on either line
    if ((infinite == false) && (xp != -1))
    {
        if ((A.x - xp) * (xp - B.x) < 0 or (C.x - xp) * (xp - D.x) < 0 or (A.y - yp) * (yp - B.y) < 0 or (C.y - yp) * (yp - D.y) < 0)
        {
            xp = -1;
            yp = -1;
        }
    }

    if (xp != -1)
    {
        xp = float(xp);
    }
    if (yp != -1)
    {
        yp = float(yp);
    }
    return cv::Point2f(xp, yp);
}

float area(std::vector<cv::Point2f> pts, bool absolute = false)
{
    // Ensure the polygon is closed by appending the first point to the end if it's not already there
    if (pts.back() != pts.front())
    {
        pts.push_back(pts.front());
    }

    float A = 0.0f;
    for (size_t i = 0; i < pts.size() - 1; ++i)
    {
        A += (pts[i + 1].x - pts[i].x) * (pts[i].y + pts[i + 1].y);
    }
    A /= 2.0f;

    if (absolute)
    {
        return std::abs(A);
    }
    else
    {
        return A;
    }
}


float sideshift_area(cv::Point2f pA, cv::Point2f pB, cv::Point2f pC, cv::Point2f pD, cv::Point2f pE, int overlap_endpt)
{
    //Returns the sideshift displacement area from ABCD to ACE.
    //    Args:
    //pA, pB, pC, pD, pE : (x, y) tuples
    //    overlap_endpt :
    //0      E is placed on AB)
    //1      E is placed on CD)
    //other  E is not placed on AB or CD)
    //Returns:
    //Total sideshift displacement area between ABCD and AED
    if (intersect(pA, pB, pC, pD)) // original line self - intersects
    {
        cv::Point2f x = intersection(pA, pB, pC, pD);
        std::vector<cv::Point2f> pts1;
        std::vector<cv::Point2f> pts2;
        pts1.push_back(x);
        pts1.push_back(pC);
        pts1.push_back(pB);

        pts2.push_back(x);
        pts2.push_back(pA);
        pts2.push_back(pE);
        pts2.push_back(pD);
        return area(pts1, true) + area(pts2, true);
    }
    else if (overlap_endpt == 0)  // E is on line through AB
    {
        if (intersect(pB, pC, pE, pD))
        {
            cv::Point2f x = intersection(pB, pC, pE, pD);

            std::vector<cv::Point2f> pts1;
            std::vector<cv::Point2f> pts2;
            pts1.push_back(pB);
            pts1.push_back(pE);
            pts1.push_back(x);

            pts2.push_back(x);
            pts2.push_back(pD);
            pts2.push_back(pC);

            return area(pts1, true) + area(pts2, true);
        }
        else
        {
            std::vector<cv::Point2f> pts1;
            pts1.push_back(pD);
            pts1.push_back(pC);
            pts1.push_back(pB);
            pts1.push_back(pE);
            return area(pts1, true);
        }
    }
    else if (overlap_endpt == 1)  // E is on line through CD
    {
        if (intersect(pC, pB, pE, pA))
        {
            cv::Point2f x = intersection(pC, pB, pE, pA);

            std::vector<cv::Point2f> pts1;
            std::vector<cv::Point2f> pts2;
            pts1.push_back(pC);
            pts1.push_back(pE);
            pts1.push_back(x);

            pts2.push_back(x);
            pts2.push_back(pA);
            pts2.push_back(pB);

            return area(pts1, true) + area(pts2, true);
        }
        else
        {
            std::vector<cv::Point2f> pts1;
            pts1.push_back(pA);
            pts1.push_back(pB);
            pts1.push_back(pC);
            pts1.push_back(pE);
            return area(pts1, true);
        }
    }
    else // possibilities : (1) AE^ BC, (2) AE^ BC& ED^ BC, (3) AE^ CD, (4) DE^ AB, (5) DE^ BC, (6) no intersections
    {
        // check if AE intersects BC or CD
        if (intersect(pA, pE, pB, pC)) // AE intersects BC
        {
            cv::Point2f x1 = intersection(pA, pE, pB, pC);
            if (intersect(pD, pE, pB, pC)) // (2)
            {
                cv::Point2f x2 = intersection(pD, pE, pB, pC);

                std::vector<cv::Point2f> pts1;
                std::vector<cv::Point2f> pts2;
                std::vector<cv::Point2f> pts3;
                pts1.push_back(pA);
                pts1.push_back(x1);
                pts1.push_back(pB);

                pts2.push_back(x1);
                pts2.push_back(pE);
                pts2.push_back(x2);

                pts3.push_back(x2);
                pts3.push_back(pD);
                pts3.push_back(pC);

                return area(pts1, true) + area(pts2, true) + area(pts3, true);
            }
            else // (1)
            {
                std::vector<cv::Point2f> pts1;
                std::vector<cv::Point2f> pts2;
                pts1.push_back(pA);
                pts1.push_back(x1);
                pts1.push_back(pB);

                pts2.push_back(x1);
                pts2.push_back(pE);
                pts2.push_back(pC);
                pts2.push_back(pD);
                return area(pts1, true) + area(pts2, true);
            }
        }
        else if (intersect(pA, pE, pC, pD)) // AE intersects CD(3)
        {
            cv::Point2f x = intersection(pA, pE, pC, pD);
            std::vector<cv::Point2f> pts1;
            std::vector<cv::Point2f> pts2;
            pts1.push_back(pA);
            pts1.push_back(x);
            pts1.push_back(pC);
            pts1.push_back(pB);

            pts2.push_back(x);
            pts2.push_back(pE);
            pts2.push_back(pD);
            return area(pts1, true) + area(pts2, true);
        }
        else if (intersect(pD, pE, pA, pB)) // DE intersects AB(4)
        {
            cv::Point2f x = intersection(pD, pE, pA, pB);

            std::vector<cv::Point2f> pts1;
            std::vector<cv::Point2f> pts2;
            pts1.push_back(pA);
            pts1.push_back(pE);
            pts1.push_back(x);

            pts2.push_back(x);
            pts2.push_back(pD);
            pts2.push_back(pC);
            pts2.push_back(pB);

            return area(pts1, true) + area(pts2, true);
        }
        else if (intersect(pD, pE, pB, pC)) // DE intersects BC(5)
        {
            cv::Point2f x = intersection(pD, pE, pB, pC);
            std::vector<cv::Point2f> pts1;
            std::vector<cv::Point2f> pts2;
            pts1.push_back(pA);
            pts1.push_back(pE);
            pts1.push_back(x);
            pts1.push_back(pB);

            pts2.push_back(x);
            pts2.push_back(pD);
            pts2.push_back(pC);
            return area(pts1, true) + area(pts2, true);
        }
        else // no intersection(6)
        {
            std::vector<cv::Point2f> pts1;
            pts1.push_back(pA);
            pts1.push_back(pE);
            pts1.push_back(pD);
            pts1.push_back(pC);
            pts1.push_back(pB);
            return area(pts1, true);
        }
    }
}

int equalAreaLine(cv::Point2f A, cv::Point2f B, cv::Point2f C, cv::Point2f D, cv::Point2f *point1, cv::Point2f *point2)
{
    // Returns two points on the line of symmetrical areal displacement.
    // Replacing points B & C with any point E on this line
    // will result in a new line feature AED whose displacement from ABCD
    // is symmetrical, i.e. equal areal displacement on either side of the original line.
    // obtain line in standard form: aX + bY + C = 0

    float a = D.y - A.y;
    float b = A.x - D.x;
    float c = -A.x * B.y + (A.y - C.y) * B.x + (B.y - D.y) * C.x + C.y * D.x;

    cv::Point2f p1;
    cv::Point2f p2;
    // create two points on line near AD
    if (a == 0)
    {
        if (b == 0)
        {
            p1.x = -1;
            p1.y = -1;
            p2.x = -1;
            p2.y = -1;
            return -1;
        }
        else // b > a so use X from A and D
        {
            float liney = -c / b;
            p1 = cv::Point2f(A.x, liney);
            p2 = cv::Point2f(D.x, liney);
        }
    }
    else if (b == 0)
    {
        float linex = -c / a;
        p1 = cv::Point2f(linex, A.y);
        p2 = cv::Point2f(linex, D.y);
    }
    else
    {
        if (abs(a) > abs(b)) // use y
        {
            p1 = cv::Point2f((-c - b * A.y) / a, A.y);
            p2 = cv::Point2f((-c - b * D.y) / a, D.y);
        }
        else // use x
        {
            p1 = cv::Point2f(A.x, (-c - a * A.x) / b);
            p2 = cv::Point2f(D.x, (-c - a * D.x) / b);
        }
    }
    
    // create line as list of points(as list of coords)
    *point1 = p1;
    *point2 = p2;

    return 0;
}

int distance_pts(cv::Point2f a, cv::Point2f b, float *distance)
{
    // computes the distance between two points.
    // inputs should be lists of two coordinates.
    *distance = sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
    return  0;
}
    

int distance_pt_line(cv::Point2f p, cv::Point2f a, cv::Point2f b, float *distance)
{
    // Computes the distance from point p to the infinite line through ab
    std::vector<cv::Point2f> pts;
    pts.push_back(a);
    pts.push_back(b);
    pts.push_back(p);
    float trianglearea = abs(area(pts));
    float line_length;
    distance_pts(a, b, &line_length);
    if (line_length == 0)
    {
        distance_pts(p, a, distance);
    }
    else
    {
        *distance = 2.0f * trianglearea / line_length;
    }

    return 0;
}

bool X_Right_of_AB(cv::Point2f A, cv::Point2f B, cv::Point2f X, float tolerance = 0.0f)
{
    //Returns True if the input point or the first point on the input line
    //is on or right of the line from A to B.

    //A: first point on baseline (tuple of two floats)
    //B: second point on baseline
    //X: either a point (tuple of two floats) or a list of points

    // get first point from X if it is a line
    //if isinstance(X[0],list):
    //    X=X[0]
    // determine area of ABC
    float ax = A.x;
    float ay = A.y;
    float bx = B.x;
    float by = B.y;
    float xx = X.x;
    float xy = X.y;
    // minx = min(ax, bx, xx)
    float minx = 0;
    minx = std::min(ax, bx);
    minx = std::min(minx, xx);

    // miny = min(ay, by, xy)
    float miny = 0;
    miny = std::min(ay, by);
    miny = std::min(miny, xy);

    ax = ax - minx;
    bx = bx - minx;
    xx = xx - minx;
    ay = ay - miny;
    by = by - miny;
    xy = xy - miny;
    float area = ax * xy + bx * ay + xx * by - ax * by - bx * xy - xx * ay;
    // check precision - THIS MIGHT REALLY SLOW THINGS DOWN!
    // tolerance=0
    if (tolerance > 0)
    {
        float d = 0.0f;
        distance_pt_line(X, A, B, &d);
        if (d < tolerance)
        {
            //A = (dec(A[0]), dec(A[1]));
            //B = (dec(B[0]), dec(B[1]));
            //X = (dec(X[0]), dec(X[1]));
            ax, ay, bx, by, xx, xy = A.x, A.y, B.x, B.y, X.x, X.y;
            //minx, miny = min(ax, bx, xx), min(ay, by, xy);

            // minx = min(ax, bx, xx)
            float minx = 0;
            minx = std::min(ax, bx);
            minx = std::min(minx, xx);

            // miny = min(ay, by, xy)
            float miny = 0;
            miny = std::min(ay, by);
            miny = std::min(miny, xy);

            ax, bx, xx = ax - minx, bx - minx, xx - minx;
            ay, by, xy = ay - miny, by - miny, xy - miny;
            area = ax * xy + bx * ay + xx * by - ax * by - bx * xy - xx * ay;
        }
    }
    
    // if area is positive, C is right of AB
    if (area >= 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

int placement_AP_EAmin(cv::Point2f pA, cv::Point2f pB, cv::Point2f pC, cv::Point2f pD, cv::Point2f *pE, float *area, int *overlap_endpt)
{
    // Determines point E to replace BC that minimizes effective area of displacement while preserving area on each side.
    // Args:
    //     pA,pB,pC,pD: xy tuples representing sequence of four vertices
    //     m_func: function to calculate displacement metric
    // Returns:
    //     pE: location of new point E
    //     displacement: (half) area of displacement between ABCD and AED
    //     overlap_endpt: 0  (E is placed on AB)
    //                    1  (E is placed on CD)
    //                    -1 (E is not placed on AB or CD)
    // Simultaneously returns error metric (1/2 of displacement area).
    if (pA == pD)
    {
        *pE = cv::Point2f(-1, -1);
        *area = -1.0f;
        *overlap_endpt = -1;
        return 0;
    }
    // get symmetry line
    cv::Point2f S1;
    cv::Point2f S2;
    equalAreaLine(pA, pB, pC, pD, &S1, &S2);

    float distance = 0.0f;
    distance_pt_line(pA, S1, S2, &distance);
    // determine situation
    if ((S1.x == -1) && (S1.y == -1) && (S2.x == -1) && (S2.y == -1))
    {
        pE->x = -1;
        pE->y = -1;
        *area = -1;
        *overlap_endpt = -1;
    }
    else if (distance == 0)
    {
        pE->x = -1;
        pE->y = -1;
        *area = -1;
        *overlap_endpt = -1;
    }
    else
    {
        // default: intersect S with CD
        cv::Point2f C2 = pC;
        cv::Point2f D2 = pD;
        *overlap_endpt = 1;

        bool pApDpB_status = X_Right_of_AB(pA, pD, pB);
        bool pApDpC_status = X_Right_of_AB(pA, pD, pC);
        if (pApDpB_status == pApDpC_status)
        {
            float dist_pBpApD = 0.0f;
            float dist_pCpApD = 0.0f;
            distance_pt_line(pB, pA, pD, &dist_pBpApD);
            distance_pt_line(pC, pA, pD, &dist_pCpApD);
            float B_further_from_AD = dist_pBpApD - dist_pCpApD;
            if (B_further_from_AD > 0) // intersect S with AB
            {
                D2 = pA;
                C2 = pB;
                *overlap_endpt = 0;
            }
        }
        else
        {
            bool pApDpB_status = X_Right_of_AB(pA, pD, pB);
            bool pApDS_status = X_Right_of_AB(pA, pD, S1);
            if (pApDpB_status == pApDS_status) // intersect S with AB
            {
                D2 = pA;
                C2 = pB;
                *overlap_endpt = 0;
            }
        }

        // calculate intersection
        *pE = intersection(D2, C2, S1, S2);
        // handle case of colinear points
        if (pE->x == -1)
        {
            *pE = pC;
        }
        // calculate displacement area
        if (*overlap_endpt == 1) // intersect S with DC
        {
            *area = sideshift_area(pA, pB, pC, pD, *pE, *overlap_endpt);
        }
        // intersect S with AB
        else
        {
            *area = sideshift_area(pD, pC, pB, pA, *pE, *overlap_endpt);
        }
    }

    return 0;
}


int addToPriorityList(int A, int B, int C, int D, std::vector<cv::Point2f> XY, std::vector<PriorityList>& priorityList)
{
    //calculates the collapse of BC to E (using the placement function p_func)
    //and the associated error metric (using the metric function m_func)
    //and adds the collapse event to the priority list pl
    //in sorted ascending order of error.
    //Args:
    //    A,B,C,D: Four vertex sequence
    //    vtree: list of vertex relations
    //    pl: priority list
    //    p_func: placement function
    //    m_func: metric funciton
    //    do_topo_check: If true, will avoid self-intersections (less efficient)
    //Adds a tuple to the priority list with the following values:
    //    (displacement,A,B,C,D,pE,overlap_endpt)

    // XY = vtree[1]
    // 
    // 增加pE限制
    float pE_xmin = FLT_MAX;
    float pE_xmax = 0.0f;
    float pE_ymin = FLT_MAX;
    float pE_ymax = 0.0f;
    for (size_t i = 0; i < XY.size(); i++)
    {
        if (XY[i].x > pE_xmax)
        {
            pE_xmax = XY[i].x;
        }
        if (XY[i].x < pE_xmin)
        {
            pE_xmin = XY[i].x;
        }
        if (XY[i].y > pE_ymax)
        {
            pE_ymax = XY[i].y;
        }
        if (XY[i].y < pE_ymin)
        {
            pE_ymin = XY[i].y;
        }
    }

    int pE_border_val = 10;
    
    // make sure all vertices are current
    if ((A != -1) && (B != -1) && (C != -1) && (D != -1))
    {
        // compute collapse point E and associated displacement error
        cv::Point2f pE = cv::Point2f(-1, -1);
        float displacement;
        int overlap_endpt;
        placement_AP_EAmin(XY[A], XY[B], XY[C], XY[D], &pE, &displacement, &overlap_endpt);

        // make sure collapse point is valid
        if ((pE.x != -1) && (pE.y != -1) 
            && (pE.x > pE_xmin - pE_border_val)
            && (pE.x < pE_xmax + pE_border_val)
            && (pE.y > pE_ymin - pE_border_val)
            && (pE.y < pE_ymax + pE_border_val))
        {
            // get index of overlapping endpoint
            if (overlap_endpt == 0)
            {
                overlap_endpt = A;
            }
            if (overlap_endpt == 1)
            {
                overlap_endpt = D;
            }

            PriorityList pl;
            pl.displacement = displacement;
            pl.A = A;
            pl.B = B;
            pl.C = C;
            pl.D = D;
            pl.pE = pE;
            pl.overlap_endpt = overlap_endpt;

            priorityList.push_back(pl);
        }
        return 0;
    }



}

// 定义比较函数
bool compareElements(const PriorityList& a, const PriorityList& b)
{
    if (a.displacement != b.displacement)
    {
        return a.displacement < b.displacement;
    }
    return a.A < b.A;
}



int init_priority_list(std::vector<cv::Point2f> XY, std::vector<PriorityList>& priorityList)
{
    for (int k = 0; k < XY.size() - 3; k++)
    {
        addToPriorityList(k, k + 1, k + 2, k + 3, XY, priorityList);
    }

    std::sort(priorityList.begin(), priorityList.end(), compareElements);

    return 0;
}

int collapse(int B, int C, int E, cv::Point2f thisxy, float thiserror, VertexRelation *vertex_tree)
{
    //records the collapse of BC to E
    //by adding a new row to the vertex relation table
    //representing the new Steiner point E
    //and adjusting the relations of B & C

    // obtain relations as lists
    //ID = vtree[0]
    //XY = vtree[1]
    //error = vtree[2]
    //parent = vtree[3]
    //LC = vtree[4]
    //RC = vtree[5]
    //LS = vtree[6]
    //RS = vtree[7]
    //current = vtree[8]

    vertex_tree->ID.push_back(E);
    vertex_tree->XY.push_back(thisxy);
    vertex_tree->error.push_back(thiserror);
    vertex_tree->parent.push_back(-1);
    vertex_tree->LC.push_back(B);
    vertex_tree->RC.push_back(C);
    vertex_tree->LS.push_back(vertex_tree->LS[B]);
    vertex_tree->RS.push_back(vertex_tree->RS[C]);
    vertex_tree->current.push_back(true);

    // add new row
    // id, xy, error, parent, LC, RC, LS, RS, current
    //ID.append(E)
    //XY.append(thisxy)
    //error.append(thiserror)
    //parent.append(-1)
    //LC.append(B)
    //RC.append(C)
    //LS.append(LS[B])
    //RS.append(RS[C])
    //current.append(True)

    // reset sibling relations to E
    int ls_index = vertex_tree->LS[B];
    if (ls_index < 0)
    {
        ls_index = vertex_tree->LS.size() + ls_index;
    }

    int rs_index = vertex_tree->RS[C];
    if (rs_index < 0)
    {
        rs_index = vertex_tree->RS.size() + rs_index;
    }

    //vertex_tree->RS[vertex_tree->LS[B]] = E;
    //vertex_tree->LS[vertex_tree->RS[C]] = E;
    // 
    vertex_tree->RS[ls_index] = E;
    vertex_tree->LS[rs_index] = E;
    // reset properties of B
    vertex_tree->LS[B] = -1;
    vertex_tree->RS[B] = C;
    vertex_tree->parent[B] = E;
    vertex_tree->current[B] = false;
    // reset properties of C
    vertex_tree->LS[C] = B;
    vertex_tree->RS[C] = -1;
    vertex_tree->parent[C] = E;
    vertex_tree->current[C] = false;

    return 0;
}



std::vector<cv::Point2f> apsc(std::vector<cv::Point2f> pts, int min_pts_num)
{
    std::vector<cv::Point2f> upsample_pts;
    // 一、simplificationTable
    VertexRelation vertex_tree;
    initVertexTree(pts, &vertex_tree);

    // build priority list of collapse events
    //priority_list
    std::vector<PriorityList> priorityList;
    init_priority_list(vertex_tree.XY, priorityList);
    int numE = (int)(pts.size() - 1);
    float maxerror = 0.0f;
    int collapsed = 0;
    int report_interval = 2000;

    while (priorityList.size() > 0)
    {
        //error, A, B, C, D, pE, overlap_endpt = priorityList.pop_back();
        PriorityList first_priority = priorityList.front();  // 获取第一个元素
        priorityList.erase(priorityList.begin());            // 移除第一个元素
        float error = first_priority.displacement;
        int A = first_priority.A;
        int B = first_priority.B;
        int C = first_priority.C;
        int D = first_priority.D;
        cv::Point2f pE = first_priority.pE;
        if (vertex_tree.current[A] && vertex_tree.current[B] && vertex_tree.current[C] && vertex_tree.current[D])
        {
            //bool topo_error = false;
            numE += 1;
            maxerror = std::max(error, maxerror);
            // perform collapse
            collapse(B, C, numE, pE, maxerror, &vertex_tree);

            collapsed = collapsed + 1;
            if (collapsed % report_interval == 0)
            {
                std::cout << "collapsed" << collapsed << "of" << pts.size() << std::endl;
            }

            int ls_index = vertex_tree.LS[A];
            if (ls_index < 0)
            {
                ls_index = vertex_tree.LS.size() + ls_index;
            }
            addToPriorityList(vertex_tree.LS[ls_index], vertex_tree.LS[A], A, numE, vertex_tree.XY, priorityList);
            addToPriorityList(vertex_tree.LS[A], A, numE, D, vertex_tree.XY, priorityList);
            addToPriorityList(A, numE, D, vertex_tree.RS[D], vertex_tree.XY, priorityList);

            int rs_index = vertex_tree.RS[D];
            if (rs_index < 0)
            {
                rs_index = vertex_tree.RS.size() + rs_index;
            }
            addToPriorityList(numE, D, vertex_tree.RS[D], vertex_tree.RS[rs_index], vertex_tree.XY, priorityList);
            std::sort(priorityList.begin(), priorityList.end(), compareElements);
        }
    }
    
    // 二、simplifiedLine
    // determine original number of points,
    // maximum ID to allow based on displacement error criterion
    int max_ID = (int)(vertex_tree.ID.size()) - 1;
    int orig_count = 0;
    int max_disp = -1;
    for (int i = 0; i < vertex_tree.ID.size(); i++)
    {
        if (vertex_tree.LC[i] == -1)
        {
            orig_count = orig_count + 1;
        }
        if (max_disp >= 0)
        {
            if (vertex_tree.error[i] <= max_disp)
            {
                max_ID = i;
            }
        }
    }

    // determine max ID to allow based on minimum number of points
    if (min_pts_num > -1)
    {
        int max_ID_by_min_pts = orig_count - 1 + orig_count - min_pts_num;
        if (max_disp > -1)
        {
            max_ID = std::max(max_ID, max_ID_by_min_pts);
        }
        else
        {
            max_ID = max_ID_by_min_pts;
        }
    }
    // initialize output as empty list
    // start at beginning
    int v = 0;
    // work to end
    while (v > -1)
    {
        if (v <= max_ID) // vertex is good
        {
            upsample_pts.push_back(vertex_tree.XY[v]); // add to result
            // if we can't move right, move up
            while ((v > -1) && (vertex_tree.RS[v] == -1))
            {
                v = vertex_tree.parent[v];
            }
            // move right if we're not at the end already
            if (v > -1)
            {
                v = vertex_tree.RS[v];
            }
        }
        else  // vertex needs to be expanded; move to left child
        {
            v = vertex_tree.LC[v];
        }
    }

    return upsample_pts;
}




